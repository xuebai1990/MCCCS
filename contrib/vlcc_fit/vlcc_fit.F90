!> \file Fitting program for vapor-liquid coexitence curve
!>
!> This program determines the critical temperature and density by
!> fitting coexistence densities--temperature data to the law of
!> rectilinear parameters and the scaling law for density. It then
!> estimates critical pressure by fitting pressure-temperature data
!> to the Antoine equation. The program also calculates the boiling
!> point temperature by fitting pressure-temperature data to the
!> Clausius-Clapeyron equation in two steps: first using data at all
!> input temperatures, second using only two closet temperatures to
!> boiling point estimated in the first step.
!>
!> Program reads from standard input, and writes to standard output
!> Run as "vlcc_fit < input | tee output" to have output display on
!> screen as well as write to the output file.
!>
!> In the input file, the first part deals with critical properties:
!> NTEMP: number of data entries for the fitting
!> Beta: critical exponent
!> NITER: max number of iterations for Antoine fit
!>
!> The program also makes xmgrace plots:
!> TMIN: lowest temperature to start generate data points of fitting result
!> NPOINT: number of points to generate
!>
!> The next part deals with boiling points:
!> NPRESS: number of data entries for the fitting
!> STP: external pressure at which to calculate the boiling point
!>
!> Sample input:
!> =================================================
!> Title
!> NDATA  Beta    NITER
!> 4      0.326   100
!>  T      RHOG      ERROR      RHOL      ERROR      PRESSURE      ERROR (standard deviation or standard error)
!> 375   0.00448    0.00006    0.6975    0.0005      369.4         4.3
!> ...
!> ...
!> TMIN NPOINT for VLCC plots
!> 450  200
!> NDATA  STP
!> 3      101.325
!> T    PRESSURE    ERROR
!> 325    54.0       0.9
!> ...
!> ...
!> TMIN TMAX NPOINT for Clausius-Clapeyron plots
!> 300  400  200
!> =================================================

PROGRAM VLCC_FIT
  IMPLICIT NONE
  INTEGER,PARAMETER::INP=5,IOUT=6,IERR=0
  REAL(KIND=8)::A,B

  CALL CRITICAL(INP,IOUT,A,B)
  CALL PLOTVLCC(INP,IERR,A,B)
  CALL BOILINGPOINT(INP,IOUT,A,B)
  CALL PLOTCC(INP,IERR,A,B)

END PROGRAM VLCC_FIT

!> \brief Make VLCC plots
SUBROUTINE PLOTVLCC(INP,IERR,A,B)
  IMPLICIT NONE
  INTEGER,PARAMETER::NDATAMAX=500
  INTEGER,INTENT(IN)::INP,IERR
  REAL(KIND=8),INTENT(IN)::A,B
  
  INTEGER::NS,I,NPOINT
  REAL(KIND=8)::T,RHOG,ERRG,RHOL,ERRL,PRESS,SIGPRESS,RHOA,ERRA&
   ,RHOS,ERRS,TC,SIGTC,RHOC,SIGRHO,BETA,TMIN,DELT,TT,RHO1,RHO2,RHO3,SUMRHO,DIFRHO

  COMMON/VLCC/T(NDATAMAX),RHOG(NDATAMAX),ERRG(NDATAMAX)&
   ,RHOL(NDATAMAX),ERRL(NDATAMAX),PRESS(NDATAMAX),SIGPRESS(NDATAMAX)&
   ,RHOA(NDATAMAX),ERRA(NDATAMAX),RHOS(NDATAMAX),ERRS(NDATAMAX),TC,SIGTC,RHOC,SIGRHO,BETA,NS

  WRITE(IERR,*) ' PLOT RESULTS '

!
! --- plot simulation results
!
  OPEN(UNIT=21,ACCESS='sequential',ACTION='write',FILE='vlcc-liq-sim.xvg',FORM='formatted',STATUS='unknown')
  OPEN(UNIT=22,ACCESS='sequential',ACTION='write',FILE='vlcc-vap-sim.xvg',FORM='formatted',STATUS='unknown')
  OPEN(UNIT=23,ACCESS='sequential',ACTION='write',FILE='vlcc-rect-sim.xvg',FORM='formatted',STATUS='unknown')
  write(21,"(A,/,A,/,A,/,A)") '@title "liq-sim"','@xaxis label "\3\xr\3"','@yaxis label "\3T"','@TYPE xydx'
  write(22,"(A,/,A,/,A,/,A)") '@title "vap-sim"','@xaxis label "\3\xr\3"','@yaxis label "\3T"','@TYPE xydx'
  write(23,"(A,/,A,/,A,/,A)") '@title "rect-sim"','@xaxis label "\3\xr\3"','@yaxis label "\3T"','@TYPE xydx'
  DO I=1,NS
     write(21,*) rhog(i),t(i),errg(i)
     write(22,*) rhol(i),t(i),errl(i)
     write(23,*) rhoa(i),t(i),erra(i)
  END DO
  CLOSE(21)
  CLOSE(22)
  CLOSE(23)

!
! === PLOT TC AND RHOC:
!
  OPEN(UNIT=24,ACCESS='sequential',ACTION='write',FILE='vlcc-crit-fit.xvg',FORM='formatted',STATUS='unknown')
  write(24,"(A,/,A,/,A,/,A)") '@title "crit-fit"','@xaxis label "\3\xr\3"','@yaxis label "\3T"','@TYPE xydxdy'
  write(24,"(4(1X,G16.9))") RHOC,TC,SIGRHO,SIGTC
  CLOSE(24)

!
! === CALCULATE SCALING CURVE AND RECLINERAR CURVE
!
  READ(INP,*)
  READ(INP,*) TMIN,NPOINT
  DELT=(TC-TMIN)/NPOINT
  OPEN(UNIT=25,ACCESS='sequential',ACTION='write',FILE='vlcc-liq-fit.xvg',FORM='formatted',STATUS='unknown')
  OPEN(UNIT=26,ACCESS='sequential',ACTION='write',FILE='vlcc-vap-fit.xvg',FORM='formatted',STATUS='unknown')
  OPEN(UNIT=27,ACCESS='sequential',ACTION='write',FILE='vlcc-rect-fit.xvg',FORM='formatted',STATUS='unknown')
  write(25,"(A,/,A,/,A,/,A)") '@title "liq-fit"','@xaxis label "\3\xr\3"','@yaxis label "\3T"','@TYPE xy'
  write(26,"(A,/,A,/,A,/,A)") '@title "vap-fit"','@xaxis label "\3\xr\3"','@yaxis label "\3T"','@TYPE xy'
  write(27,"(A,/,A,/,A,/,A)") '@title "rect-fit"','@xaxis label "\3\xr\3"','@yaxis label "\3T"','@TYPE xy'
  DO I=0,NPOINT-1
     TT=TMIN+DELT*I
     RHO3=RHOC+A*(TT-TC)
     SUMRHO=2*RHO3
     DIFRHO=B*ABS(TT-TC)**BETA
     RHO2=0.5*(SUMRHO+DIFRHO)
     RHO1=0.5*(SUMRHO-DIFRHO)
     write(25,*)RHO1,TT
     write(26,*)RHO2,TT
     write(27,*)RHO3,TT
  END DO
  CLOSE(25)
  CLOSE(26)
  CLOSE(27)

  RETURN
END SUBROUTINE PLOTVLCC

!> \brief Make Clausius-Clapeyron plots
SUBROUTINE PLOTCC(INP,IERR,A,B)
  IMPLICIT NONE
  INTEGER,PARAMETER::NDATAMAX=500
  INTEGER,INTENT(IN)::INP,IERR
  REAL(KIND=8),INTENT(IN)::A,B
  
  INTEGER::NS,I,NPOINT
  REAL(KIND=8)::INVT,LOGP,SIGLOGP,TMIN,TMAX,DELT,TT,PP

  COMMON/CC/INVT(NDATAMAX),LOGP(NDATAMAX),SIGLOGP(NDATAMAX),NS

  WRITE(IERR,*) ' PLOT RESULTS '

!
! --- plot simulation results
!
  OPEN(UNIT=21,ACCESS='sequential',ACTION='write',FILE='cc-sim.xvg',FORM='formatted',STATUS='unknown')
  write(21,"(A,/,A,/,A,/,A)") '@title "ln \3p\2 - 1/\3T"','@xaxis label "\21/\3T"','@yaxis label "\2ln \3p"','@TYPE xydy'
  DO I=1,NS
     write(21,*) INVT(i),LOGP(i),SIGLOGP(i)
  END DO
  CLOSE(21)

!
! === CALCULATE CLAUSIUS-CLAPEYRON CURVE
!
  READ(INP,*)
  READ(INP,*) TMIN,TMAX,NPOINT
  DELT=(TMAX-TMIN)/NPOINT
  OPEN(UNIT=22,ACCESS='sequential',ACTION='write',FILE='cc-fit.xvg',FORM='formatted',STATUS='unknown')
  write(22,"(A,/,A,/,A,/,A)") '@title "ln \3p\2 - 1/\3T"','@xaxis label "\21/\3T"','@yaxis label "\2ln \3p"','@TYPE xy'
  DO I=0,NPOINT-1
     TT=1.0d0/(TMIN+DELT*I)
     PP=A+B*TT
     write(22,*) TT,PP
  END DO
  CLOSE(22)

  RETURN
END SUBROUTINE PLOTCC

!> \brief Calculate critical temperature, density, and pressure
SUBROUTINE CRITICAL(INP,IOUT,A,B)
  IMPLICIT NONE
  INTEGER,PARAMETER::NDATAMAX=500,NCA=3
  REAL(KIND=8),PARAMETER::EPS=1D-9,BETASTD=0.326D0,TLOWPERC=0.9D0
  INTEGER,INTENT(IN)::INP,IOUT
  REAL(KIND=8),INTENT(OUT)::A,B

  INTEGER::NS,NITER,I
  REAL(KIND=8)::T,RHOG,ERRG,RHOL,ERRL,PRESS,SIGPRESS,RHOA,ERRA&
   ,RHOS,ERRS,TC,SIGTC,RHOC,SIGRHO,BETA,PC,SIGPC,ANTOINECOEFF(NCA)&
   ,AR,BR,SIGAR,SIGBR,CHI2R,QR,COVR,AS,BS,SIGAS,SIGBS,CHI2S,QS&
   ,COVS,SIGA,SIGB,TLOW,SUMRHO,DIFRHO,GAS,RLIQ,P,DGAS,DRLIQ,DPRESS&
   ,SDGAS,SDRLIQ,SDPRESS,DUM(NCA),DELRHO

  COMMON/VLCC/T(NDATAMAX),RHOG(NDATAMAX),ERRG(NDATAMAX)&
   ,RHOL(NDATAMAX),ERRL(NDATAMAX),PRESS(NDATAMAX),SIGPRESS(NDATAMAX)&
   ,RHOA(NDATAMAX),ERRA(NDATAMAX),RHOS(NDATAMAX),ERRS(NDATAMAX),TC,SIGTC,RHOC,SIGRHO,BETA,NS

!
! === READ INPUT DATA:
!
  READ(INP,*)
  READ(INP,*)
  READ(INP,*) NS,BETA,NITER
  READ(INP,*)
  WRITE(IOUT,"('--- GEMC: CRITICAL PROPERTIES ---',/,&
   ' CRITICAL EXPONENT (BETA):',F6.3)") BETA
  IF (abs(BETA-BETASTD).gt.EPS) WRITE(IOUT,*)&
   ' *NON-STANDARD CRITICAL EXPONENT DETECTED (DEFAULT=0.326)*'
  WRITE(IOUT,"(/,/,' SIMULATION DATA: ',/,3X,'T',19X,'RHO_G',21X,'RHO_L',18X,'RHO_G+RHO_L',11X,&
   '(RHO_L-RHO_G)**BETA',12X,'PRESSURE')")

  DO I=1,NS
     READ(INP,*) T(I),RHOG(I),ERRG(I),RHOL(I),ERRL(I)&
      ,PRESS(I),SIGPRESS(I)
     rhoa(i) = ( rhog(i) + rhol(i) ) / 2.0D0
!!??     erra(i) = ( errg(i) + errl(i) ) / 1.41
     erra(i) = ( errg(i) + errl(i) ) / 2.0D0
     DELRHO=(RHOL(I)-RHOG(I))
     RHOS(I)=DELRHO**(1.0D0/BETA)
     ERRS(I)=ABS((RHOS(I)/DELRHO/BETA)*(errg(i) + errl(i)))
     WRITE(IOUT,"(G12.5,5(1X,G10.3,' +/- ',G10.3))") T(I),RHOG(I)&
      ,ERRG(I),RHOL(I),ERRL(I),RHOA(I),ERRA(I),RHOS(I),ERRS(I)&
      ,PRESS(I),SIGPRESS(I)
  END DO

  IF (NS.lt.NCA) RETURN

!
! === FIT RESULTS RECTANGULAR RULE:
!
  CALL FIT(T,RHOA,NS,ERRA,1,AR,BR,SIGAR,SIGBR,CHI2R,QR,COVR)
!  CHECK
  WRITE(IOUT,"(/,'  CHECK FITTING RESULTS RECLINEAR LAW:',/,&
   '     NUMBER OF SAMPLES: ',I9  ,/,&
   '               CHI2/NS: ',F9.6,/,&
   '                     Q: ',F9.6,'  > 0.001 ? ')")&
   NS,CHI2R/NS,QR
!
! === FIT RESULTS TO THE SCALING LAW:
!
  CALL FIT(T,RHOS,NS,ERRS,1,AS,BS,SIGAS,SIGBS,CHI2S,QS,COVS)
!  CHECK
  WRITE(IOUT,"(/,'  CHECK FITTING RESULTS SCALING LAW:',/,&
   '     NUMBER OF SAMPLES: ',I9  ,/,&
   '               CHI2/NS: ',F9.6,/,&
   '                     Q: ',F9.6,'  > 0.001 ? ')")&
   NS,CHI2S/NS,QS
!
! === CALCULATE TC, RHOC ,A B:
!
  A=BR
  SIGA=SIGBR

  B=(-BS)**BETA
  SIGB=(BETA*(-BS)**(BETA-1))*SIGBS !!?? This is correct. The book has typos.

  TC=-AS/BS
  SIGTC=SQRT((SIGAS/BS)**2-2.0D0*COVS*AS/BS**3+(AS*SIGBS/(BS*BS))**2)

  RHOC=AR+BR*TC
  SIGRHO=SQRT(SIGAR**2+2*TC*COVR+(BR*SIGTC)**2+(TC*SIGBR)**2)

  CALL PRESSURE(IOUT,T,PRESS,SIGPRESS,NS,TC,SIGTC&
   ,PC,SIGPC,ANTOINECOEFF,NITER)

  WRITE(IOUT,"('  =========== RESULTS:',/,&
   '     NUMBER OF SAMPLES: ',I13,/,&
   '                     A: ',G13.6,' +/- ',G13.6,/,&
   '                     B: ',G13.6,' +/- ',G13.6,/,&
   '   Antoine Coefficient: (',2(G13.6,', '),G13.6,')',/,&
   '                    TC: ',G13.6,' +/- ',G13.6,/,&
   '                  RHOC: ',G13.6,' +/- ',G13.6,/,&
   '                    PC: ',G13.6,' +/- ',G13.6,/)")&
   NS,A,SIGA,B,SIGB,ANTOINECOEFF,TC,SIGTC,RHOC,SIGRHO,PC,SIGPC

  TLOW=TC*TLOWPERC
  DO I = 1,NS
     IF (T(I).lt.TLOW) THEN
        WRITE(IOUT,"('***Temperature range too large (potentially &
         &outside of the validity of the scaling law)***',/,&
         ' Consider refitting with temperatures closer to Tc.',/)")
        EXIT
     END IF
  END DO

  WRITE(IOUT,"(1X,'ISAMP',3X,'TEMP',10X,'RHOG',7X,'RHOG-p',7X,'DIF',8X,'RHOL',7X,'RHOL-p',8X,'DIF',9X,'PRESS',6X,'PRESS-p',6X,'DIF')")

  SDGAS=0.0D0
  SDRLIQ=0.0D0
  SDPRESS=0.0D0
  DO I=1,NS
     SUMRHO=2*(RHOC+A*(T(I)-TC))
     DIFRHO=B*ABS(T(I)-TC)**BETA
     RLIQ=0.5D0*(SUMRHO+DIFRHO)
     GAS=0.5D0*(SUMRHO-DIFRHO)
     CALL ANTOINE(T(I),ANTOINECOEFF,P,DUM,NCA)
     P=10.0D0**P
     DGAS=((GAS-RHOG(I))/RHOG(I))**2
     DRLIQ=((RLIQ-RHOL(I))/RHOL(I))**2
     DPRESS=((P-PRESS(I))/PRESS(I))**2
     SDGAS=SDGAS+DGAS
     SDRLIQ=SDRLIQ+DRLIQ
     SDPRESS=SDPRESS+DPRESS
     WRITE(IOUT,"(I4,2X,G12.5,9(2X,G10.3))") I,T(I),RHOG(I),GAS&
      ,DGAS,RHOL(I),RLIQ,DRLIQ,PRESS(I),P,DPRESS
  END DO

  SDGAS=SDGAS/NS
  SDRLIQ=SDRLIQ/NS
  SDPRESS=SDPRESS/NS
  WRITE(IOUT,"(/,'  MEAN SQUARE ERROR:',/,&
   '      GAS:',G16.9,/,&
   '  RLIQUID:',G16.9,/,&
   ' PRESSURE:',G16.9,/)") SDGAS,SDRLIQ,SDPRESS

  RETURN
END SUBROUTINE CRITICAL

!> \brief Calculate critical pressure from Antoine fit
SUBROUTINE PRESSURE(IOUT,TEMP,PRESS,SIGPRESS,NS&
 ,TC,SIGTC,PC,SIGPC,A,NITER)
  IMPLICIT NONE
  INTEGER,PARAMETER::ITSTMAX=4,NCA=3
  REAL(KIND=8),PARAMETER::CHI2EPS=0.1D0,INVLN10=1.0D0/LOG(10.0D0)
  INTEGER,INTENT(IN)::IOUT,NS,NITER
  REAL(KIND=8),INTENT(IN)::TEMP(NS),PRESS(NS),SIGPRESS(NS),TC,SIGTC
  REAL(KIND=8),INTENT(OUT)::PC,SIGPC,A(NCA)

  INTEGER::I,ITST,IA(NCA)
  REAL(KIND=8)::ALAMDA,CHI2,Q,COV,CHI2O,INVT(NS),LOGP(NS)&
   ,SIGLOGP(NS),COVAR(NCA,NCA),ALPHA(NCA,NCA),DUM(NCA)
  EXTERNAL ANTOINE

  DO I=1,NS
     INVT(I)=-1.0D0/TEMP(I)
     LOGP(I)=LOG10(PRESS(I))
     SIGLOGP(I)=INVLN10*(SIGPRESS(I)/PRESS(I))
  END DO

! Fit to August equation using weighted least squares,
! which are later used as initial guess for Antoine fit
  CALL FIT(INVT,LOGP,NS,SIGLOGP,1,A(1),A(2)&
   ,COVAR(1,1),COVAR(2,2),CHI2,Q,COV)

! Fit data to Antoine equation using Lavenberg-Marquardt method
! with Clausius-Clapeyron parameters as initial guess
  ALAMDA=-1.0D0
  A(3)=0.0D0
  IA=(/1,1,1/)
  ITST=0
  CHI2O=CHI2
  DO I=1,NITER
     CALL MRQMIN(TEMP,LOGP,SIGLOGP,NS,A,IA,NCA&
      ,COVAR,ALPHA,NCA,CHI2,ANTOINE,ALAMDA)
     IF (CHI2.gt.CHI2O) THEN
        ITST=0
     ELSE IF (abs(CHI2-CHI2O).lt.CHI2EPS) THEN
        ITST=ITST+1
     END IF
     IF (ITST.ge.ITSTMAX) EXIT
     CHI2O=CHI2
  END DO
  ALAMDA=0.0D0
  CALL MRQMIN(TEMP,LOGP,SIGLOGP,NS,A,IA,NCA&
   ,COVAR,ALPHA,NCA,CHI2,ANTOINE,ALAMDA)

!  CHECK
  WRITE(IOUT,"(/,'  CHECK FITTING RESULTS ANTOINE:',/,&
   '     NUMBER OF SAMPLES: ',I9  ,/,&
   '               CHI2/NS: ',F9.6,/,/)") NS,CHI2/NS

  CALL ANTOINE(TC,A,PC,DUM,NCA)
  PC=10.0D0**PC
  SIGPC=SQRT(COVAR(1,1)+COVAR(2,2)/(TC+A(3))**2&
   +COVAR(3,3)*(A(2)/(TC+A(3))**2)**2+(SIGTC*A(2)/(TC+A(3))**2)**2&
   +2*(-COVAR(1,2)/(TC+A(3))+COVAR(1,3)*A(2)/(TC+A(3))**2&
   -COVAR(2,3)*A(2)/(TC+A(3))**3))
  SIGPC=SIGPC*PC/INVLN10

  RETURN
END SUBROUTINE PRESSURE

!> \brief Calculates Antoine equation and derivatives
!>
!> Return LOGP = A(1) - A(2) / ( T + A(3) )
SUBROUTINE ANTOINE(T,A,LOGP,DLOGPDT,NP)
  IMPLICIT NONE
  INTEGER,PARAMETER::NCA=3
  INTEGER,INTENT(IN)::NP
  REAL(KIND=8),INTENT(IN)::T,A(NP)
  REAL(KIND=8),INTENT(OUT)::LOGP,DLOGPDT(NP)

  IF (NP.ne.NCA) STOP 'Antoine equation must have 3 parameters'
  LOGP=A(1)-A(2)/(T+A(3))
  DLOGPDT(1)=1.0D0
  DLOGPDT(2)=-1.0D0/(T+A(3))
  DLOGPDT(3)=A(2)/(T+A(3))**2

  RETURN
END SUBROUTINE ANTOINE

!> \brief Calculate boiling point temperature at given pressure
SUBROUTINE BOILINGPOINT(INP,IOUT,A,B)
  IMPLICIT NONE
  INTEGER,PARAMETER::NDATAMAX=500,NCA=2
  INTEGER,INTENT(IN)::INP,IOUT
  REAL(KIND=8),INTENT(OUT)::A,B

  INTEGER::NS,I,J
  REAL(KIND=8)::INVT,LOGP,SIGLOGP,BP,STP,TEMP(NDATAMAX)&
   ,PRESS(NDATAMAX),SIGPRESS,SIGA,SIGB,CHI2,Q,COV,LOGSTP&
   ,SDPRESS,P,SIGBP

  COMMON/CC/INVT(NDATAMAX),LOGP(NDATAMAX),SIGLOGP(NDATAMAX),NS

  READ(INP,*)
  READ(INP,*) NS,STP
  READ(INP,*)
  WRITE(IOUT,"(/,'--- GEMC: BOILING POINT ---',/,&
   ' STANDARD STATE PRESSURE:',G13.6,/,/,/,&
   ' SIMULATION DATA: ',/,3X,'T',12X,'INVT',12X,'PRESSURE',21X,'LOGP')") STP

  do i = 1,NS
     read(INP,*) TEMP(i),PRESS(I),SIGPRESS
     INVT(i) = 1/TEMP(i)
     LOGP(i) = LOG(PRESS(I))
!!??     SIGPRESS = SIGPRESS/press
     SIGLOGP(i) = SIGPRESS/PRESS(I)
!!??     SIGLOGP(i) = SIGPRESS
     write(IOUT,"(G12.5,G10.3,2(2X,G10.3,' +/- ',G10.3))") TEMP(i),INVT(i)&
      ,PRESS(I),SIGPRESS,LOGP(i),SIGLOGP(i)
  end do

  IF (NS.lt.NCA) RETURN

! Fit to Clausius-Clapeyron equation using weighted least squares
  CALL FIT(INVT,LOGP,NS,SIGLOGP,1,A,B,SIGA,SIGB,CHI2,Q,COV)

! ESTIMATED BOILING POINT USING ALL DATA
  LOGSTP=LOG(STP)
  BP = B/(LOGSTP-A)
  SIGBP=SQRT((SIGB/(A-LOGSTP))**2+(B*SIGA/(A-LOGSTP)**2)**2&
   -2*COV*B/(A-LOGSTP)**3)

  write(IOUT,"(/,' Clausius-Clapeyron Results:',/,&
   '     NUMBER OF SAMPLES: ',I13  ,/,&
   '               CHI2/NS: ',F13.9,/,&
   '                     Q: ',F13.9,'  > 0.001 ? ',/,&
   '                     A: ',G13.6,' +/- ',G13.6,/,&
   '                     B: ',G13.6,' +/- ',G13.6,/,&
   '                    BP: ',G13.6,' +/- ',G13.6,/)")&
   NS,CHI2/NS,Q,A,SIGA,B,SIGB,BP,SIGBP

  WRITE(IOUT,"(1X,'ISAMP',2X,'TEMP',4X,'PRESSURE',11X,'PRESSURE-p',10X,'DIF')")

  SDPRESS=0.0D0
  DO I=1,NS
     P=exp(A+B*INVT(I))
     SIGPRESS=((P-PRESS(I))/PRESS(I))**2
     SDPRESS=SDPRESS+SIGPRESS
     WRITE(IOUT,"(I4,2X,F6.1,3(2X,G16.9))") I,TEMP(I)&
      ,PRESS(I),P,SIGPRESS
  END DO

  SDPRESS=SDPRESS/NS
  WRITE(IOUT,"(/,'  MEAN SQUARE ERROR:',/,&
   ' PRESSURE:',G16.9,/)") SDPRESS

  WRITE(IOUT,"(1X,'NDATA',3X,'CHI2/NS',8X,'Q',23X,'A',31X,'B',31X,'BP')")
  write(IOUT,"(I4,1X,G13.6,1X,G13.6,3(1X,G13.6,' +/- ',G13.6))")&
   NS,CHI2/J,Q,A,SIGA,B,SIGB,BP,SIGBP

  DO J=NS-1,2,-1
     DO I=1,J
        IF (ABS(TEMP(I)-BP)>ABS(TEMP(J+1)-BP)) THEN
           CALL SWAP(TEMP(I),TEMP(J+1))
           CALL SWAP(INVT(I),INVT(J+1))
           CALL SWAP(LOGP(I),LOGP(J+1))
           CALL SWAP(SIGLOGP(I),SIGLOGP(J+1))
        END IF
     END DO

     CALL FIT(INVT,LOGP,J,SIGLOGP,1,A,B,SIGA,SIGB,CHI2,Q,COV)
      
     ! BOILING POINT AT GIVEN PRESSURE
     BP = B/(LOGSTP-A)
     SIGBP=SQRT((SIGB/(A-LOGSTP))**2+(B*SIGA/(A-LOGSTP)**2)**2&
      -2*COV*B/(A-LOGSTP)**3)

     write(IOUT,"(I4,1X,G13.6,1X,G13.6,3(1X,G13.6,' +/- ',G13.6))")&
      J,CHI2/J,Q,A,SIGA,B,SIGB,BP,SIGBP
     
  END DO
  
  RETURN
END SUBROUTINE BOILINGPOINT

SUBROUTINE SWAP(X,Y)
  IMPLICIT NONE
  REAL(KIND=8),INTENT(INOUT)::X,Y

  REAL(KIND=8)::TMP

  TMP=X
  X=Y
  Y=TMP

  RETURN
END SUBROUTINE SWAP

!> FROM NUMERICAL RECICPE, ADDED COV(A,B)
      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q,cov)
      IMPLICIT NONE
      INTEGER mwt,ndata
      REAL*8 a,b,chi2,q,cov,siga,sigb,sig(ndata),x(ndata),y(ndata)
!U    USES gammq
      INTEGER i
      REAL*8 sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      sx=0.
      sy=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
        ss=0.
        do 11 i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      COV=-SX/(SS*ST2)
      chi2=0.
      q=1.
      if(mwt.eq.0) then
        do 15 i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do 16 i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
16      continue
        if(ndata.gt.2) q=gammq(0.5d0*(ndata-2),0.5*chi2)
      endif
      return
      END

      FUNCTION gammq(a,x)
      IMPLICIT NONE
      REAL*8 a,gammq,x
!U    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END

      SUBROUTINE gser(gamser,a,x,gln)
      IMPLICIT NONE
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
!U    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

      SUBROUTINE gcf(gammcf,a,x,gln)
      IMPLICIT NONE
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!U    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

      FUNCTION gammln(xx)
      IMPLICIT NONE
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
       -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END

      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,&
       funcs,alamda)
      IMPLICIT NONE
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL*8 alamda,chisq,funcs,a(ma),alpha(nca,nca),covar(nca,nca),&
       sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20)
!U    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,mfit
      REAL*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      SAVE ochisq,atry,beta,da,mfit
      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=0.001
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
        ochisq=chisq
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      do 14 j=1,mfit
        do 13 k=1,mfit
          covar(j,k)=alpha(j,k)
13      continue
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        call covsrt(alpha,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)
        endif
15    continue
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
      if(chisq.lt.ochisq)then
        alamda=0.1*alamda
        ochisq=chisq
        do 17 j=1,mfit
          do 16 k=1,mfit
            alpha(j,k)=covar(j,k)
16        continue
          beta(j)=da(j)
17      continue
        do 18 l=1,ma
          a(l)=atry(l)
18      continue
      else
        alamda=10.*alamda
        chisq=ochisq
      endif
      return
      END

      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,&
       funcs)
      IMPLICIT NONE
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL*8 chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),&
       y(ndata)
      EXTERNAL funcs
      PARAMETER (MMAX=20)
      INTEGER mfit,i,j,k,l,m
      REAL*8 dy,sig2i,wt,ymod,dyda(MMAX)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      do 16 i=1,ndata
        call funcs(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
17      continue
18    continue
      return
      END

      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      IMPLICIT NONE
      INTEGER ma,mfit,npc,ia(ma)
      REAL*8 covar(npc,npc)
      INTEGER i,j,k
      REAL*8 swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      IMPLICIT NONE
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL*8 big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
