      module global_data

      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      use util_memory
      use util_search
      implicit none
      save

! CONTROL.INC
      logical::lgibbs,lgrand,lnpt,lanes,llj,lmipsw,lexpee
      logical::lpbc,lpbcx,lpbcy,lpbcz,lfold,lijall,lewald
      logical::lchgall,ldielect
      logical::lexpsix,lmmff
      logical::lninesix,lgaro,lionic
      logical::lgenlj
      logical::ltailc,lshift,ltailcZeo
      logical::lneigh,lvirial
      logical::lcutcm,ldual
      integer(KIND=normal_int)::nmax,numax,ntmax,nntype,nvmax,nbxmax ,npamax,npabmax,maxbin,cmax,cmaxa,tor_bin_max,nbinmax_ete,smax
      integer(KIND=normal_int)::vectormax,maxvir,maxntemp,nchmax ,nchtor_max,nchbn_max
      integer(KIND=normal_int)::ntdifmx,nbinmx,angle_max,ang_bin_max ,tor_max,iou
!      common /iounit/ iou
! --------------------------------------------------------------
! *******************************
! *** PARAMETERS FOR ENSEMBLE ***
! *******************************
! if LNPT=.TRUE..
!    then a NPT volume move is used to equilibrate with a pressure bath
!    (implies cubic simulation boxes)
!    else an NVT simulation is performed
!      parameter (lnpt = .true.)
! if LGIBBS=.TRUE..
!    then a Gibbs-ensemble simulation is performed
!    (implies cubic simulation boxes)
!      parameter (lgibbs = .true.)
! if LGRAND=.TRUE.
!    then simulation is performed in the grand-canonical ensemble
!      parameter (lgrand = .false.)
! if LANES=.TRUE.
!    then simulation is performed in the adiabatic nuclear and electronic
!    sampling technique for polarizable force fields
!      parameter (lanes = .false.)
! if LVIRIAL=.TRUE.
!    then one chain will be simulated in each box independently and
!    the second viral coefficient will be calculated for their
!    interactions at a series of distances along the x-axis
!      parameter (lvirial = .false.)
! if lmipsw is true, then thermodynamic integration is performed
! for the phases, fort.35 must be supplied if so
!      parameter (lmipsw = .false.)
! if lexpee is true, then expanded esnemble is performed
! for the phases, fort.44 must be supplied if so
!      parameter (lexpee = .false.)
! ******************************************
! *** PARAMETERS FOR BOUNDARY CONDITIONS ***
! ******************************************
! if LPBC = .TRUE.
!    then periodic boundaries are used
      parameter (lpbc = .true.)
! if LPBCX = .TRUE.
!    then periodic boundary in x-directions is used
      parameter (lpbcx = .true.)
! if LPBCY = .TRUE.
!    then periodic boundary in y-directions is used
      parameter (lpbcy = .true.)
! if LPBCZ = .TRUE.
!    then periodic boundary in z-directions is used
      parameter (lpbcz = .true.)
! if LFOLD = .TRUE.
!    then coordinates are always folded into central box
      parameter (lfold = .true.)
! ********************************************
! *** PARAMETERS FOR INTERACTION CUT-OFFS  ***
! ********************************************
! if LIJALL = .TRUE.
!    then all i-j-interactions are considered (no potential cut-off)
!    must set lcutcm and ltailc to false for lijall = true
!    * check top of sumup.f if lijall = true and lchgall is false
!      parameter (lijall = .false.)
! if LCHGALL= .TRUE.
!    then all the electrostatic interaction are considered
!      parameter (lchgall = .false.)
! if LCUTCM=.TRUE. then a cutoff of the centers of mass will be used
! with a value of rcmu as calculated in ctrmas
!      parameter (lcutcm = .true.)
! if LEWALD=.TRUE. then ewald-sum will be used to calculate the
! electrostatic interactions.
!      parameter (lewald = .true.)
! if LDIELECT=.TRUE. then dielectric constant will be calculated and
! LEWALD must be .TRUE.
! Correct only in NVT ensemble
!      parameter (ldielect = .false.)
! ********************************************
! *** PARAMETERS FOR CBMC-PSEUDO-POTENTIAL ***
! ********************************************
! if LDUAL=.TRUE. then the external potential during a CBMC growth will
! only go out to a radius of rcutin and then will be corrected to the
! full rcut at the end.  This is Dual Cutoff Configurational-bias Monte
! Carlo (DC-CBMC)
!      parameter (ldual = .true.)

! ********************************************
! ***   PARAMETER FOR LENNARD-JONES 12-6   ***
!*********************************************
! This is only used if you are using regular lj potential ie without shifting or modifying
! it in any other way.
!      parameter (llj = .false.)
! ********************************************
! *** PARAMETERS FOR EXP-6 POTENTIAL       ***
! ********************************************
! if LEXPSIX=.TRUE. exp-6 potential will be used instead of lennard jones
!      parameter (lexpsix = .false.)
! ***************************************************
! *** PARAMETERS FOR BUFFERED 14-7 POTENTIAL      ***
! ***************************************************
! if LMMFF=.TRUE. buffered 14-7 potential will be used.
!      parameter (lmmff = .false.)
! ************************************
! *** PARAMETERS FOR 9-6 POTENTIAL ***
! ************************************
! if LNINESIX=.TRUE. then the 9-6 potential will be used.
!      parameter (lninesix = .false.)
! *****************************************************
! *** PARAMETERS FOR GENERALIZED LJ POTENTIAL      ***
! *****************************************************
! if LGENLJ=.TRUE. then the the Generalized lennard jones
! potential will be used.
!      parameter (lgenlj=.false.)
! ********************************************
! *** PARAMETERS FOR GAROFALINI POTENTIAL  ***
! ********************************************
! if LGARO=.TRUE. garofalini potential will be used instead of lennard jones
!      parameter (lgaro = .false.)
! ********************************************
! *** PARAMETER FOR IONIC SYSTEMS  ***
! ********************************************
! if LIONIC=.TRUE. System contains charged species, so system may not neutral
!      parameter (lionic = .false.)
! ***************************************
! *** PARAMETERS FOR TAIL CORRECTIONS ***
! ***************************************
! if LTAILC=.TRUE. tail corrections are added
!    (WARNING:  .lsami. in external.inc switches an intrinsic
!               tailcorrection on)
!      parameter (ltailc = .true.)
! truncated and shifted potentials
!      parameter (lshift = .false.)
! ***************************************
! *** PARAMETERS OF NEIGHBOR LIST     ***
! ***************************************
! if LNIEGH=.TRUE. the nearest neighbor list will be used with a
! value of rcutnn specified by fort.4
!      parameter (lneigh = .false.)
! ***************************************
! ** DIMENSIONS FOR ARRAYS             **
! ***************************************
! - nmax = maximum number of chains +2
      parameter (nmax = 5000)
! - numax = maximum number of units
      parameter (numax = 35)
! - ntmax = maximum number of types of chains
      parameter (ntmax = 35)
! - smax = max no. of mstates
      parameter (smax = 35)
! - nbxmax = maximum number of boxes
      parameter (nbxmax = 3)
! - nntype = number of types of beads
      parameter (nntype = 450)
! - nvmax = maximum number of bond vibration and bending constants
      parameter (nvmax = 300)
! - npamax = maximum number of pairs to switch or swatch
      parameter (npamax = 40)
! - npabmax = maximum number of box pairs (for swatch and swap)
      parameter (npabmax = 3)
! - maxvir = maximum number of bins for the 2nd virial coefficient
      parameter (maxvir = 1,maxntemp=3)
! - nchmax = maximum number of choices for trial sites in CBMC growth
      parameter (nchmax=100)
! - nchbn_max = maxium number of choices for bending CBMC growth
      parameter (nchbn_max=2000)
! - nchtor_max = maximum number of choices for torsion CBMC growth
      parameter (nchtor_max = 200)
! - vectormax = the maximum number of reciprocal vectors for Ewald sum
      parameter (vectormax=100000)

! ****************************************
! *** SAFECBMC max number of bins        ***
! ****************************************
      parameter (maxbin = 201)

! *************************************************************
! *** cmax is the max number of cells for the linkcell list ***
! *** cmaxa is the max number of molecules per cell         ***
! *************************************************************
      parameter (cmax = 1000)
      parameter (cmaxa = 100)

!*************************************************************
!***
!***      FOR ANALYSIS PART
!***
!*************************************************************

! NTDIFMX = max number of diff beads in simulation

	parameter(ntdifmx=2)

! NBINMX = max number of bins for histograms

	parameter(nbinmx=2)

! THESE ARE FOR ANGLE AND TORSIONAL DISTRIBUTION

        parameter(angle_max=2)
	parameter(ang_bin_max=2)
	parameter(tor_bin_max=2)
        parameter(tor_max=2)

!  For end to end vector calculation

       parameter(nbinmax_ete=2)

! EXTERNAL.INC
      real(KIND=double_precision)::alpha1,alpha2,beta1,beta2
      real(KIND=double_precision)::a1,delta,rsol
      integer(KIND=normal_int)::ntsubst
      integer(KIND=normal_int)::tau1,tau2
      logical::ljoe,lsami,lexzeo,lpsurf,lzgrid,lmuir,lgraphite
      logical::lslit
      logical::lcorreg,lelect_field
!  external electric field
!      parameter (lelect_field=.false.)
! surfaces:
!      parameter (ljoe=.false.,lsami=.false.,lmuir=.false.)
!      parameter (lslit = .false.)
!      parameter (lpsurf=.false.)
!      parameter (lgraphite=.false.)
!      parameter (lcorreg = .false.)
! Graphite surface
! -- a1, delta are in Angstroms
! -- rsol [=] 1/A^3
      parameter (a1 = 2.460d0)
      parameter (delta = 3.40d0)
      parameter (ntsubst = 190)
      parameter (rsol = 0.114d0)
! zeolites:
!      parameter (lexzeo=.true.,lzgrid=.true.)
! Joe Hautman's parameters are read from standard input:
!               to be used with ljoe = .true.
! Sami's parameters: to be used with lsami = .true. AND lmuir = .true.
      parameter ( alpha1=21.162d0, alpha2=-21.162d0, beta1=661.6d0, beta2=6616.0d0, tau1=-32, tau2=-16 )
! parameters for zeolites are read in subroutine SUZEO
! AT PRESENT no parameters for polymeric surfactants (see LJPSUR)

! OPENMP.inc
      integer(KIND=normal_int)::thread_id,thread_num,thread_num_max,thread_num_proc

! MPI.INC
      integer(KIND=normal_int)::myid,numprocs,ierr
      integer(KIND=normal_int),parameter::numprocmax=32
!      common /mpivar/ myid,numprocs,ierr

! CBMC.INC
      logical::lovr(nchmax),lexist(numax),lexshed(numax),llplace(ntmax) ,lpnow(numax),llrig,lsave(numax),lrplc(ntmax)

! *** adding this for iswatch move ***
      logical::liswinc(numax,ntmax), liswatch
      integer(KIND=normal_int)::other,iparty,iplace(numax,numax) ,pfrom(numax),pnum(numax),pprev(numax),nplace
! *** end add ***

      real(KIND=double_precision)::bncb(ntmax,numax),bscb(ntmax,2,numax) ,rxp(numax,nchmax),ryp(numax,nchmax),rzp(numax,nchmax) ,rxnew(numax),rynew(numax),rznew(numax) ,bnregr(ntmax) ,bsregr(ntmax,2),bfac(nchmax),bsswat(npamax ,npabmax) ,bsiswat(npamax,npabmax),fbncb(ntmax,numax) ,fbscb(ntmax ,2 ,numax),vtry(nchmax),vtrext(nchmax),vtrintra(nchmax) ,vtrinter(nchmax) ,vtrelect(nchmax),vtrewald(nchmax) ,vtrorient(nchmax),vtrelect_intra(nchmax),vtrelect_inter(nchmax)
      integer(KIND=normal_int)::nchoi1(ntmax),nchoi(ntmax),icbdir(ntmax) ,icbsta(ntmax),nmaxcbmc(ntmax),nchoih(ntmax),nchtor(ntmax) ,nchbna(ntmax),nchbnb(ntmax),nrigi,nchoir(ntmax),rlist(numax ,numax),rfrom(numax),rprev(numax) ,rnum(numax)
      real(KIND=double_precision)::acchem(nbxmax,ntmax), bnchem(nbxmax ,ntmax)
!      common /pbias/  bncb, bscb
!     &   ,bnregr,bsregr
!     &   ,bsswat,bsiswat
!     &   ,fbncb, fbscb
!      common /cbpara/ nchoi1,nchoi,icbdir
!     &   ,icbsta,nmaxcbmc,nchoih
!     &   ,nchtor,nchbna,nchbnb
!     &   ,nchoir,lrplc
!      common /bbias/  bfac, vtry,vtrinter
!     &   ,vtrext,vtrintra
!     &   ,vtrelect,vtrewald,lovr
!     &   ,vtrorient,vtrelect_intra,
!     &    vtrelect_inter
!      common /xyzp/   lexist,lexshed,rxp
!     &   ,ryp, rzp
!     &   ,liswinc,liswatch,iparty,other,llplace
!     &   ,llrig,lpnow,lsave
!      common /xyznew/ rxnew, rynew, rznew
!      common /achepot/ acchem,bnchem
!      common /arig/ rlist,rfrom,rprev
!     &   ,rnum,nrigi
!      common /aplace/ nplace,iplace,pfrom
!     &   ,pnum,pprev

! POTEN.INC
      logical::lspecial,lqchg,lij,lpl,lq14scale
      logical::L_Coul_CBMC
      logical::L_tor_table, L_spline, L_linear
      integer(KIND=normal_int)::a15type(ntmax,numax,numax)
      real(KIND=double_precision)::sig2ij,epsij,epsilon,sigma,epsi,sigi ,ecut,extc12, extc3, extz0,aspecial,bspecial,q1
      real(KIND=double_precision)::xiq,jayq,jayself,fqegp(ntmax),a15(2) ,qscale,qelect
      real(KIND=double_precision)::ljscale,qscale2
      real(KIND=double_precision)::n0,n1
      dimension extc12(nntype), extc3(nntype), extz0(nntype) ,xiq(0:nntype),jayself(nntype),lqchg(0:nntype),lij(0:nntype)
      dimension sig2ij(nntype*nntype), epsij(nntype*nntype), ecut(nntype*nntype), jayq(nntype*nntype),epsilon(2,numax),sigma(2,numax),epsi(0:nntype), sigi(0:nntype),q1(nntype),qelect(0:nntype)
      dimension aspecial(nntype*nntype),bspecial(nntype*nntype) ,lspecial(nntype*nntype),lpl(nntype)
      dimension ljscale(ntmax,numax,numax),qscale2(ntmax,numax,numax)
      dimension lq14scale(ntmax),qscale(ntmax)

! ***********************************************************
! *** parameters for Generalized Lennard Jones Potential  ***
! repulsive part
      parameter(n0=12.0d0)
!   attractive part
      parameter(n1=6.0d0)
! ****  Ref:  J. Chem. Phys. 120, 4994 (2004)         ****
! ***********************************************************

!      common /ljpot/ sig2ij,epsij,ecut,epsilon,sigma,epsi,sigi
!     &   ,q1,qscale,a15,a15type
!      common /ljepot/ ljscale,qscale2,lq14scale,qelect
!      common /extpot/ extc12, extc3, extz0
!      common /specail/ aspecial,bspecial,lspecial,lpl
!      common /electr/ qqfact,xiq,jayq,jayself,fqegp,lqchg,lij
!      common /etects/ L_Coul_CBMC,L_tor_table, L_spline, L_linear

! COORD.INC
      logical::lelect,lflucq,lqtrans,lexpand,lpresim,lring,lrigid,lplace ,lrig,lrigi,lchiral,licell
      integer(KIND=normal_int)::nbox,nchain,nmolty,moltyp,nunit,nugrow ,ntype,nchbox,ncmt,ncmt2,iring,nrig,irig,frig,counthist,counttot ,boxlink
      integer(KIND=normal_int)::nboxi, eetype, parbox,iurot,iupdatefix ,nrigmin,nrigmax
      integer(KIND=normal_int)::parall,temtyp,riutry,rindex,prior ,maxgrow,nrotbd,irotbd
      real(KIND=double_precision)::rxu, ryu, rzu, mass, masst, beta,xcm ,ycm,zcm,exp_c,eta,eta2,sxcm,sycm,szcm
      real(KIND=double_precision)::pmrotbd
      real(KIND=double_precision)::qqu,rcmu,rintramax,hist,probf ,Elect_field
      character(LEN=18)::chname
      character(LEN=4)::chemid
      dimension lelect(ntmax),lflucq(ntmax),lqtrans(ntmax) ,lexpand(ntmax),lplace(ntmax,nntype)
      dimension maxgrow(ntmax),probf(30,30,maxbin) ,hist(30,30,maxbin),iring(ntmax),nrig(ntmax),irig(ntmax,6) ,prior(ntmax,numax),frig(ntmax,6)
      dimension moltyp(nmax),nrigmin(ntmax),nrigmax(ntmax)
      dimension nunit(ntmax),nugrow(ntmax)
      dimension ntype(ntmax,numax)
      dimension nboxi(nmax),eetype(ntmax),rcmu(nmax)
      dimension iurot(ntmax)
      dimension nrotbd(ntmax),irotbd(numax,ntmax),pmrotbd(numax,ntmax)
      dimension rxu(nmax,numax), ryu(nmax,numax), rzu(nmax,numax)
      dimension qqu(nmax,numax), Elect_field(nbxmax)
      dimension nchbox(nbxmax),ncmt(nbxmax,ntmax) ,ncmt2(nbxmax,ntmax,20),eta(nbxmax,ntmax,20)
      dimension parall(ntmax,nmax),temtyp(ntmax)
      dimension mass(0:nntype),masst(ntmax),parbox(nmax,nbxmax,ntmax)
      dimension xcm(nmax),ycm(nmax),zcm(nmax),exp_c(nmax) ,eta2(nbxmax,ntmax),sxcm(nmax),sycm(nmax),szcm(nmax)
      dimension riutry(ntmax,3),rindex(ntmax),lrigid(ntmax),lrig(ntmax) ,lrigi(ntmax,numax),lchiral(ntmax,numax),lring(ntmax)
      dimension chname(0:nntype), chemid(nntype)
!      common /ncnunt/ moltyp,nunit,nugrow,ntype,parbox,boxlink
!      common /ccord1/ probf, hist,rxu, ryu, rzu, qqu, rcmu, Elect_field
!      common /ccord2/ nboxi,eetype,iurot,nchbox,ncmt,ncmt2
!     &  ,parall,temtyp,licell,rintramax
!      common /ccord3/ nbox,nchain,nmolty,nrigmin,nrigmax
!      common /electa/ lelect,lflucq,lqtrans,lexpand,lring,lpresim
!     &  ,lplace,lrigid,lrig,lrigi,lchiral
!      common /masses/ mass, masst, xcm, ycm, zcm, exp_c, beta
!     &  ,eta,eta2,sxcm,sycm,szcm
!      common /histfix/ counthist,counttot,iupdatefix,maxgrow,iring
!      common /rigidstuff/ riutry,rindex,nrig,irig,frig
!      common /verbose/ chname, chemid
!      common /ccord4/ pmrotbd
!      common /ccord5/ irotbd,nrotbd

! COORD2.INC
      integer(KIND=normal_int)::moltion(2)
      real(KIND=double_precision)::rxuion(numax,2),ryuion(numax,2) ,rzuion(numax,2),qquion(numax,2),exp_cion(2)
!      common /coord2/ rxuion,ryuion,rzuion,
!     &                qquion,moltion,exp_cion

! CONNECT.INC
      logical::linclu(ntmax,numax,numax),lqinclu(ntmax,numax,numax) ,lainclu(ntmax,numax,numax),lexclu(ntmax,numax,ntmax,numax) ,wschvib(numax,6),wschben(numax,12),wschtor(numax,12)
      integer(KIND=normal_int)::invib(ntmax,numax),itvib(ntmax,numax,6) ,ijvib(ntmax,numax,6),inben(ntmax,numax),itben(ntmax,numax,12) ,ijben2(ntmax,numax,12),ijben3(ntmax,numax,12),intor(ntmax,numax) ,ittor(ntmax,numax,12),ijtor2(ntmax,numax,12),ijtor3(ntmax,numax ,12),ijtor4(ntmax,numax,12),wsched(numax)
      real(KIND=double_precision)::brvib(nvmax),brvibk(nvmax) ,brvibmin(60),brvibmax(60),brben(nvmax),brbenk(nvmax)
!      common /connev/ brvib,brvibk,brvibmin
!     &   ,brvibmax
!     &   ,itvib,ijvib
!     &   ,invib
!      common /conneb/ brben,brbenk
!     &   ,itben,ijben2
!     &   ,ijben3
!     &   ,inben
!      common /connet/ ittor,ijtor2
!     &   ,ijtor3,ijtor4
!     &   ,intor
!      common /conneg/ wsched,wschvib
!     &    ,wschben,wschtor
!      common /includ/ linclu
!     &    ,lqinclu,lainclu
!      common /exclud/ lexclu

! EWALDSUM.INC
      logical::L_Ewald_Auto
      integer(KIND=normal_int)::kmax,numvect
! - kalp, a parameter to control the real::space sum
! - kmax, a parameter to control the total number of reciprocal vectors
! - numvect, the total number of reciprocal vectors
      dimension numvect(nbxmax)
      real(KIND=double_precision)::kx,ky,kz,prefact,ssumr,ssumi,ssumrn ,ssumin,ssumro,ssumio
      real(KIND=double_precision),target::sself,correct,kalp,calp,dipolex ,dipoley,dipolez,dipolexo,dipoleyo,dipolezo
! - calp = kalp / boxlen
      dimension kx(vectormax,nbxmax),ky(vectormax,nbxmax), kz(vectormax,nbxmax),prefact(vectormax,nbxmax), ssumr(vectormax,nbxmax),ssumi(vectormax,nbxmax), ssumrn(vectormax,nbxmax),ssumin(vectormax,nbxmax), calp(nbxmax),kmax(nbxmax),kalp(nbxmax), ssumro(vectormax,nbxmax),ssumio(vectormax,nbxmax)
      dimension dipolex(nbxmax),dipoley(nbxmax),dipolez(nbxmax)
!      common /ewaldl/ L_Ewald_Auto
!      common /numvect/ numvect,kmax
!      common /vector/ kx,ky,kz
!      common /ssum/ calp,kalp,ssumr,ssumi,prefact,ssumrn,
!     &                  ssumin,ssumro,ssumio
!      common /selfterm/ sself,correct
!      common /dipoterm/ dipolex,dipoley,dipolez,dipolexo,dipoleyo,
!     &                  dipolezo

! FEPSI.INC
      logical::lfepsi
      real(KIND=double_precision)::aslope,bslope,ashift,bshift,a0,b0

! ----- If lfepsi is true, the fluctuation of epsilon is used instead of
! ----- the fluctuation of the sigma.

        parameter (lfepsi = .false.)
! *** slope=0.3 a=3.05
!	parameter (a0 = 0.2003d0)
!	parameter (b0 = 1.3946d0)
!	parameter (aslope = 8.85d5)
!	parameter (bslope = 158.25d0)
!	parameter (ashift = 7.227d5)
!	parameter (bshift = 505.97d0)
! *** slope=0.3 a=2.90
!	parameter (a0 = 0.15561d0)
!	parameter (b0 = 1.2960d0)
!	parameter (aslope = 5.5475d5)
!	parameter (bslope = 131.23d0)
!	parameter (ashift = 4.121d5)
!	parameter (bshift = 381.96d0)
! *** slope=0.3 a=2.95
!	parameter (a0 = 0.16833d0)
!	parameter (b0 = 1.3299d0)
!	parameter (aslope = 6.5125d5)
!	parameter (bslope = 139.75d0)
!	parameter (ashift = 4.9975d5)
!	parameter (bshift = 419.81d0)
! *** slope=0.3 a=2.82
!	parameter (a0 = 0.13150d0)
!	parameter (b0 = 1.2574d0)
!	parameter (aslope = 4.2863d5)
!	parameter (bslope = 117.5d0)
!	parameter (ashift = 3.0208d5)
!	parameter (bshift = 323.6d0)
! *** slope=0.3 a=2.75
!	parameter (a0 = 0.11037d0)
!	parameter (b0 = 1.1942d0)
!	parameter (aslope = 3.4013d5)
!	parameter (bslope = 108d0)
!	parameter (ashift = 2.2868d5)
!	parameter (bshift = 285.04d0)
! *** slope=0.3 a=2.85
!	parameter (a0 = 0.14035d0)
!	parameter (b0 = 1.2596d0)
!	parameter (aslope = 4.7263d5)
!	parameter (bslope = 123.25d0)
!	parameter (ashift = 3.3978d5)
!	parameter (bshift = 347.63d0)
! *** slope=0.3 a=2.88
!	parameter (a0 = 0.14820d0)
!	parameter (b0 = 1.2689d0)
!	parameter (aslope = 5.2125d5)
!	parameter (bslope = 128.75d0)
!	parameter (ashift = 3.8220d5)
!	parameter (bshift = 371.29d0)
! *** slope=0.5 a=2.71
      parameter (a0 = -0.22818d0)
      parameter (b0 = 0.44662d0)
      parameter (aslope = 14.0738d5)
      parameter (bslope = 351.25d0)
      parameter (ashift = 3.5852d5)
      parameter (bshift = 358.94d0)
! *** slope=0.5 a=2.68
!	parameter (a0 = -0.23379d0)
!	parameter (b0 = 0.43630d0)
!	parameter (aslope = 12.78125d5)
!	parameter (bslope = 337.5d0)
!	parameter (ashift = 3.1904d5)
!	parameter (bshift = 337.86d0)
! *** slope=0.5 a=2.655
!	parameter (a0 = -0.23808d0)
!	parameter (b0 = 0.425507d0)
!	parameter (aslope = 11.7775d5)
!	parameter (bslope = 322.04d0)
!	parameter (ashift = 2.8902d5)
!	parameter (bshift = 326.875d0)
! *** slope=0.5 a=2.665
!	parameter (a0 = -0.23641d0)
!	parameter (b0 = 0.43409d0)
!	parameter (aslope = 12.171d5)
!	parameter (bslope = 330d0)
!	parameter (ashift = 3.0075d5)
!	parameter (bshift = 326.52d0)
! *** slope=0.5 a=2.65
!	parameter (a0 = -0.23903d0)
!	parameter (b0 = 0.42231d0)
!	parameter (aslope = 11.5875d5)
!	parameter (bslope = 325d0)
!	parameter (ashift = 2.8339d5)
!	parameter (bshift = 319.44d0)
! *** slope=0.5 a=2.60
!	parameter (a0 = -0.24796d0)
!	parameter (b0 = 0.40345d0)
!	parameter (aslope = 9.8238d5)
!	parameter (bslope = 304d0)
!	parameter (ashift = 2.3206d5)
!	parameter (bshift = 288.72d0)
! *** slope=0.5 a=2.55
!	parameter (a0 = -0.25703d0)
!	parameter (b0 = 0.38299d0)
!	parameter (aslope = 8.3075d5)
!	parameter (bslope = 284.38d0)
!	parameter (ashift = 1.8945d5)
!	parameter (bshift = 261.09d0)
! *** slope = 0.7 a=2.52
!	parameter (a0 = -0.39106d0)
!	parameter (b0 = 0.08431d0)
!	parameter (aslope = 25.22875d5)
!	parameter (bslope = 662.75d0)
!	parameter (ashift = 3.0689d5)
!	parameter (bshift = 335.43d0)
! *** slope = 0.7 a=2.51
!	parameter (a0 = -0.39233d0)
!	parameter (b0 = 0.08203d0)
!	parameter (aslope = 24.425d5)
!	parameter (bslope = 653.75d0)
!	parameter (ashift = 2.9499d5)
!	parameter (bshift = 328.60d0)
! *** slope = 0.7 a=2.50
!	parameter (a0 = -0.39357d0)
!	parameter (b0 = 0.07946d0)
!	parameter (aslope = 23.641d5)
!	parameter (bslope = 645d0)
!	parameter (ashift = 2.835d5)
!	parameter (bshift = 322.13d0)
! *** slope = 0.7 a=2.49
!	parameter (a0 = -0.394906d0)
!	parameter (b0 = 0.075678d0)
!	parameter (aslope = 22.88625d5)
!	parameter (bslope = 637.735d0)
!	parameter (ashift = 2.72469d5)
!	parameter (bshift = 316.22d0)
! *** slope = 0.7 a=2.48
!	parameter (a0 = -0.396175d0)
!	parameter (b0 = 0.072973d0)
!	parameter (aslope = 22.149125d5)
!	parameter (bslope = 629d0)
!	parameter (ashift = 2.61798d5)
!	parameter (bshift = 309.95d0)
! *** slope = 0.7 a=2.47
!	parameter (a0 = -0.39745d0)
!	parameter (b0 = 0.070157d0)
!	parameter (aslope = 21.4335d5)
!	parameter (bslope = 620.75d0)
!	parameter (ashift = 2.5151d5)
!	parameter (bshift = 306.89d0)
! *** slope = 0.7 a=2.45
!	parameter (a0 = -0.399983d0)
!	parameter (b0 = 0.064159d0)
!	parameter (aslope = 20.06425d5)
!	parameter (bslope = 604.75d0)
!	parameter (ashift = 2.32037d5)
!	parameter (bshift = 292.09d0)
! *** slope = 0.7 a=2.40
!	parameter (a0 = -0.40629d0)
!	parameter (b0 = 0.050088d0)
!	parameter (aslope = 16.9775d5)
!	parameter (bslope = 565.50d0)
!	parameter (ashift = 1.89213d5)
!	parameter (bshift = 263.93d0)
! *** slope = 0.9 a=2.30
!	parameter (a0 = -0.48332d0)
!	parameter (b0 = -0.12334d0)
!	parameter (aslope = 34.64125d5)
!	parameter (bslope = 1014.25d0)
!	parameter (ashift = 2.2815d5)
!	parameter (bshift = 294.26d0)
! *** slope = 0.9 a=2.32
!	parameter (a0 = -0.48137d0)
!	parameter (b0 = -0.11887d0)
!	parameter (aslope = 36.99125d5)
!	parameter (bslope = 1041.25d0)
!	parameter (ashift = 2.4744d5)
!	parameter (bshift = 306.23d0)
! *** slope = 0.3 a=2.85
!	parameter (a0 = 0.0d0)
!	parameter (b0 = 0.0d0)
!	parameter (aslope = 3.0d5)
!	parameter (bslope = 0.0d0)
!	parameter (ashift = 8.0d5)
!	parameter (bshift = 1200.0d0)
      real(KIND=double_precision)::favor(nmax),favor2(nmax)
!      common /fdata/ favor,favor2

! QQLIST.INC
      integer(KIND=normal_int)::leaderq
      dimension leaderq(ntmax,numax)
!      common /qleader/ leaderq

! CLUSTERBIAS.INC
      logical::lbias(ntmax),lavbmc1(ntmax),lavbmc2(ntmax),lavbmc3(ntmax)
      real(KIND=double_precision)::rbsmax,rbsmin,vol_eff,pmbias(ntmax) ,pmbsmt(ntmax),pmbias2(ntmax)
!      common /cbiasa/ lbias,lavbmc1,lavbmc2,lavbmc3
!      common /cbiasb/ rbsmax,rbsmin,vol_eff,pmbias,pmbias2,pmbsmt

! NSIX.INC
! ***************************************************
! * stores interaction parameters for 9-6 potential *
! ***************************************************
      integer(KIND=normal_int)::nxatom
      parameter (nxatom=5)
      real(KIND=double_precision)::rzero(nxatom*nxatom),epsnx(nxatom *nxatom)
      real(KIND=double_precision)::shiftnsix(nxatom*nxatom)
!      common /formic/ rzero,epsnx,shiftnsix

! PEBOCO.INC
      real(KIND=double_precision)::bx,by,bz,hbx,hby,hbz,bxi,byi,bzi
!      common /peboco/ bx,by,bz,hbx,hby,hbz,bxi,byi,bzi

! CELL.INC
      real(KIND=double_precision),target::hmat(nbxmax,9),hmati(nbxmax,9),cell_length(nbxmax,3),min_width(nbxmax,3),cell_vol(nbxmax),cell_ang(nbxmax,3)
!      common /cella/ hmat,hmati
!      common /cellb/ cell_length,min_width,cell_vol,cell_ang

! CONVER.INC
!      real(KIND=double_precision)::raddeg, degrad ,twopi, onepi,eXV_to_K
!      common /convrs/ raddeg, degrad, twopi, onepi,eXV_to_K

! TABULATED.INC
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccc  for use with tabulated potentials
!cccc  if you are not using a tabulated potential,
!cccc  reduce the size of the arrays!!!
!cccc  to pass the test suite, arrays for tabulated electrostatics
!cccc  need to be sufficiently large
!cccc  num_int_elect(150,150), tabelect(1500,150,150), etc.
!cccc  KM 2009
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer(KIND=normal_int)::vibsplits(1), bendsplits(1),ntabvib ,ntabbend
      integer(KIND=normal_int)::vdWsplits(1,1), ntabvdW
      integer(KIND=normal_int)::electsplits(10,10), ntabelect
      integer(KIND=normal_int)::num_int_vib(1), num_int_bend(1), 	num_int_vdW(1,1), num_int_elect(1,1)
      real(KIND=double_precision)::vib(10,1), tabvib(10,1)
      real(KIND=double_precision)::bend(10,1), tabbend(10,1)
      real(KIND=double_precision)::rvdW(10,1,1), tabvdW(10,1 ,1)
      real(KIND=double_precision)::relect(10,1,1), tabelect(10 ,1,1)
      real(KIND=double_precision)::vibdiff(10,1), benddiff(10,1) ,vdWdiff(10,1,1), electdiff(10,1,1)
      logical::L_vib_table, L_bend_table, L_vdW_table, L_elect_table
!      common /tabulated/ vib, tabvib, bend, tabbend, rvdW,
!     & 	 tabvdW, relect, tabelect, vibdiff, benddiff, vdWdiff,
!     &   electdiff, electsplits, ntabelect, num_int_elect,
!     & 	 vibsplits, bendsplits, vdWsplits, ntabvib, ntabbend,
!     &   ntabvdW, num_int_vib, num_int_bend, num_int_vdW
!      common /ltab/ L_vib_table, L_bend_table, L_vdW_table,
!     &   L_elect_table

! TORSION.INC
	integer(KIND=normal_int)::splpnts(10),nttor
	real(KIND=double_precision)::deg(10,2),tabtorso(10,2) ,torderiv2(10,2)
	real(KIND=double_precision)::tordif(10,2)
!	common /torsion/ deg,tabtorso,torderiv2,tordif,splpnts,nttor

! SYSTEM.INC
      logical,target::lsolid(nbxmax),lrect(nbxmax),ltwice(nbxmax) ,lideal(nbxmax)
      integer(KIND=normal_int)::icell,nicell,iucell,solcount,numberDimensionIsIsotropic(nbxmax)
      real(KIND=double_precision),target::boxlx(nbxmax), boxly(nbxmax), boxlz(nbxmax), rcut(nbxmax)
      real(KIND=double_precision)::rmin,rcutnn(nbxmax), softcut,softlog ,rcutin,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc
      dimension icell(nmax),nicell(cmax),iucell(cmax,cmaxa) ,solcount(nbxmax,ntmax)
      dimension avsolinter(nbxmax,ntmax),avsolintra(nbxmax,ntmax) ,avsolbend(nbxmax,ntmax),avsoltor(nbxmax,ntmax) ,avsolelc(nbxmax,ntmax)
!      common /boxpha/ lsolid,lrect,ltwice
!      common /boxdim/ boxlx,boxly,boxlz,
!     &                 rcut,rcutnn,lideal
!      common /solav/ avsolinter,avsolintra,avsolbend,avsoltor,avsolelc
!      common /input1/ rmin, softcut,softlog
!     &  ,rcutin,solcount
!      common /linkc/ icell,nicell,iucell

! NEIGH.INC
! lnn(1,1):replace the 1s with nmax to use neighbor list
      logical::lnn(1,1),lneighbor
      integer(KIND=normal_int)::maxneigh
      parameter (maxneigh = 20)
      integer(KIND=normal_int)::neighbor(maxneigh,nmax),neigh_cnt(nmax) ,neigh_icnt,neighi(maxneigh),neighboro(maxneigh,nmax) ,neigh_o(nmax)
      real(KIND=double_precision)::ndij(maxneigh,nmax),nxij(maxneigh ,nmax),nyij(maxneigh,nmax),nxijo(maxneigh,nmax),nzij(maxneigh ,nmax),ndiji(nmax),nxiji(maxneigh),nyiji(maxneigh) ,nziji(maxneigh),ndijo(maxneigh,nmax),nyijo(maxneigh,nmax) ,nzijo(maxneigh,nmax)
!      common /nn2/ lnn,lneighbor
!      common /nn4/ neighbor
!      common /nn3/ neigh_cnt,neigh_icnt,neighi
!      common /nn5/ ndij,nxij,nyij,nzij
!      common /nn6/ ndiji,nxiji,nyiji,nziji
!      common /nn7/ ndijo,nxijo,nyijo,nzijo,neigh_o,neighboro

! NEIGH2.INC
! disvec(2,1,3):replace the 1 with nmax to use neighbor lists
      real(KIND=double_precision)::disvec(2,1,3)
!      common /nn1/ disvec

! ENSEMBLE.INC
      real(KIND=double_precision)::vbox(nbxmax),wbox(nbxmax) ,vinterb(nbxmax),vtailb(nbxmax),vintrab(nbxmax),vvibb(nbxmax) ,vbendb(nbxmax),vtgb(nbxmax),vextb(nbxmax), velectb(nbxmax) ,vflucqb(nbxmax),v3garob(nbxmax)
!      real(KIND=double_precision)::velectb_intra(nbxmax), velectb_inter(nbxmax)
!      common /vwbox/  vbox,wbox,vinterb
!     &   ,vtailb,vintrab,vvibb
!     &   ,vbendb,vtgb,vextb
!     &   ,velectb,vflucqb,v3garob
!      common /vabox/ velectb_intra,velectb_inter

! INPUTDATA.INC
      logical::lstop, ldie
      logical::lrdf,lintra,lstretch,lgvst,lbend,lete, lrhoz
      logical::L_movie_xyz
      character(LEN=1)::suffix
      integer(KIND=normal_int)::run_num
      logical::L_add, L_sub
      integer(KIND=normal_int)::N_add,N_box2add,N_moltyp2add
      integer(KIND=normal_int)::N_sub,N_box2sub,N_moltyp2sub
      integer(KIND=normal_int)::ghost_particles
      integer(KIND=normal_int)::ianalyze,nbin,tmcc
      integer(KIND=normal_int)::nequil,ninstf,ninsth,ndumph
      integer(KIND=normal_int)::nstep, iprint, imv, iratio, iblock, idiele,isolute,iratv,iratp,irsave,nnstep,nchoiq,numcoeff,rmexpc ,nswapb,box1,box2,nvolb,box5,box6,ntemp, iheatcapacity
      real(KIND=double_precision)::Armtrax,Armtray,Armtraz
      real(KIND=double_precision)::Num_cell_a, Num_cell_b, Num_cell_c
      real(KIND=double_precision)::pmcb, pmtra, temp, express, ewald_precision ,rmtrax,rmtray, rmtraz, rmrotx, rmroty, rmrotz ,rmflcq,taflcq,fqbeta,tatra, tarot,pmvolx,pmvoly ,pmvol, pmswap ,rmvol, rmhmat, tavol,pmfix ,upnn, upnnsq, upnndg, B, pmswmt ,pmiswat,pmisatc ,pmswat,pmsatc,pmcbmt,pmtrmt,pmromt ,pmvlmt ,pmall,pmflcq, pmfqmt,pmexpc,pmeemt ,pmswapb,pmvolb,virtemp ,bin_width,pm_atom_tra ,pmexpc1
      dimension pmswmt(ntmax),pmcbmt(ntmax) ,pmtrmt(ntmax),pmromt(ntmax) ,pmall(ntmax),pmfqmt(ntmax),pmeemt(ntmax) ,rmexpc(ntmax) ,pmswapb(ntmax,npabmax) ,box1(ntmax,npabmax),box2(ntmax,npabmax) ,box5(nbxmax),box6(nbxmax),pmvolb(nbxmax)
      dimension nchoiq(nbxmax),nswapb(ntmax)
      dimension pmsatc(npamax),pmisatc(npamax)
      dimension pmvlmt(nbxmax),virtemp(maxntemp)
      dimension numcoeff(ntmax),pmfix(ntmax)
      dimension B(ntmax),isolute(ntmax),ewald_precision(nbxmax)
      dimension express(nbxmax)
      dimension ghost_particles(nbxmax)
      dimension rmvol(nbxmax),rmhmat(nbxmax,9)
      dimension rmtrax(ntmax,nbxmax),rmtray(ntmax,nbxmax) ,rmtraz(ntmax,nbxmax),rmrotx(ntmax,nbxmax) ,rmroty(ntmax,nbxmax),rmrotz(ntmax,nbxmax)
      dimension rmflcq(ntmax,nbxmax)
!      common /indata/ nstep,iprint,imv,iratio,iblock, idiele, isolute
!     &   ,iratv,iratp,irsave,nnstep,lstop,nchoiq,numcoeff,ldie
!     &   ,rmexpc,nswapb,box1,box2,nvolb,box5,box6,ntemp,tmcc
!      common /indatc/ pmswmt,pmcbmt,pmtrmt,pmromt
!     &   ,pmfqmt,pmexpc,pmeemt,pmfix,pmexpc1
!     &   ,pmall,pmsatc,pmisatc,pmswat,pmiswat,pmcb,pmflcq
!     &   ,pmvlmt,pmtra,pmvol,pmswap,pmswapb,pmvolb,virtemp
!     &   ,ewald_precision
!      common /indatb/ rmvol,rmtrax,rmtray,rmtraz,rmrotx,rmroty,rmrotz
!     &         ,rmflcq,temp,express,tatra, tarot,tavol,taflcq,fqbeta
!     &         ,upnn, upnnsq, upnndg, B,rmhmat,pmvolx,pmvoly,bin_width,
!     &          pm_atom_tra
!      common /indatd/ ianalyze,lrdf,lintra,lstretch,lgvst,lbend
!     &         ,nbin,lete,lrhoz,L_movie_xyz
!      common /indatd2/ nequil,ninstf,ninsth,ndumph,run_num,suffix
!      common /addsol/  N_add,N_box2add,N_moltyp2add,L_add
!      common /subsol/  N_sub,N_box2sub,N_moltyp2sub,L_sub
!      common /atomdisp/ Armtrax,Armtray,Armtraz,Num_cell_a,
!     &                  Num_cell_b, Num_cell_c
!      common /ghost/ ghost_particles
!      common /heatcap/ iheatcapacity

! BNBSMA.INC
      real(KIND=double_precision)::bsvol,bnvol,bntrax,bntray,bntraz ,bshmat,bstrax,bstray,bstraz,bnrotx,bnroty,bnrotz,bnhmat,bsrotx ,bsroty,bsrotz,bnflcq,bsflcq,bnexpc,bsexpc,bnflcq2,bsflcq2
      real(KIND=double_precision)::Abntrax,Abntray,Abntraz, Abstrax,Abstray,Abstraz
      dimension bsvol(nbxmax), bnvol(nbxmax),bshmat(nbxmax,9), bnhmat(nbxmax,9)
      dimension bntrax(ntmax,nbxmax),bntray(ntmax,nbxmax), bntraz(ntmax,nbxmax)
      dimension bstrax(ntmax,nbxmax),bstray(ntmax,nbxmax), bstraz(ntmax,nbxmax)
      dimension bnrotx(ntmax,nbxmax),bnroty(ntmax,nbxmax), bnrotz(ntmax,nbxmax)
      dimension bsrotx(ntmax,nbxmax),bsroty(ntmax,nbxmax), bsrotz(ntmax,nbxmax)
      dimension bnexpc(ntmax,nbxmax),bsexpc(ntmax,nbxmax)
      dimension bnflcq(ntmax,nbxmax),bsflcq(ntmax,nbxmax), bnflcq2(ntmax,nbxmax),bsflcq2(ntmax,nbxmax)
!      common /bvolum/ bnvol,bsvol,bshmat,bnhmat
!      common /btrans/ bntrax,bntray,bntraz,bstrax,bstray,bstraz
!      common /brotat/ bnrotx,bnroty,bnrotz,bsrotx,bsroty,bsrotz
!      common /bflucq/ bnflcq,bsflcq,bnexpc,bsexpc,bnflcq2,bsflcq2
!      common /atomtra/ Abntrax,Abntray,Abntraz,Abstrax,
!     &          Abstray,Abstraz

! SWTCMOVE.INC
      integer(KIND=normal_int)::nswaty,nswatb,gswatc,nswach
      integer(KIND=normal_int)::nsampos,splist,nswtcb,box3,box4,ncut
      real(KIND=double_precision)::bnswat,bniswat,pmswtcb
      real(KIND=double_precision)::bnswat_empty
      dimension nswach(npamax),nsampos(npamax)
      dimension nswatb(npamax,npabmax)
      dimension bnswat(npamax,npabmax),ncut(npamax,2)
      dimension gswatc(npamax,2,2*npamax),splist(npamax,numax,2)
      dimension nswtcb(npamax)
      dimension pmswtcb(npamax,npabmax)
      dimension box3(npamax,npabmax)
      dimension box4(npamax,npabmax)
      dimension bnswat_empty(npamax,npabmax)
!      common /swtcma/ nswaty,nswtcb
!      common /swtcmb/ nswach,nsampos
!      common /swtcmc/ box3,box4,bnswat,bnswat_empty,nswatb
!      common /swtcmd/ gswatc,splist,ncut
!      common /swtcme/ pmswtcb

! ROSEN.INC
      integer(KIND=normal_int)::growfrom,growprev,growlist,grownum
      dimension growfrom(numax),growprev(numax),grownum(numax)
      dimension growlist(numax,numax)
      real(KIND=double_precision)::weight,  weiold, voldt,voldbb,voldtg ,voldext,voldintra,voldinter,voldbvib,voldelect,voldewald,vnewt ,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib,vnewinter,vnewelect ,vnewewald,vneworient,vnewinterr,vnewextr,vnewelectr,voldinterr ,voldextr,voldelectr,voldorient
!      common /rosena/ growfrom,growprev,growlist,grownum
!      common /rosenb/ weight,  weiold, voldt,voldbb,voldtg,voldext
!     &     ,voldintra,voldinter,voldbvib,voldelect,voldewald
!     &     ,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib
!     &     ,vnewinter,vnewelect,vnewewald,voldelectr
!     &     ,vnewinterr,vnewextr,vnewelectr,voldinterr
!     &     ,voldextr,vneworient,voldorient

! IPSWPAR.INC
      integer(KIND=normal_int)::nw,nwell,iratipsw
      parameter (nw = 4000)
      real(KIND=double_precision)::dvdl,acdvdl,acipsw,vipsw,pipsw ,vwellipsw,pwellipsw,vwellipswb,vipswb,etais,lambdais,rxwell ,rywell,rzwell,awell,bwell,vipswo,vipswn,vwellipswo,vwellipswn ,lena,lenc,pwellips,pips,dhmat,sxwell,sywell,szwell,vwellipswot ,vwellipswnt,vipswnt,vipswot
      logical::lwell,lstagea,lstageb,lstagec
      dimension vipswb(nbxmax),vwellipswb(nbxmax),rxwell(nw,ntmax) ,rywell(nw,ntmax) ,rzwell(nw,ntmax),pwellips(3,3),pips(3,3),nwell(ntmax) ,lwell(ntmax),awell(numax,numax,ntmax),dhmat(3,3),sxwell(nw,ntmax),sywell(nw,ntmax),szwell(nw,ntmax) ,vwellipswot(nchmax),vwellipswnt(nchmax) ,vipswot(nchmax),vipswnt(nchmax)
!      common /paripsw1/ dvdl,acdvdl,acipsw,vipsw,pipsw,vwellipsw
!     &   ,pwellips,pips,dhmat,sxwell,sywell,szwell
!     &   ,vwellipswot,vwellipswnt,vipswnt,vipswot
!     &   ,pwellipsw,vwellipswb,vipswb,etais,lambdais,rxwell,rywell
!      common /paripsw2/ rzwell,nwell
!      common /paripsw5/ awell,bwell
!      common /paripsw3/ vipswo,vipswn,vwellipswo,vwellipswn,
!     &                lena,lenc,iratipsw
!      common /paripsw4/ lwell,lstagea,lstageb,lstagec

! EEPAR.INC
      logical::leemove,lmstate, leeacc
      integer(KIND=normal_int)::fmstate,sstate1,sstate2,mstate,nstate ,ee_moltyp,box_state,eepointp,eeirem,boxrem1,boxins1,ee_prob ,nmolty1
      real(KIND=double_precision)::wee_ratio,psi,ee_qqu,um_markov ,eeratio,rminee
      dimension psi(smax),ee_qqu(numax,smax),box_state(smax) ,ee_moltyp(smax),um_markov(smax,smax),ee_prob(smax) ,rminee(nntype*nntype)
!      common /expepar/ rminee,ee_qqu,psi,wee_ratio,lmstate,leemove
!     &      ,fmstate,sstate1,sstate2,mstate,nstate,ee_moltyp,box_state
!     &      ,eepointp,eeirem,boxrem1,boxins1,leeacc,um_markov,eeratio
!     &      ,ee_prob,nmolty1

! BLKAVG.INC
      integer(KIND=normal_int)::nener,nprop,blockm
      real(KIND=double_precision)::naccu,nccold,accum,bccold,aver,baver ,acsvol,acnvol,acshmat,acnhmat
! ** Neeraj adding for solubility parameter and heat of vaporization
      real(KIND=double_precision)::naccu1,nccold1,accum1,bccold1,aver1 ,baver1
      integer(KIND=normal_int)::nprop1
      parameter (nener=18,nprop=4+nener+(4*ntmax)+3,blockm=100)
      parameter (nprop1=11)
      dimension naccu(nprop,nbxmax),nccold(nprop,nbxmax)
      dimension accum(nprop,nbxmax),bccold(nprop,nbxmax)
      dimension aver(nprop,nbxmax), baver(nprop,nbxmax,blockm)
      dimension acsvol(nbxmax),acnvol(nbxmax)
      dimension acshmat(nbxmax,9),acnhmat(nbxmax,9)
      dimension naccu1(nprop1,nbxmax,nbxmax), nccold1(nprop1,nbxmax,nbxmax), aver1(nprop1,nbxmax,nbxmax)
      dimension accum1(nprop1,nbxmax,nbxmax), bccold1(nprop1,nbxmax,nbxmax), baver1(nprop1,nbxmax,nbxmax,blockm)
      double precision, dimension (ntmax,nbxmax) :: acntray,acntrax, acntraz,acnrotx,acnroty,acnrotz,acstrax,acstray,acstraz, acsrotx,acsroty,acsrotz
!      common /bokavg/ baver,naccu,nccold,accum,bccold,aver,acsvol
!     &                ,acnvol,acshmat,acnhmat
!      common /blkavg/ baver1,naccu1,nccold1,accum1,bccold1,aver1
!      common /trarotavg/ acntrax,acntray,acntraz,acnrotx,acnroty,
!     &     acnrotz,acstrax,acstray,acstraz,acsrotx,acsroty,
!     &     acsrotz

! FIX.INC
      integer(KIND=normal_int)::  iend,ipast,endnum,pastnum,fclose ,fcount,iwbef,ibef,befnum,wbefnum,nextnum,inext
      logical::lcrank
      real(KIND=double_precision)::xx,yy,zz,xvec,yvec,zvec,distij,vtgtr ,vtbend,bsum_tor,vtvib
      dimension iend(numax),ipast(numax,numax) ,pastnum(numax) ,fclose(numax,numax) ,fcount(numax),iwbef(numax),ibef(numax ,numax) ,befnum(numax),xx(numax),yy(numax),zz(numax) ,xvec(numax ,numax),yvec(numax,numax) ,zvec(numax,numax),distij(numax,numax) ,nextnum(numax),inext(numax,numax),vtvib(nchmax) ,vtgtr(nchmax) ,vtbend(nchmax),bsum_tor(nchmax)
!      common /fecba/ iend,ipast,endnum,pastnum,fclose,fcount
!     &   ,befnum,wbefnum,iwbef,ibef,nextnum
!     &   ,inext,lcrank
!      common /fecbb/ xx,yy,zz,xvec,yvec,zvec,distij,vtgtr,vtbend
!     &   ,bsum_tor,vtvib

! BOLTZMANN.INC
      real(KIND=double_precision)::rxi1,ryi1,rzi1,vi1,wi1,vext1,velect1
!      common /boltzm/ rxi1,ryi1,rzi1,vi1,wi1,vext1,velect1

! EXPSIX.INC
      integer(KIND=normal_int)::natom
      parameter (natom = 3)
      real(KIND=double_precision)::aexsix,bexsix,cexsix,sexsix,consp ,consu
      dimension aexsix(natom),bexsix(natom),cexsix(natom) ,sexsix(natom),consu(natom),consp(natom)
!      common /stuff/ aexsix,bexsix,cexsix,sexsix,consu,consp

! MERCK.INC
      integer(KIND=normal_int)::natomtyp
      parameter (natomtyp = 3)
      real(KIND=double_precision)::epsimmff,sigimmff,smmff,sigisq,ammff ,nmmff
      real(KIND=double_precision)::gmmff,alphammff,coru_cons,corp_cons
      dimension epsimmff(natomtyp),sigimmff(natomtyp),smmff(natomtyp)
      dimension ammff(natomtyp),nmmff(natomtyp),gmmff(natomtyp)
      dimension sigisq(natomtyp),alphammff(natomtyp)
      dimension coru_cons(natomtyp),corp_cons(natomtyp)
!      common /merk/ epsimmff,sigimmff,smmff,sigisq,ammff,nmmff
!     &   ,alphammff,gmmff,coru_cons,corp_cons

! INPAR.INC
      logical::lbranch
      integer(KIND=normal_int)::ininch,inix, iniy, iniz, inirot, inimix
      real(KIND=double_precision)::zshift, dshift
      dimension inix(nbxmax),iniy(nbxmax),iniz(nbxmax), inirot(nbxmax),inimix(nbxmax)
      dimension zshift(nbxmax),dshift(nbxmax)
      dimension ininch(ntmax,nbxmax)
      dimension lbranch(ntmax)
!      common /input2/ zshift,dshift,inix,iniy,iniz,inirot,inimix
!     &               ,ininch,lbranch

! GAROFALINI.INC
      integer(KIND=normal_int)::pair_max
      parameter(pair_max=50000)
      real(KIND=double_precision)::galpha(6),grho(6),gbeta(6),ga(6,3), 	 gb(6,3),gc(6,3),glambda(4),grij(4,2),ggamma(4,2), gtheta(4),grijsq(4,2)
      real(KIND=double_precision)::dxij(pair_max),dyij(pair_max), 	 dzij(pair_max),dij(pair_max),dik(pair_max),dxik(pair_max) ,dyik(pair_max),dzik(pair_max)
      integer(KIND=normal_int)::itr1(pair_max),itr2(pair_max) ,itr3(pair_max),ntr,tagged
      real(KIND=double_precision)::v3garo
!      common /garo/ ga,gb,gc,ggamma,grij,galpha,grho,gbeta,
!     &              glambda,gtheta,grijsq,tagged
!      common /Ddist/ dij,dxij,dyij,dzij,dik,dxik,dyik,dzik
!      common /triads/ itr1,itr2,itr3,ntr
!      common /garoenrgy/ v3garo

! EXTERNALMUIR.INC
      real(KIND=double_precision)::sigpri,c9ch2,c3ch2,c9ch3,c3ch3 ,zprmin,v2prmin,v3prmin,betac2,betac3
!      common /constmuir/ sigpri,c9ch2,c3ch2,c9ch3,c3ch3
!     &                ,zprmin,v2prmin,v3prmin,betac2,betac3

! EXPAND.INC
      real(KIND=double_precision)::epsil(ntmax,numax,100),sigm(ntmax ,numax,100)
      real(KIND=double_precision)::qcharge(ntmax,numax,100)
!      common /expand2/ epsil,sigm,qcharge

! LJSAMIPARA.INC
      real(KIND=double_precision)::sij(9), eij(9), vsh(9), vsha(9)
!      common /ljspar/ sij,eij,vsh,vsha

! NRTAB.INC
      integer(KIND=normal_int)::nrtcon,nrtmax,nrtbin
! - parameters currently set low to minimize storage
      parameter(nrtcon=50,nrtmax=10,nrtbin=10)
      real(KIND=double_precision)::dmrtab(nrtmax),dbrtab(nrtmax), wnrtab(nrtmax,nrtbin)
!      common /nrtabc/ dmrtab,dbrtab,wnrtab

! CONTORSION.INC
      integer(KIND=normal_int)::ntormax
      parameter (ntormax=700)
      real(KIND=double_precision)::vtt0(ntormax),vtt1(ntormax) ,vtt2(ntormax),vtt3(ntormax),vtt4(ntormax),vtt5(ntormax) ,vtt6(ntormax),vtt7(ntormax),vtt8(ntormax),vtt9(ntormax)
!      common /contor/ vtt0,vtt1,vtt2, vtt3,vtt4,vtt5,vtt6,vtt7,vtt8, vtt9

      type(LookupTable)::atoms !,bonds,angles,dihedrals

    CONTAINS
      subroutine checkAtom()
        if (.not.allocated(atoms%list)) call cleanup(TRIM(__FILE__)//integer_to_string(__LINE__)//": ATOMS section has not been defined!")
      end subroutine checkAtom
    end module global_data
