      subroutine pressure ( press, surf, ibox )

!    *****************************************************
!    ** calculates the pressure for a configuration.    **
!    *****************************************************
 
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
 
!$$$      include 'mpi.inc'
!$$$      include 'mpif.h'
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'system.inc'
!$$$      include 'poten.inc' 
!$$$      include 'expsix.inc'
!$$$      include 'merck.inc' 
!$$$      include 'connect.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'conver.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'qqlist.inc'
!$$$      include 'nsix.inc'
!$$$      include 'ipswpar.inc'
!$$$      include 'cell.inc'

      logical::lcoulo,lexplt,lqimol,lqjmol,lij2
      integer(KIND=normal_int)::ibox
      integer(KIND=normal_int)::i,ii,j,jj,ntii,ntjj,ntij,imolty,jmolty
     & ,iii,jjj,k
 
      real(KIND=double_precision)::press,repress,erfunc

      real(KIND=double_precision)::rxui,ryui,rzui,rxuij,ryuij,rzuij
     & ,rijsq,rcutsq,sr2,sr6,rhosq,corp,rij,rs1,sr1,rs2,sr7,rs7,rs6
      real(KIND=double_precision)::srij
      real(KIND=double_precision)::surf,pxx,pyy,pzz,rpxx,rpyy,rpzz,rpxy
     & ,rpyx,rpxz,rpzx,rpyz,rpzy

      real(KIND=double_precision)::fxcmi,fycmi,fzcmi,fij,xcmi,ycmi,zcmi
     & ,flj
      real(KIND=double_precision)::rcm,rcmsq,rcmi
      real(KIND=double_precision)::rbcut,qave, volsq,epsilon2,sigma2
     & ,pwell,vol 

      dimension lcoulo(numax,numax)

! RP added for MPI
      real(KIND=double_precision)::pips12,pips13,pips21,pips23,pips31
     & ,pips32
      real(KIND=double_precision)::diff_pips12,diff_pips13,diff_pips21
     & ,diff_pips23,diff_pips31,diff_pips32,diff_pxx,diff_pyy,diff_pzz
     & ,sum_press
! --------------------------------------------------------------------
      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut  = rcut(ibox)

      press = 0.0d0
      pxx = 0.0d0
      pyy = 0.0d0
      pzz = 0.0d0

      do i = 1, 3
         do j = 1, 3
            pips(i,j) = 0.0d0
         end do
      end do

! RP added for MPI
      pips12 = 0.0d0
      pips13 = 0.0d0
      pips21 = 0.0d0
      pips23 = 0.0d0
      pips31 = 0.0d0
      pips32 = 0.0d0
! KM for MPI
      diff_pips12 = 0.0d0
      diff_pips13 = 0.0d0
      diff_pips21 = 0.0d0
      diff_pips23 = 0.0d0
      diff_pips31 = 0.0d0
      diff_pips32 = 0.0d0
      diff_pxx = 0.0d0
      diff_pyy = 0.0d0
      diff_pzz = 0.0d0
      sum_press = 0.0d0
! ----------------------------------


! *******************************
! *** INTERCHAIN INTERACTIONS ***
! *******************************

      

!      if(LSOLPAR.and.(ibox.eq.2)) then
!         press= 1.380662d4 * ( ( nchbox(ibox) / beta) -
!     +     ( press/3.0d0 ) ) /
!     +     ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
!         surf = 0.0d0
!         return
!      end if


! --- loop over all chains i 
! RP added for MPI
      do 100 i = myid+1, nchain - 1,numprocs
!      do 100 i = 1, nchain - 1 
! ### check if i is in relevant box ###
         if ( nboxi(i) .eq. ibox ) then

            imolty = moltyp(i)
            lqimol = lelect(imolty)
            if ( nugrow(imolty) .eq. nunit(imolty) ) then
               lexplt = .false.
            else
               lexplt = .true.
            end if
            xcmi = xcm(i)
            ycmi = ycm(i)
            zcmi = zcm(i)
            if (lcutcm) then
               rcmi = rcmu(i)
            else
               lij2 = .true.
            end if

! --- loop over all chains j with j>i 
            do 99 j = i + 1, nchain
               
! ### check for simulation box ###
               if ( nboxi(j) .eq. ibox ) then
                  
                  jmolty = moltyp(j)
                  lqjmol = lelect(jmolty)
                  fxcmi = 0.0d0
                  fycmi = 0.0d0
                  fzcmi = 0.0d0
                  if ( lcutcm ) then
!                    --- check if ctrmas within rcmsq
                     rxuij = xcmi - xcm(j)
                     ryuij = ycmi - ycm(j)
                     rzuij = zcmi - zcm(j)
!                    --- minimum image the ctrmas pair separations
                     if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     rcm = rbcut + rcmi + rcmu(j)
                     rcmsq = rcm*rcm
                     if ( rijsq .gt. rcmsq ) then
                        if (lqimol .and. lqjmol .and. lchgall) then
                           lij2 = .false.
                           goto 108
                        else
                           goto 99
                        end if
                     else
                        lij2 = .true.
                     end if
                  end if


! --- loop over all beads ii of chain i 
 108              do 98 ii = 1, nunit(imolty)

                     ntii = ntype(imolty,ii)

                     rxui = rxu(i,ii)
                     ryui = ryu(i,ii)
                     rzui = rzu(i,ii)

! --- loop over all beads jj of chain j 
                     do 97 jj = 1, nunit(jmolty)
                        
                        ntjj = ntype(jmolty,jj)
                        if ( lij2 ) then
                           if ( (.not. (lij(ntii) .and. lij(ntjj))) 
     &                          .and. 
     &                          (.not. (lqchg(ntii) .and. lqchg(ntjj)))) 
     &                          goto 97
                        else
                           if (.not. (lqchg(ntii) .and. lqchg(ntjj)))
     &                          goto 97
                        end if
                        if ( lexpsix .or. lmmff ) then
                           ntij = (ntii+ntjj)/2
                        elseif (lninesix) then
                           ntij = (ntii-1)*nxatom + ntjj
                        elseif (lgenlj) then
                           ntij = (ntii-1)*nntype + ntjj
                        else
                           ntij = (ntii-1)*nntype + ntjj
                        end if

                        rxuij = rxui - rxu(j,jj)
                        ryuij = ryui - ryu(j,jj)
                        rzuij = rzui - rzu(j,jj)

! *** minimum image the pair separations ***
                        if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)

!     write(2,*) 'bead ruij',rxuij,ryuij,rzuij
                        
                        rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                        
! --- compute whether the charged groups will interact & fij
                        fij = 0.0d0
                        if ( lchgall ) then
                           if ( lewald ) then
                              fij = fij - (2.0d0*calp(ibox)*
     &                             dexp(-calp(ibox)*calp(ibox)
     &                             *rijsq)/dsqrt(onepi)
     &                             +erfunc(calp(ibox)*dsqrt(rijsq))/
     &                             dsqrt(rijsq))*qqu(i,ii)*qqu(j,jj)*
     &                             qqfact/rijsq
                           else
! -- KM 06/09/09
                              fij = -(qqfact*qqu(i,ii)*qqu(j,jj)
     &                             /dsqrt(rijsq))/rijsq
                           end if

                        elseif ( lqimol .and. lqjmol .and. 
     &                          lqchg(ntii) .and. lqchg(ntjj) ) then
                           
                           iii = leaderq(imolty,ii)
                           jjj = leaderq(jmolty,jj)
                            
                           if ( ii .eq. iii .and. jj .eq. jjj ) then
! --- set up the charge-interaction table
                              if ( rijsq .lt. rcutsq ) then
                                 lcoulo(ii,jj) = .true.
                              else
                                 lcoulo(ii,jj) = .false.
                              end if
                           end if
                           
                           if ( lcoulo(iii,jjj) ) then
                              if ( lewald ) then
                                 fij = fij - (2.0d0*calp(ibox)*
     &                                dexp(-calp(ibox)*calp(ibox)
     &                                *rijsq)/dsqrt(onepi)
     &                                +erfunc(calp(ibox)*dsqrt(rijsq))/
     &                                dsqrt(rijsq))*qqu(i,ii)*qqu(j,jj)*
     &                                qqfact/rijsq
                              else
                                 fij = -(qqfact*qqu(i,ii)*qqu(j,jj)
     &                                /dsqrt(rijsq))/rijsq
                              end if
                           end if
                        end if

                        if ( rijsq .lt. rcutsq .or. lijall) then
                           if ( lexpsix ) then
                              rij = dsqrt(rijsq)
                              fij = fij + (-6.0d0*aexsix(ntij)/(rijsq
     &                             *rijsq*rijsq) + bexsix(ntij)
     &                             *cexsix(ntij)*rij*dexp( cexsix(ntij)
     &				   *rij ))/rijsq
!                              write(2,*) 'rij,fij,ntij',rij,fij,ntij
                           elseif ( lmmff ) then
                              rs2 = rijsq/(sigisq(ntij))
                              rs1 = dsqrt(rs2)
                              rs6 = rs2*rs2*rs2
                              rs7 = rs1*rs6
                              sr1 = 1.07d0/(rs1+0.07d0)
                              sr7 = sr1**7.0d0
                              fij = fij - 7.0d0*epsimmff(ntij)*rs1*sr7*(
     &                           sr1*(1.12d0/(rs7+0.12d0)-2.0d0)/1.07d0+
     &                           1.12d0*rs6/((rs7+0.12d0)*(rs7+0.12d0)))
     &                             /rijsq
                           elseif (lninesix) then
                              rij = dsqrt(rijsq)
                              fij = fij + ( 72.0d0*epsnx(ntij)/
     &                            (rij*rzero(ntij)) ) * 
     &                            (rzero(ntij)/rij)**7 *
     &                            (1.0d0-(rzero(ntij)/rij)**3)
                           elseif (lgenlj) then
                              if ( lexpand(imolty)
     &                             .and. lexpand(jmolty) ) then
                                 sigma2=(sigma(imolty,ii)+
     &                                sigma(jmolty,jj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsilon(imolty,ii)
     &                                *epsilon(jmolty,jj))
                              elseif ( lexpand(imolty) ) then
                                 sigma2=(sigma(imolty,ii)+
     &                                sigi(ntjj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsilon(imolty,ii)
     &                                *epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 sigma2=(sigma(jmolty,jj)+
     &                                sigi(ntii))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsilon(jmolty,jj)
     &                                *epsi(ntii))
                              else
                                 rij = dsqrt(rijsq)
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2 = epsij(ntij)
                                 srij = dsqrt (sr2)
                              end if
                              if ( (rij) .le. (rij*srij)*
     &                             2.0d0**(2.0d0/n0) ) then
                                 flj=-4.0d0*epsilon2*((n0*((srij)**n0))-
     &                                ((n0/2.0d0)*((srij)**(n0/2.0d0))))
     &                                /rijsq
                              else
                                 flj =-epsilon2* ( ((2.0d0*n1) *
     &                                ((srij) **(2.0d0*n1))*
     &                                ( 2.0d0** ((4.0d0*n1)/n0)))-
     &                                (n1*((srij)**(n1))*
     &                                (2.0d0 ** ((2.0d0*n1/n0)+1.0d0))))
     &                                /rijsq
                              end if
                              fij = fij + flj
                           else
                              if ( lfepsi ) then
                                 sr2 = 1.0d0/rijsq
                                 sr6 = sr2 * sr2 * sr2
                                 if ( (.not. lqchg(ntii)) .and. 
     &                                (.not. lqchg(ntjj)) ) then
                                    if ( nunit(imolty) .eq. 4 ) then
! *** TIP-4P structure (temperary use ???)
                                       qave = (qqu(i,4)+qqu(j,4))/2.0d0
                                    else
                                       qave = (qqu(i,4)+qqu(i,5)+
     &                                      qqu(j,4)+qqu(j,5))*0.85d0
                                    end if
                                 else
                                    qave = (qqu(i,ii)+qqu(j,jj))/2.0d0
                                 end if

                                 if ( lexpand(imolty) 
     &                                .and. lexpand(jmolty)) then
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsilon(jmolty,jj))
                                 elseif ( lexpand(imolty)) then
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsi(ntjj))
                                 elseif ( lexpand(jmolty) ) then
                                    epsilon2=dsqrt(epsilon(jmolty,jj)
     &                                   *epsi(ntii))
                                 else
                                    epsilon2=epsij(ntij)
                                 end if
                                 flj = 48.0d0*epsilon2*
     &                                sr6*(-sr6*(aslope*(qave-a0)*
     &                                (qave-a0)+ashift)+0.5d0*
     &                                (bslope*(qave-b0)*(qave-b0)+
     &                                bshift)) / rijsq
                              else
                                 if ( lexpand(imolty) 
     &                                .and. lexpand(jmolty) ) then
                                    sigma2=(sigma(imolty,ii)+
     &                                   sigma(jmolty,jj))/2.0d0
                                    sr2 = sigma2*sigma2/rijsq
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsilon(jmolty,jj))
                                 elseif ( lexpand(imolty) ) then
                                    sigma2=(sigma(imolty,ii)+
     &                                   sigi(ntjj))/2.0d0
                                    sr2 = sigma2*sigma2/rijsq
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsi(ntjj))
                                 elseif ( lexpand(jmolty) ) then
                                    sigma2=(sigma(jmolty,jj)+
     &                                   sigi(ntii))/2.0d0
                                    sr2 = sigma2*sigma2/rijsq
                                    epsilon2=dsqrt(epsilon(jmolty,jj)
     &                                   *epsi(ntii))
                                 else
                                    sr2 = sig2ij(ntij) / rijsq
                                    epsilon2 = epsij(ntij)
                                 end if
                                 sr6 = sr2 * sr2 * sr2
                                 flj = 48.0d0*sr6*(-sr6+0.5d0)
     &                                *epsilon2 / rijsq
                              end if
                              fij = fij + flj
                           end if
                        end if
                        fxcmi = fxcmi + fij * rxuij 
                        fycmi = fycmi + fij * ryuij 
                        fzcmi = fzcmi + fij * rzuij 
 97                  continue
                      
 98               continue


! --- calculate distance between c-o-m ---
                  rxuij = xcmi - xcm(j)
                  ryuij = ycmi - ycm(j)
                  rzuij = zcmi - zcm(j)

! *** minimum image the pair separations ***
                  if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)

!                  write(2,*) 'COM  ruij',rxuij,ryuij,rzuij

                  press = press + 
     &                 fxcmi*rxuij + fycmi*ryuij + fzcmi*rzuij

! * for surface tension                  
! * this is correct for the coulombic part and for LJ.  Note sign difference!
                  pxx = pxx - fxcmi*rxuij
                  pyy = pyy - fycmi*ryuij
                  pzz = pzz - fzcmi*rzuij
                  pips12 = pips12 - rxuij*fycmi
                  pips13 = pips13 - rxuij*fzcmi
                  pips21 = pips21 - ryuij*fxcmi
                  pips23 = pips23 - ryuij*fzcmi
                  pips31 = pips31 - rzuij*fxcmi
                  pips32 = pips32 - rzuij*fycmi
                  
               end if
 99         continue
         end if
 100  continue

! RP for MPI

       CALL MPI_ALLREDUCE(press,sum_press,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pxx,diff_pxx,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pyy,diff_pyy,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pzz,diff_pzz,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pips12,diff_pips12,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pips13,diff_pips13,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pips21,diff_pips21,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pips23,diff_pips23,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pips31,diff_pips31,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pips32,diff_pips32,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

     
      press = sum_press
      pxx = diff_pxx
      pyy = diff_pyy
      pzz = diff_pzz
      pips(1,2) = diff_pips12
      pips(1,3) = diff_pips13
      pips(2,1) = diff_pips21
      pips(2,3) = diff_pips23
      pips(3,1) = diff_pips31
      pips(3,2) = diff_pips32
 
!      if(myid .eq. 0)then
!        write(6,*)'press=',press,'pxx=',pxx,'pyy=',pyy
!     &   ,'pzz=',pzz,'pips(1,2)',pips(1,2),'pips(1,3)=',pips(1,3)
!     &   ,'pips(2,1)=',pips(2,1),'pips(2,3)=',pips(2,3),'pips(3,1)=',
!     &   pips(3,1),'pips(3,2)=',pips(3,2)
!      end if

! ################################################################

      if ( lewald ) then

! *** Compute the reciprocal space contribution
! *** by using the thermodynamic definition
! *** 
         call recippress(ibox,repress,rpxx,rpyy,rpzz,rpxy,rpyx,rpxz,
     &                  rpzx,rpyz,rpzy)
         press = press - repress

         pxx = pxx + rpxx
         pyy = pyy + rpyy
         pzz = pzz + rpzz
         pips(1,2) = pips(1,2) + qqfact*rpxy
         pips(1,3) = pips(1,3) + qqfact*rpxz
         pips(2,1) = pips(2,1) + qqfact*rpyx
         pips(2,3) = pips(2,3) + qqfact*rpyz
         pips(3,1) = pips(3,1) + qqfact*rpzx
         pips(3,2) = pips(3,2) + qqfact*rpzy

      end if

      pips(1,1) = pxx
      pips(2,2) = pyy
      pips(3,3) = pzz

! Compute the Gaussian well contribution to the pressure
      pwell = 0.0d0
      do i = 1,3
         do j = 1,3
             pwellips(i,j) = 0.0d0
         end do
      end do


! KM for MPI - comment:
! this could likely be parallelized in the future
      do i = 1,nchain
         imolty = moltyp(i)
         if (lwell(imolty)) then
            rxui = xcm(i)
            ryui = ycm(i)
            rzui = zcm(i)
            do j = 1, nwell(imolty)*nunit(imolty)
               k = j - int(j/nunit(imolty))*nunit(imolty)
               if (k.eq.0) k = nunit(imolty)
               rxuij = rxui-rxwell(j,imolty)
               ryuij = ryui-rywell(j,imolty)
               rzuij = rzui-rzwell(j,imolty)
               call mimage (rxuij,ryuij,rzuij,ibox)
               rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
               rcm = rcut(ibox)+rcmu(i)
               rcmsq = rcm*rcm
               if (rijsq.lt.rcmsq) then
                  do ii = 1, nunit(imolty)
                     if (awell(ii,k,imolty).lt.1.0d-6) goto 666
                     rxui = rxu(i,ii)
                     ryui = ryu(i,ii)
                     rzui = rzu(i,ii)
                     rxuij = rxui-rxwell(j,imolty)
                     ryuij = ryui-rywell(j,imolty)
                     rzuij = rzui-rzwell(j,imolty)
                     call mimage (rxuij,ryuij,rzuij,ibox)
                     rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                     fij = 2.0d0*awell(ii,k,imolty)*bwell*
     &                    dexp(-bwell*rijsq)
                     pwell = pwell+fij*rijsq
                     pwellips(1,1) = pwellips(1,1)+fij*rxuij*rxuij
                     pwellips(2,2) = pwellips(2,2)+fij*ryuij*ryuij
                     pwellips(3,3) = pwellips(3,3)+fij*rzuij*rzuij
                     pwellips(1,2) = pwellips(1,2)+fij*rxuij*ryuij
                     pwellips(1,3) = pwellips(1,3)+fij*rxuij*rzuij
                     pwellips(2,3) = pwellips(2,3)+fij*ryuij*rzuij
 666              end do
               end if
            end do
         end if
      end do
      pwellips(2,1) = pwellips(1,2)
      pwellips(3,1) = pwellips(1,3)
      pwellips(3,2) = pwellips(2,3)

      if (lsolid(ibox).and.(.not.lrect(ibox))) then
         vol = cell_vol(ibox) 
      else
         vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
      end if
     
      pipsw = -(press/3.0d0)/vol
      pwellipsw = -(pwell/3.0d0)/vol

      do i = 1, 3
         do j = 1, 3
            pips(i,j) = pips(i,j)/vol
            pwellips(i,j) = -pwellips(i,j)/vol
         end do
      end do

! Chainging the order of calculation of tail correction
! So that it could be included in the scaled version for
! thermodynamic integration

!$$$      if (ltailc) then
!$$$c --     add tail corrections for the Lennard-Jones energy
!$$$c --     Not adding tail correction for the ghost particles
!$$$c --     as they are ideal (no interaction) Neeraj.
!$$$         volsq = ( vol )**2
!$$$         do imolty=1, nmolty
!$$$            do jmolty=1, nmolty
!$$$               rhosq = ncmt(ibox,imolty)*ncmt(ibox,jmolty)
!$$$     +              / volsq
!$$$               press=press + 1.380662d4 * corp(imolty,jmolty,rhosq,ibox)
!$$$            end do
!$$$        end do
!$$$      end if

      if (lstagea) then
         press = (1.0d0-lambdais*(1.0d0-etais))*press
      elseif (lstageb) then
         press = etais*press+lambdais*pwell
      elseif (lstagec) then
         press = (etais+(1.0d0-etais)*lambdais)*press
     &           +(1.0d0-lambdais)*pwell
      end if

      press = 1.380662d4 * ( ( (nchbox(ibox)+ ghost_particles(ibox))
     &      / beta) -
     &     ( press/3.0d0 ) ) / 
     &     ( vol )

      surf = pzz - 0.5d0*(pxx + pyy)
! * divide by surface area and convert from K to put surf in mN/m 
      surf = 1.380658d0*surf / (2.0d0*boxlx(ibox)*boxly(ibox))

!----check pressure tail correction

      if (ltailc) then
! --     add tail corrections for the Lennard-Jones energy
! --     Not adding tail correction for the ghost particles
! --     as they are ideal (no interaction) Neeraj.
         volsq = ( vol )**2
         do imolty=1, nmolty        
            do jmolty=1, nmolty  
               rhosq = ncmt(ibox,imolty)*ncmt(ibox,jmolty)
     &              / volsq
               press=press + 1.380662d4 * corp(imolty,jmolty,rhosq,ibox)
            end do
        end do
      end if
!      write (iou,*) ' press tail' ,  press

      return
      end











