      subroutine ee_energy ( i,imolty, v, vintra, vinter,vext
     &     ,velect,vewald,vtail,flagon,ibox, istart, iend,lljii,ovrlap
     &     ,ltors,vtors,lcharge_table,lfavor)
 
!    *******************************************************************
!    ** calculates the total potential energy for a configuration.    **
!    *******************************************************************
 
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

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'system.inc'
!$$$      include 'neigh.inc'
!$$$      include 'poten.inc'
!$$$      include 'coord2.inc' 
!$$$      include 'cell.inc'
!$$$      include 'external.inc'
!$$$      include 'connect.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'qqlist.inc'
!$$$      include 'clusterbias.inc'
!$$$      include 'nsix.inc'
!$$$      include 'eepar.inc'

      logical::lqimol,lqjmol,lexplt,lcoulo,lfavor,lij2,liji,lqchgi
      logical::lljii,ovrlap,ltors,lcharge_table,lt,lfound

      integer(KIND=int)::growii,growjj,k,cellinc,jcell,ic,nmole
      integer(KIND=int)::i,ibox, istart, iend,ii,ntii,flagon,jjj,iii
     &       ,j,jj,ntjj,ntij,ntj,imolty,jmolty,ncell
      integer(KIND=int)::iivib,jjtor,ip1,ip2,ip3,it,nchp2,acellinc,kmolty

      real(KIND=double_precision)::ljsami,ljpsur,ljmuir,v,vintra, vinter,vext 
     &                ,rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq
     &                ,ryuij,rzuij,sr2,sr6,rij,rijsq,dzui,dz3,dz12
     &                ,exgrph,exsami,exmuir,exzeo,vtors,exsix,velect
     &                ,vewald,mmff,rbcut,ninesix,vharo,genlj
      real(KIND=double_precision)::erfunc,qave,rho,vol,vtail
      real(KIND=double_precision)::xvec,yvec,zvec,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2
     &     ,daa1,da1a2,dot,thetac,vtorso,coru
      real(KIND=double_precision)::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,epsilon2,sigma2

      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)
      dimension lcoulo(numax,numax),cellinc(cmax),jcell(nmax)
      dimension acellinc(numax,27)

! --------------------------------------------------------------------

!      write(iou,*) 'start ENERGY'
      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut  = rcut(ibox)
 
      if (ldual) rcinsq = rcutin*rcutin

!      rminsq = rmin * rmin
      v = 0.0d0
      vinter = 0.0d0
      vintra = 0.0d0
      vext = 0.0d0
      velect = 0.0d0
      vewald = 0.0d0
      vtail = 0.0d0
      vtors = 0.0d0
      vharo = 0.0d0
      ovrlap = .false.
      if ( istart .eq. 1 .and. flagon .eq. 2) then
         neigh_icnt = 0
      end if

! *******************************
! *** INTERCHAIN INTERACTIONS ***
! *******************************
      lqimol = lelect(imolty)

      if (nugrow(imolty) .eq. nunit(imolty)) then
         lexplt = .false.
      else
         lexplt = .true.
         growii = nugrow(imolty)
      end if

      if ( lcutcm .or. lfavor ) then
! --- calculate the center of mass of chain i and give it a dummy #
         nchp2 = nchain + 2
         do ii = 1,nunit(imolty)
            rxu(nchp2,ii) = rxuion(ii,flagon) 
            ryu(nchp2,ii) = ryuion(ii,flagon) 
            rzu(nchp2,ii) = rzuion(ii,flagon) 
!	write(iou,*) ii,flagon,rxu(nchp2,ii),ryu(nchp2,ii),rzu(nchp2,ii)
         end do
         nboxi(nchp2) = ibox
         moltyp(nchp2) = imolty
         call ctrmas(.false.,ibox,nchp2,9)
         xcmi = xcm(nchp2)
         ycmi = ycm(nchp2)
         zcmi = zcm(nchp2)
         rcmi = rcmu(nchp2)
!         write(iou,*) 'rcmi:',rcmi
      else
         lij2 = .true.
      end if

      if (licell.and.(ibox.eq.boxlink)) then
         do ii = istart, iend

            rxui = rxuion(ii,flagon)
            ryui = ryuion(ii,flagon)
            rzui = rzuion(ii,flagon)

!     --- check perodic boundaries
            if (rxui.gt.boxlx(ibox)) then
               rxui = rxui - boxlx(ibox)
            elseif (rxui.lt.0) then
               rxui = rxui + boxlx(ibox)
            end if

            if (ryui.gt.boxly(ibox)) then
               ryui = ryui - boxly(ibox)
            elseif (ryui.lt.0) then
               ryui = ryui + boxly(ibox)
            end if

            if (rzui.gt.boxlz(ibox)) then
               rzui = rzui - boxlz(ibox)
            elseif (rzui.lt.0) then
               rzui = rzui + boxlz(ibox)
            end if
                        
            call linkcell(3,i,rxui,ryui,rzui,cellinc)
            
            do j = 1, 27
               acellinc(ii,j) = cellinc(j)
            end do
         end do
                  
         ncell = 0
         
         do j = 1, 27
            ncell = ncell + 1
            cellinc(j) = acellinc(1,j)
         end do

         if (abs(iend-istart).gt.0) then
            do ii = istart+1, iend
               
               do j = 1, 27
                  
                  ic = acellinc(ii,j)

                  lfound = .false.
                  do jj = 1, ncell

                     if (ic.eq.cellinc(jj)) then
                        lfound = .true.
                        goto 92
                     end if
                                          
                  end do
 92               continue
                  
                  if (.not.lfound) then
                     ncell = ncell + 1
                     cellinc(ncell) = ic
                  end if
               end do
            end do
         end if
         
!         lt = .true.
         nmole = 0
         do j = 1, ncell
            ic = cellinc(j)
            
            do k = 1, nicell(ic)

               nmole = nmole + 1
               jcell(nmole) = iucell(ic,k)

! *** what is this??? always compute interactions with chain number 1?
! solute? lt = solute? removing...
!               if (jcell(nmole).eq.1) then
!                  lt = .false.
!               end if
                  
            end do
         end do

!         if (lt) then
!            nmole = nmole + 1
!            jcell(nmole) = 1
!         end if
      else
         nmole = nchain
      end if

!     --- loop over all chains except i - not for grand can. with ibox=2 !
      if (.not.(lgrand.and.(ibox.eq.2))) then    
         do 98 k = 1, nmole

            if (licell.and.(ibox.eq.boxlink)) then
               j = jcell(k)
            else
               j = k
            end if
            
            jmolty = moltyp(j)
            lqjmol = lelect(jmolty)
            growjj = nugrow(jmolty)

! ### check for simulation box ###
            if ( ( ibox .eq. nboxi(j) ) .and. (i .ne. j )) then

               if ( lneigh ) then
                  if ( .not. lnn(j,i) ) goto 98
               end if

               if (lcutcm .or. lfavor) then
!              --- check if ctrmas within rcmsq
                  rxuij = xcmi - xcm(j)
                  ryuij = ycmi - ycm(j)
                  rzuij = zcmi - zcm(j)
!              --- minimum image the ctrmas pair separations ***
                  if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  rcm = rbcut + rcmi + rcmu(j)
                  rcmsq = rcm*rcm
                  if ( lfavor ) then
                     favor(j) = (rminsq/rijsq)**2*5.0d0
                     favor2(j) = rminsq/rijsq
                  end if
                  if ( rijsq .gt. rcmsq .and. lcutcm) then 
                     if ( lqimol .and. lqjmol .and. lchgall ) then
                        lij2 = .false.
                        goto 108
                     else
                        goto 98
                     end if
                  else
                     lij2 = .true.
                  end if
               end if

               if ( lcharge_table .and. (.not. lchgall) ) then
! --- called from CBMC and must set up charge-interaction table ---
                  do ii = 1,nugrow(imolty)
                     do jj = 1,nugrow(jmolty)
                        iii = leaderq(imolty,ii)
                        jjj = leaderq(jmolty,jj)
                        if ( iii .eq. ii .and. jjj .eq. jj ) then
                           rxuij = rxuion(ii,flagon) - rxu(j,jj)
                           ryuij = ryuion(ii,flagon) - ryu(j,jj)
                           rzuij = rzuion(ii,flagon) - rzu(j,jj)
                           if ( lpbc ) 
     &                          call mimage(rxuij,ryuij,rzuij,ibox)
                           rijsq = rxuij*rxuij + ryuij*ryuij 
     &                          + rzuij*rzuij
                           if ((rijsq .lt. rcutsq) .or. lijall) then
                              lcoulo(ii,jj) = .true.
                           else
                              lcoulo(ii,jj) = .false.
                           end if
                        end if
                     end do
                  end do
               end if

!              --- loop over all beads ii of chain i 
 108           do ii = istart, iend
                  ntii = ntype(imolty,ii)
                  liji = lij(ntii)
                  lqchgi = lqchg(ntii)
                  rxui = rxuion(ii,flagon)
                  ryui = ryuion(ii,flagon)
                  rzui = rzuion(ii,flagon)
                  
!                 --- loop over all beads jj of chain j 
                  do 97 jj = 1, nunit(jmolty) 
!                    --- check exclusion table
                     if ( lexclu(imolty,ii,jmolty,jj) ) goto 97
                     
                     ntjj = ntype(jmolty,jj)
                     if ( lij2 ) then
                        if ( (.not. (liji .and. lij(ntjj))) 
     &                       .and. 
     &                       (.not. (lqchgi .and. lqchg(ntjj)))) 
     &                       goto 97
                     else
                        if (.not. (lqchgi .and. lqchg(ntjj)))
     &                       goto 97
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
                     rminsq = rminee(ntij)*rminee(ntij)
                     
                     rxuij = rxui - rxu(j,jj)
                     ryuij = ryui - ryu(j,jj)
                     rzuij = rzui - rzu(j,jj)
                   
! *** minimum image the pair separations ***
                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox)
                   
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   
                     if ( rijsq .lt. rminsq .and. .not. 
     &                    (lexpand(imolty) .or. lexpand(jmolty))) then
                        ovrlap = .true.
!	if (leemove) then
!                        write(iou,*) 'inter ovrlap:',i,j
!                        write(iou,*) 'i xyz',rxui,ryui,rzui
!                        write(iou,*) 'j xyz',rxu(j,jj),ryu(j,jj),rzu(j,jj) 
!                        write(iou,*) 'ii:',ii,'jj:',jj
!                        write(iou,*) 'distance', dsqrt(rijsq)
!	end if
                        return
                     end if
                     if ( (rijsq .lt. rcutsq) .or. lijall) then
                        if ( lsami ) then
                           vinter = vinter + ljsami(rijsq,ntij)
                        elseif (lexpsix) then
                           vinter = vinter + exsix(rijsq,ntij)
                        elseif (lmmff) then
                           vinter = vinter + mmff(rijsq,ntij)
                        elseif (lninesix) then
                           vinter = vinter + ninesix(rijsq,ntij)
                        elseif (lgenlj) then
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2=epsij(ntij)
                           vinter = vinter + genlj(rijsq,sr2,epsilon2)
                        elseif ( lmuir ) then
                           vinter = vinter + ljmuir(rijsq,ntij)
                        elseif ( lpsurf ) then
                           vinter = vinter + ljpsur(rijsq,ntij)
                        else if (lshift) then
                           sr2 = sig2ij(ntij) / rijsq
                           sr6 = sr2 * sr2 * sr2
                           vinter = vinter + sr6*(sr6-1.0d0)
     &                          *epsij(ntij)-ecut(ntij) 
                        elseif ( lij(ntii) .and. lij(ntjj) ) then
                           if ( lfepsi ) then
                              sr6 = rijsq*rijsq*rijsq
                              if ( (.not. lqchg(ntii)) .and. 
     &                             (.not. lqchg(ntjj)) ) then
                                 if ( nunit(imolty) .eq. 4 ) then
! *** TIP-4P structure (temperary use ???)
                                    qave = (qquion(4,flagon)
     &                                   +qqu(j,4))/2.0d0
                                 else
                                    qave=(qquion(4,flagon)
     &                                   +qquion(5,flagon)
     &                                   +qqu(j,4)+qqu(j,5))*0.85d0
                                 end if
                              else
                                 qave = (qquion(ii,flagon)
     &                                +qqu(j,jj))/2.0d0
                              end if
                              if ( lexpand(imolty) 
     &                             .and. lexpand(jmolty)) then
                                 epsilon2=dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif (lexpand(imolty)) then
                                 epsilon2=dsqrt(epsilon(imolty,ii)
     &                                *epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 epsilon2=dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 epsilon2 = epsij(ntij)
                              end if
                              vinter = vinter + 
     &                             ((aslope*(qave-a0)*(qave-a0)
     &                             +ashift)/sr6 - (bslope*(qave-
     &                             b0)*(qave-b0)+bshift))/
     &                             sr6*epsilon2
                           else
                              if ( lexpand(imolty) 
     &                             .and. lexpand(jmolty)) then
                                 sigma2=(sigma(imolty,ii)+
     &                                sigma(jmolty,jj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif ( lexpand(imolty) ) then
                                 sigma2=(sigma(imolty,ii)+
     &                                sigi(ntjj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsilon(imolty,ii)*
     &                                epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 sigma2=(sigma(jmolty,jj)+
     &                                sigi(ntii))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2 = epsij(ntij)
                              end if
                              sr6 = sr2 * sr2 * sr2
                              vinter = vinter 
     &                             + sr6*(sr6-1.0d0)*epsilon2
                           end if
                           
                        end if
                     end if

!                    --- electrostatics
                     if ( lchgall .and. lqchg(ntii) 
     &                    .and. lqchg(ntjj) ) then
                        if ( lewald ) then
                           velect = velect + qquion(ii,flagon)*qqu(j,jj)
     &                          *erfunc(calp(ibox)*dsqrt(rijsq))
     &                          /dsqrt(rijsq)
                        else
                           velect = velect + qquion(ii,flagon)
     &                          *qqu(j,jj)/dsqrt(rijsq)
                        end if
                     elseif ( lqimol .and. lqjmol .and. lqchg(ntii) 
     &                       .and. lqchg(ntjj) ) then

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
                              velect = velect + qquion(ii,flagon)*
     &                             qqu(j,jj)*erfunc(calp(ibox)*
     &                             dsqrt(rijsq))/dsqrt(rijsq)
                           else
                              velect = velect + qquion(ii,flagon)
     &                             *qqu(j,jj)/dsqrt(rijsq)
                           end if
                        end if
                     end if

                     if ( lneighbor .and. ii .eq. 1 .and. 
     &                    jj .eq. 1 .and. flagon .eq. 2
     &                    .and. rijsq .lt. rbsmax**2 
     &                    .and. rijsq .gt. rbsmin**2) then
!                           neigh_icnt1(jmolty)=neigh_icnt1(jmolty)+1
!                           neighi1(neigh_icnt1(jmolty),jmolty)=j
                        neigh_icnt=neigh_icnt+1
                        neighi(neigh_icnt)=j
                     end if

! *** additional repulsion between negative charge on the aromatic
! *** ring and the positive H in e.g. phenol

 97               continue
               end do
           end if
 98      continue
      end if

      if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff
     & .and. .not. lgenlj .and. .not. lninesix ) then
         vinter = 4.0d0 * vinter
      end if
      vinter = vinter+vharo

! ################################################################

! * the intramolecular van der waals and ewald terms have to be calculated 
! for the explicit atom placement models
! *******************************
! *** INTRACHAIN INTERACTIONS ***
! *******************************
      
! *** for expanded ensemble

   
! --- calculate intramolecular energy correction for chain i 
      do ii = istart, iend

         ntii = ntype(imolty,ii)
         rxui = rxuion(ii,flagon)
         ryui = ryuion(ii,flagon)
         rzui = rzuion(ii,flagon)
            
         do jj = 1,ii-1
               
            ntjj = ntype(imolty,jj)
            if ( lexpsix .or. lmmff ) then
               ntij = (ntii+ntjj)/2
            elseif (lninesix) then
               ntij = (ntii-1)*nxatom + ntjj
            elseif (lgenlj) then
               ntij = (ntii-1)*nntype + ntjj
            else
               ntij = (ntii-1)*nntype + ntjj
            end if
            rminsq = rminee(ntij)*rminee(ntij)
            
            rxuij = rxuion(ii,flagon) - rxuion(jj,flagon)
            ryuij = ryuion(ii,flagon) - ryuion(jj,flagon)
            rzuij = rzuion(ii,flagon) - rzuion(jj,flagon)
               
            rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  
! * calculation of intramolecular electrostatics
            if ( lqinclu(imolty,ii,jj) ) then
                  
               if ( lchgall .and. lqchg(ntii)
     &              .and. lqchg(ntjj) ) then
                  if ( lewald ) then

!                        write(iou,*) 'energy including ii,jj',ii,jj,
!     &                       qqfact*qquion(ii,flagon)*qquion(jj,flagon)
!     &                       *erfunc(calp(ibox)*dsqrt(rijsq))
!     &                       /dsqrt(rijsq)

!                    * 1,5 and beyond which gets full interaction
                     if (linclu(imolty,ii,jj)) then
                        velect = velect + 
     &                       qquion(ii,flagon)*qquion(jj,flagon)
     &                       *erfunc(calp(ibox)*dsqrt(rijsq))
     &                       /dsqrt(rijsq)
!                    * scale 1,4 interactions by qscale2(imolty,ii,jj)
                     else
                        velect = velect + qscale2(imolty,ii,jj)*
     &                       qquion(ii,flagon)
     &                       *qquion(jj,flagon)*erfunc(calp(ibox)
     &                       *dsqrt(rijsq))/dsqrt(rijsq)
                     end if
                  else
!                    * 1,5 and beyond which gets full interaction
                     if (linclu(imolty,ii,jj)) then
                        velect = velect + qquion(ii,flagon)
     &                       *qquion(jj,flagon)/dsqrt(rijsq)
!                    * scale 1,4 interactions by qscale2(imolty,ii,jj)
                     else
                        velect = velect + qscale2(imolty,ii,jj)*
     &                       qquion(ii,flagon)
     &                       *qquion(jj,flagon)/dsqrt(rijsq)
                     end if
                  end if

               elseif ( lqimol .and. lqchg(ntii)
     &                 .and. lqchg(ntjj) ) then

                  iii = leaderq(imolty,ii)
                  jjj = leaderq(imolty,jj)

                  if ( ii .eq. iii .and. jj .eq. jjj ) then
! --- set up the charge-interaction table
                     if ( rijsq .lt. rcutsq ) then
                        lcoulo(ii,jj) = .true.
                     else
                        lcoulo(ii,jj) = .false.
                     end if
                  end if

! *** set up table for neighboring groups- make sure they interact when 
! *** leaderqs are only 2 bonds apart.
                  if (.not. lqinclu(imolty,iii,jjj)) then
                     lcoulo(iii,jjj)  = .true.
                  end if

                  if ( lcoulo(iii,jjj) ) then
!                    * 1,5 and beyond which gets full interaction
                     if (linclu(imolty,ii,jj)) then

                        if ( lewald ) then
                           velect = velect + qquion(ii,flagon)*
     &                          qquion(jj,flagon)*erfunc(calp(ibox)*
     &                          dsqrt(rijsq))/dsqrt(rijsq)
                        else
                           velect = velect + qquion(ii,flagon)
     &                          *qquion(jj,flagon)/dsqrt(rijsq)
                        end if

!                    * scale 1,4 interactions by qscale2(imolty,ii,jj)
                     else

                        if ( lewald ) then
                           velect = velect + qscale2(imolty,ii,jj)*
     &                          qquion(ii,flagon)*
     &                          qquion(jj,flagon)*erfunc(calp(ibox)*
     &                          dsqrt(rijsq))/dsqrt(rijsq)
                        else
                           velect = velect + qscale2(imolty,ii,jj)*
     &                          qquion(ii,flagon)
     &                          *qquion(jj,flagon)/dsqrt(rijsq)
                        end if

                     end if
                  end if
               end if
               
            end if

! * calculation of other non-bonded interactions
            if ( linclu(imolty,ii,jj) ) then

               if (lljii) then

                  if ( rijsq .lt. rminsq .and. .not. 
     &                 lexpand(imolty)) then
                     ovrlap = .true.
                
!                     write(iou,*) 'intra ovrlap:',ii,jj
                     return
                  elseif ( rijsq .lt. rcutsq .or. lijall) then
                     if ( lsami ) then
                        vintra = vintra + ljsami(rijsq,ntij)
                     elseif (lexpsix) then
                        vintra = vintra + exsix(rijsq,ntij)
                     elseif (lmmff) then
                        vintra = vintra + mmff(rijsq,ntij)
                     elseif (lninesix) then
                        vintra = vintra + ninesix(rijsq,ntij)
                     elseif (lgenlj) then
                        sr2 = sig2ij(ntij) / rijsq
                        epsilon2=epsij(ntij)
                        vintra = vintra + genlj(rijsq,sr2,epsilon2)
                     elseif ( lmuir ) then
                        vintra = vintra + ljmuir(rijsq,ntij)
                     elseif ( lpsurf ) then
                        vintra = vintra + ljpsur(rijsq,ntij)
                     else if (lshift) then
                        sr2 = sig2ij(ntij) / rijsq
                        sr6 = sr2 * sr2 * sr2
                        vintra = vintra + 
     &                       sr6*(sr6-1.0d0)*epsij(ntij)-ecut(ntij) 
                     else
                        if ( lexpand(imolty) ) then
                           sigma2=(sigma(imolty,ii)+
     &                          sigma(imolty,jj))/2.0d0
                           sr2 = sigma2*sigma2/rijsq
                           epsilon2 = dsqrt(epsilon(imolty,ii)
     &                          *epsilon(imolty,jj))
                        else
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2 = epsij(ntij)
                        end if
                        sr6 = sr2 * sr2 * sr2
                        vintra = vintra + sr6*(sr6-1.0d0)
     &                       *epsilon2

! * OH 1-5 interaction
                             if (lainclu(imolty,ii,jj)) then
                                vintra = vintra + 0.25d0 * 
     &                               a15(a15type(imolty,ii,jj)) /
     &                               ((rijsq**2)*(rijsq**2)*(rijsq**2))
                             end if

                     end if
                  end if
               end if

            end if

            if (lewald ) then
!     compute the ewald intramolecular (self and correction) terms for 
!     the interactions of the placed atoms with themselves, and with the
!     rest of their own molecule, if there's no interaction
               rij = dsqrt(rijsq)
!              * these are 1,2 and 1,3
               if (.not. lqinclu(imolty,ii,jj)) then
                  vewald = vewald + qquion(ii,flagon)*qquion(jj,flagon)*
     &                 (erfunc(calp(ibox)*rij)-1.0d0)/rij
!              * 1,4 interaction which we scale by qscale2(imolty,ii,jj)
               elseif (lqinclu(imolty,ii,jj) .and. 
     &                 (.not. linclu(imolty,ii,jj))) then
                  vewald = vewald + (1.0d0-qscale2(imolty,ii,jj))*
     &                 qquion(ii,flagon)*
     &                 qquion(jj,flagon)*
     &                 (erfunc(calp(ibox)*rij)-1.0d0)/rij
               end if
            end if
         end do
         if ( lewald ) then
!           -- self interaction term, 1.772 is sqrt of pi
            vewald = vewald - qquion(ii,flagon)*qquion(ii,flagon)*
     &           calp(ibox)/1.772453851d0
         end if
      end do
      if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     & .and. .not. lgenlj  .and. .not. lninesix ) then 
           vintra = 4.0d0 * vintra 
      end if

      if (ltailc) then
!--     add tail corrections for the Lennard-Jones energy
          if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
             vol = (hmat(ibox,1) * (hmat(ibox,5) * hmat(ibox,9) -
     &            hmat(ibox,8) * hmat(ibox,6)) + hmat(ibox,4)
     &            * (hmat(ibox,8) * hmat(ibox,3) - hmat(ibox,2)
     &            * hmat(ibox,9)) + hmat(ibox,7) * (hmat(ibox,2)
     &            * hmat(ibox,6) - hmat(ibox,5)*hmat(ibox,3)))
          else
             vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
          end if
          do kmolty = 1, nmolty
             do jmolty = 1, nmolty
!                rho = ncmt(ibox,jmolty) /
!     &               ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
                if (flagon.eq.1) then
                   rho = ncmt(ibox,jmolty) / vol
                   vtail = vtail +
     &               ncmt(ibox,kmolty) * coru(kmolty,jmolty,rho,ibox)
!                write(iou,*) 'vtail',vtail
                else
                   if (jmolty.eq.ee_moltyp(mstate)) then
                      rho = (ncmt(ibox,jmolty)-1) / vol
                   elseif (jmolty.eq.imolty) then
                      rho = (ncmt(ibox,jmolty)+1) / vol
                   else
                      rho = ncmt(ibox,jmolty) / vol
                   end if
                   if (kmolty.eq.imolty) then
                      vtail = vtail +
     &                  (ncmt(ibox,kmolty)+1) * coru(kmolty,jmolty,rho
     &                   ,ibox)
                   elseif (kmolty.eq.ee_moltyp(mstate)) then
                      vtail = vtail +
     &                  (ncmt(ibox,kmolty)-1) * coru(kmolty,jmolty,rho
     &                    ,ibox)
                   else
                      vtail = vtail +
     &                  ncmt(ibox,kmolty) * coru(kmolty,jmolty,rho,ibox)
                   end if
                end if
             end do
          end do
!-----
          vinter = vinter + vtail
!----
       end if

! ################################################################

! ***************************************************************
! *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************
      
      if ( ljoe .or. lsami .or. lmuir .or. lexzeo
     &				      .or. lgraphite) then
! ---  not for grand can. with ibox=2 !
         if (.not.(lgrand.and.(ibox.eq.2))) then    
            do 399 j = istart,iend

               ntj = ntype(imolty,j)
 
               if ( ljoe ) then
                  if ( extc12(ntj) .gt. 0.1d0 ) then
                     dzui = rzuion(j,flagon) - extz0(ntj)
                     dz3  = dzui * dzui * dzui
                     dz12 = dz3**4
                     vext = vext + 
     &                    (extc12(ntj)/dz12) - (extc3(ntj)/dz3)
                  end if
               end if
	       
 	       if( lgraphite ) then
	       		ntij = (ntj-1)*nntype + ntsubst
			vext = vext + exgrph(rxuion(j,flagon),
     &		    	       ryuion(j,flagon),rzuion(j,flagon),ntij)
	       end if
	       
               if ( lsami ) vext = vext + exsami(rzuion(j,flagon),ntj)
               if ( lmuir ) vext = vext + exmuir(rzuion(j,flagon),ntj)

               if ( lexzeo ) vext = vext + exzeo(rxuion(j,flagon)
     &              ,ryuion(j,flagon),rzuion(j,flagon),ntj)
               
 399        continue

         end if
      end if

! *********************************************************************
! *** calculation of torsion energy for explicit atom methy groups ****
! *********************************************************************

      if ( ltors ) then
         do ii = 1, nunit(imolty)
            rxui=rxuion(ii,flagon)
            ryui=ryuion(ii,flagon)
            rzui=rzuion(ii,flagon)
            do iivib = 1, invib(imolty,ii)
               jj = ijvib(imolty,ii,iivib)
               xvec(ii,jj) = rxuion(jj,flagon) - rxui
               yvec(ii,jj) = ryuion(jj,flagon) - ryui
               zvec(ii,jj) = rzuion(jj,flagon) - rzui
!              --- account for explicit atoms in opposite direction
               xvec(jj,ii)   = -xvec(ii,jj)
               yvec(jj,ii)   = -yvec(ii,jj)
               zvec(jj,ii)   = -zvec(ii,jj)
            end do
         end do

         do j = nugrow(imolty)+1, nunit(imolty)
            do jjtor = 1, intor(imolty,j)
               ip3 = ijtor4(imolty,j,jjtor)
               if ( ip3 .lt. j ) then
                  ip1 = ijtor2(imolty,j,jjtor)
                  ip2 = ijtor3(imolty,j,jjtor)
                  it  = ittor(imolty,j,jjtor)
!*** calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                  xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     &                 zvec(ip1,j) * yvec(ip1,ip2)
                  yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     &                 xvec(ip1,j) * zvec(ip1,ip2)
                  zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     &                 yvec(ip1,j) * xvec(ip1,ip2)
                  xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     &                 zvec(ip1,ip2) * yvec(ip3,ip2)
                  ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     &                 xvec(ip1,ip2) * zvec(ip3,ip2)
                  za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     &                 yvec(ip1,ip2) * xvec(ip3,ip2)
! *** calculate lengths of cross products ***
                  daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                  da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
! *** calculate dot product of cross products ***
                  dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                  thetac = - dot / ( daa1 * da1a2 )
                  vtors = vtors + vtorso( thetac, it )
               end if
            end do
         end do
      end if

! ----------------------------------------------------------------------------
      
      velect = velect * qqfact
      vewald = vewald * qqfact

!     note that vintra is only computed when the flag lljii is true
      v = vinter + vext + vintra + velect + vewald 

!      write(iou,*) 'end ENERGY'

      return
      end





