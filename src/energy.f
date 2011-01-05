      subroutine energy ( i,imolty, v, vintra, vinter,vext
     &     ,velect,vewald,flagon,ibox, istart,iuend,lljii,ovrlap
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

!$$$      include 'mpi.inc'
!$$$      include 'mpif.h'
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'system.inc'
!$$$      include 'neigh.inc'
!$$$      include 'poten.inc'
!$$$      include 'coord2.inc' 
!$$$      include 'external.inc'
!$$$      include 'connect.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'qqlist.inc'
!$$$      include 'clusterbias.inc'
!$$$      include 'nsix.inc'
!$$$      include 'peboco.inc'
!$$$      include 'cell.inc'
!$$$      include 'ipswpar.inc'
!$$$      include 'eepar.inc'
!$$$      include 'conver.inc'
!$$$!kea include for garofalini 3 body term
!$$$      include 'garofalini.inc'
!$$$      include 'tabulated.inc'
      logical::lqimol,lqjmol,lexplt,lcoulo,lfavor,lij2,liji,lqchgi
      logical::lljii,ovrlap,ltors,lcharge_table,lfound

      integer(KIND=normal_int)::growii,growjj,k,cellinc,jcell,ic,nmole
      integer(KIND=normal_int)::i,ibox, istart, iuend,ii,ntii,flagon,jjj
     & ,iii,mmm,j,jj,ntjj,ntij,ntj,imolty,jmolty,ncell
      integer(KIND=normal_int)::iivib,jjtor,ip1,ip2,ip3,it,nchp2
     & ,acellinc

      real(KIND=double_precision)::ljsami,ljpsur,ljmuir,v,vintra, vinter
     & ,vext,rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq,ryuij,rzuij,sr2
     & ,sr6,rij,rijsq,dzui,dz3,dz12,exgrph,exsami,exmuir,exzeo,vtors
     & ,exsix,velect,vewald,mmff,rbcut,ninesix, genlj,garofalini
      real(KIND=double_precision)::erfunc,qave
      real(KIND=double_precision)::rxvec,ryvec,rzvec,xaa1,yaa1,zaa1
     & ,xa1a2,ya1a2,za1a2,daa1,da1a2,dot,thetac,vtorso,vwell
      real(KIND=double_precision)::xcc,ycc,zcc,tcc,theta,spltor
      real(KIND=double_precision)::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq
     & ,epsilon2,sigma2
      real(KIND=double_precision)::slitpore,v_elect_field, field

      dimension rxvec(numax,numax),ryvec(numax,numax),rzvec(numax,numax)
      dimension lcoulo(numax,numax),cellinc(cmax),jcell(nmax)
      dimension acellinc(numax,27)
! KEA
      integer(KIND=normal_int)::neigh_j,neighj(maxneigh)
      real(KIND=double_precision)::ndijj(maxneigh),nxijj(maxneigh),
     & nyijj(maxneigh),nzijj(maxneigh)
! KM
      real(KIND=double_precision)::tabulated_vdW, tabulated_bend,
     & tabulated_elect
! Neeraj & RP added for MPI
      real(KIND=double_precision)::sum_velect, sum_vinter
      logical::all_ovrlap
! --------------------------------------------------------------------

!      write(iou,*) 'start ENERGY'
      if ( lpbc ) call setpbc (ibox)

!      lljii = .false.

      rcutsq = rcut(ibox) * rcut(ibox)

      rbcut = rcut(ibox)

      field = Elect_field(ibox)

      if (ldual) rcinsq = rcutin*rcutin

      rminsq = rmin * rmin

! KM for MPI
      all_ovrlap = .false.
      sum_velect = 0.0d0
      sum_vinter = 0.0d0

      v = 0.0d0
      vinter = 0.0d0
      vintra = 0.0d0
      vext = 0.0d0
      velect = 0.0d0
      vewald = 0.0d0
      vtors = 0.0d0
      ovrlap = .false.
      all_ovrlap = .false.
      sself  = 0.0d0
      correct = 0.0d0     
!kea
      v3garo = 0.0d0

      if ( istart .eq. 1 .and. flagon .eq. 2) then
         neigh_icnt = 0
      elseif(lgaro.and.flagon.eq.1) then
         neigh_j = 0
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
         do ii = istart, iuend

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

         if (abs(iuend-istart).gt.0) then
            do ii = istart+1, iuend
               
               do j = 1, 27
                  
                  ic = acellinc(ii,j)

                  lfound = .false.
                  do jj = 1, ncell

                     if (ic.eq.cellinc(jj)) then
                        lfound = .true.
                        exit
                     end if
                                          
                  end do
                  
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
! --- JLR 11-24-09  also don't loop if box is ideal gas
!      if (.not.(lgrand.and.(ibox.eq.2))) then
      if (.not.(lgrand.and.(ibox.eq.2)) .and. 
     &     .not.(lideal(ibox)) ) then    
! --- END JLR 11-24-09
! RP added for MPI
      do 98 k = myid+1, nmole,numprocs
!         do 98 k = 1, nmole

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
                  rij  = dsqrt(rijsq)
                  rcm = rbcut + rcmi + rcmu(j)
                  rcmsq = rcm*rcm
!                  write(6,*) rcm,rcmi,rcmu(j)
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
 108           do ii = istart, iuend
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
!kea
                     elseif (lgaro) then
                        if(ntii.eq.ntjj) then
                           ntij = ntii
                        else
                           ntij = ntii+ntjj+1
                        end if
                     else
                        ntij = (ntii-1)*nntype + ntjj
                     end if
                     if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
                     
                     rxuij = rxui - rxu(j,jj)
                     ryuij = ryui - ryu(j,jj)
                     rzuij = rzui - rzu(j,jj)
                   
! *** minimum image the pair separations ***
                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox)
                     
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     rij  = dsqrt(rijsq)
                     if ( rijsq .lt. rminsq .and. .not. 
     &                    (lexpand(imolty) .or. lexpand(jmolty))) then
                        ovrlap = .true.
!                        write(iou,*) 'inter ovrlap:',i,j, myid
!                        write(iou,*) 'i xyz',rxui,ryui,rzui
!                        write(iou,*) 'j xyz',rxu(j,jj),ryu(j,jj),rzu(j,jj) 
!                        write(iou,*) 'ii:',ii,'jj:',jj
!                        write(iou,*) 'distance', dsqrt(rijsq)
! RP added for MPI
!                        return
                        goto 99
                     end if
                     if ( (rijsq .lt. rcutsq) .or. lijall) then
                        if (L_vdW_table.and.(.not.(lexpand(imolty)
     &                       .or.lexpand(jmolty)))) then
                           call lininter_vdW(rij, tabulated_vdW,
     &                          ntii, ntjj) 
                           vinter = vinter + tabulated_vdW

                        elseif (llj.and.(.not.(lexpand(imolty).or.
     &                       lexpand(jmolty)))) then
                           if ( lij(ntii) .and. lij(ntjj) ) then
                              sr2 = sig2ij(ntij) / rijsq
                              epsilon2=epsij(ntij)
                              sr6 = sr2 * sr2 * sr2
                              vinter = vinter 
     &                             + sr6*(sr6-1.0d0)*epsilon2
                           end if
                        elseif ( lsami ) then
                           vinter = vinter + ljsami(rijsq,ntij)
                        elseif (lexpsix) then
                           vinter = vinter + exsix(rijsq,ntij)
                        elseif (lmmff) then
                           vinter = vinter + mmff(rijsq,ntij)
                        elseif (lninesix) then
                           vinter = vinter + ninesix(rijsq,ntij)
! generalized lennard jones potential
                        elseif ( lgenlj ) then
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2=epsij(ntij)
                           vinter = vinter + genlj(rijsq,sr2,epsilon2)
                        elseif ( lmuir ) then
                           vinter = vinter + ljmuir(rijsq,ntij)
                        elseif ( lpsurf ) then
                           vinter = vinter + ljpsur(rijsq,ntij)
!kea
                        elseif (lgaro) then
                           vinter = vinter + garofalini(rijsq,
     &                          ntij,qquion(ii,flagon),qqu(j,jj),i,j)
                           if(lshift) then
                              vinter = vinter-ecut(ntij)
                           end if
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
                     if(lgaro) then
!kea --- skip for garofalini; included in vinter

                     elseif ( lchgall .and. lqchg(ntii) 
     &                    .and. lqchg(ntjj) ) then
                        if ( lewald ) then
                           velect = velect + qquion(ii,flagon)*qqu(j,jj)
     &                          *erfunc(calp(ibox)*rij)
     &                          /rij
                        else
                           velect = velect + qquion(ii,flagon)
     &                          *qqu(j,jj)/rij
                        end if
                     elseif ( lqimol .and. lqjmol .and. lqchg(ntii) 
     &                       .and. lqchg(ntjj) ) then

                        if (lewald) then
                           if (rijsq.lt.rcutsq) then
                               velect = velect + qquion(ii,flagon)*
     &                         qqu(j,jj)*erfunc(calp(ibox)*
     &                          rij)/rij
                           end if  
                        else
                           iii = leaderq(imolty,ii)
                           jjj = leaderq(jmolty,jj)

                           if (ii .eq. iii .and. jj .eq. jjj) then
! --- set up the charge-interaction table
                              if ( rijsq .lt. rcutsq ) then
                                 lcoulo(ii,jj) = .true.
                              else
                                 lcoulo(ii,jj) = .false.
                              end if
                           end if

                           if ( lcoulo(iii,jjj) ) then
                              if (L_elect_table) then
                                 call lininter_elect(rij, 
     &                                tabulated_elect, ntii, ntjj)
                                 velect = velect + qquion(ii,flagon)*
     &                                qqu(j,jj)*tabulated_elect
                              else
                                 velect = velect + qquion(ii,flagon)
     &                                *qqu(j,jj)/rij
                              end if
                           end if
                        end if
                     end if
! KM lneighbor and lgaro does not work in parallel
                     if ( lneighbor .and. ii .eq. 1 .and. 
     &                    jj .eq. 1 .and. flagon .eq. 2
     &                    .and. rijsq .lt. rbsmax**2 
     &                    .and. rijsq .gt. rbsmin**2) then
!                           neigh_icnt1(jmolty)=neigh_icnt1(jmolty)+1
!                           neighi1(neigh_icnt1(jmolty),jmolty)=j
                        neigh_icnt=neigh_icnt+1
                        neighi(neigh_icnt)=j
                     elseif(lgaro) then
                        if((ntij.eq.4.and.rijsq.lt.grijsq(2,1)).or.
     &                       (ntij.eq.6.and.rijsq.lt.grijsq(3,1))) then
                           if(flagon.eq.2) then
                              neigh_icnt=neigh_icnt+1
                              neighi(neigh_icnt)=j
                              ndiji(neigh_icnt) = rij
                              nxiji(neigh_icnt) = rxuij
                              nyiji(neigh_icnt) = ryuij
                              nziji(neigh_icnt) = rzuij
                           elseif(flagon.eq.1) then
                              neigh_j = neigh_j+1
                              neighj(neigh_j) = j
                              ndijj(neigh_j) = rij
                              nxijj(neigh_j) = rxuij
                              nyijj(neigh_j) = ryuij
                              nzijj(neigh_j) = rzuij
                           end if
                        end if
                     end if
   
 97               continue
               end do
           end if
 98      continue
      end if

! RP added for MPI
!----- Returning from ovrlap--------------

 99   continue

! KM don't check overlap until after allreduce
!      if(ovrlap .eq. .true.)then
!         write(iou,*)'576: in energy ovrlap=',ovrlap,'myid=',myid
!      end if
! -----------------------------------------

      CALL MPI_ALLREDUCE(ovrlap,all_ovrlap,1,MPI_LOGICAL,MPI_LOR,
     &          MPI_COMM_WORLD,ierr)

       ovrlap = all_ovrlap
      if(ovrlap)then
!            write(iou,*)'630 in energy ovrlap=',ovrlap,'myid=',myid
          return 
      end if

      CALL MPI_ALLREDUCE(vinter, sum_vinter,1,MPI_DOUBLE_PRECISION,
     &          MPI_SUM,MPI_COMM_WORLD,ierr)

      CALL MPI_ALLREDUCE(velect, sum_velect,1,MPI_DOUBLE_PRECISION,
     &          MPI_SUM,MPI_COMM_WORLD,ierr)

      velect = sum_velect
      vinter = sum_vinter
! -----------------------------------------------

      if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff
     & .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro
     &     .and..not.L_vdW_table) then
         vinter = 4.0d0 * vinter
      end if

!kea - garo: add three body loop for intermolecular interactions
      if (lgaro) then
         tagged = i
         if(flagon.eq.2) then
            call triad_en(i,v3garo,neigh_icnt,neighi,ndiji,nxiji,
     &           nyiji,nziji,.true.)
         elseif(flagon.eq.1) then
            call triad_en(i,v3garo,neigh_j,neighj,ndijj,nxijj,nyijj,
     &           nzijj,.false.)
         end if
      end if

! ################################################################

! * the intramolecular van der waals and ewald terms have to be calculated 
! for the explicit atom placement models
! *******************************
! *** INTRACHAIN INTERACTIONS ***
! *******************************

! --- JLR 11-19-09 commenting this out, alway do mimage for intrachain      
! *** for expanded ensemble
!      lmim = .false. 
!      nchp2=nchain+2
!      mlen2 = rcmu(nchp2)*2d0 
!      if ( mlen2>boxlx(ibox) .or. mlen2>boxly(ibox) .or.
!     $     mlen2>boxlz(ibox)) lmim = .true.
! --- END JLR 11-19-09 
  
! --- calculate intramolecular energy correction for chain i 
      do ii = istart, iuend
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
            if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
            rxuij = rxuion(ii,flagon) - rxuion(jj,flagon)
            ryuij = ryuion(ii,flagon) - ryuion(jj,flagon)
            rzuij = rzuion(ii,flagon) - rzuion(jj,flagon)
! --- JLR 11-19-09 always do mimage for intrachain
!            if (lpbc .and. lmim) call mimage( rxuij,ryuij,rzuij,ibox)
            if (lpbc) call mimage( rxuij,ryuij,rzuij,ibox)
! --- END JLR 11-19-09 
            rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
            rij = dsqrt(rijsq)
! * calculation of intramolecular electrostatics
            if ( lqinclu(imolty,ii,jj) ) then
               if ( lchgall .and. lqchg(ntii)
     &              .and. lqchg(ntjj) ) then
                  if ( lewald ) then
                        velect = velect + qscale2(imolty,ii,jj)* 
     &                       qquion(ii,flagon)*qquion(jj,flagon)
     &                       *erfunc(calp(ibox)*rij)
     &                       /rij
                   else
                        velect = velect +qscale2(imolty,ii,jj)*
     &                       qquion(ii,flagon)
     &                       *qquion(jj,flagon)/rij
                   end if
               elseif ( lqimol .and. lqchg(ntii)
     &                 .and. lqchg(ntjj) ) then
                  if(lewald) then
                     velect = velect + qscale2(imolty,ii,jj)*
     &                           qquion(ii,flagon)*
     &                          qquion(jj,flagon)*erfunc(calp(ibox)*
     &                          rij)/rij
                  else
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
                       if (L_elect_table) then
                          call lininter_elect(rij,tabulated_elect,
     &                         ntii, ntjj)
                          velect = velect + qscale2(imolty,ii,jj)*
     &                         qquion(ii,flagon)*qquion(jj,flagon)*
     &                         tabulated_elect
                       else
                          velect = velect + qscale2(imolty,ii,jj)*
     &                         qquion(ii,flagon)
     &                         *qquion(jj,flagon)/rij
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
!     write(iou,*) 'intra ovrlap:',ii,jj
                     return
                  elseif ( rijsq .lt. rcutsq .or. lijall) then
                     if (L_vdW_table.or.L_bend_table.and.
     &                    (.not.(lexpand(imolty)))) then

                        do mmm=1,inben(imolty,ii)
                           if (ijben3(imolty,ii,mmm).eq.jj) then
                              
                              call lininter_bend(rij,
     &                             tabulated_bend, 
     &                             itben(imolty,ii,mmm))
                              vintra = vintra + tabulated_bend
                              
                              goto 96
                           end if
                        end do


                        call lininter_vdW(rij,tabulated_vdW,ntii,ntjj)
                        vintra = vintra + tabulated_vdW
                     elseif (llj.and.(.not.(lexpand(imolty)
     &                    ))) then
                        sr2 = sig2ij(ntij) / rijsq
                        epsilon2=epsij(ntij)
                        sr6 = sr2 * sr2 * sr2
                        vintra = vintra 
     &                       + sr6*(sr6-1.0d0)*epsilon2
     &				  *ljscale(imolty,ii,jj)
                        
!     * OH 1-5 interaction
                        if (lainclu(imolty,ii,jj)) then
                           vintra = vintra + 0.25d0 * 
     &                          a15(a15type(imolty,ii,jj)) /
     &                          ((rijsq**2)*(rijsq**2)*(rijsq**2))
                        end if
                     elseif ( lsami ) then
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
                        vintra = vintra + genlj (rijsq,sr2,epsilon2)
                     elseif ( lmuir ) then
                        vintra = vintra + ljmuir(rijsq,ntij)
                     elseif ( lpsurf ) then
                        vintra = vintra + ljpsur(rijsq,ntij)
                     else if (lshift) then
                        sr2 = sig2ij(ntij) / rijsq
                        sr6 = sr2 * sr2 * sr2
                        vintra = vintra + 
     &                       (sr6*(sr6-1.0d0)*epsij(ntij)-ecut(ntij)) 
     &			      *ljscale(imolty,ii,jj)
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
     &                       *epsilon2*ljscale(imolty,ii,jj)

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

 96         if (lewald ) then
!     compute the ewald intramolecular (self and correction) terms for 
!     the interactions of the placed atoms with themselves, and with the
!     rest of their own molecule, if there's no interaction
               rij = dsqrt(rijsq)
!              * these are 1,2 and 1,3
               if (.not. lqinclu(imolty,ii,jj)) then
                  correct=correct+qquion(ii,flagon)*qquion(jj,flagon)*
     &                 (erfunc(calp(ibox)*rij)-1.0d0)/rij
!              * 1,4 interaction which we scale by qscale
               else
                  correct=correct+(1.0d0-qscale2(imolty,ii,jj))*
     &                        qquion(ii,flagon)*
     &                 qquion(jj,flagon)*
     &                 (erfunc(calp(ibox)*rij)-1.0d0)/rij
               end if
            end if
         end do
         if ( lewald ) then
            sself = sself + qquion(ii,flagon)*qquion(ii,flagon)
         end if
      end do
      if (lewald) then
         sself = -sself * calp(ibox)/dsqrt(onepi)
         vewald = sself + correct  
      end if
      if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     & .and. .not. lgenlj  .and. .not. lninesix
     & .and..not.L_vdW_table) then 
           vintra = 4.0d0 * vintra 
      end if

! ################################################################

! ***************************************************************
! *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************
      
      if ( ljoe .or. lsami .or. lmuir .or. lexzeo
     &			   .or. lgraphite .or. lslit) then
! ---  not for grand can. with ibox=2 !
!         if (.not.(lgrand.and.(ibox.eq.2))) then    
         if (ibox .eq. 1) then
            do j = istart,iuend

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
! -- Carbon slitpore	       
	       if (lslit) then
	          ntij = (ntj-1)*nntype + ntsubst
! -- calculate interaction with surface at the bottom of the box		  		  
		  vext = vext + slitpore(rzuion(j,flagon),ntij)
! -- calculate interaction with the surface at the top of the box
		  dzui = boxlz(ibox)-rzuion(j,flagon)
		  vext = vext +slitpore(dzui,ntij)
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
               
            end do

         end if
      end if

! ********************************************************************
! *** calculate of interaction energy with external electric field ***
! *** added 06/24/07 by KM
! ********************************************************************
      if(lelect_field) then
        if(lelect(moltyp(i))) then
           if (nboxi(i).eq.ibox) then
              do j = 1,nunit(moltyp(i))
                 vext = vext + v_elect_field(i,j,rzuion(j,flagon),
     &                field)
              end do
           end if
        end if
        vext = vext * eXV_to_K
      end if

! *********************************************************************
! *** calculation of torsion energy for explicit atom methyl groups ****
! *********************************************************************

      if ( ltors ) then
         do ii = 1, nunit(imolty)
            rxui=rxuion(ii,flagon)
            ryui=ryuion(ii,flagon)
            rzui=rzuion(ii,flagon)
            do iivib = 1, invib(imolty,ii)
               jj = ijvib(imolty,ii,iivib)
               rxvec(ii,jj) = rxuion(jj,flagon) - rxui
               ryvec(ii,jj) = ryuion(jj,flagon) - ryui
               rzvec(ii,jj) = rzuion(jj,flagon) - rzui
!              --- account for explicit atoms in opposite direction
               rxvec(jj,ii)   = -rxvec(ii,jj)
               ryvec(jj,ii)   = -ryvec(ii,jj)
               rzvec(jj,ii)   = -rzvec(ii,jj)
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
                  xaa1 = ryvec(ip1,j) * rzvec(ip2,ip1) +
     &                 rzvec(ip1,j) * ryvec(ip1,ip2)
                  yaa1 = rzvec(ip1,j) * rxvec(ip2,ip1) +
     &                 rxvec(ip1,j) * rzvec(ip1,ip2)
                  zaa1 = rxvec(ip1,j) * ryvec(ip2,ip1) +
     &                 ryvec(ip1,j) * rxvec(ip1,ip2)
                  xa1a2 = ryvec(ip1,ip2) * rzvec(ip2,ip3) +
     &                 rzvec(ip1,ip2) * ryvec(ip3,ip2)
                  ya1a2 = rzvec(ip1,ip2) * rxvec(ip2,ip3) +
     &                 rxvec(ip1,ip2) * rzvec(ip3,ip2)
                  za1a2 = rxvec(ip1,ip2) * ryvec(ip2,ip3) +
     &                 ryvec(ip1,ip2) * rxvec(ip3,ip2)
! *** calculate lengths of cross products ***
                  daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                  da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
! *** calculate dot product of cross products ***
                  dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                  thetac = - dot / ( daa1 * da1a2 )
!     KEA -- added for extending range to +/- 180
!     additional definitions for torsions
                  if (L_tor_table) then
!     *** calculate cross product of cross products ***
                     xcc = yaa1*za1a2 - zaa1*ya1a2
                     ycc = zaa1*xa1a2 - xaa1*za1a2
                     zcc = xaa1*ya1a2 - yaa1*xa1a2
!     *** calculate scalar triple product ***
                     tcc = xcc*rxvec(ip1,ip2) + ycc*ryvec(ip1,ip2)
     &                    + zcc*rzvec(ip1,ip2)
                     theta = dacos(thetac)
                     if (tcc .lt. 0.0d0) theta = -theta

                     if (L_spline) then
                        call splint(theta,spltor,it)
                     elseif(L_linear) then
                        call lininter(theta,spltor,it)
                     end if

                     vtors = vtors + spltor
                  else
                     vtors = vtors + vtorso( thetac, it )
                  end if
               end if
            end do
         end do
      end if

! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
      vwell = 0.0d0
      if (lwell(imolty).and.lmipsw) then
         rxui = xcmi
         ryui = ycmi
         rzui = zcmi
         do j = 1, nwell(imolty)*nunit(imolty)
            k = j - int(j/nunit(imolty))*nunit(imolty)
            if (k.eq.0) k = nunit(imolty)
            rxuij = rxui-rxwell(j,imolty)
            ryuij = ryui-rywell(j,imolty)
            rzuij = rzui-rzwell(j,imolty)
            call mimage(rxuij,ryuij,rzuij,ibox)
            rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
            rcm = rcut(ibox)+rcmi
            rcmsq = rcm*rcm
            if (rijsq.lt.rcmsq) then
            do ii = 1, nunit(imolty)
               if (awell(ii,k,imolty).lt.1.0d-6) cycle
               rxui = rxuion(ii,flagon)
               ryui = ryuion(ii,flagon)
               rzui = rzuion(ii,flagon)
               rxuij = rxui-rxwell(j,imolty)
               ryuij = ryui-rywell(j,imolty)
               rzuij = rzui-rzwell(j,imolty)
               call mimage(rxuij,ryuij,rzuij,ibox)
               rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
               vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
            end do
            end if
         end do
      end if

! ----------------------------------------------------------------------------
     
      if (.not.L_elect_table) then
         velect = velect*qqfact
         vewald = vewald*qqfact
      end if

 
!     note that vintra is only computed when the flag lljii is true
      v = vinter + vext + vintra + velect + vewald + v3garo


!      write(iou,*) 'vinter:',vinter,'vext:',vext,'vintra:',vintra,'velect'
!     & ,velect,'vewald:',vewald,'v'


      if (flagon.eq.1) then
         vipswo = v
         vwellipswo = vwell
      else
         vipswn = v
         vwellipswn = vwell
      end if

      if (lmipsw) then                                                
         if (lstagea) then
             v = (1.0d0-lambdais*(1.0d0-etais))*v
         elseif (lstageb) then
             v = etais*v+lambdais*vwell
         elseif (lstagec) then
             v = (etais+(1.0d0-etais)*lambdais)*v+(1.0d0-lambdais)*vwell
         end if 
      end if 

!      write(iou,*) 'v :', v

!      write(iou,*) 'end ENERGY'

      return
      end





