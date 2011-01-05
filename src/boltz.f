      subroutine boltz( lnew,lfirst,ovrlap,i,icharge,imolty,ibox
     &     ,ichoi,iufrom,ntogrow,glist,maxlen)

!    *******************************************************************
!    ** calculates the potential energy and the boltzmann factor      **
!    ** for ichoi trial positions.   
!    *******************************************************************
!     lnew: true for new configurations
!     lfirst: true for insertion of the first bead in swap moves
!     ovrlap: logical variable, true for walk termination
!     i:
!     icharge:
!     imolty:
!     ibox:
!     ichoi: number of trial positions
!     iufrom:
!     ntogrow:
!     glist:
!     maxlen:
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
!$$$      include 'poten.inc'
!$$$      include 'neigh.inc'
!$$$      include 'cbmc.inc' 
!$$$      include 'external.inc'
!$$$      include 'connect.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'qqlist.inc'
!$$$      include 'nsix.inc'
!$$$      include 'cell.inc'
!$$$      include 'peboco.inc'
!$$$      include 'ipswpar.inc'
!$$$      include 'eepar.inc'
!$$$      include 'conver.inc'
!$$$      include 'tabulated.inc'
!$$$      include 'mpif.h'
!$$$      include 'mpi.inc'

      logical::lnew,ovrlap,lcmno(nmax),lfirst,lcompute
      logical::lqimol,lqjmol,liji,lqchgi
      integer(KIND=normal_int)::ichoi,growjj,igrow,count,glist(numax)
     & ,icharge,cnt,jcell(nmax),ic
      integer(KIND=normal_int)::i,imolty,ibox,ntogrow,itrial,ntii,j,jj
     & ,ntjj,ntij,iu,jmolty,iufrom,ii,cellinc(27),k,nmole
!      integer(KIND=normal_int)::NRtype 
      real(KIND=double_precision)::ljsami,rminsq,rxui,sr6,ryui,rzui
     & ,rxuij,ryuij,rzuij,rij,rijsq,sr2,dzui,dz3,dz12 ,exzeo,exsami
     & ,exmuir,exgrph,ljpsur,ljmuir,exsix ,mmff,maxlen,rcm,rcmsq ,corr
     & ,erfunc,rcutmax,ninesix, genlj
      real(KIND=double_precision)::vinter,vintra,vext,velect,vewald,qave
     & ,epsilon2,sigma2,vwell,v,rcutsq,rcinsq
      real(KIND=double_precision)::v_elect_field, field
      real(KIND=double_precision)::slitpore
      real(KIND=double_precision)::tabulated_bend, tabulated_vdW,
     & tabulated_elect
      integer(KIND=normal_int)::mmm

!------------- RP added for MPI
      integer(KIND=normal_int)::my_start,my_end,loops_per_proc
     & ,my_itrial
      real(KIND=double_precision)::my_vtry(nchmax),my_vtrintra(nchmax),
     & my_vtrext(nchmax),my_vtrinter(nchmax),my_vtrelect(nchmax)
     & ,my_vtrewald(nchmax),my_bfac(nchmax),my_vipswot(nchmax)
     & ,my_vwellipswot(nchmax),my_vipswnt(nchmax),my_vwellipswnt(nchmax)
      logical::my_lovr(nchmax)
      integer(KIND=normal_int)::ncount_arr(numprocmax+1)
     & ,ncount_displs(numprocmax+1)
! ------------------------------------------
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
!      write(iou,*) 'start BOLTZ'

!      lcompute = .false.

      if ( lpbc ) call setpbc(ibox)

!     --- determine the potential cutoffs

      rcutsq = rcut(ibox)*rcut(ibox)
      field = Elect_field(ibox)

! KM
! initialize variables
      do j=1,ichoi
         my_lovr(j) = .false.
         lovr(j) = .false.
         my_vtry(j) = 0.0d0
         my_vtrintra(j) = 0.0d0
         my_vtrelect(j) = 0.0d0
         my_vtrext(j) = 0.0d0
         my_vtrinter(j) = 0.0d0
         my_vtrewald(j) = 0.0d0
         my_bfac(j) = 0.0d0
         my_vipswot(j) = 0.0d0
         my_vwellipswot(j) = 0.0d0
         my_vipswnt(j) = 0.0d0
         my_vwellipswnt(j) = 0.0d0
      end do
      do j=1,numprocs
         ncount_arr(j) = 0
         ncount_displs(j) = 0
      end do
      

      if ( ldual ) then
!        --- use rcutin for both types of interactions (except intra)
         rcinsq = rcutin*rcutin
      else
!        --- compute the cutoffs squared for each interaction
         rcinsq =  rcutsq
         if ( lcutcm ) then
!           --- not needed when ldual is true since will use rcutin then
               rcutmax = rcut(ibox)
         end if
      end if

!     --- compute minimum cutoff squared
      rminsq = rmin * rmin

      lqimol = lelect(imolty)
      igrow = nugrow(imolty)

      if ( lcutcm .and. (.not. lfirst) ) then
!     --- check previous bead (iufrom) COM for each molecule in the box
         if ( lnew ) then
!           ### for trial chain ###
            rxui  = rxnew(iufrom)
            ryui  = rynew(iufrom)
            rzui  = rznew(iufrom)
         else
!           ### for old chain ###
            rxui  = rxu(i,iufrom)
            ryui  = ryu(i,iufrom)
            rzui  = rzu(i,iufrom)
         end if
         
         do j = 1,nchain
            lcmno(j) = .false.
            if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
               rxuij = rxui-xcm(j)
               ryuij = ryui-ycm(j)
               rzuij = rzui-zcm(j)
!                --- minimum image the pseudo-ctrmas pair separation
               if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
               rij  = dsqrt(rijsq)                

               if ( ldual ) then
                  rcm = rcutin + rcmu(j) + maxlen
                  rcmsq = rcm*rcm
               else
                  rcm = rcutmax + rcmu(j) + maxlen
                  rcmsq = rcm*rcm
               end if
               if (rijsq .gt. rcmsq ) lcmno(j) = .true.
            end if
         end do
      end if

! RP added for MPI

      loops_per_proc = ichoi/numprocs

      my_start = (myid*loops_per_proc)+1
      if(myid .eq. (numprocs-1))then
        my_end = ichoi
      else
         my_end   = (myid + 1)*loops_per_proc
      end if

! RP added for MPI
!      write(6,*)'170: boltz my_start=',my_start,'my_end=',my_end
!     &    ,'ichoi=',ichoi,'loops_per_proc=',loops_per_proc,'myid=',myid
!
      my_itrial  = 0
!      do itrial = 1, ichoi
      do itrial = my_start,my_end
         my_itrial  = my_itrial + 1
         my_lovr(my_itrial) = .false.
       
!         lovr(itrial) = .false.
         vinter = 0.0d0
         vintra = 0.0d0
         vext = 0.0d0
         velect = 0.0d0
         vewald = 0.0d0 

! -- if L_Coul_CBMC is true  only then compute electrostatic interactions/corrections
         if(L_Coul_CBMC) then
            do count = 1,ntogrow
               ii = glist(count)

               if (lewald) then
!              -- This part does not change for fixed charge moves, but is
!              -- used in the swap rosenbluth weight. - ewald self term
!              -- 1.772 is sqrt of pi
                  vewald = vewald - qqu(icharge,ii)*qqu(icharge,ii)
     &             *calp(ibox)/1.772453851d0
               end if
            
            end do
         end if

!        --- no intramolecular interactions if this is the first bead
         if ( .not. lfirst ) then

! *****************************************
! *** INTRACHAIN BEAD-BEAD INTERACTIONS ***
! *****************************************

!        --- cycle through molecule and check bead by bead
         do iu = 1, igrow

!           --- see if iu exists in the new chain yet
            if (.not. lexist(iu)) cycle
            
!           --- loop over all the grown beads
            do count = 1,ntogrow
               ii = glist(count)
!           --- see if iu has nonbonded intramolecular interaction with ii
               if ( linclu(imolty,ii,iu) .or. lewald ) then
!                 --- assign bead type for ii,iu, and the cross term
                  ntii = ntype(imolty,ii)
                  ntjj = ntype(imolty,iu)
                  if (lexpsix .or. lmmff ) then
                     ntij = (ntii+ntjj)/2
                  elseif (lninesix) then
                     ntij = (ntii-1)*nxatom + ntjj
                  elseif (lgenlj) then
                     ntij = (ntii-1)*nntype + ntjj
                  else
                     ntij = (ntii-1)*nntype + ntjj
                  end if
                  if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
!                 --- determine distances
                  if ( lnew ) then
!                    --- use new trial chain coordinates
                     rxuij  = rxnew(iu) - rxp(count,itrial)
                     ryuij  = rynew(iu) - ryp(count,itrial)
                     rzuij  = rznew(iu) - rzp(count,itrial)
                  else
!                    --- use old chain coordinates
                     rxuij  = rxu(i,iu) - rxp(count,itrial)
                     ryuij  = ryu(i,iu) - ryp(count,itrial)
                     rzuij  = rzu(i,iu) - rzp(count,itrial)
                  end if
                  if (lpbc) call mimage ( rxuij,ryuij,rzuij,ibox )
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  rij = dsqrt(rijsq)
               end if
               if ( linclu(imolty,ii,iu) .or.
     &              lqinclu(imolty,ii,iu)) then
                  if ( linclu(imolty,ii,iu) ) then
                     if ( rijsq .lt. rminsq .and. 
     &                    .not. lexpand(imolty) ) then
! RP added for MPI
                        my_lovr(my_itrial) = .true.

!                     write(iou,*) 'intra overlap'
                        goto 19
                     elseif ( rijsq .lt. rcutsq .or. lijall) then
                        if (L_vdW_table.or. L_bend_table.and.
     &                       (.not.lexpand(imolty))) then

                           do mmm=1,inben(imolty,ii)
                              if (ijben3(imolty,ii,mmm).eq.jj) then
                                 
                                 call lininter_bend(rij,
     &                                tabulated_bend, 
     &                                itben(imolty,ii,mmm))
                                 vintra = vintra + tabulated_bend
                                 
                                 goto 96
                              end if
                           end do

                           call lininter_vdW(rij, 
     &                          tabulated_vdW, ntii, ntjj)
                           vintra = vintra + tabulated_vdW

                        elseif (llj.and.(.not.(lexpand(imolty)
     &                       ))) then
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2=epsij(ntij)
                           sr6 = sr2 * sr2 * sr2
                           vintra = vintra 
     &                          + sr6*(sr6-1.0d0)*epsilon2
     &				* ljscale(imolty,ii,iu)
! * OH 1-5 interaction
                           if (lainclu(imolty,ii,iu)) then
                              vintra = vintra + 0.25d0 * 
     &                             a15(a15type(imolty,ii,iu)) /
     &                             ((rijsq**2)*(rijsq**2)*(rijsq**2))
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
                           vintra = vintra + genlj(rijsq,sr2,epsilon2)
                        elseif ( lmuir ) then
                           vintra = vintra + ljmuir(rijsq,ntij)
                        elseif ( lpsurf ) then
                           vintra = vintra + ljpsur(rijsq,ntij)
                        else if (lshift) then
                           sr2 = sig2ij(ntij) / rijsq
                           sr6 = sr2 * sr2 * sr2
                           vintra = vintra + sr6 
     &                          *((sr6-1.0d0)*epsij(ntij)-ecut(ntij))
     &				*ljscale(imolty,ii,iu)
                        else
                           if ( lexpand(imolty) ) then
                              sigma2 = (sigma(imolty,ii)
     &                             +sigma(imolty,iu))/2.0d0
                              sr2 = sigma2*sigma2 / rijsq
                              epsilon2 = dsqrt(epsilon(imolty,ii)
     &                             *epsilon(imolty,iu))
                           else
                              sr2 = sig2ij(ntij) / rijsq
                              epsilon2 = epsij(ntij)
                           end if
                           sr6 = sr2 * sr2 * sr2
                           vintra = vintra + sr6*(sr6-1.0d0)*epsilon2
     &					*ljscale(imolty,ii,iu)		   
! * OH 1-5 interaction
                           if (lainclu(imolty,ii,iu)) then
                              vintra = vintra + 0.25d0 * 
     &                             a15(a15type(imolty,ii,iu)) /
     &                             ((rijsq**2)*(rijsq**2)*(rijsq**2))
                           end if

                        end if
                     end if
                  end if

! intramolecular charge interaction 
!                    --- compute velect (coulomb and ewald)
 96               if (L_Coul_CBMC) then
                  if ( lelect(imolty) .and.
     &                 lqinclu(imolty,ii,iu)) then

! *** boltz.f has problem to compute the electrostatic interactions
! *** in a group-based way because the leader q might not be grown at
! *** present, so it calculates electrostatic interaction not based on 
! *** group but on its own distance in SC, but should be corrected 
! *** later by calling energy subroutine.

                     lcompute = .false.
                     if ( rijsq .lt. rcinsq )
     &                    lcompute = .true.
!     --- if lcompute then compute the electrostatic energy 
                        if ( lcompute ) then
                           if (L_elect_table) then
                              call lininter_elect(rij, 
     &                             tabulated_elect, ntii, ntjj)
                              velect = velect + qscale2(imolty,ii,iu)*
     &                             qqu(icharge,ii)*
     &                             qqu(icharge,iu)*tabulated_elect
                           elseif (lewald) then
!                          --- compute real::space term of vewald
                              velect = velect + 
     &                          qscale2(imolty,ii,iu)*qqu(icharge,ii)
     &                             *qqu(icharge,iu)*
     &                             erfunc(calp(ibox)*rij)/
     &                             rij
!                 --- ewald sum correction term
                              corr = (1.0d0 - qscale2(imolty,ii,iu))*
     &                           qqu(icharge,ii)
     &                             *qqu(icharge,iu)*(erfunc(calp(ibox)
     &                             * rij)-1.0d0) /rij
                              vewald = vewald + corr
                           else
                              velect = velect + qscale2(imolty,ii,iu)
     &                             *qqu(icharge,ii)
     &                             *qqu(i,iu)/rij
                            end if
                        end if
                     end if
                     end if
! end charge calculation 

! will only add correction if lqinclu is false.
                  elseif ( lewald .and.L_Coul_CBMC) then
!                 --- ewald sum correction term
                     corr = qqu(icharge,ii)*qqu(icharge,iu)
     &                 *(erfunc(calp(ibox)
     &                 * rij)-1.0d0) /rij
                     vewald = vewald + corr
                  end if
            end do

         end do

         if (L_Coul_CBMC) then 
         if ( lewald .and. ntogrow .gt. 1) then
!          --- ewald sum correction term for interactions of the 
!          --- growing beads with each other
! * this is 1-3, so don't need to consult lqinclu
! * should change this since it corrects for all currently grown beads in this
! * step, which could at somepoint be further than 1-3 apart!!! (say for rigrot...)
            do cnt = 1,ntogrow-1
               iu = glist(cnt)
               do count = cnt+1,ntogrow
                  ii = glist(count)
                  ntii = ntype(imolty,ii)
                  
!     --- assign bead type for ii,iu, and the cross term
                  ntjj = ntype(imolty,iu)
                  ntij = (ntii-1)*nntype + ntjj
!     --- determine distances - use trial chain coordinates
                  rxuij  = rxp(cnt,itrial) - rxp(count,itrial)
                  ryuij  = ryp(cnt,itrial) - ryp(count,itrial)
                  rzuij  = rzp(cnt,itrial) - rzp(count,itrial)
                  
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  rij   = dsqrt(rijsq)                  
!     --- ewald sum correction term
                  corr = qqu(icharge,ii)*qqu(icharge,iu)
     &                 *(erfunc(calp(ibox)
     &                 * rij)-1.0d0) /rij
                  vewald = vewald + corr
               end do
            end do
         end if
         end if
 
         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     & .and. .not. lgenlj .and. .not. lninesix
     & .and..not.L_vdW_table.and..not.L_bend_table) then 
                 vintra = 4.0d0 * vintra
         end if
 
! 
        end if

!     grand-canonical: if ibox = 2 (ideal gas box) only intra-chain
! --- JLR 11-24-09 don't compute if lideal          
!         if ( .not. ((lgrand .and. ibox .eq. 2))) then 
         if ( .not. (lgrand .and. ibox .eq. 2) 
     &       .and. .not.lideal(ibox) ) then
! --- END JLR 11-24-09
            if (licell.and.(ibox.eq.boxlink)) then
!     --- we only use count = 1, the rest should be taken care of
!     --- in rintramax
               count = 1
               rxui = rxp(count,itrial)
               ryui = ryp(count,itrial)
               rzui = rzp(count,itrial)

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
                  
               nmole = 0
               do j = 1, 27
                  ic = cellinc(j)
                     
                  do k = 1, nicell(ic)

                     nmole = nmole + 1
                     jcell(nmole) = iucell(ic,k)
                             
                  end do
               end do
            else
               nmole = nchain                  
            end if

! *******************************
! *** INTERCHAIN INTERACTIONS ***
! *******************************
!        ---    loop over all chains except i
         do 98 k = 1, nmole

            if (licell.and.(ibox.eq.boxlink)) then
               j = jcell(k)
            else
               j = k
            end if

!     --- check for simulation box 
            if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
               if ( lneigh ) then
!                 --- check neighbor list
                  if ( .not. lnew ) then
                     if ( .not. lnn(j,i) ) goto 98
                  end if
               end if

               jmolty = moltyp(j)
               lqjmol = lelect(jmolty)
               growjj = nugrow(jmolty)

               if ( .not. lfirst ) then
!                 --- check COM table calculated above
                  if (lcmno(j) .and. lcutcm) 
     &                 goto 98
               end if
               
!              --- loop over all beads of molecule i grown this step
 108           do count = 1,ntogrow
!                 --- assign bead type for ii
                  ii = glist(count)
                  ntii = ntype(imolty,ii)
                  liji = lij(ntii)
                  lqchgi = lqchg(ntii)
                  
!                 --- assign positions to r*ui
                  rxui = rxp(count,itrial)
                  ryui = ryp(count,itrial)
                  rzui = rzp(count,itrial)
            
                  if ( lfirst .and. lcutcm ) then
!                    --- check if ctrmas within rcmsq
    
                     rxuij = rxui-xcm(j)
                     ryuij = ryui-ycm(j)
                     rzuij = rzui-zcm(j)

!                    --- minimum image the ctrmas pair separations
                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     rij   = dsqrt(rijsq)
!                    --- determine cutoff
                     if ( ldual ) then
!                       --- must be lfirst so no previous bead
                        rcm = rcutin + rcmu(j)
                     else
!                       --- standard lcutcm cutoff
                        rcm = rcutmax + rcmu(j)
                     end if

!                    --- check if interaction distance is greater than cutoff
                     rcmsq = rcm*rcm
                     if ( rijsq .gt. rcmsq )
     &                    goto 98
                  end if

!                 --- loop over all beads jj of chain j
                  do 97 jj = 1,nunit(jmolty)
               
!                    --- check exclusion table
                     if ( lexclu(imolty,ii,jmolty,jj) ) goto 97
                     
! *** start iswatch add-on ***
! *** is there a way to pull this out of the loops? ***
                     if ( liswatch ) then
                        if (j .eq. other) then
                              if ( .not. liswinc(jj,jmolty) ) then
!                                 write(iou,*) 'iSwatch-skipping:',jj
                                 goto 97
                              end if
                        end if
                     end if
! *** end iswatch add-on ***

                     ntjj = ntype(jmolty,jj)
                     if ( (.not. (liji .and. lij(ntjj))) 
     &                    .and. 
     &                    (.not. (lqchgi .and. lqchg(ntjj)))
     &                    .and..not.L_vdW_table) 
     &                    goto 97
                     if (lexpsix .or. lmmff) then
                        ntij = (ntii+ntjj)/2
                     elseif (lninesix) then
                        ntij = (ntii-1)*nxatom + ntjj
                     elseif (lgenlj) then
                        ntij = (ntii-1)*nntype + ntjj
                     else
                        ntij = (ntii-1)*nntype + ntjj
                     end if
                     if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

                     rxuij = rxui - rxu(j,jj)
                     ryuij = ryui - ryu(j,jj)
                     rzuij = rzui - rzu(j,jj)

!                    --- minimum image the pair separations ***
                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     rij   = dsqrt(rijsq)
!                    --- compute vinter (ex. lennard-jones)
                     if ( rijsq .lt. rminsq .and. .not. 
     &                    (lexpand(imolty) .or. lexpand(jmolty))) then
                        my_lovr(my_itrial) = .true.
!                        write(iou,*) 'j:',j,jj
!                        write(iou,*) 'rjsq:',rijsq,rminsq
                        goto 19
                     elseif (rijsq .lt. rcinsq .or. lijall) then
                        if (L_vdW_table.and.(.not.
     &                       (lexpand(imolty).or.lexpand(jmolty))))
     &                       then
                           call lininter_vdW(rij, 
     &                          tabulated_vdW, ntii, ntjj)
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
                        elseif ( lsami) then
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
                           vinter = vinter + sr6*(sr6-1.0d0)*epsij(ntij) 
     &                          - ecut(ntij) 
                        else if ( lij(ntii) .and. lij(ntjj) ) then
                           if (lfepsi) then
                              sr6 = rijsq*rijsq*rijsq
                              if ( (.not. lqchg(ntii)) .and. 
     &                             (.not. lqchg(ntjj)) ) then
                                 if ( nunit(imolty) .eq. 4 ) then
! *** TIP-4P structure (temperary use ???)
                                    qave = (qqu(i,4)+qqu(j,4))/2.0d0
                                 else
                                    qave = (qqu(i,4)+qqu(i,5)+
     &                                   qqu(j,4)+qqu(j,5))*0.85d0
                                 end if
                              else
                                 qave = (qqu(i,ii)+qqu(j,jj))/2.0d0
                              end if
                              if ( lexpand(imolty) 
     &                             .and. lexpand(jmolty)) then
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif (lexpand(imolty)) then
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 epsilon2 = dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 epsilon2 = epsij(ntij)
                              end if
                              vinter = vinter + ((aslope*(qave-a0)*
     &                             (qave-a0)+ashift)/sr6 - (bslope*
     &                             (qave-b0)*(qave-b0)+bshift))/
     &                             sr6*epsilon2
                           else
                              if ( lexpand(imolty) .and. 
     &                             lexpand(jmolty)) then
                                 sigma2 = (sigma(imolty,ii)+
     &                                sigma(jmolty,jj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif ( lexpand(imolty) ) then
                                 sigma2 = (sigma(imolty,ii)+
     &                                sigi(ntjj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 sigma2 = (sigma(jmolty,jj)+
     &                                sigi(ntii))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2 = dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2 = epsij(ntij)
                              end if
                              sr6 = sr2 * sr2 * sr2
                              vinter = vinter + 
     &                             sr6*(sr6-1.0d0)*epsilon2

                           end if
                        end if
                     end if

!                    --- compute velect (coulomb and ewald)
! * lcompute has not been set yet, this was wrong
!                     if ( lcompute.and.lqimol .and. lqjmol .and. 
!     &                    lqchg(ntii) .and. lqchg(ntjj) ) then
                     if (L_Coul_CBMC) then
                     if ( lqimol .and. lqjmol .and. 
     &                    lqchg(ntii) .and. lqchg(ntjj) ) then
! *** boltz.f has problem to compute the electrostatic interactions
! *** in a group-based way because the leader q might not be grown at
! *** present, so it calculates electrostatic interaction not based on 
! *** group but on its own distance in SC, but should be corrected 
! *** later by calling energy subroutine.
                        lcompute = .false.
                        if ( rijsq .lt. rcinsq ) 
     &                       lcompute = .true.
!     --- if lcompute then compute the electrostatic energy 
                        if ( lcompute ) then
                           if (L_elect_table) then
                              call lininter_elect(rij, 
     &                             tabulated_elect, ntii, ntjj)
                              velect = velect + qqu(icharge,ii)*
     &                             qqu(j,jj)*tabulated_elect
                           elseif (lewald) then
!                             --- compute real::space term of velect
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(j,jj)*
     &                             erfunc(calp(ibox)*rij)/
     &                             rij
                           else
!                             --- compute all electrostatic interactions
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(j,jj)/
     &                             rij
                           end if
                        end if
                     end if
                     end if 
 97               continue 
               end do
            end if
            
 98      continue
         
         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &   .and. .not. lgenlj   .and. .not. lninesix
     &   .and..not.L_vdW_table ) then
                 vinter = 4.0d0 * vinter
         end if
      end if         
! ################################################################

! **************************************************************
! *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************
 
! ---  not for grand can. with ibox=2 !
! -- required for histogram reweighting to work for monolayer 
! -- phase diagrams.
! -- not used for adsorption isotherms
!         if (.not. (lslit .and. ibox.eq.2)) then
      if (ibox .eq. 1) then
         if ( ljoe .or. lsami .or. lmuir .or. lexzeo
     &         .or. lgraphite .or. lslit ) then
            do count = 1,ntogrow
!              --- assign bead type for ii
               ii = glist(count)
               ntii = ntype(imolty,ii)
               rxui = rxp(count,itrial)
               ryui = ryp(count,itrial)
               rzui = rzp(count,itrial)
               if ( ljoe ) then
                  if ( extc12(ntii) .gt. 0.1d0 ) then
                     dzui = rzui - extz0(ntii)
                     dz3  = dzui * dzui * dzui
                     dz12 = dz3**4
                     vext = vext +
     &                    (extc12(ntii)/dz12) - (extc3(ntii)/dz3)    
                  end if
               end if
	        
	       if (lslit) then
	          ntij = (ntii-1)*nntype + ntsubst		  
! -- calculate interaction with surface at the bottom of the box		  
		  vext = vext+slitpore(rzui,ntij)
! -- calculate interaction with the surface at the top of the box
		  dzui = boxlz(ibox)-rzui
		  vext = vext +slitpore(dzui,ntij)
	       end if
	       	  
               if( lgraphite ) then
                  ntij = (ntii-1)*nntype + ntsubst
                  vext = vext + exgrph(rxui,ryui,rzui,ntij)
               end if

               if ( lsami ) vext = vext + exsami(rzui,ntii)
               if ( lmuir ) vext = vext + exmuir(rzui,ntii)
               if ( lexzeo ) vext = vext + exzeo(rxui,ryui,rzui,ntii)
            end do
         end if
	 end if


         if (.not.ldual) then
            if (lelect_field) then
               if(lelect(moltyp(i))) then
                  if (nboxi(i).eq.ibox) then
                     do count = 1,ntogrow
                        rzui = rzp(count,itrial)
                        vext = vext + v_elect_field(i,count,rzui,field)
                     end do
                  end if
               end if
               vext = vext * eXV_to_K
            end if 
         end if
         
! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
                                                                                
      vwell = 0.0d0
      if (lwell(imolty).and.lmipsw) then
                                                                                
         rxui = xcm(i)
         ryui = ycm(i)
         rzui = zcm(i)
         do j = 1, nwell(imolty)*nunit(imolty)
            k = j - int(j/nunit(imolty))*nunit(imolty)
            if (k.eq.0) k = nunit(imolty)
            rxuij = rxui-rxwell(j,imolty)
            ryuij = ryui-rywell(j,imolty)
            rzuij = rzui-rzwell(j,imolty)
            call mimage(rxuij,ryuij,rzuij,ibox)
            rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
            rcm = rcut(ibox)+rcmu(i)+maxlen
            rcmsq = rcm*rcm
            if (rijsq.lt.rcmsq) then
            do count = 1, ntogrow
               ii = glist(count)
               if (awell(ii,k,imolty).lt.1.0d-6) goto 666
               rxui = rxp(count,itrial)
               ryui = ryp(count,itrial)
               rzui = rzp(count,itrial)
               rxuij = rxui-rxwell(j,imolty)
               ryuij = ryui-rywell(j,imolty)
               rzuij = rzui-rzwell(j,imolty)
               call mimage(rxuij,ryuij,rzuij,ibox)
               rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
               vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
 666        end do
            end if
         end do
                                                                                
      end if
! ----------------------------------------------------------------------------
 
! *********************************************
! *** CALCULATION OF TOTAL POTENTIAL ENERGY ***
! *********************************************
!         write(23,*) 
!         write(23,*) 'wirting out total energy'
!         write(23,*) 'Lnew', lnew 
!         if (NRtype.eq.1) then
!             write(23,*) 'called from swap'
!         else
!             write(23,*) 'called from rigrot'
!         end if 


 19      if ( my_lovr(my_itrial) ) then
            my_bfac(my_itrial) = 0.0d0
         else
            if (.not.L_elect_table) then
               velect = velect*qqfact
               vewald = vewald*qqfact
            end if
            v = vinter+vintra+vext+velect+vewald
            if (.not.lnew) then
               my_vipswot(my_itrial) = v
               my_vwellipswot(my_itrial) = vwell
            else
               my_vipswnt(my_itrial) = v
               my_vwellipswnt(my_itrial) = vwell
            end if
                                                                                
            if (lstagea) then
               v = (1.0d0-lambdais*(1.0d0-etais))*v
            elseif (lstageb) then
               v = etais*v+lambdais*vwell
            elseif (lstagec) then
               v = (etais+(1.0d0-etais)*lambdais)*v+
     &                         (1.0d0-lambdais)*vwell
            end if

            my_vtry(my_itrial) = v
            my_vtrintra(my_itrial) = vintra
            my_vtrext(my_itrial)   = vext
            my_vtrinter(my_itrial) = vinter
            my_vtrelect(my_itrial) = velect
            my_vtrewald(my_itrial) = vewald
!            write(23,*) 'itrial' ,itrial
!            write(23,*) vtry(itrial), vtrintra(itrial), vtrext(itrial),
!     &       vtrinter(itrial),vtrelect(itrial), vtrewald(itrial)


            if ((my_vtry(my_itrial)*beta).gt.(2.3d0*softcut))then
!               write(iou,*) 'caught by softcut',vtry(itrial)*beta
               my_lovr(my_itrial) = .true.
               my_bfac(my_itrial) = 0.0d0
            elseif((my_vtry(my_itrial)*beta).lt.-2.303d0*308)then
!               write(iou,*) '### warning: weight too big out of range'
               my_lovr(my_itrial) = .true.
               my_bfac(my_itrial) = 0.0d0
            else
               my_bfac(my_itrial) = dexp ( -(my_vtry(my_itrial)*beta) )
            end if
         end if
      end do

!      scount = loops_per_proc
       loops_per_proc = (my_end-my_start) + 1

       CALL MPI_ALLGATHER(loops_per_proc,1,MPI_INTEGER,ncount_arr,
     &       1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       ncount_displs(1) = 0

! KM for MPI
! cannot loop over i, must use j
! i is an argument to the boltz subroutine!
       do j = 2,numprocs
           ncount_displs(j) = ncount_displs(j-1) + ncount_arr(j-1)
       end do

      call MPI_ALLGATHERV(my_lovr,loops_per_proc,MPI_LOGICAL
     &     ,lovr,ncount_arr,ncount_displs,
     &     MPI_LOGICAL,MPI_COMM_WORLD,ierr)

      ovrlap = .true.
      do itrial = 1, ichoi
         if ( .not. lovr(itrial)) ovrlap = .false.
      end do

      call MPI_ALLGATHERV(my_vtry, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtry, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrintra, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrintra, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrext,loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrext, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrinter, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrinter, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrelect,loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrelect, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrewald, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrewald, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_bfac, loops_per_proc, 
     &     MPI_DOUBLE_PRECISION, bfac, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

!       do my_itrial=my_start,my_end
!        write(6,*)"938:boltz my_bfac(",my_itrial,")=",my_bfac(my_itrial)
!     &      ,'myid=',myid
!       end do

      call MPI_ALLGATHERV(my_vipswot, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vipswot, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vwellipswot, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vwellipswot, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vipswnt, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vipswnt, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vwellipswnt, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vwellipswnt, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

! ----------------------------------------------------------------------------


!      write(iou,*) 'end BOLTZ'

      return

      end


