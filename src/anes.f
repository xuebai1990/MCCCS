!    *********************************************************************
!    ** optimize the electronic configuration for trans, rot, and swap  **
!    ** (config, swatch in the future) moves and accept/reject the      **
!    ** combined move.                                                  **
!    ** written on June 25/99 by Bin Chen.                              **
!    *********************************************************************
subroutine anes(i,ibox,boxrem,mtype,laccept,deltv,vintern,vintran, vextn,velectn,vintero,vintrao,vexto,velecto,vinsta,vremta, vnewflucq,voldflucq,lswapinter)
  use util_random,only:random
  use sim_system
  use energy_kspace,only:recip
  implicit none

      logical::laccept,lswapinter
      integer::i,ibox,boxrem,mtype,imolty,iunit,ichoiq ,ip,ibox2,j
      real::deltv,vintern,vintran,vextn,velectn ,vintero,vintrao,vexto,velecto,vinsta,vremta,vnewflucq,voldflucq ,vboxo(nbxmax),vinterbo(nbxmax),vintrabo(nbxmax),vextbo(nbxmax) ,velectbo(nbxmax),vflucqbo(nbxmax),vtailbo(nbxmax),vvibbo(nbxmax) ,vtgbo(nbxmax),vbendbo(nbxmax),rxuo(numax),ryuo(numax) ,rzuo(numax),xcmo,ycmo,zcmo,vdum,wratio,volins,volrem,deltvb,vnewt2,voldt2
      real::qquo(nmax,numax) 

!      write(io_output,*) 'START the optimization of the charge configuration'

      imolty = moltyp(i)
      iunit = nunit(imolty)

! *** store the old energy, old coordinates and ewald sum
      do ibox2 = 1,nbox
         vboxo(ibox2) = vbox(ibox2)
         vinterbo(ibox2) = vinterb(ibox2)
         vintrabo(ibox2) = vintrab(ibox2)
         vextbo(ibox2) = vextb(ibox2)
         vflucqbo(ibox2) = vflucqb(ibox2)
         velectbo(ibox2) = velectb(ibox2)

         if ( mtype .eq. 3 ) then
            vtailbo(ibox2) = vtailb(ibox2)
            vvibbo(ibox2) = vvibb(ibox2)
            vtgbo(ibox2) = vtgb(ibox2)
            vbendbo(ibox2) = vbendb(ibox2)
         end if
         
      end do
      do j = 1,iunit
         rxuo(j) = rxu(i,j)
         ryuo(j) = ryu(i,j)
         rzuo(j) = rzu(i,j)
      end do
      xcmo = xcm(i)
      ycmo = ycm(i)
      zcmo = zcm(i)
! *** store the old charges
      do ip = 1,nchain
         do j = 1,nunit(moltyp(ip))
            qquo(ip,j) = qqu(ip,j)
         end do
      end do
      if (lewald) then
! *** store the reciprocal-space sum
         do ibox2 = 1,nbox
            call recip(ibox2,vdum,vdum,3)
         end do
         if ( ldielect ) then
            call dipole(ibox,2)
         end if
      end if

! *** on the new coordinates, continue to use the fluctuating charge
! *** algorithm to optimize the charge configurations, update the
! *** energy, coordinates and the ewald sum
      
      if ( mtype .eq. 3 ) then
! *** for swap move
         vbox(ibox)     = vbox(ibox) + vnewt
         vinterb(ibox)  = vinterb(ibox) + vnewinter
         vintrab(ibox)  = vintrab(ibox) + vnewintra
         vextb(ibox)    = vextb(ibox)   + vnewext
         velectb(ibox)   = velectb(ibox)  + vnewelect+vnewewald
         vtailb(ibox)   = vtailb(ibox)   + vinsta
         vvibb(ibox)    =  vvibb(ibox)   + vnewbvib
         vtgb(ibox)     = vtgb(ibox)     + vnewtg
	 vbendb(ibox)   = vbendb(ibox)   + vnewbb
         vflucqb(ibox)  = vflucqb(ibox)  + vnewflucq

         vbox(boxrem)     = vbox(boxrem)     - voldt 
         vinterb(boxrem)  = vinterb(boxrem)  - voldinter
         vtailb(boxrem)   = vtailb(boxrem)   - vremta
         vintrab(boxrem)  = vintrab(boxrem)  - voldintra
         vvibb(boxrem)    = vvibb(boxrem)    - voldbvib
         vtgb(boxrem)     = vtgb(boxrem)     - voldtg
         vextb(boxrem)    = vextb(boxrem)    - voldext
         vbendb(boxrem)   = vbendb(boxrem)   - voldbb
         velectb(boxrem)  = velectb(boxrem)  -  (voldelect+voldewald)
         vflucqb(boxrem)  = vflucqb(boxrem)  - voldflucq
      else

         vbox(ibox)     = vbox(ibox) + deltv
         vinterb(ibox)  = vinterb(ibox) + (vintern - vintero)
         vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
         vextb(ibox)    = vextb(ibox)   + (vextn   - vexto)
         velectb(ibox)   = velectb(ibox)  + (velectn - velecto)
      end if
      
      do j = 1,iunit
         rxu(i,j) = rxuion(j,2)
         ryu(i,j) = ryuion(j,2)
         rzu(i,j) = rzuion(j,2)
      end do
! *** update chain center of mass
      call ctrmas(.false.,ibox,i,mtype)

      if (lewald) then
! *** update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
         if (mtype .eq. 3) call recip(boxrem,vdum,vdum,2)
         if ( ldielect ) then
! *** update the dipole term
            call dipole(ibox,1)
         end if
      end if            

! *** begin to optimize the charge configuration

      if ( mtype .eq. 3 ) then
         do ichoiq = 1,500
            call flucq(-1,0)
         end do
!         do ichoiq = 1,0
!            call flucq(-2,0)
!         end do
         do ichoiq = 1,500
            call flucq(2,0) 
         end do
         deltv = vbox(ibox) - vboxo(ibox)
         weight = dexp(-deltv*beta)
         deltv = vboxo(boxrem) - vbox(boxrem)
         weiold = dexp(-deltv*beta)
         volins=boxlx(ibox)*boxly(ibox)*boxlz(ibox)
         volrem=boxlx(boxrem)*boxly(boxrem)*boxlz(boxrem)

         if ( lswapinter ) then
            if (lgibbs) then
!     --- Note: acceptance based on only molecules of type imolty
               wratio = ( weight / weiold ) * ( volins * dble( ncmt(boxrem,imolty)+1 ) /  ( volrem * dble( ncmt(ibox,imolty) ) ) )
            else if (lgrand) then
               if (ibox.eq.1) then
!           --- molecule added to box 1
                  wratio = (weight /  weiold ) *  volins * B(imolty) / (ncmt(ibox,imolty)) 
               else
!            --- molecule removed from box 1
                  wratio = (weight /  weiold ) *  (ncmt(boxrem,imolty)+1)/ (B(imolty)*volrem) 
               end if
            end if
         else
               wratio = weight / weiold
         end if
         if ( wratio .gt. random() ) then
            laccept = .true.
         else
            laccept = .false.
         end if

      else

         do ichoiq = 1,nchoiq(ibox)
            call flucq(0,ibox)
         end do
         vnewt2 = 0.0d0
         voldt2 = 0.0d0
         do ibox2 = 1, nbox
            vnewt2 = vnewt2 + vbox(ibox2)
            voldt2 = voldt2 + vboxo(ibox2)
         end do
         deltv = vnewt2 - voldt2

         deltvb = beta * deltv
         if ( deltvb .gt. (2.3d0*softcut) ) then
            laccept = .false.
         else if ( deltv .le. 0.0d0 ) then
            laccept = .true.
         else if ( dexp(-deltvb) .gt. random() ) then
            laccept = .true.
         else
            laccept = .false.
         end if

      end if

      if ( laccept ) then

! *** combined move can be accepted now !!!

      else
! *** restore the old energy and old coordinates and ewald sum
         do ibox2 = 1,nbox
            vbox(ibox2) = vboxo(ibox2)
            vinterb(ibox2) = vinterbo(ibox2)
            vintrab(ibox2) = vintrabo(ibox2)
            vextb(ibox2) = vextbo(ibox2)
            velectb(ibox2) = velectbo(ibox2)
            vflucqb(ibox2) = vflucqbo(ibox2)
            if ( mtype .eq. 3 ) then
               vtailb(ibox2) = vtailbo(ibox2)
               vvibb(ibox2) = vvibbo(ibox2)
               vtgb(ibox2) = vtgbo(ibox2)
               vbendb(ibox2) = vbendbo(ibox2)
            end if
         end do
         do j = 1,iunit
            rxu(i,j) = rxuo(j)
            ryu(i,j) = ryuo(j)
            rzu(i,j) = rzuo(j)
         end do
         do ip = 1,nchain
            do j = 1,nunit(moltyp(ip))
               qqu(ip,j) = qquo(ip,j)
            end do
         end do
         xcm(i) = xcmo
         ycm(i) = ycmo
         zcm(i) = zcmo
         if (lewald) then
! *** restore the reciprocal-space sum
            do ibox2 = 1,nbox
               call recip(ibox2,vdum,vdum,4)
            end do
            if ( ldielect ) then
! *** restore old dipole moment
               call dipole(ibox,3)
            end if
         end if
      end if    

!      write(io_output,*) 'END the optimization of the charge configuration'

      return         
      end





