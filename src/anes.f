!    *********************************************************************
! optimize the electronic configuration for trans, rot, and swap  **
!    ** (config, swatch in the future) moves and accept/reject the      **
! combined move.                                                  **
! written on June 25/99 by Bin Chen.                              **
!    *********************************************************************
subroutine anes(i,ibox,boxrem,mtype,laccept,deltv,vn,vo,vinsta,vremta,vnewflucq,voldflucq,lswapinter)
  use util_random,only:random
  use sim_system
  use energy_kspace,only:recip
  implicit none

  logical::laccept,lswapinter
  integer::i,ibox,boxrem,mtype,imolty,iunit,ichoiq,ip,ibox2,j
  real::deltv,vn(nEnergy),vo(nEnergy),vinsta,vremta,vnewflucq,voldflucq,vboxo(nbxmax)&
   ,vinterbo(nbxmax),vintrabo(nbxmax),vextbo(nbxmax),velectbo(nbxmax),vflucqbo(nbxmax),vtailbo(nbxmax),vvibbo(nbxmax)&
   ,vtgbo(nbxmax),vbendbo(nbxmax),rxuo(numax),ryuo(numax),rzuo(numax),xcmo,ycmo,zcmo,vdum,wratio,volins,volrem,deltvb,vnewt2,voldt2
  real::qquo(nmax,numax)

#ifdef __DEBUG__
      write(io_output,*) 'START the optimization of the charge configuration in ',myid
#endif

      imolty = moltyp(i)
      iunit = nunit(imolty)

! store the old energy, old coordinates and ewald sum
      do ibox2 = 1,nbox
         vboxo(ibox2) = vbox(1,ibox2)
         vinterbo(ibox2) = vbox(2,ibox2)
         vintrabo(ibox2) = vbox(4,ibox2)
         vextbo(ibox2) = vbox(9,ibox2)
         vflucqbo(ibox2) = vbox(11,ibox2)
         velectbo(ibox2) = vbox(8,ibox2)

         if ( mtype .eq. 3 ) then
            vtailbo(ibox2) = vbox(3,ibox2)
            vvibbo(ibox2) = vbox(5,ibox2)
            vtgbo(ibox2) = vbox(7,ibox2)
            vbendbo(ibox2) = vbox(6,ibox2)
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
! store the old charges
      do ip = 1,nchain
         do j = 1,nunit(moltyp(ip))
            qquo(ip,j) = qqu(ip,j)
         end do
      end do
      if (lewald) then
! store the reciprocal-space sum
         do ibox2 = 1,nbox
            call recip(ibox2,vdum,vdum,3)
         end do
         if ( ldielect ) then
            call dipole(ibox,2)
         end if
      end if

! on the new coordinates, continue to use the fluctuating charge
! algorithm to optimize the charge configurations, update the
! energy, coordinates and the ewald sum

      if ( mtype .eq. 3 ) then
! for swap move
         vbox(1,ibox)     = vbox(1,ibox) + vnew(1)
         vbox(2,ibox)  = vbox(2,ibox) + vnew(2)
         vbox(4,ibox)  = vbox(4,ibox) + vnew(4)
         vbox(9,ibox)    = vbox(9,ibox)   + vnew(9)
         vbox(8,ibox)   = vbox(8,ibox)  + vnew(8)+vnew(14)
         vbox(3,ibox)   = vbox(3,ibox)   + vinsta
         vbox(5,ibox)    =  vbox(5,ibox)   + vnew(5)
         vbox(7,ibox)     = vbox(7,ibox)     + vnew(7)
	 vbox(6,ibox)   = vbox(6,ibox)   + vnew(6)
         vbox(11,ibox)  = vbox(11,ibox)  + vnewflucq

         vbox(1,boxrem)     = vbox(1,boxrem)     - vold(1)
         vbox(2,boxrem)  = vbox(2,boxrem)  - vold(2)
         vbox(3,boxrem)   = vbox(3,boxrem)   - vremta
         vbox(4,boxrem)  = vbox(4,boxrem)  - vold(4)
         vbox(5,boxrem)    = vbox(5,boxrem)    - vold(5)
         vbox(7,boxrem)     = vbox(7,boxrem)     - vold(7)
         vbox(9,boxrem)    = vbox(9,boxrem)    - vold(9)
         vbox(6,boxrem)   = vbox(6,boxrem)   - vold(6)
         vbox(8,boxrem)  = vbox(8,boxrem)  -  (vold(8)+vold(14))
         vbox(11,boxrem)  = vbox(11,boxrem)  - voldflucq
      else

         vbox(1,ibox)     = vbox(1,ibox) + deltv
         vbox(2,ibox)  = vbox(2,ibox) + (vn(2) - vo(2))
         vbox(4,ibox)  = vbox(4,ibox) + (vn(4) - vo(4))
         vbox(9,ibox)    = vbox(9,ibox)   + (vn(9)   - vo(9))
         vbox(8,ibox)   = vbox(8,ibox)  + (vn(8) - vo(8))
      end if

      do j = 1,iunit
         rxu(i,j) = rxuion(j,2)
         ryu(i,j) = ryuion(j,2)
         rzu(i,j) = rzuion(j,2)
      end do
! update chain center of mass
      call ctrmas(.false.,ibox,i,mtype)

      if (lewald) then
! update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
         if (mtype .eq. 3) call recip(boxrem,vdum,vdum,2)
         if ( ldielect ) then
! update the dipole term
            call dipole(ibox,1)
         end if
      end if

! begin to optimize the charge configuration

      if ( mtype .eq. 3 ) then
         do ichoiq = 1,500
            call flucq(-1,0)
         end do
! do ichoiq = 1,0
! call flucq(-2,0)
! end do
         do ichoiq = 1,500
            call flucq(2,0)
         end do
         deltv = vbox(1,ibox) - vboxo(ibox)
         weight = exp(-deltv*beta)
         deltv = vboxo(boxrem) - vbox(1,boxrem)
         weiold = exp(-deltv*beta)
         volins=boxlx(ibox)*boxly(ibox)*boxlz(ibox)
         volrem=boxlx(boxrem)*boxly(boxrem)*boxlz(boxrem)

         if ( lswapinter ) then
            if (lgibbs) then
! Note: acceptance based on only molecules of type imolty
               wratio = ( weight / weiold ) * ( volins * dble( ncmt(boxrem,imolty)+1 ) /  ( volrem * dble( ncmt(ibox,imolty) ) ) )
            else if (lgrand) then
               if (ibox.eq.1) then
! molecule added to box 1
                  wratio = (weight /  weiold ) *  volins * B(imolty) / (ncmt(ibox,imolty))
               else
! molecule removed from box 1
                  wratio = (weight /  weiold ) *  (ncmt(boxrem,imolty)+1)/ (B(imolty)*volrem)
               end if
            end if
         else
               wratio = weight / weiold
         end if
         if ( wratio .gt. random(-1) ) then
            laccept = .true.
         else
            laccept = .false.
         end if

      else

         do ichoiq = 1,nchoiq(ibox)
            call flucq(0,ibox)
         end do
         vnewt2 = 0.0E0_dp
         voldt2 = 0.0E0_dp
         do ibox2 = 1, nbox
            vnewt2 = vnewt2 + vbox(1,ibox2)
            voldt2 = voldt2 + vboxo(ibox2)
         end do
         deltv = vnewt2 - voldt2

         deltvb = beta * deltv
         if ( deltvb .gt. (2.3E0_dp*softcut) ) then
            laccept = .false.
         else if ( deltv .le. 0.0E0_dp ) then
            laccept = .true.
         else if ( exp(-deltvb) .gt. random(-1) ) then
            laccept = .true.
         else
            laccept = .false.
         end if

      end if

      if ( laccept ) then

! combined move can be accepted now !!!

      else
! restore the old energy and old coordinates and ewald sum
         do ibox2 = 1,nbox
            vbox(1,ibox2) = vboxo(ibox2)
            vbox(2,ibox2) = vinterbo(ibox2)
            vbox(4,ibox2) = vintrabo(ibox2)
            vbox(9,ibox2) = vextbo(ibox2)
            vbox(8,ibox2) = velectbo(ibox2)
            vbox(11,ibox2) = vflucqbo(ibox2)
            if ( mtype .eq. 3 ) then
               vbox(3,ibox2) = vtailbo(ibox2)
               vbox(5,ibox2) = vvibbo(ibox2)
               vbox(7,ibox2) = vtgbo(ibox2)
               vbox(6,ibox2) = vbendbo(ibox2)
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
! restore the reciprocal-space sum
            do ibox2 = 1,nbox
               call recip(ibox2,vdum,vdum,4)
            end do
            if ( ldielect ) then
! restore old dipole moment
               call dipole(ibox,3)
            end if
         end if
      end if

#ifdef __DEBUG__
      write(io_output,*) 'END the optimization of the charge configuration in ',myid
#endif

      return
    end subroutine anes





