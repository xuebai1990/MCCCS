MODULE moves_simple
  use util_random,only:random
  use util_runtime,only:err_exit
  use sim_system
  use energy_kspace,only:recip
  use energy_pairwise,only:energy
  implicit none
  private
  save
  public::traxyz,rotxyz,Atom_traxyz,volume,prvolume,init_moves_simple,output_translation_rotation_stats,output_volume_stats

  real,public::acntrax(ntmax,nbxmax),acntray(ntmax,nbxmax),acntraz(ntmax,nbxmax),acnrotx(ntmax,nbxmax),acnroty(ntmax,nbxmax),acnrotz(ntmax,nbxmax),acstrax(ntmax,nbxmax),acstray(ntmax,nbxmax),acstraz(ntmax,nbxmax),acsrotx(ntmax,nbxmax),acsroty(ntmax,nbxmax),acsrotz(ntmax,nbxmax),bntrax(ntmax,nbxmax),bntray(ntmax,nbxmax),bntraz(ntmax,nbxmax),bstrax(ntmax,nbxmax),bstray(ntmax,nbxmax),bstraz(ntmax,nbxmax),bnrotx(ntmax,nbxmax),bnroty(ntmax,nbxmax),bnrotz(ntmax,nbxmax),bsrotx(ntmax,nbxmax),bsroty(ntmax,nbxmax),bsrotz(ntmax,nbxmax),acsvol(nbxmax),acnvol(nbxmax),acshmat(nbxmax,9),acnhmat(nbxmax,9),bsvol(nbxmax),bnvol(nbxmax),bshmat(nbxmax,9),bnhmat(nbxmax,9)

contains
!    *******************************************************************
!    ** makes a translational movement in x,y,or z-direction.         **
!    ** the maximum displacement is controlled by rmtrax(yz) and the  **
!    ** number of successful trial moves is stored in bstrax(yz).     **
!    ** The attempts are stored in bntrax(yz)                         **
!    *******************************************************************
  subroutine traxyz (lx,ly,lz)
    use sim_particle,only:update_neighbor_list
    use sim_cell,only:update_linked_cell
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'neigh2.inc'
!$$$      include 'system.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'bnbsma.inc'
!$$$      include 'neigh.inc'
!$$$      include 'ipswpar.inc'
!$$$      include 'eepar.inc'
!$$$!kea
!$$$      include 'garofalini.inc'

      logical::lx,ly,lz,ovrlap
      logical::lneighij

      integer::i,ibox,flagon,iunit,j,imolty,icbu,ic,ip
      real::rx,ry,rz,dchain,ddx,ddy,ddz,vnew,vold,vintran,vintrao,deltv,deltvb,vintern,vintero ,vextn,vexto,rchain,velectn,velecto,vdum,vrecipo,vrecipn,v3n,v3o

      logical::laccept

! --------------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'start TRAXYZ',lx,ly,lz
#endif
      ovrlap = .false.
!     ***    select a chain at random ***
      rchain  = random()
      do icbu = 1,nmolty
         if ( rchain .lt. pmtrmt(icbu) ) then
            imolty = icbu
            rchain = 2.0d0
         end if
      end do

      if ((lexpee).and.(imolty.ge.nmolty1)) imolty = ee_moltyp(mstate)

      if (temtyp(imolty).eq.0) return

      if (lgrand) then
! ---    select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) bntrax(imolty,ibox) = bntrax(imolty,ibox) + 1.0d0
            if(ly) bntray(imolty,ibox) = bntray(imolty,ibox) + 1.0d0
            if(lz) bntraz(imolty,ibox) = bntraz(imolty,ibox) + 1.0d0
            return
         end if
         i = idint( dble(ncmt(1,imolty))*random() ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(io_output,*) 'screwup traxyz'


      else

         dchain = dble(temtyp(imolty))
         i = int( dchain*random() + 1 )
         i = parall(imolty,i)
         ibox = nboxi(i)

      end if

! *** store number of units of i in iunit ***

      iunit = nunit(imolty)

      do  j = 1, iunit
         rxuion(j,1) = rxu(i,j)
         ryuion(j,1) = ryu(i,j)
         rzuion(j,1) = rzu(i,j)
         qquion(j,1) = qqu(i,j)
      end do
      moltion(1) = imolty

! *** move i ***
      if (lx) then
         rx =  ( 2.0*random() - 1.0d0 ) * rmtrax(imolty,ibox)
         bntrax(imolty,ibox) = bntrax(imolty,ibox) + 1.0d0
      else
         rx=0
      end if
      if (ly) then
         ry =  ( 2.0*random() - 1.0d0 ) * rmtray(imolty,ibox)
         bntray(imolty,ibox) = bntray(imolty,ibox) + 1.0d0
      else
         ry=0
      end if
      if (lz) then
         rz =  ( 2.0*random() - 1.0d0 ) * rmtraz(imolty,ibox)
         bntraz(imolty,ibox) = bntraz(imolty,ibox) + 1.0d0
      else
         rz=0
      end if
      do j = 1, iunit
         ddx=0
         ddy=0
         ddz=0
         if ( .not. (lsolid(ibox) .and. .not. lrect(ibox)) ) then
            if (lx) then
               if ( lfold .and. lpbcx ) then
                  if ((xcm(i) + rx) .lt. 0.0d0) ddx = boxlx(ibox)
                  if ((xcm(i) + rx) .gt. boxlx(ibox)) ddx = -boxlx(ibox)
               end if
            end if
            if (ly) then
               if ( lfold .and. lpbcy ) then
                  if ((ycm(i) + ry).lt.0.0d0)  ddy=boxly(ibox)
                  if ((ycm(i) + ry).gt.boxly(ibox))  ddy=-boxly(ibox)
               end if
            end if
            if (lz) then
               if ( lfold .and. lpbcz ) then
                  if ((zcm(i) + rz).lt.0.0d0)  ddz=boxlz(ibox)
                  if ((zcm(i) + rz).gt.boxlz(ibox))  ddz=-boxlz(ibox)
               end if
            end if
         end if
         rxuion(j,2) = rxuion(j,1) + rx + ddx
         ryuion(j,2) = ryuion(j,1) + ry + ddy
         rzuion(j,2) = rzuion(j,1) + rz + ddz
         qquion(j,2) = qquion(j,1)
      end do
      moltion(2) = imolty

!     *** calculate the energy of i in the new configuration ***
      flagon = 2
      call energy(i,imolty, vnew,vintran, vintern,vextn,velectn ,vdum,flagon, ibox,1, iunit,.false.,ovrlap,.false. ,vdum,.false.,.false.,.false.)
      v3n=v3garo
      if (ovrlap) return

! *** calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto ,vdum,flagon,ibox,1, iunit,.false.,ovrlap,.false. ,vdum,.false.,.false.,.false.)
      v3o = v3garo

      if (ovrlap) then
         call err_exit('disaster ovrlap in old conf of TRAXYZ')
      end if

      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
            vrecipn    =   (1.0d0-(1.0d0-etais)*lambdais)*vrecipn
            vrecipo    =   (1.0d0-(1.0d0-etais)*lambdais)*vrecipo
         else if (lstageb) then
            vrecipn =  etais*vrecipn
            vrecipo =  etais*vrecipo
         else if (lstagec) then
            vrecipn =  (etais+(1.0d0-etais)*lambdais)*vrecipn
            vrecipo =  (etais+(1.0d0-etais)*lambdais)*vrecipo
         end if
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if
! *** check for acceptance ***

      deltv  = vnew - vold
      deltvb = beta * deltv


! *** For ANES algorithm, do the Fluctuating charge moves.

      if ( lanes ) then
         call anes(i,ibox,ibox,1,laccept,deltv,vintern,vintran,vextn, velectn,vintero,vintrao,vexto,velecto,vdum,vdum, vdum,vdum,.false.)
         if ( laccept ) then
            if (lx) bstrax(imolty,ibox) = bstrax(imolty,ibox) + 1.0d0
            if (ly) bstray(imolty,ibox) = bstray(imolty,ibox) + 1.0d0
            if (lz) bstraz(imolty,ibox) = bstraz(imolty,ibox) + 1.0d0
         end if
         return
      end if

      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
!        --- accept move
      else if ( dexp(-deltvb) .gt. random() ) then
!        --- accept move
      else
!        --- move rejected
         return
      end if

!      write(io_output,*) 'TRAXYZ accepted i',i

      vbox(ibox)     = vbox(ibox) + deltv
      vinterb(ibox)  = vinterb(ibox) + (vintern - vintero)
      vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
      vextb(ibox)    = vextb(ibox)   + (vextn   - vexto)
      velectb(ibox)   = velectb(ibox)  + (velectn - velecto)
      vipswb(ibox) = vipswb(ibox) + (vipswn-vipswo)
      vwellipswb(ibox) = vwellipswb(ibox) + (vwellipswn-vwellipswo)
      vipsw = vipswb(ibox)
      vwellipsw = vwellipswb(ibox)
      v3garob(ibox) = v3garob(ibox) + (v3n-v3o)

      do j = 1,iunit
         if (lx) rxu(i,j) = rxuion(j,2)
         if (ly) ryu(i,j) = ryuion(j,2)
         if (lz) rzu(i,j) = rzuion(j,2)
      end do

      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
! *** update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if ( ldielect ) then
! *** update the dipole term
         call dipole(ibox,1)
      end if

! *** update chain center of mass

      call ctrmas(.false.,ibox,i,1)

      if (lx) bstrax(imolty,ibox) = bstrax(imolty,ibox) + 1.0d0
      if (ly) bstray(imolty,ibox) = bstray(imolty,ibox) + 1.0d0
      if (lz) bstraz(imolty,ibox) = bstraz(imolty,ibox) + 1.0d0

      if ( licell .and. (ibox.eq.boxlink)) then
!     --- update linkcell list
         call update_linked_cell(i)
      end if

      if ( lneigh ) then
         call update_neighbor_list(i,rx,ry,rz,.false.)
      end if

      if ( lneighbor .or. lgaro) then

         do ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
!            write(io_output,*) ic,i,'j:',j
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                  ndij(ip,j) = ndij(neigh_cnt(j),j)
                  nxij(ip,j) = nxij(neigh_cnt(j),j)
                  nyij(ip,j) = nyij(neigh_cnt(j),j)
                  nzij(ip,j) = nzij(neigh_cnt(j),j)
                  neigh_cnt(j) = neigh_cnt(j)-1
                  cycle
               end if
            end do
         end do
         neigh_cnt(i) = neigh_icnt
         do ic = 1,neigh_icnt
            j = neighi(ic)
            neighbor(ic,i)=j
            ndij(ic,i) = ndiji(ic)
            nxij(ic,i) = nxiji(ic)
            nyij(ic,i) = nyiji(ic)
            nzij(ic,i) = nziji(ic)
            lneighij = .false.
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  lneighij = .true.
               end if
            end do
            if ( .not. lneighij ) then
               neigh_cnt(j) = neigh_cnt(j)+1
               neighbor(neigh_cnt(j),j) = i
               ndij(neigh_cnt(j),j) = ndiji(ic)
               nxij(neigh_cnt(j),j) = -nxiji(ic)
               nyij(neigh_cnt(j),j) = -nyiji(ic)
               nzij(neigh_cnt(j),j) = -nziji(ic)
            end if
         end do
      end if

#ifdef __DEBUG__
      write(io_output,*) 'end TRAXYZ',i
#endif

      return
    end subroutine traxyz

!    *******************************************************************
!    ** makes a rotational movement around "x" space-fixed axis.      **
!    ** the maximum displacement is controlled by rmrotx and the      **
!    ** number of successful rotation is given by bsrotx.             **
!    **                                                               **
!    ** rotxyz chooses one of the three space-fixed axes at random    **
!    ** and rotates the molecule around this axis by dgamma radians.  **
!    ** the maximum angular displacement is dgamax.                   **
!    *******************************************************************
    subroutine rotxyz (lx,ly,lz )
      use sim_particle,only:update_neighbor_list
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'neigh2.inc'
!$$$      include 'system.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'bnbsma.inc'
!$$$      include 'neigh.inc'
!$$$      include 'ipswpar.inc'
!$$$      include 'eepar.inc'

      logical::lx,ly,lz,ovrlap,lneighij
      integer::i,ibox,flagon,iunit,j,imolty,iuroty,icbu ,ic,ip
      real::rx,ry,rz,dchain,rchain,vnew,vold,vintrao,dgamma,rxorig,ryorig,rzorig,rxnew2,rynew2,rznew2,vintran,deltv,deltvb,vintern,vintero,vextn,vexto,vdum ,velectn,velecto
! *** further variable definitions
      real::cosdg, sindg, rmrot
      real::vrecipn,vrecipo


      logical::laccept

! --------------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'start ROTXYZ',lx,ly,lz
#endif
      ovrlap = .false.
      if (lgrand) then
! ---    select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         rchain  = random()
         do icbu = 1,nmolty
            if ( rchain .lt. pmromt(icbu) ) then
               imolty = icbu
               rchain = 2.0d0
            end if
         end do
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) bnrotx(imolty,ibox) = bnrotx(imolty,ibox) + 1.0d0
            if(ly) bnroty(imolty,ibox) = bnroty(imolty,ibox) + 1.0d0
            if(lz) bnrotz(imolty,ibox) = bnrotz(imolty,ibox) + 1.0d0
            return
         end if

 10      dchain = dble(temtyp(imolty))
         i = int( dchain*random() + 1 )
         i = parall(imolty,i)
         if (nboxi(i) .ne. 1) goto 10
         iuroty = iurot(imolty)
      else
! ***    select a chain type at random ***
         rchain  = random()
         do icbu = 1,nmolty
            if ( rchain .lt. pmromt(icbu) ) then
               imolty = icbu
               rchain = 2.0d0
            end if
         end do

         if ((lexpee).and.(imolty.ge.nmolty1)) imolty = ee_moltyp(mstate)

         if (temtyp(imolty).eq.0) return

         dchain = dble(temtyp(imolty))
         i = int( dchain*random() + 1 )
         i = parall(imolty,i)

         ibox = nboxi(i)
         iuroty = iurot(imolty)
      end if

!kea 6/4/09 -- for multiple rotation centers
      if(iuroty.lt.0) then
         if(nrotbd(imolty).gt.1) then
            rchain = random()
            do icbu = 1,nrotbd(imolty)
               if( rchain .lt. pmrotbd(icbu,imolty)) then
                  iuroty = irotbd(icbu,imolty)
               end if
            end do
         else
            iuroty = irotbd(1,imolty)
         end if
      end if

! *** store number of units of i in iunit ***
      iunit = nunit(imolty)

! *** store current positions in old-new array #1# ***
      do j = 1, iunit
         rxuion(j,1) = rxu(i,j)
         ryuion(j,1) = ryu(i,j)
         rzuion(j,1) = rzu(i,j)
         qquion(j,1) = qqu(i,j)
      end do
      moltion(1) = imolty

! *** choose a random angular displacement ***

      if (lx) then
         rmrot = rmrotx(imolty,ibox)
         bnrotx(imolty,ibox) = bnrotx(imolty,ibox) + 1.0d0
      else if (ly) then
         rmrot = rmroty(imolty,ibox)
         bnroty(imolty,ibox) = bnroty(imolty,ibox) + 1.0d0
      else if (lz) then
         rmrot = rmrotz(imolty,ibox)
         bnrotz(imolty,ibox) = bnrotz(imolty,ibox) + 1.0d0
      end if
      dgamma = ( 2.0d0*random() - 1.0d0 ) * rmrot

! *** set up the rotation marix ***

      cosdg = dcos( dgamma )
      sindg = dsin( dgamma )

! *** Determine the rotation coordinates ***
      if (iuroty .eq. 0) then
! *** Use the center of mass for rotation
         rxorig = xcm(i)
         ryorig = ycm(i)
         rzorig = zcm(i)
      else
! *** Use iurot for rotation
         rxorig = rxuion(iuroty,1)
         ryorig = ryuion(iuroty,1)
         rzorig = rzuion(iuroty,1)
      end if

!      write(io_output,*) 'before rotating'
!      write(io_output,*) xcm(i),ycm(i),zcm(i)
!      write(io_output,*) rxu(i,1),ryu(i,1),rzu(i,1)


      if (lx) then

! ***    ROTATE UNITS OF I AROUND X-AXIS ***

         do  j = 1, iunit
            ry = ryuion(j,1) - ryorig
            rz = rzuion(j,1) - rzorig
            rynew2 = cosdg * ry + sindg * rz
            rznew2 = cosdg * rz - sindg * ry

            rxuion(j,2) = rxuion(j,1)
            ryuion(j,2) = ryorig + rynew2
            rzuion(j,2) = rzorig + rznew2
            qquion(j,2) = qquion(j,1)
         end do
      end if
      if (ly) then

! ***    ROTATE UNITS OF I AROUND y-AXIS ***

         do  j = 1, iunit
            rx = rxuion(j,1) - rxorig
            rz = rzuion(j,1) - rzorig
            rxnew2 = cosdg * rx - sindg * rz
            rznew2 = cosdg * rz + sindg * rx

            rxuion(j,2) = rxorig + rxnew2
            ryuion(j,2) = ryuion(j,1)
            rzuion(j,2) = rzorig + rznew2
            qquion(j,2) = qquion(j,1)
         end do
      end if
      if (lz) then

! ***    ROTATE UNITS OF I AROUND z-AXIS ***

         do  j = 1, iunit
            rx = rxuion(j,1) - rxorig
            ry = ryuion(j,1) - ryorig
            rxnew2 = cosdg * rx + sindg * ry
            rynew2 = cosdg * ry - sindg * rx

            rxuion(j,2) = rxorig + rxnew2
            ryuion(j,2) = ryorig + rynew2
            rzuion(j,2) = rzuion(j,1)
            qquion(j,2) = qquion(j,1)
         end do
      end if
      moltion(2) = imolty

!  *** calculate the energy of i in the new configuration ***

      flagon = 2
      call energy(i,imolty, vnew,vintran,vintern,vextn,velectn,vdum ,flagon, ibox,1,iunit,.false.,ovrlap,.false.,vdum, .false.,.false.,.false.)
      if (ovrlap) return

! *** calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty, vold,vintrao,vintero,vexto,velecto,vdum ,flagon,ibox, 1, iunit,.false.,ovrlap,.false.,vdum, .false.,.false.,.false.)

      if (ovrlap) call err_exit('disaster- overlap for old conf in ROTXYZ')
      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
           vrecipn  =  (1.0d0-(1.0d0-etais)*lambdais)*vrecipn
           vrecipo  =  (1.0d0-(1.0d0-etais)*lambdais)*vrecipo
         else if (lstageb) then
           vrecipn  =  etais*vrecipn
           vrecipo  =  etais*vrecipo
         else if (lstagec) then
           vrecipn  =  (etais+(1.0d0-etais)*lambdais)*vrecipn
           vrecipo  =  (etais+(1.0d0-etais)*lambdais)*vrecipo
         end if
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if

! *** check for acceptance ***

      deltv  = vnew - vold
      deltvb = beta * deltv

! *** For ANES algorithm, do the Fluctuating charge moves.

      if ( lanes ) then
         call anes(i,ibox,ibox,2,laccept,deltv,vintern,vintran,vextn, velectn,vintero,vintrao,vexto,velecto,vdum,vdum,vdum, vdum,.false.)
         if (laccept) then
            if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0d0
            if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0d0
            if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0d0
         end if
         return
      end if

      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
!        --- move accepted
      else if ( dexp(-deltvb) .gt. random() ) then
!        --- move accepted
      else
!        --- move rejected
         return
      end if

!      write(io_output,*) 'ROTXYZ accepted',i
      vbox(ibox) = vbox(ibox) + deltv
      vinterb(ibox)  = vinterb(ibox) + (vintern -  vintero)
      vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
      vextb(ibox)    = vextb(ibox)   + (vextn-vexto)
      velectb(ibox)  = velectb(ibox) + (velectn-velecto)
      vipswb(ibox) = vipswb(ibox) + (vipswn-vipswo)
      vwellipswb(ibox) = vwellipswb(ibox) + (vwellipswn-vwellipswo)
      vipsw = vipswb(ibox)
      vwellipsw = vwellipswb(ibox)

      do j = 1, iunit
         rxu(i,j) = rxuion(j,2)
         ryu(i,j) = ryuion(j,2)
         rzu(i,j) = rzuion(j,2)
      end do

      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
! *** update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if ( ldielect ) then
           call dipole(ibox,1)
      end if


! *** update center of mass
      call ctrmas(.false.,ibox,i,2)

      if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0d0
      if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0d0
      if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0d0

      if ( lneigh ) then
         call update_neighbor_list(i,rxuion(1,2) - rxuion(1,1),ryuion(1,2) - ryuion(1,1),rzuion(1,2) - rzuion(1,1),.false.)
         call update_neighbor_list(i,rxuion(iunit,2) - rxuion(iunit,1),ryuion(iunit,2) - ryuion(iunit,1),rzuion(iunit,2) - rzuion(iunit,1),.false.)
      end if

      if ( lneighbor ) then
!         write(io_output,*) 'in rotxyz:',i,neigh_cnt(i)
         do 11 ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                  neigh_cnt(j) = neigh_cnt(j)-1
                  goto 11
               end if
            end do
 11      continue
         neigh_cnt(i) = neigh_icnt
         do ic = 1,neigh_icnt
            j = neighi(ic)
            neighbor(ic,i)=j
            lneighij = .false.
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  lneighij = .true.
               end if
            end do
            if ( .not. lneighij ) then
               neigh_cnt(j) = neigh_cnt(j)+1
               neighbor(neigh_cnt(j),j) = i
            end if
         end do
      end if

#ifdef __DEBUG__
      write(io_output,*) 'end ROTXYZ'
#endif

      return
    end subroutine rotxyz

!    *******************************************************************
!    ** makes a translational movement in x,y,or z-direction.         **
!    ** the maximum displacement is controlled by rAtrax(yz) and the  **
!    ** number of successful trial moves is stored in Abstrax(yz).     **
!    ** The attempts are stored in Abntrax(yz)                         **
!    *******************************************************************
    subroutine Atom_traxyz (lx,ly,lz )
      use sim_particle,only:update_neighbor_list
      use sim_cell,only:update_linked_cell
      use energy_kspace,only:recip_atom
      use energy_intramolecular,only:U_bonded
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'neigh2.inc'
!$$$      include 'system.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'bnbsma.inc'
!$$$      include 'neigh.inc'

      logical::lx,ly,lz,ovrlap

      logical::lneighij

      integer::i,ibox,flagon,iunit,j,imolty,icbu ,ic,ip
      integer::pick_unit, pick_chain
      real::rx,ry,rz,dchain,vnew,vold,vintran,vintrao,deltv,deltvb,vintern,vintero ,vextn,vexto,rchain,velectn,velecto,vdum,vrecipo,vrecipn

      real::vvibn,vbendn,vtgn,vvibo,vbendo,vtgo

      real::vewaldn, vewaldo

      logical::laccept

! --------------------------------------------------------------------

!      write(io_output,*) 'start TRAXYZ'
      ovrlap = .false.
!     ***    select a chain at random ***
      rchain  = random()
      do icbu = 1,nmolty
         if ( rchain .lt. pmtrmt(icbu) ) then
            imolty = icbu
            rchain = 2.0d0
         end if
      end do

      if (lgrand) then
! ---    select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) Abntrax = Abntrax + 1.0d0
            if(ly) Abntray = Abntray + 1.0d0
            if(lz) Abntraz = Abntraz + 1.0d0
            return
         end if
         pick_chain = idint( dble(ncmt(1,imolty))*random() ) + 1
         pick_chain = parbox(pick_chain,1,imolty)
         if ( moltyp(pick_chain) .ne. imolty )  write(io_output,*) 'screwup traxyz'


      else

         dchain = dble(temtyp(imolty))
         pick_chain = int( dchain*random() + 1 )
         pick_chain = parall(imolty,pick_chain)
         ibox = nboxi(pick_chain)

      end if

! *** store number of units of i in iunit ***

      iunit = nunit(imolty)

      pick_unit = int(dble(iunit*random()) + 1 )

!      write(io_output,*) pick_unit, imolty, pick_chain

      i = pick_chain

      do j = 1,iunit
        rxuion(j,1) = rxu(i,j)
        ryuion(j,1) = ryu(i,j)
        rzuion(j,1) = rzu(i,j)
        qquion(j,1) = qqu(i,j)
      end do

      moltion(1) = imolty

! *** move i ***
      if (lx) then
         rx =  ( 2.0*random() - 1.0d0 ) * Armtrax
         Abntrax = Abntrax + 1.0d0
      else
         rx=0
      end if
      if (ly) then
         ry =  ( 2.0*random() - 1.0d0 ) * Armtray
         Abntray = Abntray + 1.0d0
      else
         ry=0
      end if
      if (lz) then
         rz =  ( 2.0*random() - 1.0d0 ) * Armtraz
         Abntraz = Abntraz + 1.0d0
      else
         rz=0
      end if

      do j = 1,iunit
        if (j .eq. pick_unit) then
           rxuion(j,2) = rxuion(j,1) + rx
           ryuion(j,2) = ryuion(j,1) + ry
           rzuion(j,2) = rzuion(j,1) + rz
           qquion(j,2) = qquion(j,1)
        else
           rxuion(j,2) = rxuion(j,1)
           ryuion(j,2) = ryuion(j,1)
           rzuion(j,2) = rzuion(j,1)
           qquion(j,2) = qquion(j,1)
        end if
      end do

      moltion(2) = imolty

!  *** calculate the energy of i in the new configuration ***
      flagon = 2
      call energy(i,imolty,vnew,vintran,vintern,vextn,velectn,vewaldn,flagon,ibox,pick_unit, pick_unit,.true.,ovrlap,.false.,vdum,.false.,.false.,.true.)
      if (ovrlap) return
      call U_bonded(i,imolty,vvibn,vbendn,vtgn)

! *** calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto,vewaldo,flagon,ibox,pick_unit, pick_unit,.true.,ovrlap,.false.,vdum,.false.,.false.,.true.)
      if (ovrlap) call err_exit('disaster ovrlap in old conf of ATOM_TRAXYZ')
      call U_bonded(i,imolty,vvibo,vbendo,vtgo)

      if ( lewald .and. lelect(imolty) ) then
         call recip_atom(ibox,vrecipn,vrecipo,1,pick_unit)
         velectn = velectn + vrecipn + vewaldn
         velecto = velecto + vrecipo + vewaldo
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if
! *** check for acceptance ***

      deltv  = vnew - vold
      deltvb = beta * deltv

! *** For ANES algorithm, do the Fluctuating charge moves.
! *** For time being it will not work for atom disp [Neeraj]****
      if ( lanes ) then
         call anes(i,ibox,ibox,1,laccept,deltv,vintern,vintran,vextn, velectn,vintero,vintrao,vexto,velecto,vdum,vdum, vdum,vdum,.false.)
         if ( laccept ) then
            if (lx) Abstrax = Abstrax + 1.0d0
            if (ly) Abstray = Abstray + 1.0d0
            if (lz) Abstraz = Abstraz + 1.0d0
         end if
         return
      end if

      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
!        --- accept move
      else if ( dexp(-deltvb) .gt. random() ) then
!        --- accept move
      else
!        --- move rejected
         return
      end if

!      write(io_output,*) 'TRAXYZ accepted i',i
      vbox(ibox)     = vbox(ibox) + deltv
      vinterb(ibox)  = vinterb(ibox) + (vintern - vintero)
      vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
      vextb(ibox)    = vextb(ibox)   + (vextn   - vexto)
      velectb(ibox)   = velectb(ibox)  + (velectn - velecto)

      vtgb(ibox) = vtgb(ibox) + (vtgn-vtgo)
      vbendb(ibox) = vbendb(ibox) + (vbendn-vbendo)
      vvibb(ibox) = vvibb(ibox) + (vvibn-vvibo)


      if (lx) rxu(i,pick_unit) = rxuion(pick_unit,2)
      if (ly) ryu(i,pick_unit) = ryuion(pick_unit,2)
      if (lz) rzu(i,pick_unit) = rzuion(pick_unit,2)

      if (lewald .and. lelect(imolty)) then
! *** update reciprocal-space sum
         call recip_atom(ibox,vdum,vdum,2,pick_unit)
      end if

      if ( ldielect ) then
! *** update the dipole term
         call dipole(ibox,1)
      end if


! *** update chain center of mass

      call ctrmas(.false.,ibox,i,10)

      if (lx) Abstrax = Abstrax + 1.0d0
      if (ly) Abstray = Abstray + 1.0d0
      if (lz) Abstraz = Abstraz + 1.0d0

      if ( licell .and. (ibox.eq.boxlink)) then
!     --- update linkcell list
         call update_linked_cell(i)
      end if

      if ( lneigh ) then
         call update_neighbor_list(i,rx,ry,rz,.false.)
      end if

      if ( lneighbor ) then

         do 10 ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
!            write(io_output,*) ic,i,'j:',j
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                  neigh_cnt(j) = neigh_cnt(j)-1
                  goto 10
               end if
            end do
 10      continue
         neigh_cnt(i) = neigh_icnt
         do ic = 1,neigh_icnt
            j = neighi(ic)
            neighbor(ic,i)=j
            lneighij = .false.
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  lneighij = .true.
               end if
            end do
            if ( .not. lneighij ) then
               neigh_cnt(j) = neigh_cnt(j)+1
               neighbor(neigh_cnt(j),j) = i
            end if
         end do
      end if

!      write(io_output,*) 'end TRAXYZ',i

      return
    end subroutine Atom_traxyz

!    *************************************************************
!    ** makes an isotropic volume change                        **
!    ** the maximum change is controlled by rmtrax and the      **
!    ** number of successful trial moves is stored in bsvol.    **
!    *************************************************************
!
! --- perform change of the volume: random walk in ln(vol)
!
    subroutine volume
      use sim_particle,only:rebuild_neighbor_list
      use sim_cell
      use energy_pairwise,only:sumup
      use energy_kspace,only:calp,save_kvector,restore_kvector
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'system.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'bnbsma.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'cell.inc'
!$$$!kea
!$$$      include 'garofalini.inc'
!$$$      include 'neigh.inc'

!      logical::ovrlap,lvol,lx(nbxmax),ly(nbxmax),lz(nbxmax),lncubic
      logical::ovrlap,lvol,lx,ly,lz,lncubic
      logical::la,lb,lc

      integer::i,j,ibox,imolty,boxa,boxb,ipair,ipairb
      real::bxo(nbxmax),byo(nbxmax),bzo(nbxmax)
      real::volo(nbxmax),vboxo(nbxmax) ,dfac(nbxmax),voln(nbxmax),vflucqo(nbxmax),vintero(nbxmax) ,vtailo(nbxmax),vexto(nbxmax),velecto(nbxmax)
      real::vboxn(nbxmax),vintern(nbxmax), vtailn(nbxmax),vextn(nbxmax),velectn(nbxmax)
      real::rxuo(nmax,numax),ryuo(nmax,numax) ,rzuo(nmax,numax),qquo(nmax,numax)
      real::volt,expdv,df,dx,dy,dz,v,dele, vinter,vtail,vext,vdum,velect
      real::vminim(nbxmax)
      real::xcmo,ycmo,zcmo
      real::rbcut(nbxmax),rbcuta,rbcutb,rpair,rm
      real::w(3),min_boxl

      dimension xcmo(nmax),ycmo(nmax),zcmo(nmax)

      integer::ichoiq
      integer::hbox,jbox,jhmat
      real::rbox,hmato(9),hmatio(9)

!kea
      real::v3n(nbxmax),v3o(nbxmax)

! --------------------------------------------------------------------

!      write(io_output,*) 'start VOLUME'

! *** select pair of boxes to do the volume move
      if ( nvolb .gt. 1 ) then
         rpair = random()
         do ipair = 1, nvolb
            if ( rpair .lt. pmvolb(ipair) ) then
               ipairb = ipair
               exit
            end if
         end do
      else
         ipairb = 1
      end if

      boxa = box5(ipairb)
      boxb = box6(ipairb)

      bnvol(ipairb) = bnvol(ipairb) + 1.0d0

      lx = .false.
      ly = .false.
      lz = .false.

      if ( lsolid(boxa) ) then
! *** volume move independently in x, y, z directions
! * changing to only move in z direction:
         rm = random()
         if ( rm .le. pmvolx ) then
!            lx(boxa) = .true.
!            ly(boxa) = .false.
!            lz(boxa) = .false.
            lx = .true.
            ly = .false.
            lz = .false.
         else if ( rm .le. pmvoly ) then
!            lx(boxa) = .false.
!            ly(boxa) = .true.
!            lz(boxa) = .false.
            lx = .false.
            ly = .true.
            lz = .false.
         else
!            lx(boxa) = .false.
!            ly(boxa) = .false.
!            lz(boxa) = .true.
            lx = .false.
            ly = .false.
            lz = .true.
         end if
         if (.not. lrect(boxa)) then
            do i = 1,9
               hmato(i) = hmat(boxa,i)
               hmatio(i) = hmati(boxa,i)
            end do
         end if
      end if
      if ( lsolid(boxb) ) then
! *** volume move independently in x, y, z directions
! * changing to only move in z direction:
         rm = random()
         if ( rm .le. pmvolx ) then
!            lx(boxb) = .true.
!            ly(boxb) = .false.
!            lz(boxb) = .false.
            lx = .true.
            ly = .false.
            lz = .false.
         else if ( rm .le. pmvoly ) then
!            lx(boxb) = .false.
!            ly(boxb) = .true.
!            lz(boxb) = .false.
            lx = .false.
            ly = .true.
            lz = .false.
         else
!            lx(boxb) = .false.
!            ly(boxb) = .false.
!            lz(boxb) = .true.
            lx = .false.
            ly = .false.
            lz = .true.
         end if
         if ( .not. lrect(boxb) ) then
            do i = 1,9
               hmato(i) = hmat(boxb,i)
               hmatio(i) = hmati(boxb,i)
            end do
         end if
      end if

      if (lsolid(boxa) .and. .not. lrect(boxa)) then
         if (lsolid(boxb) .and. .not. lrect(boxb)) then
            call err_exit('can not perform volume move between two non-re ctangular boxes')
         end if
      end if

! --- store old box lengths, energy, configuration etc
      lncubic = .false.

      do ibox = 1, 2
         if (ibox .eq. 1) i = boxa
         if (ibox .eq. 2) i = boxb

         bxo(i) = boxlx(i)
         byo(i) = boxly(i)

         if ( lpbcz ) then
            bzo(i) = boxlz(i)
         end if

         if (lsolid(i) .and. .not. lrect(i)) then
             volo(i) = cell_vol(i)
             lncubic = .true.
         else
            if ( lpbcz ) then
               volo(i)   = bxo(i)*byo(i)*bzo(i)
            else
               volo(i)   = bxo(i)*byo(i)
            end if
         end if
         vboxo(i)    = vbox(i)
         vintero(i)  = vinterb(i)
         vtailo(i)   = vtailb(i)
         vexto(i)    = vextb(i)
         velecto(i)  = velectb(i)
         vflucqo(i) = vflucqb(i)
         v3o(i)      = v3garob(i)

! --- store neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_o(i) = neigh_cnt(i)
            do j=1,neigh_o(i)
               neighboro(j,i) = neighbor(j,i)
               ndijo(j,i) = ndij(j,i)
               nxijo(j,i) = nxij(j,i)
               nyijo(j,i) = nyij(j,i)
               nzijo(j,i) = nzij(j,i)
            end do
         end do
      end if

         if ( lewald ) then
            call save_kvector(i)
            call recip(i,vdum,vdum,3)
         end if
      end do

      do i = 1, nchain
         ibox = nboxi(i)
         if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
            imolty = moltyp(i)
               xcmo(i) = xcm(i)
               ycmo(i) = ycm(i)
               if (lpbcz) zcmo(i) = zcm(i)
               do j = 1, nunit(imolty)
                  rxuo(i,j) = rxu(i,j)
                  ryuo(i,j) = ryu(i,j)
                  if ( lpbcz ) rzuo(i,j) = rzu(i,j)
                  qquo(i,j) = qqu(i,j)
               end do
         end if

      end do


!      write(io_output,*) 'before lncubic', lx,ly,lz

! --- calculate total volume
      volt = volo(boxa) + volo(boxb)

      if ( lncubic ) then

! *** select one of the cell edge
         rbox = 3.0d0*random()
         if ( lx ) then
            if ( rbox .le. 1.0d0 ) then
               la = .true.
               lb = .false.
               lc = .false.
            else if (rbox .le. 2.0d0 ) then
               la = .false.
               lb = .true.
               lc = .false.
            else
               la = .false.
               lb = .false.
               lc = .true.
            end if
         else if ( ly ) then
            la = .false.
            if ( rbox .le. 1.5d0 ) then
               lb = .true.
               lc = .false.
            else
               lb = .false.
               lc = .true.
            end if
         else
            la = .false.
            lb = .false.
            lc = .true.
         end if

         if ( la ) then
            if ( lx ) jhmat = 1
            if ( ly ) jhmat = 2
            if ( lz ) jhmat = 3
         else if ( lb ) then
            if ( lx ) jhmat = 4
            if ( ly ) jhmat = 5
            if ( lz ) jhmat = 6
         else if ( lc ) then
            if ( lx ) jhmat = 7
            if ( ly ) jhmat = 8
            if ( lz ) jhmat = 9
         end if

         if ( lsolid(boxa) .and. .not. lsolid(boxb) ) then
            hbox = boxa
            jbox = boxb

            hmat(boxa,jhmat) = hmat(boxa,jhmat) + rmhmat(boxa,jhmat)* ( 2.0d0*random() - 1.0d0 )
            bnhmat(boxa,jhmat) = bnhmat(boxa,jhmat) + 1.0d0

         else
            hbox = boxb
            jbox = boxa

            hmat(boxb,jhmat) = hmat(boxb,jhmat) + rmhmat(boxb,jhmat)* ( 2.0d0*random() - 1.0d0 )
            bnhmat(boxb,jhmat) = bnhmat(boxb,jhmat) + 1.0d0

         end if

         call matops(hbox)

         voln(hbox) = cell_vol(hbox)

         rbcut(hbox) = rcut(hbox)

         w(1) = min_width(hbox,1)
         w(2) = min_width(hbox,2)
         w(3) = min_width(hbox,3)


         if (rbcut(hbox)/w(1) .gt. 0.5d0 .or. rbcut(hbox)/w(2) .gt. 0.5d0 .or. rbcut(hbox)/w(3) .gt. 0.5d0) then
            write(io_output,*) 'Problem with line 381 in volume.f'
            write(io_output,*) 'non-rectangular volume move rejected-', ' box width below cutoff size'
            write(io_output,*) 'w1:',w(1),'w2:',w(2),'w3:',w(3)
            hmat(hbox,jhmat) = hmato(jhmat)
            call dump
            call err_exit('')
!            goto 500
         end if

         voln(jbox) = volt-voln(hbox)
!         if ( la ) boxlx(hbox) = hmat(hbox,1)
!         if ( lb ) boxly(hbox) = dsqrt(hmat(hbox,4)*hmat(hbox,4)+
!     &                 hmat(hbox,5)*hmat(hbox,5))
!         if ( lc ) boxlz(hbox) = dsqrt(hmat(hbox,7)*hmat(hbox,7)
!     &        +hmat(hbox,8)*hmat(hbox,8) +
!     &        hmat(hbox,9)*hmat(hbox,9))

! *** determine the displacement of the COM

         df=(voln(jbox)/volo(jbox))**(1.0d0/3.0d0)
         boxlx(jbox) = boxlx(jbox)*df
         boxly(jbox) = boxly(jbox)*df
         boxlz(jbox) = boxlz(jbox)*df
         df = df - 1.0d0

         do i = 1,nchain
            imolty = moltyp(i)
            if (nboxi(i) .eq. hbox) then
               if ( lx ) then
                  dx = sxcm(i)*(hmat(hbox,1)-hmato(1))+ sycm(i)*(hmat(hbox,4)-hmato(4))+ szcm(i)*(hmat(hbox,7)-hmato(7))
                  xcm(i) = xcm(i) + dx
                  do j = 1, nunit(imolty)
                     rxu(i,j) = rxu(i,j) + dx
                  end do
               else if ( ly ) then
                  dy = sxcm(i)*(hmat(hbox,2)-hmato(2))+ sycm(i)*(hmat(hbox,5)-hmato(5))+ szcm(i)*(hmat(hbox,8)-hmato(8))
                  ycm(i) = ycm(i) + dy
                  do j = 1, nunit(imolty)
                     ryu(i,j) = ryu(i,j) + dy
                  end do
               else
                  dz = sxcm(i)*(hmat(hbox,3)-hmato(3))+ sycm(i)*(hmat(hbox,6)-hmato(6))+ szcm(i)*(hmat(hbox,9)-hmato(9))
                  zcm(i) = zcm(i) + dz
                  do j = 1, nunit(imolty)
                     rzu(i,j) = rzu(i,j) + dz
                  end do
               end if
            else if (nboxi(i) .eq. jbox) then
               dx = xcm(i) * df
               dy = ycm(i) * df
               if ( lpbcz ) dz = zcm(i) * df

               xcm(i) = xcm(i) + dx
               ycm(i) = ycm(i) + dy
               if ( lpbcz ) zcm(i) = zcm(i) + dz
               do j = 1, nunit(imolty)
                  rxu(i,j) = rxu(i,j) + dx
                  ryu(i,j) = ryu(i,j) + dy
                  if ( lpbcz ) rzu(i,j) = rzu(i,j) + dz
               end do
            end if

         end do

      else


! --- calculate new volume
         expdv = dexp(dlog(volo(boxa)/volo(boxb)) + rmvol(ipairb)*(2.0d0*random()-1.0d0))
         voln(boxa)= expdv*volt/(1+expdv)
         voln(boxb)= volt-voln(boxa)
         rbcut(boxa) = rcut(boxa)
         rbcut(boxb) = rcut(boxb)

         if ( lpbcz ) then

            if ( lsolid(boxa).and.lrect(boxa) ) then
! *** volume move independently in x, y, z directions
               dfac(boxa)=voln(boxa)/volo(boxa)
            else
               dfac(boxa)= (voln(boxa)/volo(boxa))**(1.0d0/3.0d0)
            end if


            if ( lsolid(boxb).and.lrect(boxb) ) then
! *** volume move independently in x, y, z directions
               dfac(boxb)=voln(boxb)/volo(boxb)
            else
               dfac(boxb)= (voln(boxb)/volo(boxb))**(1.0d0/3.0d0)
            end if
            vminim(boxa) = (2.0d0*rbcut(boxa))**(3.0d0)
         else
            dfac(boxa)= dsqrt(voln(boxa)/volo(boxa))
            dfac(boxb)= dsqrt(voln(boxb)/volo(boxb))
            vminim(boxb) = (2.0d0*rbcut(boxb))**(2.0d0)
         end if

         if ( lsolid(boxa).and. lrect(boxa) ) then
            if (lx) boxlx(boxa) = boxlx(boxa) * dfac(boxa)
            if (ly) boxly(boxa) = boxly(boxa) * dfac(boxa)
            if (lz) boxlz(boxa) = boxlz(boxa) * dfac(boxa)
         else

            boxlx(boxa) = boxlx(boxa) * dfac(boxa)
            boxly(boxa) = boxly(boxa) * dfac(boxa)
            if ( lpbcz ) then
               boxlz(boxa) = boxlz(boxa) * dfac(boxa)
            end if
         end if

         if ( lsolid(boxb).and.lrect(boxb) ) then
            if (lx) boxlx(boxb) = boxlx(boxb) * dfac(boxb)
            if (ly) boxly(boxb) = boxly(boxb) * dfac(boxb)
            if (lz) boxlz(boxb) = boxlz(boxb) * dfac(boxb)
         else
            boxlx(boxb) = boxlx(boxb) * dfac(boxb)
            boxly(boxb) = boxly(boxb) * dfac(boxb)
            if ( lpbcz ) then
               boxlz(boxb) = boxlz(boxb) * dfac(boxb)
            end if
         end if

         rbcuta = 2.0d0*rbcut(boxa)
         rbcutb = 2.0d0*rbcut(boxb)
         if ( boxlx(boxa) .lt. rbcuta .or.  boxly(boxa) .lt. rbcuta .or.  (lpbcz .and. boxlz(boxa) .lt. rbcuta) .or. boxlx(boxb) .lt. rbcutb .or.  boxly(boxb) .lt. rbcutb .or.  (lpbcz .and. boxlz(boxb) .lt. rbcutb) ) then

            write(io_output,*) 'Problem in line 552 of subroutine volume.f'
            write(io_output,*) 'A move was attempted that would lead to a  boxlength less than twice rcut'

            boxlx(boxa) = bxo(boxa)
            boxlx(boxb) = bxo(boxb)
            boxly(boxa) = byo(boxa)
            boxly(boxb) = byo(boxb)
            if ( lpbcz ) then
               boxlz(boxa) = bzo(boxa)
               boxlz(boxb) = bzo(boxb)
            end if
            call dump
            call err_exit('')
            return
         end if


! *** determine new positions of the molecules
! *** calculate centre of mass and its displacement

! - WARNING
         if ( .not. lfold )  call err_exit('volume move only correct with folded coordinates')

         do i = 1, nchain

            ibox = nboxi(i)

            if ( ibox .eq. boxa .or. ibox .eq. boxb ) then

               imolty = moltyp(i)
               df = dfac(ibox) - 1.0d0

!     if ( lsolid(ibox) ) then
               if ( lsolid(ibox).and.lrect(ibox) ) then
!     if ( lx(ibox) ) then
                  if ( lx ) then
                     dx = xcm(i) * df
                     xcm(i) = xcm(i) + dx
                     do j = 1, nunit(imolty)
                        rxu(i,j) = rxu(i,j) + dx
                     end do
                  end if
!     if ( ly(ibox) ) then
                  if ( ly ) then
                     dy = ycm(i) * df
                     ycm(i) = ycm(i) + dy
                     do j = 1, nunit(imolty)
                        ryu(i,j) = ryu(i,j) + dy
                     end do
                  end if
!     if ( lz(ibox) ) then
                  if ( lz ) then
                     dz = zcm(i) * df
                     zcm(i) = zcm(i) + dz
                     do j = 1, nunit(imolty)
                        rzu(i,j) = rzu(i,j) + dz
                     end do
                  end if
               else

                  dx = xcm(i) * df
                  dy = ycm(i) * df
                  if ( lpbcz ) dz = zcm(i) * df

                  xcm(i) = xcm(i) + dx
                  ycm(i) = ycm(i) + dy
                  if ( lpbcz ) zcm(i) = zcm(i) + dz
                  do j = 1, nunit(imolty)
                     rxu(i,j) = rxu(i,j) + dx
                     ryu(i,j) = ryu(i,j) + dy
                     if ( lpbcz ) rzu(i,j) = rzu(i,j) + dz
                  end do
               end if
            end if
         end do

      end if

      lvol = .true.
      if ( lchgall ) then
         if (lsolid(boxa).and.(.not.lrect(boxa))) then
             min_boxl = min(min_width(boxa,1),min_width(boxa,2), min_width(boxa,3))
         else
              min_boxl = min(boxlx(boxa),boxly(boxa),boxlz(boxa))
         end if
         calp(boxa) = kalp(boxa)/min_boxl
         if (lsolid(boxb).and.(.not.lrect(boxb))) then
             min_boxl = min(min_width(boxb,1),min_width(boxb,2), min_width(boxb,3))
         else
              min_boxl = min(boxlx(boxb),boxly(boxb),boxlz(boxb))
         end if
         calp(boxb) = kalp(boxb)/min_boxl
      end if

      do i = 1,2
         if ( i .eq. 1 ) ibox = boxa
         if ( i .eq. 2 ) ibox = boxb
         call sumup( ovrlap, v, vinter,vtail, vdum,vdum, vdum,vdum,vext,velect,vdum, ibox, lvol)
         if ( ovrlap ) goto 500
         vintern(ibox) = vinter
         vtailn(ibox)  = vtail
         vextn(ibox)   = vext
         velectn(ibox) = velect
         v3n(ibox) = v3garo
         vboxn(ibox)   = vboxo(ibox) + (vintern(ibox)-vintero(ibox)) + (vextn(ibox)-vexto(ibox)) + (velectn(ibox)-velecto(ibox)) + (v3n(ibox)-v3o(ibox))
!kea
      end do

      if ( lanes ) then
! *** for ANES algorithm, optimize the charge configuration
         do i = 1,2
            if ( i .eq. 1 ) ibox = boxa
            if ( i .eq. 2 ) ibox = boxb
! *** on the new coordinates, continue to use the fluctuating charge
! *** algorithm to optimize the charge configurations, update the
! *** energy, coordinates and the ewald sum

            vbox(ibox) = vbox(ibox) + (vboxn(ibox) - vboxo(ibox))
            vinterb(ibox)  = vinterb(ibox) +  (vintern(ibox) - vintero(ibox))
            vtailb(ibox) = vtailb(ibox) + (vtailn(ibox) - vtailo(ibox))
            vextb(ibox) = vextb(ibox) + (vextn(ibox) - vexto(ibox))
            velectb(ibox) = velectb(ibox) +  (velectn(ibox) - velecto(ibox))
            do ichoiq = 1,nchoiq(ibox)
               call flucq(0,ibox)
            end do
         end do

         dele = (vbox(boxa) - vboxo(boxa))+( vbox(boxb)- vboxo(boxb)) - ((nchbox(boxa)+1+ghost_particles(boxa)) *dlog(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+1+ghost_particles(boxb)) *dlog(voln(boxb)/volo(boxb))/beta)

      else if (lncubic) then

         dele = (vboxn(boxa)-vboxo(boxa)) + (vboxn(boxb)-vboxo(boxb)) - ((nchbox(boxa)+ghost_particles(boxa)) *dlog(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+ghost_particles(boxb)) *dlog(voln(boxb)/volo(boxb))/beta)

      else

         dele = (vboxn(boxa)-vboxo(boxa)) + (vboxn(boxb)-vboxo(boxb)) - ((nchbox(boxa)+1+ghost_particles(boxa)) *dlog(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+1+ghost_particles(boxb)) *dlog(voln(boxb)/volo(boxb))/beta)
      end if

! --- acceptance test

!      write(io_output,*) 'dele',dele
      if (random() .lt. dexp(-beta*dele) ) then
! --      accepted
         if ( lncubic ) then
            bshmat(hbox,jhmat) = bshmat(hbox,jhmat) + 1.0d0
         else
            bsvol(ipairb) = bsvol(ipairb) + 1.0d0
         end if

          if ( .not. lanes ) then
             do ibox = 1,2
                if ( ibox .eq. 1 ) i = boxa
                if ( ibox .eq. 2 ) i = boxb
                vbox(i)    = vbox(i) + (vboxn(i) - vboxo(i))
                vinterb(i) = vinterb(i) + (vintern(i) - vintero(i))
                vtailb(i)  = vtailb(i) + (vtailn(i) - vtailo(i))
                vextb(i)   = vextb(i) + (vextn(i) - vexto(i))
                velectb(i) = velectb(i) + (velectn(i) - velecto(i))
                v3garob(i)     = v3garob(i) + (v3n(i)-v3o(i))
             end do
          end if

! --      update centers of mass
          call ctrmas(.true.,boxa,0,5)
          call ctrmas(.true.,boxb,0,5)
! *** update linkcell, if applicable
          if (licell .and. (boxa .eq. boxlink .or. boxb .eq. boxlink)) then
             call build_linked_cell()
          end if
          if (lneigh) then
             call rebuild_neighbor_list(boxa)
             call rebuild_neighbor_list(boxb)
          end if
          return
      end if
! ---     rejected
! --- restore old energy, box lengths
 500  do ibox = 1, 2
         if ( ibox .eq. 1 ) i = boxa
         if ( ibox .eq. 2 ) i = boxb
         vbox(i)    = vboxo(i)
         vinterb(i)  = vintero(i)
         vtailb(i)   = vtailo(i)
         vextb(i)    = vexto(i)
         velectb(i)  = velecto(i)
         vflucqb(i) = vflucqo(i)
         v3garob(i)      = v3o(i)

         if (lsolid(i) .and. .not. lrect(i)) then
            do j = 1,9
               hmat(ibox,j) = hmato(j)
            end do
            call matops(i)
         end if

         boxlx(i)   = bxo(i)
         boxly(i)   = byo(i)
         if ( lpbcz ) boxlz(i)   = bzo(i)

         if ( lewald ) then
            call restore_kvector(i)
            call recip(i,vdum,vdum,4)
         end if
      end do

! --- restore old neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_cnt(i) = neigh_o(i)
            do j=1,neigh_cnt(i)
               neighbor(j,i) = neighboro(j,i)
               ndij(j,i) = ndijo(j,i)
               nxij(j,i) = nxijo(j,i)
               nyij(j,i) = nyijo(j,i)
               nzij(j,i) = nzijo(j,i)
            end do
         end do
      end if

      do i = 1, nchain
         ibox = nboxi(i)
         if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
            imolty = moltyp(i)
            xcm(i) = xcmo(i)
            ycm(i) = ycmo(i)
            if ( lpbcz ) zcm(i) = zcmo(i)
            do j = 1, nunit(imolty)
               rxu(i,j) = rxuo(i,j)
               ryu(i,j) = ryuo(i,j)
               if ( lpbcz ) rzu(i,j) = rzuo(i,j)
               qqu(i,j) = qquo(i,j)
            end do
         end if
      end do
!      write(io_output,*) 'end VOLUME'
      return
    end subroutine volume

!    *************************************************************
!    ** makes an isotropic volume change under const. pressure  **
!    ** the maximum change is controlled by rmtrax and the      **
!    ** number of successful trial moves is stored in bsvol.    **
!    *************************************************************
!
! --- perform change of the volume: random walk in ln(vol)
!
    subroutine prvolume
      use sim_particle,only:rebuild_neighbor_list
      use sim_cell
      use energy_pairwise,only:sumup
      use energy_kspace,only:calp,save_kvector,restore_kvector
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'system.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'bnbsma.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'neigh.inc'
!$$$      include 'cell.inc'
!$$$!kea
!$$$      include 'garofalini.inc'

      logical::ovrlap,lvol,lx,ly,lz,la,lb,lc

      integer::i,j,imolty,jhmat
      real::bxo,byo,bzo
      real::volo,vboxo,dfac, vintero,vtailo,vexto ,velecto,vflucqo
      real::voln,vboxn, vintern,vtailn,vextn ,velectn
      real::rxuo(nmax,numax),ryuo(nmax,numax) ,rzuo(nmax,numax),qquo(nmax,numax)
      real::df,dx,dy,dz,v,dele ,vdum,vinter,vtail,vext,velect
      real::xcmo(nmax),ycmo(nmax),zcmo(nmax)
      real::rbcut,rbox
      real::w(3)
      integer::boxvch,ichoiq,ibox
      integer::neigho_cnt(nmax),neigho(100,nmax)
      real::hmato(9),hmatio(9)
      real::min_boxl
! KEA
      real::v3o,v3n

! --------------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'start VOLUME of LNPTGIBBS'
#endif

!     Select a box at  random to change the volume of box


      lx = .false.
      ly = .false.
      lz = .false.


      rbox = random()
      do ibox = 1,nbox
         if (rbox .lt. pmvlmt(ibox) ) then
            boxvch=ibox
            goto 99
         end if
      end do


 99   bnvol(boxvch) = bnvol(boxvch) + 1.0d0

      if ( lsolid(boxvch) ) then
! *** volume move independently in x, y, z directions
         rbox = random()
         if ( rbox .le. pmvolx ) then
            lx = .true.
            ly = .false.
            lz = .false.
         else if ( rbox .le. pmvoly ) then
            lx = .false.
            ly = .true.
            lz = .false.
         else
            lx = .false.
            ly = .false.
            lz = .true.
         end if
         if ( .not. lrect(boxvch) ) then
            do i = 1,9
               hmato(i) = hmat(boxvch,i)
               hmatio(i) = hmati(boxvch,i)
            end do
         end if
      end if

! --- store old box lengths, energy, configuration etc

!---- in the box of volume change ----
      bxo    = boxlx(boxvch)
      byo    = boxly(boxvch)
      if ( lpbcz ) bzo = boxlz(boxvch)

      if (lsolid(boxvch) .and. .not. lrect(boxvch)) then
         volo = cell_vol(boxvch)
      else
         if ( lpbcz ) then
            volo   = bxo*byo*bzo
         else
            volo   = bxo*byo
         end if
      end if

      vboxo    = vbox(boxvch)
      vintero  = vinterb(boxvch)
      vtailo   = vtailb(boxvch)
      vexto    = vextb(boxvch)
      velecto  = velectb(boxvch)
      vflucqo = vflucqb(boxvch)
!kea
      v3o = v3garob(boxvch)
! --- store neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_o(i) = neigh_cnt(i)
            do j=1,neigh_o(i)
               neighboro(j,i) = neighbor(j,i)
               ndijo(j,i) = ndij(j,i)
               nxijo(j,i) = nxij(j,i)
               nyijo(j,i) = nyij(j,i)
               nzijo(j,i) = nzij(j,i)
            end do
         end do
      end if

! --- store old k vectors and reciprocal sum
      if ( lewald ) then
         call save_kvector(boxvch)
         call recip(boxvch,vdum,vdum,3)
      end if
      do i = 1, nchain
! ----- Check if the chain i is in the correct box
         if (nboxi(i).eq. boxvch) then

            imolty = moltyp(i)
            xcmo(i) = xcm(i)
            ycmo(i) = ycm(i)
            if (lpbcz) zcmo(i) = zcm(i)
            do j = 1, nunit(imolty)
               rxuo(i,j) = rxu(i,j)
               ryuo(i,j) = ryu(i,j)
               if ( lpbcz ) rzuo(i,j) = rzu(i,j)
               qquo(i,j) = qqu(i,j)
            end do
            if (lneighbor) then
               neigho_cnt(i) = neigh_cnt(i)
               do j = 1,neigho_cnt(i)
                  neigho(j,i)=neighbor(j,i)
               end do
            end if
         end if
      end do

      if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then
! *** select one of the cell edge
         rbox = 3.0d0*random()
         if ( lx ) then
            if ( rbox .le. 1.0d0 ) then
               la = .true.
               lb = .false.
               lc = .false.
            else if (rbox .le. 2.0d0 ) then
               la = .false.
               lb = .true.
               lc = .false.
            else
               la = .false.
               lb = .false.
               lc = .true.
            end if
         else if ( ly ) then
            la = .false.
            if ( rbox .le. 1.5d0 ) then
               lb = .true.
               lc = .false.
            else
               lb = .false.
               lc = .true.
            end if
         else
            la = .false.
            lb = .false.
            lc = .true.
         end if

         if ( la ) then
            if ( lx ) jhmat = 1
            if ( ly ) jhmat = 2
            if ( lz ) jhmat = 3
         else if ( lb ) then
            if ( lx ) jhmat = 4
            if ( ly ) jhmat = 5
            if ( lz ) jhmat = 6
         else if ( lc ) then
            if ( lx ) jhmat = 7
            if ( ly ) jhmat = 8
            if ( lz ) jhmat = 9
         end if

         hmat(boxvch,jhmat) = hmat(boxvch,jhmat) + rmhmat(boxvch,jhmat)* ( 2.0d0*random() - 1.0d0 )
         bnhmat(boxvch,jhmat) = bnhmat(boxvch,jhmat) + 1.0d0

         call matops(boxvch)

         voln = cell_vol(boxvch)

         w(1) = min_width(boxvch,1)
         w(2) = min_width(boxvch,2)
         w(3) = min_width(boxvch,3)

         rbcut = rcut(boxvch)

         if (rbcut/w(1) .gt. 0.5d0 .or. rbcut/w(2) .gt. 0.5d0 .or. rbcut/w(3) .gt. 0.5d0) then
            write(io_output,*) 'Problem with line 275 prvolume.f'
            write(io_output,*) 'non-rectangular prvolume move rejected-', ' box width below cutoff size'
            write(io_output,*) 'w1:',w(1),'w2:',w(2),'w3:',w(3)
            hmat(boxvch,jhmat) = hmato(jhmat)
            call dump
            call err_exit('')
!            goto 500
         end if


! *** determine the displacement of the COM
         do i = 1,nchain
            imolty = moltyp(i)
            if (nboxi(i) .eq. boxvch) then
               if ( lx ) then
                  dx = sxcm(i)*(hmat(boxvch,1)-hmato(1))+ sycm(i)*(hmat(boxvch,4)-hmato(4))+ szcm(i)*(hmat(boxvch,7)-hmato(7))
                  xcm(i) = xcm(i) + dx
                  do j = 1, nunit(imolty)
                     rxu(i,j) = rxu(i,j) + dx
                  end do
               else if ( ly ) then
                  dy = sxcm(i)*(hmat(boxvch,2)-hmato(2))+ sycm(i)*(hmat(boxvch,5)-hmato(5))+ szcm(i)*(hmat(boxvch,8)-hmato(8))
                  ycm(i) = ycm(i) + dy
                  do j = 1, nunit(imolty)
                     ryu(i,j) = ryu(i,j) + dy
                  end do
               else
                  dz = sxcm(i)*(hmat(boxvch,3)-hmato(3))+ sycm(i)*(hmat(boxvch,6)-hmato(6))+ szcm(i)*(hmat(boxvch,9)-hmato(9))
                  zcm(i) = zcm(i) + dz
                  do j = 1, nunit(imolty)
                     rzu(i,j) = rzu(i,j) + dz
                  end do
               end if
            end if
         end do

      else

! --- calculate new volume
         voln = volo + rmvol(boxvch) * ( 2.0d0*random() - 1.0d0 )
         rbcut = rcut(boxvch)
         if ( lpbcz ) then
            if (lsolid(boxvch).and.lrect(boxvch)) then
! *** volume move independently in x, y, z directions
               dfac=voln/volo
            else
               dfac = (voln/volo)**(1.0d0/3.0d0)
            end if
!            vminim = rbcut**(3.0d0)
         else
            dfac= dsqrt(voln/volo)
!            vminim = rbcut**(2.0d0)
         end if

!         if ( voln .lt. vminim ) then
!            write(io_output,*) 'prvolume move rejected - below cut-off size'
!            return
!         end if

         if ( lsolid(boxvch).and.lrect(boxvch) ) then
            if (lx) boxlx(boxvch) = boxlx(boxvch) * dfac
            if (ly) boxly(boxvch) = boxly(boxvch) * dfac
            if (lz) boxlz(boxvch) = boxlz(boxvch) * dfac
         else
            boxlx(boxvch) = boxlx(boxvch) * dfac
            boxly(boxvch) = boxly(boxvch) * dfac
            if ( lpbcz ) then
               boxlz(boxvch) = boxlz(boxvch) * dfac
            end if
         end if


         rbcut = 2.0d0*rbcut

         if (boxlx(boxvch) .lt. rbcut .or. boxly(boxvch) .lt. rbcut .or. (lpbcz .and. boxlz(boxvch) .lt. rbcut) ) then
            boxlx(boxvch) = bxo
            boxly(boxvch) = byo
            if ( lpbcz ) then
               boxlz(boxvch) = bzo
            end if
            write(io_output,*) 'boxvch',boxvch
            write(io_output,*) 'Problem in line 381 of subroutine prvolume.f'
            write(io_output,*) 'A move was attempted that would lead to a  boxlength less than twice rcut'
            call dump
            call err_exit('')
            return
         end if

! *** determine new positions of the molecules
! *** calculate centre of mass and its displacement

! - WARNING
         if ( .not. lfold ) call err_exit('volume move only correct with folded coordinates')

         df = dfac - 1.0d0

         do i = 1, nchain
            imolty = moltyp(i)

! ----- Check if the chain i is in the correct box
            if (nboxi(i) .eq. boxvch) then

               if (lsolid(boxvch).and.lrect(boxvch)) then
                  if ( lx ) then
                     dx = xcm(i) * df
                     xcm(i) = xcm(i) + dx
                     do j = 1, nunit(imolty)
                        rxu(i,j) = rxu(i,j) + dx
                     end do
                  end if
                  if ( ly ) then
                     dy = ycm(i) * df
                     ycm(i) = ycm(i) + dy
                     do j = 1, nunit(imolty)
                        ryu(i,j) = ryu(i,j) + dy
                     end do
                  end if
                  if ( lz ) then
                     dz = zcm(i) * df
                     zcm(i) = zcm(i) + dz
                     do j = 1, nunit(imolty)
                        rzu(i,j) = rzu(i,j) + dz
                     end do
                  end if
               else
                  dx = xcm(i) * df
                  dy = ycm(i) * df
                  if ( lpbcz ) dz = zcm(i) * df
                  xcm(i) = xcm(i) + dx
                  ycm(i) = ycm(i) + dy
                  if ( lpbcz ) zcm(i) = zcm(i) + dz
                  do j = 1, nunit(imolty)
                     rxu(i,j) = rxu(i,j) + dx
                     ryu(i,j) = ryu(i,j) + dy
                     if ( lpbcz ) rzu(i,j) = rzu(i,j) + dz
                  end do
               end if
            end if
         end do

      end if

      lvol = .true.
      if ( lchgall ) then
           if (lsolid(boxvch).and.(.not.lrect(boxvch))) then
             min_boxl = min(min_width(boxvch,1),min_width(boxvch,2), min_width(boxvch,3))
           else
              min_boxl = min(boxlx(boxvch),boxly(boxvch),boxlz(boxvch))
           end if
           calp(boxvch) = kalp(boxvch)/boxlx(boxvch)
      end if
      call sumup( ovrlap, v, vinter, vtail, vdum,vdum, vdum,vdum,vext,velect,vdum, boxvch, lvol)
      if ( ovrlap ) then
!         write(io_output,*) 'move rejected due to overlap in PRVOLUME'
         goto 500
      end if

      vintern  = vinter
      vtailn   = vtail
      vextn    = vext
      velectn  = velect
      v3n = v3garo
      vboxn    = vboxo + (vintern-vintero) + (vextn-vexto)  + (velectn-velecto) + (v3n-v3o)
!      write(io_output,*) 'new  energy',  vboxn

      if ( lanes ) then
! *** for ANES algorithm, optimize the charge configuration
! *** on the new coordinates, continue to use the fluctuating charge
! *** algorithm to optimize the charge configurations, update the
! *** energy, coordinates and the ewald sum

         vbox(boxvch)    = vbox(boxvch) + (vboxn - vboxo)
         vinterb(boxvch) = vinterb(boxvch) + (vintern-vintero)
         vtailb(boxvch)  = vtailb(boxvch) + (vtailn-vtailo)
         vextb(boxvch)   = vextb(boxvch) + (vextn-vexto)
         velectb(boxvch) = velectb(boxvch) + (velectn-velecto)
         do ichoiq = 1,nchoiq(boxvch)
            call flucq(0,boxvch)
         end do

         dele = (vbox(boxvch) - vboxo)+ express(boxvch)*(voln-volo)  - ((nchbox(boxvch)+ghost_particles(boxvch)) * dlog(voln/volo) / beta )

      else

         dele = ( vboxn - vboxo ) +  express(boxvch)*(voln-volo)  - ((nchbox(boxvch)+ghost_particles(boxvch))  * dlog(voln/volo)/beta )
      end if

! --- acceptance test

!--- check problem--
!      write(io_output,*) 'length of box',      boxlx(boxvch)
!      write(io_output,*) 'change of  vol',     voln - volo
!      write(io_output,*) 'change of  energy',  vboxn, vboxo
!      write(io_output,*) 'dele',     dele


      if (random() .lt. dexp(-(beta*dele)) ) then
! --      accepted
          if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then
             bshmat(boxvch,jhmat) = bshmat(boxvch,jhmat) + 1.0d0
          else
             bsvol(boxvch) = bsvol(boxvch) + 1.0d0
          end if
          if ( .not. lanes ) then
             vbox(boxvch)    = vbox(boxvch) + (vboxn - vboxo)
             vinterb(boxvch) = vinterb(boxvch) + (vintern-vintero)
             vtailb(boxvch)  = vtailb(boxvch) + (vtailn-vtailo)
             vextb(boxvch)   = vextb(boxvch) + (vextn-vexto)
             velectb(boxvch) = velectb(boxvch) + (velectn-velecto)
             v3garob(boxvch) = v3garob(boxvch) + (v3n-v3o)
          end if
! ---     update centers of mass
          call ctrmas(.true.,boxvch,0,5)
! *** update linkcell, if applicable
          if (licell .and. (boxvch .eq. boxlink)) then
             call build_linked_cell()
          end if
          if (lneigh) then
             call rebuild_neighbor_list(boxvch)
          end if
          return
      end if

! ---     rejected
! --- re store old box lengths, energy, configuration etc
 500  continue
      if ( lsolid(boxvch) .and. .not. lrect( boxvch)) then
         do i = 1,9
            hmat(boxvch,i) = hmato(i)
         end do
         call matops(boxvch)
      else
         boxlx(boxvch)   = bxo
         boxly(boxvch)   = byo
         if ( lpbcz ) boxlz(boxvch)   = bzo
      end if
      vbox(boxvch) = vboxo
      vinterb(boxvch) = vintero
      vtailb(boxvch) = vtailo
      vextb(boxvch) = vexto
      velectb(boxvch) = velecto
      vflucqb(boxvch) = vflucqo
      v3garob(boxvch) = v3o

! --- restore old neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_cnt(i) = neigh_o(i)
            do j=1,neigh_cnt(i)
               neighbor(j,i) = neighboro(j,i)
               ndij(j,i) = ndijo(j,i)
               nxij(j,i) = nxijo(j,i)
               nyij(j,i) = nyijo(j,i)
               nzij(j,i) = nzijo(j,i)
            end do
         end do
      end if

      if ( lewald ) then
! --- restore old k vectors and reciprocal sum and calp
         call restore_kvector(boxvch)
         call recip(boxvch,vdum,vdum,4)
      end if
      do i = 1, nchain
         imolty = moltyp(i)
! ----- Check if the chain i is in the correct box
         if (nboxi(i) .eq. boxvch) then
            xcm(i) = xcmo(i)
            ycm(i) = ycmo(i)
            if ( lpbcz ) zcm(i) = zcmo(i)
            do j = 1, nunit(imolty)
               rxu(i,j) = rxuo(i,j)
               ryu(i,j) = ryuo(i,j)
               if ( lpbcz ) rzu(i,j) = rzuo(i,j)
               qqu(i,j) = qquo(i,j)
            end do
            if (lneighbor) then
               neigh_cnt(i) = neigho_cnt(i)
               do j = 1,neigh_cnt(i)
                  neighbor(j,i)=neigho(j,i)
               end do
            end if
         end if
      end do

#ifdef __DEBUG__
      write(io_output,*) 'end VOLUME of LNPTGIBBS'
#endif

      return
  end subroutine prvolume

  subroutine init_moves_simple
    acntrax = 0.d0
    acntray = 0.d0
    acntraz = 0.d0
    acnrotx = 0.d0
    acnroty = 0.d0
    acnrotz = 0.d0
    acstrax = 0.d0
    acstray = 0.d0
    acstraz = 0.d0
    acsrotx = 0.d0
    acsroty = 0.d0
    acsrotz = 0.d0
    bstrax = 0.0d0
    bstray = 0.0d0
    bstraz = 0.0d0
    bsrotx = 0.0d0
    bsroty = 0.0d0
    bsrotz = 0.0d0
    bntrax = 0.0d0
    bntray = 0.0d0
    bntraz = 0.0d0
    bnrotx = 0.0d0
    bnroty = 0.0d0
    bnrotz = 0.0d0

    acsvol = 0.d0
    acnvol = 0.d0
    acshmat = 0.0d0
    acnhmat = 0.0d0
    bsvol = 0.0d0
    bnvol = 0.0d0
    bshmat = 0.0d0
    bnhmat = 0.0d0
  end subroutine init_moves_simple

! *** write some information about translations and rotations
  subroutine output_translation_rotation_stats(io_output)
    integer,intent(in)::io_output
    integer::ibox,i
    real::ratvol
    character(LEN=default_path_length)::fmt(3)=(/"(' x-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)","(' y-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)","(' z-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"/)

    write(io_output,*)
    write(io_output,*) '### Translations ###'
    write(io_output,*)
    do ibox = 1,nbox
       do i=1,nmolty
          write(io_output,*) 'molecule typ =',i,' in box',ibox
          acntrax(i,ibox) = acntrax(i,ibox) + bntrax(i,ibox)
          acstrax(i,ibox) = acstrax(i,ibox) + bstrax(i,ibox)
          if ( acntrax(i,ibox) .ne. 0.0d0 ) then
             ratvol = acstrax(i,ibox) / acntrax(i,ibox)
          else
             ratvol = 0.0d0
          end if
          write(io_output,fmt(1)) acntrax(i,ibox),ratvol,rmtrax(i,ibox)

          acntray(i,ibox) = acntray(i,ibox) + bntray(i,ibox)
          acstray(i,ibox) = acstray(i,ibox) + bstray(i,ibox)
          if ( acntray(i,ibox) .ne. 0.0d0 ) then
             ratvol = acstray(i,ibox) / acntray(i,ibox)
          else
             ratvol = 0.0d0
          end if
          write(io_output,fmt(2)) acntray(i,ibox),ratvol,rmtray(i,ibox)

          acntraz(i,ibox) = acntraz(i,ibox) + bntraz(i,ibox)
          acstraz(i,ibox) = acstraz(i,ibox) + bstraz(i,ibox)
          if ( acntraz(i,ibox) .ne. 0.0d0 ) then
             ratvol = acstraz(i,ibox) / acntraz(i,ibox)
          else
             ratvol = 0.0d0
          end if
          write(io_output,fmt(3)) acntraz(i,ibox),ratvol,rmtraz(i,ibox)
          write(io_output,*)

       end do
    end do

    write(io_output,*) '### Rotations ###'
    write(io_output,*)
    do ibox = 1,nbox
       do i=1,nmolty
          write(io_output,*) 'molecule typ =',i,' in box',ibox
          acnrotx(i,ibox) = acnrotx(i,ibox) + bnrotx(i,ibox)
          acsrotx(i,ibox) = acsrotx(i,ibox) + bsrotx(i,ibox)
          if ( acnrotx(i,ibox) .ne. 0.0d0 ) then
             ratvol = acsrotx(i,ibox) / acnrotx(i,ibox)
          else
             ratvol = 0.0d0
          end if
          write(io_output,fmt(1)) acnrotx(i,ibox),ratvol,rmrotx(i,ibox)

          acnroty(i,ibox) = acnroty(i,ibox) + bnroty(i,ibox)
          acsroty(i,ibox) = acsroty(i,ibox) + bsroty(i,ibox)
          if ( acnroty(i,ibox) .ne. 0.0d0 ) then
             ratvol = acsroty(i,ibox) / acnroty(i,ibox)
          else
             ratvol = 0.0d0
          end if
          write(io_output,fmt(2)) acnroty(i,ibox),ratvol,rmroty(i,ibox)

          acnrotz(i,ibox) = acnrotz(i,ibox) + bnrotz(i,ibox)
          acsrotz(i,ibox) = acsrotz(i,ibox) + bsrotz(i,ibox)
          if ( acnrotz(i,ibox) .ne. 0.0d0 ) then
             ratvol = acsrotz(i,ibox) / acnrotz(i,ibox)
          else
             ratvol = 0.0d0
          end if
          write(io_output,fmt(3)) acnrotz(i,ibox),ratvol,rmrotz(i,ibox)
          write(io_output,*)
       end do
    end do
  end subroutine output_translation_rotation_stats

  subroutine output_volume_stats(io_output)
    integer,intent(in)::io_output
    integer::ibox,j
    real::ratvol
    character(LEN=default_path_length)::fmt="(' h-matrix attempts =',f8.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"

    write(io_output,*)
    write(io_output,*) '### Volume change       ###'
    do ibox = 1,nbox
       if (lsolid(ibox) .and. .not. lrect(ibox)) then
          do j = 1,9
             acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
             acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
             if ( acshmat(ibox,j) .gt. 0.5d0) then
                write(io_output,fmt) acnhmat(ibox,j), acshmat(ibox,j)/acnhmat(ibox,j),rmhmat(ibox,j)
             else
                write(io_output,fmt) acnhmat(ibox,j),0.0d0,rmhmat(ibox,j)
             end if
          end do
       else
          acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
          acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
          if ( acnvol(ibox) .ne. 0.0d0 ) then
             ratvol = acsvol(ibox) / acnvol(ibox)
          else
             ratvol = 0.0d0
          end if
          write(io_output,"(' attempts =',f8.1,'   ratio =',f6.3, '   max.displ. =',e11.4)") acnvol(ibox),ratvol,rmvol(ibox)
       end if
    end do
  end subroutine output_volume_stats
end MODULE moves_simple
