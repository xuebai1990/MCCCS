MODULE moves_simple
  use util_random,only:random
  use util_runtime,only:err_exit
  use sim_system
  use energy_kspace,only:recip
  use energy_pairwise,only:energy
  implicit none
  private
  save
  public::traxyz,rotxyz,Atom_traxyz,init_moves_simple,output_translation_rotation_stats

  real,allocatable,public::acntrax(:,:),acntray(:,:),acntraz(:,:),acnrotx(:,:),acnroty(:,:),acnrotz(:,:),acstrax(:,:)&
   ,acstray(:,:),acstraz(:,:),acsrotx(:,:),acsroty(:,:),acsrotz(:,:),bntrax(:,:),bntray(:,:),bntraz(:,:),bstrax(:,:)&
   ,bstray(:,:),bstraz(:,:),bnrotx(:,:),bnroty(:,:),bnrotz(:,:),bsrotx(:,:),bsroty(:,:),bsrotz(:,:)

contains
!    *******************************************************************
! makes a translational movement in x,y,or z-direction.         **
! the maximum displacement is controlled by rmtrax(yz) and the  **
! number of successful trial moves is stored in bstrax(yz).     **
! The attempts are stored in bntrax(yz)                         **
!    *******************************************************************
  subroutine traxyz (lx,ly,lz)
    use sim_particle,only:update_neighbor_list
    use sim_cell,only:update_linked_cell

      logical::lx,ly,lz,ovrlap
      logical::lneighij

      integer::i,ibox,flagon,iunit,j,imolty,icbu,ic,ip
      real::rx,ry,rz,dchain,ddx,ddy,ddz,vnew,vold,vintran,vintrao,deltv,deltvb,vintern,vintero ,vextn,vexto,rchain,velectn,velecto,vdum,vrecipo,vrecipn,v3n,v3o

      logical::laccept

! --------------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'start TRAXYZ in ',myid,lx,ly,lz
#endif
      ovrlap = .false.
! select a chain at random ***
      rchain  = random(-1)
      do icbu = 1,nmolty
         if ( rchain .lt. pmtrmt(icbu) ) then
            imolty = icbu
            rchain = 2.0E0_dp
         end if
      end do

      if ((lexpee).and.(imolty.ge.nmolty1)) imolty = ee_moltyp(mstate)

      if (temtyp(imolty).eq.0) return

      if (lgrand) then
! select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) bntrax(imolty,ibox) = bntrax(imolty,ibox) + 1.0E0_dp
            if(ly) bntray(imolty,ibox) = bntray(imolty,ibox) + 1.0E0_dp
            if(lz) bntraz(imolty,ibox) = bntraz(imolty,ibox) + 1.0E0_dp
            return
         end if
         i = int( real(ncmt(1,imolty),dp)*random(-1) ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(io_output,*) 'screwup traxyz'


      else

         dchain = real(temtyp(imolty),dp)
         i = int( dchain*random(-1) + 1 )
         i = parall(imolty,i)
         ibox = nboxi(i)

      end if

! store number of units of i in iunit ***

      iunit = nunit(imolty)

      do  j = 1, iunit
         rxuion(j,1) = rxu(i,j)
         ryuion(j,1) = ryu(i,j)
         rzuion(j,1) = rzu(i,j)
         qquion(j,1) = qqu(i,j)
      end do
      moltion(1) = imolty

! move i ***
      if (lx) then
         rx =  ( 2.0*random(-1) - 1.0E0_dp ) * rmtrax(imolty,ibox)
         bntrax(imolty,ibox) = bntrax(imolty,ibox) + 1.0E0_dp
      else
         rx=0
      end if
      if (ly) then
         ry =  ( 2.0*random(-1) - 1.0E0_dp ) * rmtray(imolty,ibox)
         bntray(imolty,ibox) = bntray(imolty,ibox) + 1.0E0_dp
      else
         ry=0
      end if
      if (lz) then
         rz =  ( 2.0*random(-1) - 1.0E0_dp ) * rmtraz(imolty,ibox)
         bntraz(imolty,ibox) = bntraz(imolty,ibox) + 1.0E0_dp
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
                  if ((xcm(i) + rx) .lt. 0.0E0_dp) ddx = boxlx(ibox)
                  if ((xcm(i) + rx) .gt. boxlx(ibox)) ddx = -boxlx(ibox)
               end if
            end if
            if (ly) then
               if ( lfold .and. lpbcy ) then
                  if ((ycm(i) + ry).lt.0.0E0_dp)  ddy=boxly(ibox)
                  if ((ycm(i) + ry).gt.boxly(ibox))  ddy=-boxly(ibox)
               end if
            end if
            if (lz) then
               if ( lfold .and. lpbcz ) then
                  if ((zcm(i) + rz).lt.0.0E0_dp)  ddz=boxlz(ibox)
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

! calculate the energy of i in the new configuration ***
      flagon = 2
      call energy(i,imolty, vnew,vintran, vintern,vextn,velectn ,vdum,flagon, ibox,1, iunit,.false.,ovrlap,.false. ,vdum,.false.,.false.,.false.)
      v3n=v3garo
      if (ovrlap) return

! calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto ,vdum,flagon,ibox,1, iunit,.false.,ovrlap,.false. ,vdum,.false.,.false.,.false.)
      v3o = v3garo

      if (ovrlap) then
         call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf of TRAXYZ',myid+1)
      end if

      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
            vrecipn    =   (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipn
            vrecipo    =   (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipo
         else if (lstageb) then
            vrecipn =  etais*vrecipn
            vrecipo =  etais*vrecipo
         else if (lstagec) then
            vrecipn =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipn
            vrecipo =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipo
         end if
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if
! check for acceptance ***

      deltv  = vnew - vold
      deltvb = beta * deltv


! For ANES algorithm, do the Fluctuating charge moves.

      if ( lanes ) then
         call anes(i,ibox,ibox,1,laccept,deltv,vintern,vintran,vextn, velectn,vintero,vintrao,vexto,velecto,vdum,vdum, vdum,vdum,.false.)
         if ( laccept ) then
            if (lx) bstrax(imolty,ibox) = bstrax(imolty,ibox) + 1.0E0_dp
            if (ly) bstray(imolty,ibox) = bstray(imolty,ibox) + 1.0E0_dp
            if (lz) bstraz(imolty,ibox) = bstraz(imolty,ibox) + 1.0E0_dp
         end if
         return
      end if

      if ( deltvb .gt. (2.3E0_dp*softcut) ) return

      if ( deltv .le. 0.0E0_dp ) then
! accept move
      else if ( exp(-deltvb) .gt. random(-1) ) then
! accept move
      else
! move rejected
         return
      end if

! write(io_output,*) 'TRAXYZ accepted i',i

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
! update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if ( ldielect ) then
! update the dipole term
         call dipole(ibox,1)
      end if

! update chain center of mass

      call ctrmas(.false.,ibox,i,1)

      if (lx) bstrax(imolty,ibox) = bstrax(imolty,ibox) + 1.0E0_dp
      if (ly) bstray(imolty,ibox) = bstray(imolty,ibox) + 1.0E0_dp
      if (lz) bstraz(imolty,ibox) = bstraz(imolty,ibox) + 1.0E0_dp

      if ( licell .and. (ibox.eq.boxlink)) then
! update linkcell list
         call update_linked_cell(i)
      end if

      if ( lneigh ) then
         call update_neighbor_list(i,rx,ry,rz,.false.)
      end if

      if ( lneighbor .or. lgaro) then

         do ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
! write(io_output,*) ic,i,'j:',j
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
      write(io_output,*) 'end TRAXYZ in ',myid,i
#endif

      return
    end subroutine traxyz

!    *******************************************************************
! makes a rotational movement around "x" space-fixed axis.      **
! the maximum displacement is controlled by rmrotx and the      **
! number of successful rotation is given by bsrotx.             **
!    **                                                               **
! rotxyz chooses one of the three space-fixed axes at random    **
! and rotates the molecule around this axis by dgamma radians.  **
! the maximum angular displacement is dgamax.                   **
!    *******************************************************************
    subroutine rotxyz (lx,ly,lz )
      use sim_particle,only:update_neighbor_list

      logical::lx,ly,lz,ovrlap,lneighij
      integer::i,ibox,flagon,iunit,j,imolty,iuroty,icbu ,ic,ip
      real::rx,ry,rz,dchain,rchain,vnew,vold,vintrao,dgamma,rxorig,ryorig,rzorig,rxnew2,rynew2,rznew2,vintran,deltv,deltvb,vintern,vintero,vextn,vexto,vdum ,velectn,velecto
! further variable definitions
      real::cosdg, sindg, rmrot
      real::vrecipn,vrecipo


      logical::laccept

! --------------------------------------------------------------------

#ifdef __DEBUG__
      write(io_output,*) 'start ROTXYZ in ',myid,lx,ly,lz
#endif
      ovrlap = .false.
      if (lgrand) then
! select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         rchain  = random(-1)
         do icbu = 1,nmolty
            if ( rchain .lt. pmromt(icbu) ) then
               imolty = icbu
               rchain = 2.0E0_dp
            end if
         end do
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) bnrotx(imolty,ibox) = bnrotx(imolty,ibox) + 1.0E0_dp
            if(ly) bnroty(imolty,ibox) = bnroty(imolty,ibox) + 1.0E0_dp
            if(lz) bnrotz(imolty,ibox) = bnrotz(imolty,ibox) + 1.0E0_dp
            return
         end if

 10      dchain = real(temtyp(imolty),dp)
         i = int( dchain*random(-1) + 1 )
         i = parall(imolty,i)
         if (nboxi(i) .ne. 1) goto 10
         iuroty = iurot(imolty)
      else
! select a chain type at random ***
         rchain  = random(-1)
         do icbu = 1,nmolty
            if ( rchain .lt. pmromt(icbu) ) then
               imolty = icbu
               rchain = 2.0E0_dp
            end if
         end do

         if ((lexpee).and.(imolty.ge.nmolty1)) imolty = ee_moltyp(mstate)

         if (temtyp(imolty).eq.0) return

         dchain = real(temtyp(imolty),dp)
         i = int( dchain*random(-1) + 1 )
         i = parall(imolty,i)

         ibox = nboxi(i)
         iuroty = iurot(imolty)
      end if

!kea 6/4/09 -- for multiple rotation centers
      if(iuroty.lt.0) then
         if(nrotbd(imolty).gt.1) then
            rchain = random(-1)
            do icbu = 1,nrotbd(imolty)
               if( rchain .lt. pmrotbd(icbu,imolty)) then
                  iuroty = irotbd(icbu,imolty)
               end if
            end do
         else
            iuroty = irotbd(1,imolty)
         end if
      end if

! store number of units of i in iunit ***
      iunit = nunit(imolty)

! store current positions in old-new array #1# ***
      do j = 1, iunit
         rxuion(j,1) = rxu(i,j)
         ryuion(j,1) = ryu(i,j)
         rzuion(j,1) = rzu(i,j)
         qquion(j,1) = qqu(i,j)
      end do
      moltion(1) = imolty

! choose a random angular displacement ***

      if (lx) then
         rmrot = rmrotx(imolty,ibox)
         bnrotx(imolty,ibox) = bnrotx(imolty,ibox) + 1.0E0_dp
      else if (ly) then
         rmrot = rmroty(imolty,ibox)
         bnroty(imolty,ibox) = bnroty(imolty,ibox) + 1.0E0_dp
      else if (lz) then
         rmrot = rmrotz(imolty,ibox)
         bnrotz(imolty,ibox) = bnrotz(imolty,ibox) + 1.0E0_dp
      end if
      dgamma = ( 2.0E0_dp*random(-1) - 1.0E0_dp ) * rmrot

! set up the rotation marix ***

      cosdg = cos( dgamma )
      sindg = sin( dgamma )

! Determine the rotation coordinates ***
      if (iuroty .eq. 0) then
! Use the center of mass for rotation
         rxorig = xcm(i)
         ryorig = ycm(i)
         rzorig = zcm(i)
      else
! Use iurot for rotation
         rxorig = rxuion(iuroty,1)
         ryorig = ryuion(iuroty,1)
         rzorig = rzuion(iuroty,1)
      end if

! write(io_output,*) 'before rotating'
! write(io_output,*) xcm(i),ycm(i),zcm(i)
! write(io_output,*) rxu(i,1),ryu(i,1),rzu(i,1)


      if (lx) then

! ROTATE UNITS OF I AROUND X-AXIS ***

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

! ROTATE UNITS OF I AROUND y-AXIS ***

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

! ROTATE UNITS OF I AROUND z-AXIS ***

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

! calculate the energy of i in the new configuration ***

      flagon = 2
      call energy(i,imolty, vnew,vintran,vintern,vextn,velectn,vdum ,flagon, ibox,1,iunit,.false.,ovrlap,.false.,vdum, .false.,.false.,.false.)
      if (ovrlap) return

! calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty, vold,vintrao,vintero,vexto,velecto,vdum ,flagon,ibox, 1, iunit,.false.,ovrlap,.false.,vdum, .false.,.false.,.false.)

      if (ovrlap) call err_exit(__FILE__,__LINE__,'disaster- overlap for old conf in ROTXYZ',myid+1)
      if (lewald.and.lelect(imolty).and..not.lideal(ibox)) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
           vrecipn  =  (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipn
           vrecipo  =  (1.0E0_dp-(1.0E0_dp-etais)*lambdais)*vrecipo
         else if (lstageb) then
           vrecipn  =  etais*vrecipn
           vrecipo  =  etais*vrecipo
         else if (lstagec) then
           vrecipn  =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipn
           vrecipo  =  (etais+(1.0E0_dp-etais)*lambdais)*vrecipo
         end if
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if

! check for acceptance ***

      deltv  = vnew - vold
      deltvb = beta * deltv

! For ANES algorithm, do the Fluctuating charge moves.

      if ( lanes ) then
         call anes(i,ibox,ibox,2,laccept,deltv,vintern,vintran,vextn, velectn,vintero,vintrao,vexto,velecto,vdum,vdum,vdum, vdum,.false.)
         if (laccept) then
            if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0E0_dp
            if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0E0_dp
            if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0E0_dp
         end if
         return
      end if

      if ( deltvb .gt. (2.3E0_dp*softcut) ) return

      if ( deltv .le. 0.0E0_dp ) then
! move accepted
      else if ( exp(-deltvb) .gt. random(-1) ) then
! move accepted
      else
! move rejected
         return
      end if

! write(io_output,*) 'ROTXYZ accepted',i
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
! update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if ( ldielect ) then
           call dipole(ibox,1)
      end if


! update center of mass
      call ctrmas(.false.,ibox,i,2)

      if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0E0_dp
      if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0E0_dp
      if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0E0_dp

      if ( lneigh ) then
         call update_neighbor_list(i,rxuion(1,2) - rxuion(1,1),ryuion(1,2) - ryuion(1,1),rzuion(1,2) - rzuion(1,1),.false.)
         call update_neighbor_list(i,rxuion(iunit,2) - rxuion(iunit,1),ryuion(iunit,2) - ryuion(iunit,1),rzuion(iunit,2) - rzuion(iunit,1),.false.)
      end if

      if ( lneighbor ) then
! write(io_output,*) 'in rotxyz:',i,neigh_cnt(i)
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
      write(io_output,*) 'end ROTXYZ in ',myid,i
#endif

      return
    end subroutine rotxyz

!    *******************************************************************
! makes a translational movement in x,y,or z-direction.         **
! the maximum displacement is controlled by rAtrax(yz) and the  **
! number of successful trial moves is stored in Abstrax(yz).     **
! The attempts are stored in Abntrax(yz)                         **
!    *******************************************************************
    subroutine Atom_traxyz (lx,ly,lz )
      use sim_particle,only:update_neighbor_list
      use sim_cell,only:update_linked_cell
      use energy_kspace,only:recip_atom
      use energy_intramolecular,only:U_bonded

      logical::lx,ly,lz,ovrlap

      logical::lneighij

      integer::i,ibox,flagon,iunit,j,imolty,icbu ,ic,ip
      integer::pick_unit, pick_chain
      real::rx,ry,rz,dchain,vnew,vold,vintran,vintrao,deltv,deltvb,vintern,vintero ,vextn,vexto,rchain,velectn,velecto,vdum,vrecipo,vrecipn

      real::vvibn,vbendn,vtgn,vvibo,vbendo,vtgo

      real::vewaldn, vewaldo

      logical::laccept

! --------------------------------------------------------------------
#ifdef __DEBUG__
      write(io_output,*) 'start ATOM_TRAXYZ in ',myid
#endif
      ovrlap = .false.
! select a chain at random ***
      rchain  = random(-1)
      do icbu = 1,nmolty
         if ( rchain .lt. pmtrmt(icbu) ) then
            imolty = icbu
            rchain = 2.0E0_dp
         end if
      end do

      if (lgrand) then
! select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) Abntrax = Abntrax + 1.0E0_dp
            if(ly) Abntray = Abntray + 1.0E0_dp
            if(lz) Abntraz = Abntraz + 1.0E0_dp
            return
         end if
         pick_chain = int( real(ncmt(1,imolty),dp)*random(-1) ) + 1
         pick_chain = parbox(pick_chain,1,imolty)
         if ( moltyp(pick_chain) .ne. imolty )  write(io_output,*) 'screwup traxyz'


      else

         dchain = real(temtyp(imolty),dp)
         pick_chain = int( dchain*random(-1) + 1 )
         pick_chain = parall(imolty,pick_chain)
         ibox = nboxi(pick_chain)

      end if

! store number of units of i in iunit ***

      iunit = nunit(imolty)

      pick_unit = int(real(iunit*random(-1),dp) + 1 )

! write(io_output,*) pick_unit, imolty, pick_chain

      i = pick_chain

      do j = 1,iunit
        rxuion(j,1) = rxu(i,j)
        ryuion(j,1) = ryu(i,j)
        rzuion(j,1) = rzu(i,j)
        qquion(j,1) = qqu(i,j)
      end do

      moltion(1) = imolty

! move i ***
      if (lx) then
         rx =  ( 2.0*random(-1) - 1.0E0_dp ) * Armtrax
         Abntrax = Abntrax + 1.0E0_dp
      else
         rx=0
      end if
      if (ly) then
         ry =  ( 2.0*random(-1) - 1.0E0_dp ) * Armtray
         Abntray = Abntray + 1.0E0_dp
      else
         ry=0
      end if
      if (lz) then
         rz =  ( 2.0*random(-1) - 1.0E0_dp ) * Armtraz
         Abntraz = Abntraz + 1.0E0_dp
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

! calculate the energy of i in the new configuration ***
      flagon = 2
      call energy(i,imolty,vnew,vintran,vintern,vextn,velectn,vewaldn,flagon,ibox,pick_unit, pick_unit,.true.,ovrlap,.false.,vdum,.false.,.false.,.true.)
      if (ovrlap) return
      call U_bonded(i,imolty,vvibn,vbendn,vtgn)

! calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto,vewaldo,flagon,ibox,pick_unit, pick_unit,.true.,ovrlap,.false.,vdum,.false.,.false.,.true.)
      if (ovrlap) call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf of ATOM_TRAXYZ',myid+1)
      call U_bonded(i,imolty,vvibo,vbendo,vtgo)

      if ( lewald .and. lelect(imolty) ) then
         call recip_atom(ibox,vrecipn,vrecipo,1,pick_unit)
         velectn = velectn + vrecipn + vewaldn
         velecto = velecto + vrecipo + vewaldo
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if
! check for acceptance ***

      deltv  = vnew - vold
      deltvb = beta * deltv

! For ANES algorithm, do the Fluctuating charge moves.
! For time being it will not work for atom disp [Neeraj]****
      if ( lanes ) then
         call anes(i,ibox,ibox,1,laccept,deltv,vintern,vintran,vextn, velectn,vintero,vintrao,vexto,velecto,vdum,vdum, vdum,vdum,.false.)
         if ( laccept ) then
            if (lx) Abstrax = Abstrax + 1.0E0_dp
            if (ly) Abstray = Abstray + 1.0E0_dp
            if (lz) Abstraz = Abstraz + 1.0E0_dp
         end if
         return
      end if

      if ( deltvb .gt. (2.3E0_dp*softcut) ) return

      if ( deltv .le. 0.0E0_dp ) then
! accept move
      else if ( exp(-deltvb) .gt. random(-1) ) then
! accept move
      else
! move rejected
         return
      end if

! write(io_output,*) 'TRAXYZ accepted i',i
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
! update reciprocal-space sum
         call recip_atom(ibox,vdum,vdum,2,pick_unit)
      end if

      if ( ldielect ) then
! update the dipole term
         call dipole(ibox,1)
      end if


! update chain center of mass

      call ctrmas(.false.,ibox,i,10)

      if (lx) Abstrax = Abstrax + 1.0E0_dp
      if (ly) Abstray = Abstray + 1.0E0_dp
      if (lz) Abstraz = Abstraz + 1.0E0_dp

      if ( licell .and. (ibox.eq.boxlink)) then
! update linkcell list
         call update_linked_cell(i)
      end if

      if ( lneigh ) then
         call update_neighbor_list(i,rx,ry,rz,.false.)
      end if

      if ( lneighbor ) then

         do 10 ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
! write(io_output,*) ic,i,'j:',j
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

#ifdef __DEBUG__
      write(io_output,*) 'end ATOM_TRAXYZ in ',myid,i
#endif
      return
    end subroutine Atom_traxyz

  subroutine init_moves_simple
    integer::jerr
    allocate(acntrax(ntmax,nbxmax),acntray(ntmax,nbxmax),acntraz(ntmax,nbxmax),acnrotx(ntmax,nbxmax),acnroty(ntmax,nbxmax)&
     ,acnrotz(ntmax,nbxmax),acstrax(ntmax,nbxmax),acstray(ntmax,nbxmax),acstraz(ntmax,nbxmax),acsrotx(ntmax,nbxmax)&
     ,acsroty(ntmax,nbxmax),acsrotz(ntmax,nbxmax),bntrax(ntmax,nbxmax),bntray(ntmax,nbxmax),bntraz(ntmax,nbxmax)&
     ,bstrax(ntmax,nbxmax),bstray(ntmax,nbxmax),bstraz(ntmax,nbxmax),bnrotx(ntmax,nbxmax),bnroty(ntmax,nbxmax)&
     ,bnrotz(ntmax,nbxmax),bsrotx(ntmax,nbxmax),bsroty(ntmax,nbxmax),bsrotz(ntmax,nbxmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'init_moves_simple: allocation failed',jerr)
    end if

    acntrax = 0.E0_dp
    acntray = 0.E0_dp
    acntraz = 0.E0_dp
    acnrotx = 0.E0_dp
    acnroty = 0.E0_dp
    acnrotz = 0.E0_dp
    acstrax = 0.E0_dp
    acstray = 0.E0_dp
    acstraz = 0.E0_dp
    acsrotx = 0.E0_dp
    acsroty = 0.E0_dp
    acsrotz = 0.E0_dp
    bstrax = 0.0E0_dp
    bstray = 0.0E0_dp
    bstraz = 0.0E0_dp
    bsrotx = 0.0E0_dp
    bsroty = 0.0E0_dp
    bsrotz = 0.0E0_dp
    bntrax = 0.0E0_dp
    bntray = 0.0E0_dp
    bntraz = 0.0E0_dp
    bnrotx = 0.0E0_dp
    bnroty = 0.0E0_dp
    bnrotz = 0.0E0_dp
  end subroutine init_moves_simple

! write some information about translations and rotations
  subroutine output_translation_rotation_stats(io_output)
    integer,intent(in)::io_output
    integer::ibox,i
    real::ratvol
    character(LEN=default_path_length)::fmt(3)=(/"(' x-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"&
     ,"(' y-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"&
     ,"(' z-dir: attempts =',F10.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"/)

    write(io_output,*)
    write(io_output,*) '### Translations ###'
    write(io_output,*)
    do ibox = 1,nbox
       do i=1,nmolty
          write(io_output,*) 'molecule typ =',i,' in box',ibox
          acntrax(i,ibox) = acntrax(i,ibox) + bntrax(i,ibox)
          acstrax(i,ibox) = acstrax(i,ibox) + bstrax(i,ibox)
          if ( acntrax(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acstrax(i,ibox) / acntrax(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(1)) acntrax(i,ibox),ratvol,rmtrax(i,ibox)

          acntray(i,ibox) = acntray(i,ibox) + bntray(i,ibox)
          acstray(i,ibox) = acstray(i,ibox) + bstray(i,ibox)
          if ( acntray(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acstray(i,ibox) / acntray(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(2)) acntray(i,ibox),ratvol,rmtray(i,ibox)

          acntraz(i,ibox) = acntraz(i,ibox) + bntraz(i,ibox)
          acstraz(i,ibox) = acstraz(i,ibox) + bstraz(i,ibox)
          if ( acntraz(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acstraz(i,ibox) / acntraz(i,ibox)
          else
             ratvol = 0.0E0_dp
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
          if ( acnrotx(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acsrotx(i,ibox) / acnrotx(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(1)) acnrotx(i,ibox),ratvol,rmrotx(i,ibox)

          acnroty(i,ibox) = acnroty(i,ibox) + bnroty(i,ibox)
          acsroty(i,ibox) = acsroty(i,ibox) + bsroty(i,ibox)
          if ( acnroty(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acsroty(i,ibox) / acnroty(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(2)) acnroty(i,ibox),ratvol,rmroty(i,ibox)

          acnrotz(i,ibox) = acnrotz(i,ibox) + bnrotz(i,ibox)
          acsrotz(i,ibox) = acsrotz(i,ibox) + bsrotz(i,ibox)
          if ( acnrotz(i,ibox) .ne. 0.0E0_dp ) then
             ratvol = acsrotz(i,ibox) / acnrotz(i,ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,fmt(3)) acnrotz(i,ibox),ratvol,rmrotz(i,ibox)
          write(io_output,*)
       end do
    end do
  end subroutine output_translation_rotation_stats
end MODULE moves_simple
