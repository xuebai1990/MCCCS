      subroutine rotxyz (lx,ly,lz )

!    *******************************************************************
!    ** makes a rotational movement around "x" space-fixed axis.      **
!    ** the maximum displacement is controlled by rmrotx and the      **
!    ** number of successful rotation is given by bsrotx.             **
!    **                                                               **
!    ** rotxyz chooses one of the three space-fixed axes at random    **
!    ** and rotates the molecule around this axis by dgamma radians.  **
!    ** the maximum angular displacement is dgamax.                   **
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
      integer(KIND=normal_int)::i,ibox,flagon,iunit,j,imolty,iuroty,icbu
     & ,ic,ip
      real(KIND=double_precision)::rx,ry,rz,dchain,rchain,random,vnew
     & ,vold,vintrao,dgamma,rxorig,ryorig,rzorig,rxnew2,rynew2,rznew2
     & ,vintran,disvsq,deltv,deltvb,vintern,vintero,vextn,vexto,vdum
     & ,velectn,velecto
! *** further variable definitions
      real(KIND=double_precision)::cosdg, sindg, rmrot
      real(KIND=double_precision)::vrecipn,vrecipo


      logical::laccept

! --------------------------------------------------------------------

!      write(iou,*) 'start ROTXYZ'
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

         if ((lexpee).and.(imolty.ge.nmolty1))
     &      imolty = ee_moltyp(mstate)
                                                                                
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

!      write(iou,*) 'before rotating'
!      write(iou,*) xcm(i),ycm(i),zcm(i)
!      write(iou,*) rxu(i,1),ryu(i,1),rzu(i,1)


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
      call energy(i,imolty, vnew,vintran,vintern,vextn,velectn,vdum
     &     ,flagon, ibox,1,iunit,.false.,ovrlap,.false.,vdum,
     &     .false.,.false.)
      if (ovrlap) return
 
! *** calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty, vold,vintrao,vintero,vexto,velecto,vdum
     &     ,flagon,ibox, 1, iunit,.false.,ovrlap,.false.,vdum,
     &     .false.,.false.)

      if (ovrlap) call cleanup('disaster- overlap for
     &                      old conf in ROTXYZ')
      if ( lewald .and. lelect(imolty) ) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
           vrecipn  =  (1.0d0-(1.0d0-etais)*lambdais)*vrecipn
           vrecipo  =  (1.0d0-(1.0d0-etais)*lambdais)*vrecipo
         elseif (lstageb) then
           vrecipn  =  etais*vrecipn
           vrecipo  =  etais*vrecipo
         elseif (lstagec) then
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
         call anes(i,ibox,ibox,2,laccept,deltv,vintern,vintran,vextn,
     &        velectn,vintero,vintrao,vexto,velecto,vdum,vdum,vdum,
     &        vdum,.false.)
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
      elseif ( dexp(-deltvb) .gt. random() ) then
!        --- move accepted
      else
!        --- move rejected
         return
      end if

!      write(iou,*) 'ROTXYZ accepted',i
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

      if (lewald .and. lelect(imolty)) then
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
! *** check for update of near neighbour bitmap ***
! *** check for headgroup ***
         disvec(1,i,1) = disvec(1,i,1) + rxuion(1,2) - rxuion(1,1)
         disvec(1,i,2) = disvec(1,i,2) + ryuion(1,2) - ryuion(1,1)
         disvec(1,i,3) = disvec(1,i,3) + rzuion(1,2) - rzuion(1,1)
         disvsq = disvec(1,i,1) * disvec(1,i,1) +
     &        disvec(1,i,2) * disvec(1,i,2) +
     &        disvec(1,i,3) * disvec(1,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
! *** check for last unit ***
         disvec(2,i,1) = disvec(2,i,1)+rxuion(iunit,2)-rxuion(iunit,1)
         disvec(2,i,2) = disvec(2,i,2)+ryuion(iunit,2)-ryuion(iunit,1)
         disvec(2,i,3) = disvec(2,i,3)+rzuion(iunit,2)-rzuion(iunit,1)
         disvsq = disvec(2,i,1) * disvec(2,i,1) +
     &        disvec(2,i,2) * disvec(2,i,2) +
     &        disvec(2,i,3) * disvec(2,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
      end if

      if ( lneighbor ) then
!         write(iou,*) 'in rotxyz:',i,neigh_cnt(i)
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
 
!      write(iou,*) 'end ROTXYZ'

      return
      end

