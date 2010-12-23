      subroutine traxyz (lx,ly,lz )

!    *******************************************************************
!    ** makes a translational movement in x,y,or z-direction.         **
!    ** the maximum displacement is controlled by rmtrax(yz) and the  **
!    ** number of successful trial moves is stored in bstrax(yz).     **
!    ** The attempts are stored in bntrax(yz)                         **
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
!$$$!kea
!$$$      include 'garofalini.inc'

      logical::lx,ly,lz,ovrlap,idum,ddum

      logical::lneighij,lclu_cmp,lexclude(nmax)

      integer(KIND=int)::i,ibox,flagon,iunit,j,imolty,icbu,ncount,ic,ip,k
      real(KIND=double_precision)::rx,ry,rz,dchain,ddx,ddy,ddz,random,vnew,vold
     &                 ,vintran,vintrao,deltv,deltvb,disvsq
     &                 ,vintern,vintero,vextn,vexto,rchain
     &                 ,velectn,velecto,vdum
     &                 ,vrecipo,vrecipn,v3n,v3o   
     & ,velectn_intra,velectn_inter,velecto_intra,velecto_inter 
      dimension ddum(27)

      logical::laccept

! --------------------------------------------------------------------

!      write(iou,*) 'start TRAXYZ'
      ovrlap = .false.
!     ***    select a chain at random ***
      rchain  = random()
      do icbu = 1,nmolty
         if ( rchain .lt. pmtrmt(icbu) ) then
            imolty = icbu
            rchain = 2.0d0
         end if
      end do

      if ((lexpee).and.(imolty.ge.nmolty1))
     &   imolty = ee_moltyp(mstate)
                                                                                
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
         if ( moltyp(i) .ne. imolty ) write(iou,*) 'screwup traxyz'


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
      call energy(i,imolty, vnew,vintran, vintern,vextn,velectn
     &     ,vdum,flagon, ibox,1, iunit,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)
      v3n=v3garo
      if (ovrlap) return

! *** calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto
     &     ,vdum,flagon,ibox,1, iunit,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)
      v3o = v3garo

      if (ovrlap) then
         call cleanup('disaster ovrlap in old conf of TRAXYZ')
      end if

      if ( lewald .and. lelect(imolty) ) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
            vrecipn    =   (1.0d0-(1.0d0-etais)*lambdais)*vrecipn
            vrecipo    =   (1.0d0-(1.0d0-etais)*lambdais)*vrecipo
         elseif (lstageb) then
            vrecipn =  etais*vrecipn
            vrecipo =  etais*vrecipo
         elseif (lstagec) then
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
         call anes(i,ibox,ibox,1,laccept,deltv,vintern,vintran,vextn,
     &        velectn,vintero,vintrao,vexto,velecto,vdum,vdum,
     &        vdum,vdum,.false.)
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
      elseif ( dexp(-deltvb) .gt. random() ) then
!        --- accept move
      else
!        --- move rejected
         return
      end if

!      write(iou,*) 'TRAXYZ accepted i',i

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

      if (lewald .and. lelect(imolty)) then
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
         call linkcell(2,i,vdum,vdum,vdum,ddum)
      end if

      if ( lneigh ) then
! *** check for update of near neighbour bitmap ***
! *** check for headgroup ***
         disvec(1,i,1) = disvec(1,i,1) + rx
         disvsq = disvec(1,i,1) * disvec(1,i,1) +
     &        disvec(1,i,2) * disvec(1,i,2) +
     &        disvec(1,i,3) * disvec(1,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
! *** check for last unit ***
         disvec(2,i,1) = disvec(2,i,1) + rx
         disvsq = disvec(2,i,1) * disvec(2,i,1) +
     &        disvec(2,i,2) * disvec(2,i,2) +
     &        disvec(2,i,3) * disvec(2,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
      end if

      if ( lneighbor .or. lgaro) then
         
         do ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
!            write(iou,*) ic,i,'j:',j
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

!      write(iou,*) 'end TRAXYZ',i

      return
      end
