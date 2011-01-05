      subroutine Atom_traxyz (lx,ly,lz )
 
!    *******************************************************************
!    ** makes a translational movement in x,y,or z-direction.         **
!    ** the maximum displacement is controlled by rAtrax(yz) and the  **
!    ** number of successful trial moves is stored in Abstrax(yz).     **
!    ** The attempts are stored in Abntrax(yz)                         **
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

      logical::lx,ly,lz,ovrlap,ddum

      logical::lneighij

      integer(KIND=normal_int)::i,ibox,flagon,iunit,j,imolty,icbu
     & ,ic,ip
      integer(KIND=normal_int)::pick_unit, pick_chain
      real(KIND=double_precision)::rx,ry,rz,dchain,random
     & ,vnew,vold,vintran,vintrao,deltv,deltvb,disvsq,vintern,vintero
     & ,vextn,vexto,rchain,velectn,velecto,vdum,vrecipo,vrecipn
 
      real(KIND=double_precision)::vvibn,vbendn,vtgn,vvibo,vbendo,vtgo

      real(KIND=double_precision)::vewaldn, vewaldo   

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
         if ( moltyp(pick_chain) .ne. imolty ) 
     &           write(iou,*) 'screwup traxyz'


      else

         dchain = dble(temtyp(imolty))
         pick_chain = int( dchain*random() + 1 )
         pick_chain = parall(imolty,pick_chain)
         ibox = nboxi(pick_chain)

      end if

! *** store number of units of i in iunit ***

      iunit = nunit(imolty)

      pick_unit = int(dble(iunit*random()) + 1 )

!      write(iou,*) pick_unit, imolty, pick_chain

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
       call Atom_energy(i,imolty, vnew,vintran, vintern,vextn,velectn
     &     ,vewaldn,flagon, ibox,pick_unit, pick_unit,.true.,ovrlap,
     &      .false.
     &     ,vdum,.false.,.false.,vvibn,vbendn,vtgn)
      if (ovrlap) return

! *** calculate the energy of i in the old configuration ***
      flagon = 1
      call Atom_energy(i,imolty,vold,vintrao,vintero,vexto,velecto
     &     ,vewaldo,flagon,ibox,pick_unit, pick_unit,.true.,ovrlap,
     &       .false.
     &     ,vdum,.false.,.false.,vvibo,vbendo,vtgo)

      if (ovrlap) call cleanup('disaster ovrlap in old conf of TRAXYZ')
      
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
         call anes(i,ibox,ibox,1,laccept,deltv,vintern,vintran,vextn,
     &        velectn,vintero,vintrao,vexto,velecto,vdum,vdum,
     &        vdum,vdum,.false.)
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

      if ( lneighbor ) then
         
         do 10 ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
!            write(iou,*) ic,i,'j:',j
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

!      write(iou,*) 'end TRAXYZ',i

      return
      end
