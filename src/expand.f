      subroutine expand

      implicit none

!    **********************************************************************
!    ** make a transition of a selected molecule from state i to state j **
!    ** in an expanded-ensemble sampling.                                **
!    ** written on Aug. 4/99 by Bin Chen.                                **
!    **********************************************************************
      
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'bnbsma.inc'
      include 'inputdata.inc'
      include 'expand.inc'
      include 'poten.inc'

      logical::ovrlap
      integer::i,ibox,iunit,flagon,itype,j,imolty,icbu,ic,imt,jmt,
     &     itype2,disp
      real(8)::dchain,random,vnew,vold
     &                 ,vintran,vintrao,deltv,deltvb,disvsq
     &                 ,vintern,vintero,vextn,vexto,rchain
     &                 ,velectn,velecto,vdum
     &                 ,vrecipo,vrecipn,vexpta,vexptb,volume,rho,coru
      logical::laccept

!      write(iou,*) 'start expand-ensemble move'
! ***    select a chain at random ***
      dchain  = random()
      do icbu = 1,nmolty
         if ( dchain .lt. pmeemt(icbu) ) then
            imolty = icbu
            dchain = 2.0d0
         end if
      end do
      if ( .not. lexpand(imolty) ) 
     &     call cleanup('wrong type of molecule for the ES-move')

      if (lgrand) then
! ---    select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         if (nchbox(ibox).eq.0) then
            bnexpc(imolty,ibox) = bnexpc(imolty,ibox) + 1.0d0
            return
         end if
         i = dint( dble(ncmt(1,imolty))*random() ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(iou,*) 'screwup'

      else

         dchain = dble(temtyp(imolty))
         i = int( dchain*random() ) + 1
         i = parall(imolty,i)
         ibox = nboxi(i)
      end if

      iunit = nunit(imolty)

! *** perform a move in the expanded coefficients

 10   disp = int( rmexpc(imolty)*(2.0d0*random()-1.0d0) )
      itype = mod(eetype(imolty)+disp+numcoeff(imolty)
     &     , numcoeff(imolty))
      if ( disp .eq. 0 ) goto 10
      if ( itype .eq. 0 ) itype = numcoeff(imolty)
      do ic = 1,2
         if ( ic .eq. 1 ) then
            itype2 = eetype(imolty)
         else
            itype2 = itype 
         end if
         do j = 1,iunit
            rxuion(j,ic) = rxu(i,j)
            ryuion(j,ic) = ryu(i,j)
            rzuion(j,ic) = rzu(i,j)
            qquion(j,ic) = qcharge(imolty,j,itype2)
         end do
         moltion(ic) = imolty
      end do

! *** calculate the energy of i in the new configuration ***

      flagon = 2
      do j = 1,iunit
         epsilon(imolty,j) = epsil(imolty,j,itype)
         sigma(imolty,j) = sigm(imolty,j,itype)
      end do
      call energy(i,imolty, vnew,vintran, vintern,vextn,velectn
     &     ,vdum,flagon, ibox,1, iunit,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)
      if (ovrlap) return

!     Start of intermolecular tail correction for new

      if ( ltailc ) then

         volume = boxlx(ibox)*boxly(ibox)*boxlz(ibox)

         vexpta = 0.0d0

         do imt = 1, nmolty
            do jmt = 1, nmolty
               rho = dble(ncmt(ibox,jmt))/volume
               vexpta = vexpta +
     &              dble( ncmt(ibox,imt) ) * coru(imt,jmt,rho,ibox)
            end do
         end do
         vnew = vnew + vexpta
         vintern = vintern + vexpta
      end if

! *** calculate the energy of i in the old configuration ***
      flagon = 1
      do j = 1, iunit
         epsilon(imolty,j) = epsil(imolty,j,eetype(imolty))
         sigma(imolty,j) = sigm(imolty,j,eetype(imolty))
      end do
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto
     &     ,vdum,flagon,ibox,1, iunit,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)


!     Start of intermolecular tail correction for old

      if (ltailc ) then
         vexptb = 0.0d0
         do imt = 1, nmolty
            do jmt = 1, nmolty
               rho = ncmt(ibox,jmt) / volume
               vexptb = vexptb + 
     &              dble(ncmt(ibox,imt)) * coru(imt,jmt,rho,ibox)
            end do
         end do
         vold = vold + vexptb
         vintero = vintero + vexptb
         
      end if

      bnexpc(imolty,ibox) = bnexpc(imolty,ibox) + 1.0d0

      if (ovrlap) call cleanup('disaster ovrlap in old conf of TRAXYZ')
      
      if ( lewald ) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if

! *** check for acceptance ***
 
      deltv  = vnew - vold + eta(ibox,imolty,itype) 
     &     - eta(ibox,imolty,eetype(imolty))
      deltvb = beta * deltv

      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
!        --- accept move
      elseif ( dexp(-deltvb) .gt. random() ) then
!        --- accept move
      else
!        --- move rejected
         return
      end if

!      write(iou,*) 'expanded move accepted i',i,exp_cion(2)
      vbox(ibox)     = vbox(ibox) + vnew - vold
      vinterb(ibox)  = vinterb(ibox) + (vintern - vintero)
      vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
      vextb(ibox)    = vextb(ibox)   + (vextn   - vexto)
      velectb(ibox)   = velectb(ibox)  + (velectn - velecto)
      vtailb(ibox) = vtailb(ibox) + vexpta - vexptb
      
      ncmt2(ibox,imolty,itype) = ncmt2(ibox,imolty,itype) + 1
      ncmt2(ibox,imolty,eetype(imolty)) = 
     &     ncmt2(ibox,imolty,eetype(imolty)) - 1
      eetype(imolty) = itype
      do j = 1,iunit
         qqu(i,j) = qquion(j,2)
         epsilon(imolty,j) = epsil(imolty,j,itype)
         sigma(imolty,j) = sigm(imolty,j,itype)
      end do


      if (lewald) then
! *** update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if (ldielect) then
          call dipole(ibox,1)
      end if

      bsexpc(imolty,ibox) = bsexpc(imolty,ibox) + 1.0d0

!      write(iou,*) 'end expand-ensemble move'

      return
      end




