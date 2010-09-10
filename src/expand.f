      subroutine expand

c expand
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
c John Stubbs, and Collin Wick and Ilja Siepmann  
c                     
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to 
c
c Free Software Foundation, Inc. 
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

c    **********************************************************************
c    ** make a transition of a selected molecule from state i to state j **
c    ** in an expanded-ensemble sampling.                                **
c    ** written on Aug. 4/99 by Bin Chen.                                **
c    **********************************************************************
      
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
     +                 ,vintran,vintrao,deltv,deltvb,disvsq
     +                 ,vintern,vintero,vextn,vexto,rchain
     &                 ,velectn,velecto,vdum
     &                 ,vrecipo,vrecipn,vexpta,vexptb,volume,rho,coru
      logical::laccept

c      write(iou,*) 'start expand-ensemble move'
c ***    select a chain at random ***
      dchain  = random()
      do icbu = 1,nmolty
         if ( dchain .lt. pmeemt(icbu) ) then
            imolty = icbu
            dchain = 2.0d0
         endif
      enddo
      if ( .not. lexpand(imolty) ) 
     &     call cleanup('select a wrong type of molecule for the ES-move')

      if (lgrand) then
c ---    select a chain at random in box 1!
c         (in box 2 is an ideal gas!)
         ibox = 1
         if (nchbox(ibox).eq.0) then
            bnexpc(imolty,ibox) = bnexpc(imolty,ibox) + 1.0d0
            return
         endif
         i = dint( dble(ncmt(1,imolty))*random() ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(iou,*) 'screwup'

      else

         dchain = dble(temtyp(imolty))
         i = int( dchain*random() ) + 1
         i = parall(imolty,i)
         ibox = nboxi(i)
      endif

      iunit = nunit(imolty)

c *** perform a move in the expanded coefficients

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
         endif
         do j = 1,iunit
            rxuion(j,ic) = rxu(i,j)
            ryuion(j,ic) = ryu(i,j)
            rzuion(j,ic) = rzu(i,j)
            qquion(j,ic) = qcharge(imolty,j,itype2)
         enddo
         moltion(ic) = imolty
      enddo

c *** calculate the energy of i in the new configuration ***

      flagon = 2
      do j = 1,iunit
         epsilon(imolty,j) = epsil(imolty,j,itype)
         sigma(imolty,j) = sigm(imolty,j,itype)
      enddo
      call energy(i,imolty, vnew,vintran, vintern,vextn,velectn
     &     ,vdum,flagon, ibox,1, iunit,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)
      if (ovrlap) return

c     Start of intermolecular tail correction for new

      if ( ltailc ) then

         volume = boxlx(ibox)*boxly(ibox)*boxlz(ibox)

         vexpta = 0.0d0

         do imt = 1, nmolty
            do jmt = 1, nmolty
               rho = dble(ncmt(ibox,jmt))/volume
               vexpta = vexpta +
     &              dble( ncmt(ibox,imt) ) * coru(imt,jmt,rho,ibox)
            enddo
         enddo
         vnew = vnew + vexpta
         vintern = vintern + vexpta
      endif

c *** calculate the energy of i in the old configuration ***
      flagon = 1
      do j = 1, iunit
         epsilon(imolty,j) = epsil(imolty,j,eetype(imolty))
         sigma(imolty,j) = sigm(imolty,j,eetype(imolty))
      enddo
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto
     &     ,vdum,flagon,ibox,1, iunit,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)


c     Start of intermolecular tail correction for old

      if (ltailc ) then
         vexptb = 0.0d0
         do imt = 1, nmolty
            do jmt = 1, nmolty
               rho = ncmt(ibox,jmt) / volume
               vexptb = vexptb + 
     &              dble(ncmt(ibox,imt)) * coru(imt,jmt,rho,ibox)
            enddo
         enddo
         vold = vold + vexptb
         vintero = vintero + vexptb
         
      endif

      bnexpc(imolty,ibox) = bnexpc(imolty,ibox) + 1.0d0

      if (ovrlap) call cleanup('disaster ovrlap in old conf of TRAXYZ')
      
      if ( lewald ) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      endif

c *** check for acceptance ***
 
      deltv  = vnew - vold + eta(ibox,imolty,itype) 
     &     - eta(ibox,imolty,eetype(imolty))
      deltvb = beta * deltv

      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
c        --- accept move
      elseif ( dexp(-deltvb) .gt. random() ) then
c        --- accept move
      else
c        --- move rejected
         return
      endif

c      write(iou,*) 'expanded move accepted i',i,exp_cion(2)
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
      enddo


      if (lewald) then
c *** update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      endif

      if (ldielect) then
          call dipole(ibox,1)
      endif

      bsexpc(imolty,ibox) = bsexpc(imolty,ibox) + 1.0d0

c      write(iou,*) 'end expand-ensemble move'

      return
      end




