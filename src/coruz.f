      function coruz(imolty,rho,ibox)
!      function coruz(iunit,rho,ibox)

! coruz
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
! John Stubbs, and Collin Wick and Ilja Siepmann  
!                     
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to 
!
! Free Software Foundation, Inc. 
! 59 Temple Place - Suite 330
! Boston, MA  02111-1307, USA.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ***************************************************
! *** tail-corrections in energy with the zeolite ***
! ***************************************************

      implicit none
      include 'zeopoten.inc'
      include 'control.inc'
      include 'coord.inc'
      include 'poten.inc'
      include 'conver.inc'
      include 'system.inc'
      include 'expsix.inc'
      include 'merck.inc'
      include 'nsix.inc'
      include 'external.inc'

      real(8)::coruz,eps,rci3,rho
      integer::iunit,ibox
      real(8)::epsilon2,sigma2
      real(8)::rci1
      integer::imolty,jmolty,ii,jj, ntii, ntjj, ntij

! --- note works only for alkanes!!!
!      if (iunit.ne.1) then
!        eps=zeps(1,4)+(iunit-2)*zeps(2,4)+zeps(3,4)
!      else
!        eps=zeps(1,4)
!      end if
!      rci3=zsig2(1,4)**(3.d0/2.d0)/rcut(ibox)**3 
!      coruz=8.*onepi*eps*rho*(rci3*rci3*rci3/9.-rci3/3.)

      coruz=0.
      ntjj=ntsubst
      do ii = 1, nunit(imolty) 
         ntii = ntype(imolty,ii)
!         do jj = 1, nunit(jmolty) 
!            ntjj = ntype(jmolty,jj)
            ntij = (ntii-1)*nntype + ntjj
            rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
            if ( lexpand(imolty) ) then
               sigma2 = (sigma(imolty,ii)+sigi(ntjj))**2/4.0d0
               epsilon2 = dsqrt(epsilon(imolty,ii)*epsi(ntjj))
            else
               sigma2 = sig2ij(ntij)
               epsilon2 = epsij(ntij)
            end if
            coruz = coruz + 
     &           8.0d0 * onepi * epsilon2 * 
     &           sigma2**(1.5d0) *rho * 
     &           (rci3 * rci3 * rci3 / 9.0d0 - rci3 / 3.0d0)
!            end do
      end do
      return
      end
