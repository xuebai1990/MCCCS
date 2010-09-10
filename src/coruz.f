      function coruz(imolty,rho,ibox)
c      function coruz(iunit,rho,ibox)

c coruz
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

c ***************************************************
c *** tail-corrections in energy with the zeolite ***
c ***************************************************

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

c --- note works only for alkanes!!!
c      if (iunit.ne.1) then
c        eps=zeps(1,4)+(iunit-2)*zeps(2,4)+zeps(3,4)
c      else
c        eps=zeps(1,4)
c      endif
c      rci3=zsig2(1,4)**(3.d0/2.d0)/rcut(ibox)**3 
c      coruz=8.*onepi*eps*rho*(rci3*rci3*rci3/9.-rci3/3.)

      coruz=0.
      ntjj=ntsubst
      do ii = 1, nunit(imolty) 
         ntii = ntype(imolty,ii)
c         do jj = 1, nunit(jmolty) 
c            ntjj = ntype(jmolty,jj)
            ntij = (ntii-1)*nntype + ntjj
            rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
            if ( lexpand(imolty) ) then
               sigma2 = (sigma(imolty,ii)+sigi(ntjj))**2/4.0d0
               epsilon2 = dsqrt(epsilon(imolty,ii)*epsi(ntjj))
            else
               sigma2 = sig2ij(ntij)
               epsilon2 = epsij(ntij)
            endif
            coruz = coruz + 
     +           8.0d0 * onepi * epsilon2 * 
     +           sigma2**(1.5d0) *rho * 
     +           (rci3 * rci3 * rci3 / 9.0d0 - rci3 / 3.0d0)
c            enddo
      enddo
      return
      end
