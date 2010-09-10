      function corp(imolty,jmolty,rhosq,ibox)

c corp 
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

c ************************************
c *** tail-corrections in pressure ***
c ************************************

      implicit none
      include 'control.inc'
      include 'coord.inc'
      include 'conver.inc'
      include 'poten.inc'
      include 'system.inc'
      include 'expsix.inc'
      include 'merck.inc'
      include 'nsix.inc'

      real(8)::corp,rci3,rhosq, epsilon2, sigma2
      real(8)::rci1
      integer::imolty,jmolty,ii,jj, ntii, ntjj, ntij,ibox
      corp = 0.0d0

      do ii = 1, nunit(imolty) 
         ntii = ntype(imolty,ii)

         do jj = 1, nunit(jmolty) 
            ntjj = ntype(jmolty,jj)

            if (lexpsix) then
               ntij = (ntii+ntjj)/2
               corp = corp + rhosq*consp(ntij)
            elseif (lmmff) then
               ntij = (ntii+ntjj)/2
               corp = corp+((-2.0d0)/3.0d0)*onepi*rhosq*epsimmff(ntij)
     +              *sigimmff(ntij)**3.0d0*corp_cons(ntij)
            elseif (lninesix) then
               ntij = (ntii-1)*nxatom + ntjj
               corp = corp + 16.0d0 * onepi * epsnx(ntij) *
     &            rhosq * rzero(ntij)**3 * 
     &            (0.5d0*(rzero(ntij)/rcut(ibox))**6 - 
     &            (rzero(ntij)/rcut(ibox))**3)
            elseif (lgenlj) then
               ntij = (ntii-1)*nntype + ntjj
               rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
               rci1 = rci3 **(1.0d0/3.0d0)
               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigma(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)
     &                 *epsilon(jmolty,jj))
               elseif ( lexpand(imolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigi(ntjj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)*epsi(ntjj))
               elseif ( lexpand(jmolty) ) then
                  sigma2 = (sigma(jmolty,jj)+sigi(ntii))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               endif
                 corp = corp
     &  + (2.0d0/3.0d0) * onepi * epsilon2 * sigma2 ** (1.50d0)
     &    * rhosq * n1 *
     &   ( (( (2.0d0**((4.0d0*n1/n0)+1.0d0))/(2.0d0*n1-3.0d0)) *
     & rci1 **(2.0d0*n1-3.0d0) ) -
     &   ( ((2.0d0**((2.0d0*n1/n0)+1.0d0))/(n1-3.0d0))*
     & rci1 **(n1-3.0d0) )  )
            else
               ntij = (ntii-1)*nntype + ntjj
               
               rci3 = sig2ij(ntij)**(3.d0/2.d0) / rcut(ibox)**3 
               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigma(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                 epsilon(jmolty,jj))
               elseif ( lexpand(imolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigi(ntjj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)*epsi(ntjj))
               elseif ( lexpand(jmolty) ) then
                  sigma2 = (sigma(jmolty,jj)+sigi(ntii))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               endif
               corp = corp + 
     +          32.0d0 * onepi * epsilon2 * sigma2**(1.5d0)*
     +          rhosq * ( rci3*rci3*rci3 / 9.0d0 - rci3 / 6.0d0 )
            endif
         enddo
      enddo

      return
      end


