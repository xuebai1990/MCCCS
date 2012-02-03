      subroutine susami

c susami
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
 
c    *******************************
c    ** Set-Up SAMI's potentials. **
c    *******************************
 
      implicit none

c *** common blocks ***
      include 'external.inc'
      include 'ljsamipara.inc'

      double precision hsig,heps,tsig,teps

      parameter ( hsig=4.220d0, heps=110.68816d0, 
     &            tsig=3.527d0, teps=79.982210d0 )

      double precision rcsami
      integer ij

C --------------------------------------------------------------------

      rcsami = 2.5d0*tsig

      sij(1)=hsig
      sij(2)=0.5d0*(hsig+tsig)
      sij(3)=sij(2)
      sij(4)=sij(2)
      sij(5)=tsig
      sij(6)=tsig
      sij(7)=sij(2)
      sij(8)=tsig
      sij(9)=tsig

      eij(1)=heps
      eij(2)=dsqrt(heps*teps)
      eij(3)=eij(2)
      eij(4)=eij(2)
      eij(5)=teps
      eij(6)=teps
      eij(7)=eij(2)
      eij(8)=teps
      eij(9)=teps

      vsh(1)  = eij(1) *
     &          ( ( 13.0d0 * (sij(1)/rcsami)**12 ) +
     &            (  4.0d0 * (sij(1)/rcsami)**3  ) )
      vsha(1) = eij(1) *
     &          ( ( 12.0d0 * sij(1)**12 / rcsami**13 ) +
     &            (  3.0d0 * sij(1)**3  / rcsami**4  ) )

      do ij = 2, 9
         vsh(ij)  = 4.0d0 * eij(ij) * 
     &          ( ( 13.0d0 * (sij(ij)/rcsami)**12 ) -
     &            (  7.0d0 * (sij(ij)/rcsami)**6  ) )
         vsha(ij) = 4.0d0 * eij(ij) *
     &          ( ( 12.0d0 * sij(ij)**12 / rcsami**13 ) -
     &            (  6.0d0 * sij(ij)**6  / rcsami**7  ) )
      enddo

c      do ij = 1,9
c         write(2,*) 'ij',ij,'vsh',(vsh(ij)/80.0d0),
c     +                      'vsha',(vsha(ij)/80.0d0),
c     +                      'eij',(eij(ij)/80.0d0)
c      enddo

      return

C ----------------------------------------------------------------------------

      end
