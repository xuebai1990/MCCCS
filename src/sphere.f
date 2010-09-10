      subroutine sphere ( x, y, z )

c sphere
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

c     ************************************************************
c     **  calculates a random vector on the unit sphere        ***
c     **  M. G. Martin   2-3-98                                ***
c     ************************************************************

      implicit none
      integer::ii
      real(8)::random, x, y, z
      real(8)::xi1,xi2,xisq

c     --- calculate random vector on the unit sphere ---
      do ii = 1,100
         xi1 = ( 2.0d0 * random() ) - 1.0d0
         xi2 = ( 2.0d0 * random() ) - 1.0d0
         xisq = xi1**2 + xi2**2
         if ( xisq .lt. 1.0d0 ) then
            x = 2.0d0 * xi1 * dsqrt( 1.0d0 - xisq )
            y = 2.0d0 * xi2 * dsqrt( 1.0d0 - xisq )
            z = ( 1.0d0 - 2.0d0 * xisq )
            return
         endif
      enddo

      stop 'exceeded 100 tries to get a vector in sphere'

      end

      
