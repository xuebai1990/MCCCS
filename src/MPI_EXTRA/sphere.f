      subroutine sphere ( x, y, z )

! sphere
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

!     ************************************************************
!     **  calculates a random vector on the unit sphere        ***
!     **  M. G. Martin   2-3-98                                ***
!     ************************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      integer(KIND=int)::ii
      real(KIND=double_precision)::random2, x, y, z
      real(KIND=double_precision)::xi1,xi2,xisq

!     --- calculate random vector on the unit sphere ---
      do ii = 1,100
! RP added for MPI ------changed random() to random2()
         xi1 = ( 2.0d0 * random2() ) - 1.0d0
         xi2 = ( 2.0d0 * random2() ) - 1.0d0
         xisq = xi1**2 + xi2**2
         if ( xisq .lt. 1.0d0 ) then
            x = 2.0d0 * xi1 * dsqrt( 1.0d0 - xisq )
            y = 2.0d0 * xi2 * dsqrt( 1.0d0 - xisq )
            z = ( 1.0d0 - 2.0d0 * xisq )
            return
         end if
      end do

      call cleanup('exceeded 100 tries to get a vector in sphere')

      end

      
