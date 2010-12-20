      subroutine sphere ( x, y, z )

!     ************************************************************
!     **  calculates a random vector on the unit sphere        ***
!     **  M. G. Martin   2-3-98                                ***
!     ************************************************************

      implicit none
      integer::ii
      real(8)::random, x, y, z
      real(8)::xi1,xi2,xisq

!     --- calculate random vector on the unit sphere ---
      do ii = 1,100
         xi1 = ( 2.0d0 * random() ) - 1.0d0
         xi2 = ( 2.0d0 * random() ) - 1.0d0
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

      
