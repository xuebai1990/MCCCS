      subroutine sphere ( x, y, z )

!     ************************************************************
!     **  calculates a random vector on the unit sphere        ***
!     **  M. G. Martin   2-3-98                                ***
!     ************************************************************

      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_runtime,only:err_exit
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
      integer(KIND=normal_int)::ii
      real(KIND=double_precision)::random, x, y, z
      real(KIND=double_precision)::xi1,xi2,xisq

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

      call err_exit('exceeded 100 tries to get a vector in sphere')

      end

      
