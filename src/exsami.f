       function exsami ( z, ntj )

!    **********************************************************
!    ** calculates the SAMI's external energy for a bead.    **
!    **********************************************************
 
      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'external.inc'

      real(KIND=double_precision)::exsami, z
      integer(KIND=normal_int)::ntj

! --------------------------------------------------------------------

      if ( ntj .eq. 1 ) then
! --- HEADgroup potential ---
         if ( z .le. alpha2 ) then
            exsami = 0.0d0
         else
            exsami = beta2 / ( 1.0d0 + ( (z/alpha2) - 1.0d0 )**tau2 )
         end if
      else
! --- TAILgroup potential ---
         if ( z .ge. alpha1 ) then
            exsami = 0.0d0
         else
            exsami = beta1 / ( 1.0d0 + ( 1.0d0 - (z/alpha1) )**tau1 )
         end if
      end if

      return

! ----------------------------------------------------------------------------

      end
