       function ljpsur ( rijsq, ntij )
 
!    *********************************************************
!    ** calculates energy for a polymeric surfactant bead.  **
!    *********************************************************
 
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none

!$$$      include 'external.inc'

      real(KIND=double_precision)::ljpsur, rijsq, sr, sr6
      integer(KIND=int)::ntij

! --------------------------------------------------------------------
! AT PRESENT: all sigma = 1.0
!             epsilon   = 1.0
! --------------------------------------------------------------------

      sr = 1.0d0 / rijsq
      sr6 = sr**3

      if ( ntij .eq. 1 ) then
! *** nonpolar-nonpolar interaction ( LJ interaction )
         ljpsur = sr6*sr6 - sr6
      else
! *** polar-polar or polar-nonpolar interaction ( repulsive LJ interaction )
         if ( rijsq .le. 1.259921d0 ) then
            ljpsur = sr6*sr6 - sr6 + 0.25d0
         else
            ljpsur = 0.0d0
         end if
      end if

      return

! ----------------------------------------------------------------------------

      end
