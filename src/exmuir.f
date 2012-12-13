       function exmuir ( z, ntj )
 
!    *********************************************************
!    ** calculates the lmuir external energy for a bead.    **
!    *********************************************************
 
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
!$$$      include 'externalmuir.inc'

      real(KIND=double_precision)::exmuir, z
      integer(KIND=normal_int)::ntj

! --------------------------------------------------------------------

      if ( ntj .eq. 1 ) then
! --- HEADgroup potential ---
         if ( z .le. alpha2 ) then
            exmuir = 0.0d0
         else
            exmuir = beta2 / ( 1.0d0 + ( (z/alpha2) - 1.0d0 )**tau2 )
         end if
      else
! --- TAILgroup potential ---
         if ( z .ge. alpha1 ) then
            exmuir = 0.0d0
         else
            if ( z .lt. zprmin ) then
               if ( ntj .eq. 2 ) then
                  exmuir = betac2 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + v2prmin
               else
                  exmuir = betac3 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + v3prmin
               end if
            else
               if ( ntj .eq. 2 ) then
                  exmuir = betac2 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + c9ch2 / z**9 -  c3ch2 / z**3
               else
                  exmuir = betac3 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + c9ch3 / z**9 -  c3ch3 / z**3
               end if
            end if
         end if
      end if
      return

! ----------------------------------------------------------------------------

      end
