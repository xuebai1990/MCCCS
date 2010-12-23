      function ninesix(rijsq,ntij)

!     ***************************************************
!     ***  calculates the energy of the 9-6 potential ***
!     ***  parameters defined in suijtab.f  JMS       ***
!     ***************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'control.inc'
!$$$      include 'nsix.inc'
      real(KIND=double_precision)::rijsq,rij,ror,ninesix
      integer(KIND=int)::ntij

      rij=dsqrt(rijsq)
      ror = rzero(ntij)/rij
      ninesix = 4.0d0*epsnx(ntij)*ror**6*(2.0d0*ror**3 - 3.0d0)
      if (lshift) then
         ninesix = ninesix - shiftnsix(ntij)
      end if

      return
      end

      
