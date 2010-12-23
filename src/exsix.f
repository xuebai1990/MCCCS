      function exsix(rijsq,ntij)

!     **************************************************************
!     ***  calculates the energy using the exp-6 potential       ***
!     ***  parameters defined in suijtab.f  M.G. Martin          ***
!     **************************************************************
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
!$$$      include 'expsix.inc'
      real(KIND=double_precision)::rijsq,rij,exsix
      integer(KIND=int)::ntij

      rij=dsqrt(rijsq)
      exsix = aexsix(ntij)/(rijsq*rijsq*rijsq)
     &     + bexsix(ntij)*dexp(cexsix(ntij)*rij)
      if (lshift) exsix = exsix-sexsix(ntij)
      return
      end

      
