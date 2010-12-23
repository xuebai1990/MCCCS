      function mmff(rijsq,ntij)

!   ******************************************************************
!   ***  calculates the energy using the Buffered 14-7 potential   ***
!   ***  parameters defined in suijtab.f          Bin Chen         *** 
!   ******************************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
!$$$      include 'control.inc'
!$$$      include 'merck.inc'
      real(KIND=double_precision)::rijsq,mmff,rs2,rs1,sr1,sr7,rs7
      integer(KIND=int)::ntij

        rs2 = rijsq / (sigisq(ntij))
        rs1 = dsqrt(rs2)
        rs7 = rs1*rs2*rs2*rs2
        sr1 = 1.07d0/(rs1+0.07d0)
        sr7 = sr1**7.0d0
        mmff = sr7*(1.12d0/(rs7+0.12d0) - 2.0d0)
     &           *epsimmff(ntij)

      if (lshift) mmff = mmff-smmff(ntij)
      return
      end
