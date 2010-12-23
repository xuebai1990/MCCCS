      subroutine lininter_vib(len, tabulated_vib, vibtyp)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  calculates vibrational potential using linear interpolation
!ccc  between two points
!ccc  must specify equilibrium bond length in suvibe, force
!ccc  constant must be zero
!ccc  requires a file (fort.41) that starts with 0.0 (not 0.5)
!ccc  fort.41: number of tabulated potentials, potential number from
!ccc  suvibe, number of points per angstrom, tabulated potential
!ccc  (repeat last three parts for each additional potential)
!ccc  KM 12/02/08
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!$$$      include 'tabulated.inc'
!$$$      include 'conver.inc'

      integer(KIND=normal_int)::vibtyp, vlow, vhigh, xa, bin, add1
      real(KIND=double_precision)::len, tabulated_vib, left, lenrem
      
      vlow=1
      vhigh=vibsplits(vibtyp)

      xa = int(len)
      left = len-xa

!     select correct bin
!     multiply by num_int_vib - number of intervals per angstrom
      add1 = int(left*num_int_vib(vibtyp))
      bin = (xa-vib(1,vibtyp))*num_int_vib(vibtyp) + 1 + add1
      vlow = bin
      vhigh = vlow+1

      if (vib(vlow,vibtyp).gt. len.or.vib(vhigh,vibtyp).lt.len) then
         write(2,*) 'problem in lininter_vib!'
         write(2,*) 'len', len, ' vibtyp', vibtyp
         write(2,*) 'vlow ', vlow, vib(vlow, vibtyp),len
         write(2,*) 'vhigh ', vhigh, vib(vhigh, vibtyp), len
         write(2,*)
      end if

      lenrem=len-vib(vlow, vibtyp)
      lenrem=lenrem*dble(num_int_vib(vibtyp))

      tabulated_vib=lenrem*vibdiff(vlow,vibtyp)
      tabulated_vib=tabulated_vib+tabvib(vlow,vibtyp)

      return
      end

      subroutine lininter_bend(r, tabulated_bend, bendtyp)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  calculates 1-3 nonbonded 'bending' potential using linear
!ccc  interpolation between two points
!ccc  must include 1-3 interactions
!ccc  must specify equilibrium angle in suvibe, force constant
!ccc  must be very small but non-zero
!ccc  requires a file (fort.42) with distances in A
!ccc  fort.42: number of tabulated potentials, potential number from
!ccc  suvibe, number of points per degree, tabulated potential
!ccc  (repeat last three parts for each additional potential,
!ccc   separated by 1000 1000)
!ccc  make sure potential does not go up to infinity!
!ccc  KM 12/03/08
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!$$$      include 'tabulated.inc'
!$$$      include 'conver.inc'

      integer(KIND=normal_int)::bendtyp, blow, bhigh, xa, bin, add1
      real(KIND=double_precision)::r, tabulated_bend, left, rem

      blow=1
      bhigh=bendsplits(bendtyp)

      xa = int(r)
      left = r-xa

!     select correct bin
!     multiply by num_int_bend - number of intervals per A
      add1 = int(left*num_int_bend(bendtyp))
      bin = (xa-bend(1,bendtyp))*num_int_bend(bendtyp) + 1 + add1
      blow = bin
      bhigh = blow+1

      if (bend(blow,bendtyp).gt. r.or.
     &     bend(bhigh,bendtyp).lt.r) then
         write(2,*) 'problem in lininter_bend!'
         write(2,*) 'r', r, ' bendtyp', bendtyp
         write(2,*) 'blow ', blow, bend(blow, bendtyp),r
         write(2,*) 'bhigh ', bhigh, bend(bhigh, bendtyp), r
         write(2,*)
      end if

      rem=r-bend(blow, bendtyp)
      rem=rem*dble(num_int_bend(bendtyp))

      tabulated_bend=rem*benddiff(blow,bendtyp)
      tabulated_bend=tabulated_bend+tabbend(blow,bendtyp)
!      write(2,*) 'tabulated_bend ', tabulated_bend

      return
      end

      subroutine lininter_vdW(r, tabulated_vdW, typi, typj)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  calculates nonbonding van der Waals potential using linear 
!ccc  interpolation between two points
!ccc  requires a file (fort.43)
!ccc  fort.43: number of tabulated potentials, potential number,
!ccc  number of points per angstrom, tabulated potential
!ccc  (repeat last three parts for each additional potential)
!ccc  for unlike interactions, list 1-2 and 2-1
!ccc  separate potentials with 1000
!ccc  make sure potential does not go up to infinity!
!ccc  bead type numbers should be defined in suijtab, but it doesn't
!ccc  matter what they are (the parameters aren't used)
!ccc  KM 12/03/08
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!$$$      include 'tabulated.inc'
!$$$      include 'conver.inc'

      integer(KIND=normal_int)::typi, typj, low, high, xa, bin, add1
      real(KIND=double_precision)::r,  tabulated_vdW, left, rem

      low=1
      high=vdWsplits(typi, typj)

      xa = int(r)
      left = r-xa

!     select correct bin
!     multiply by num_int_vdW - number of intervals per angstrom
      add1 = int(left*num_int_vdW(typi,typj))
      bin = (xa-rvdW(1,typi,typj))*num_int_vdW(typi,typj) + 
     &     1 + add1
      low = bin
      high = low+1

      if (rvdW(low,typi,typj).gt. r.or.rvdW(high,typi,typj)
     &     .lt.r) then
         write(2,*) 'problem in lininter_vdW!'
         write(2,*) 'r', r, ' typi', typi, ' typj ', typj
         write(2,*) 'low ', low, rvdW(low, typi, typj)
         write(2,*) 'high ', high, rvdW(high, typi, typj)
         write(2,*)
      end if

      rem=r-rvdW(low, typi, typj)
      rem=rem*dble(num_int_vdW(typi,typj))

      tabulated_vdW=rem*vdWdiff(low,typi, typj)
      tabulated_vdW=tabulated_vdW+tabvdW(low,typi, typj)

      return
      end


      subroutine lininter_elect(r, tabulated_elect, typi, typj)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  calculates electrostatic potential using linear interpolation
!ccc  between two points
!ccc  requires a file (fort.44)
!ccc  fort.44: number of tabulated potentials, potential number,
!ccc  number of points per angstrom, tabulated potential
!ccc  (repeat last three parts for each additional potential)
!ccc  for unlike interactions, list 1-2 and 2-1
!ccc  separate potentials with 1000
!ccc  KM 04/23/09
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
!$$$      include 'tabulated.inc'
!$$$      include 'conver.inc'

      integer(KIND=normal_int)::typi, typj, low, high, xa, bin, add1
      real(KIND=double_precision)::r,  tabulated_elect, left, rem

      low=1
      high=electsplits(typi, typj)

      xa = int(r)
      left = r-xa

!     select correct bin
!     multiply by num_int_nonbond - number of intervals per angstrom
      add1 = int(left*num_int_elect(typi,typj))
      bin = (xa-relect(1,typi,typj))*num_int_elect(typi,typj) + 
     &     1 + add1
      low = bin
      high = low+1

      if (relect(low,typi,typj).gt. r.or.relect(high,typi,typj)
     &     .lt.r) then
         write(2,*) 'problem in lininter_elect!'
         write(2,*) 'r', r, ' typi', typi, ' typj ', typj
         write(2,*) 'low ', low, relect(low, typi, typj)
         write(2,*) 'high ', high, relect(high, typi, typj)
         write(2,*)
      end if

      rem=r-relect(low, typi, typj)
      rem=rem*dble(num_int_elect(typi,typj))

      tabulated_elect=rem*electdiff(low,typi, typj)
      tabulated_elect=tabulated_elect+tabelect(low,typi, typj)

      return
      end

