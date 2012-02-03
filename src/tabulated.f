      subroutine lininter_vib(len, tabulated_vib, vibtyp)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  calculates vibrational potential using linear interpolation
cccc  between two points
cccc  must specify equilibrium bond length in suvibe, force
cccc  constant must be zero
cccc  requires a file (fort.41) that starts with 0.0 (not 0.5)
cccc  fort.41: number of tabulated potentials, potential number from
cccc  suvibe, number of points per angstrom, tabulated potential
cccc  (repeat last three parts for each additional potential)
cccc  KM 12/02/08
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'tabulated.inc'
      include 'conver.inc'

      integer vibtyp, vlow, vhigh, xa, bin, add1
      double precision len, tabulated_vib, left, lenrem
      
      vlow=1
      vhigh=vibsplits(vibtyp)

      xa = int(len)
      left = len-xa

c     select correct bin
c     multiply by num_int_vib - number of intervals per angstrom
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
      endif

      lenrem=len-vib(vlow, vibtyp)
      lenrem=lenrem*dble(num_int_vib(vibtyp))

      tabulated_vib=lenrem*vibdiff(vlow,vibtyp)
      tabulated_vib=tabulated_vib+tabvib(vlow,vibtyp)

      return
      end

      subroutine lininter_bend(r, tabulated_bend, bendtyp)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  calculates 1-3 nonbonded 'bending' potential using linear
cccc  interpolation between two points
cccc  must include 1-3 interactions
cccc  must specify equilibrium angle in suvibe, force constant
cccc  must be very small but non-zero
cccc  requires a file (fort.42) with distances in A
cccc  fort.42: number of tabulated potentials, potential number from
cccc  suvibe, number of points per degree, tabulated potential
cccc  (repeat last three parts for each additional potential,
cccc   separated by 1000 1000)
cccc  make sure potential does not go up to infinity!
cccc  KM 12/03/08
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'tabulated.inc'
      include 'conver.inc'

      integer bendtyp, blow, bhigh, xa, bin, add1
      double precision r, tabulated_bend, left, rem

      blow=1
      bhigh=bendsplits(bendtyp)

      xa = int(r)
      left = r-xa

c     select correct bin
c     multiply by num_int_bend - number of intervals per A
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
      endif

      rem=r-bend(blow, bendtyp)
      rem=rem*dble(num_int_bend(bendtyp))

      tabulated_bend=rem*benddiff(blow,bendtyp)
      tabulated_bend=tabulated_bend+tabbend(blow,bendtyp)
c      write(2,*) 'tabulated_bend ', tabulated_bend

      return
      end

      subroutine lininter_vdW(r, tabulated_vdW, typi, typj)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  calculates nonbonding van der Waals potential using linear 
cccc  interpolation between two points
cccc  requires a file (fort.43)
cccc  fort.43: number of tabulated potentials, potential number,
cccc  number of points per angstrom, tabulated potential
cccc  (repeat last three parts for each additional potential)
cccc  for unlike interactions, list 1-2 and 2-1
cccc  separate potentials with 1000
cccc  make sure potential does not go up to infinity!
cccc  bead type numbers should be defined in suijtab, but it doesn't
cccc  matter what they are (the parameters aren't used)
cccc  KM 12/03/08
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'tabulated.inc'
      include 'conver.inc'

      integer typi, typj, low, high, xa, bin, add1
      double precision r,  tabulated_vdW, left, rem

      low=1
      high=vdWsplits(typi, typj)

      xa = int(r)
      left = r-xa

c     select correct bin
c     multiply by num_int_vdW - number of intervals per angstrom
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
      endif

      rem=r-rvdW(low, typi, typj)
      rem=rem*dble(num_int_vdW(typi,typj))

      tabulated_vdW=rem*vdWdiff(low,typi, typj)
      tabulated_vdW=tabulated_vdW+tabvdW(low,typi, typj)

      return
      end


      subroutine lininter_elect(r, tabulated_elect, typi, typj)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  calculates electrostatic potential using linear interpolation
cccc  between two points
cccc  requires a file (fort.44)
cccc  fort.44: number of tabulated potentials, potential number,
cccc  number of points per angstrom, tabulated potential
cccc  (repeat last three parts for each additional potential)
cccc  for unlike interactions, list 1-2 and 2-1
cccc  separate potentials with 1000
cccc  KM 04/23/09
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'tabulated.inc'
      include 'conver.inc'

      integer typi, typj, low, high, xa, bin, add1
      double precision r,  tabulated_elect, left, rem

      low=1
      high=electsplits(typi, typj)

      xa = int(r)
      left = r-xa

c     select correct bin
c     multiply by num_int_nonbond - number of intervals per angstrom
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
      endif

      rem=r-relect(low, typi, typj)
      rem=rem*dble(num_int_elect(typi,typj))

      tabulated_elect=rem*electdiff(low,typi, typj)
      tabulated_elect=tabulated_elect+tabelect(low,typi, typj)

      return
      end

