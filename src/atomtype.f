      function atomtype(ntype,atom)

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
!$$$      include 'zeopoten.inc'
      character(LEN=*)::atom
      integer(KIND=int)::ntype,atomtype

!     Assigning id's to atomtypes

!      if (ntype.ne.2) call cleanup('** atomtype: ntype ne 2 not allowed **')

! List of conversions
      atom=trim(atom)
      if (atom.eq."  Si") then
         atomtype = ztype(1)
         znum(1)=znum(1)+1
      elseif (atom.eq."   O") then
         atomtype = ztype(2)
         znum(2)=znum(2)+1
      else
         print*,len(atom),atom,atom.eq."Si",atom.eq."O"
         call cleanup('** atomtype: unknown atomtype **')
      end if

      return
      end
