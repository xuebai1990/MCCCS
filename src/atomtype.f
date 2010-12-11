      function atomtype(ntype,atom)

      implicit none
      include 'zeopoten.inc'
      character(len=*)::atom
      integer::ntype,atomtype

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
