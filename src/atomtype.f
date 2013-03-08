      function atomtype(ntype,atom)

! atomtype 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
! John Stubbs, and Collin Wick and Ilja Siepmann  
!                     
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to 
!
! Free Software Foundation, Inc. 
! 59 Temple Place - Suite 330
! Boston, MA  02111-1307, USA.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include "zeopoten.inc"
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
