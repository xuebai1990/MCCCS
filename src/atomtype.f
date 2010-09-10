      function atomtype(ntype,atom)

c atomtype 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
c John Stubbs, and Collin Wick and Ilja Siepmann  
c                     
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to 
c
c Free Software Foundation, Inc. 
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      character(len=*)::atom
      integer::ntype,atomtype

C     Assigning id's to atomtypes

c      if (ntype.ne.2) call cleanup('** atomtype: ntype ne 2 not allowed **')

C List of conversions
      atom=trim(atom)
      if (atom.eq."O") then
         atomtype = 190
      elseif (atom.eq."Si") then
         atomtype = 191
      else
         call cleanup('** atomtype: unknown atomtype **')
      endif

      return
      end
