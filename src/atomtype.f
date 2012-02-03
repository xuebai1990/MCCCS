      subroutine atomtype(ntype,atom,id)

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
      character atom*4
      integer ntype,id
      character at1*1,at2*2,at3*3

C     Assigning id's to atomtypes

      at1 = atom 
      at2 = atom 
      at3 = atom 

c      if (ntype.ne.2) stop '** atomtype: ntype ne 2 not allowed **'

C List of conversions

      if (at1.eq."O".or.at1.eq."o") then
         id =4
      else 
         stop "** atomtype: unknown atomtype **"
      endif

      return
      end

