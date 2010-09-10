      subroutine suzeo(rczeo,newztab)

c suzeo
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
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'external.inc'
      double precision rczeo
      logical newztab

c === load force field
 
c      call forcefield(rczeo)
 
c === load positions of zeolite atoms
 
      call zeocoord()
      zeorxi=1./zeorx
      zeoryi=1./zeory
      zeorzi=1./zeorz
c
c === tabulation of the zeolite potential
c
      if (lzgrid) then
         if (newztab) then
            call ztab
         else
c --        read table from file
            call rwztab(0)
         endif
         call ztest
      endif
      return

      end  
