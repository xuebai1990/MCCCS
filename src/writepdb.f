      subroutine writepdb(nmol,natom,ibox)

c writepdb
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

      include 'zeolite.inc'
      include 'control.inc'
      include 'coord.inc'

      integer::im,ia,nmol,natom,iguest,imm,ibox
      character::atom*4

      open(unit=48, file='system.pdb', form='formatted')
      rewind(48)
      write(iou,*) 'most likely does not work as before 9-3-96'
      write(iou,100)

      atom = 'C   '
      iguest=0 
      do im=1,nmol
         do ia=1,natom
            imm=parbox(im,ibox,1)
            iguest=iguest+1
            write(48,99) "ATOM  ",iguest,atom,
     +                 rxu(imm,ia),ryu(imm,ia),rzu(imm,ia)
         enddo
      enddo

   99 format(a6,i5,1x,a4,14x,3f8.3)
  100 format(/,' WRITING GUEST ATOMS TO FILE system.pdb')

      return
      end

