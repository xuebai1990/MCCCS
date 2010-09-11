      subroutine writepdb(nmol,natom,ibox)

! writepdb
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
     &                 rxu(imm,ia),ryu(imm,ia),rzu(imm,ia)
         end do
      end do

   99 format(a6,i5,1x,a4,14x,3f8.3)
  100 format(/,' WRITING GUEST ATOMS TO FILE system.pdb')

      return
      end

