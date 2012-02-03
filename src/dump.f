      subroutine dump  

c dump
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

C -----------------------------------------------------------------
C subroutine dump
C dumps the final configuration before stopping the program (Neeraj).
C -----------------------------------------------------------------

      implicit none

      include 'control.inc'
      include 'coord.inc'      
      include 'poten.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'cbmc.inc'
      include 'conver.inc'
      include 'inpar.inc' 
      include 'external.inc'
      include 'zeolite.inc'
      include 'inputdata.inc'
      include 'blkavg.inc'
      include 'bnbsma.inc'
      include 'swtcmove.inc'
      include 'ewaldsum.inc'
      include 'neigh.inc'
      include 'clusterbias.inc'
      include 'cell.inc'

      integer i,j,im,imolty,ibox

      open (unit=8, file="final-config",status="unknown")
      write(8,*) tmcc
      if ( tmcc .gt. 0 ) then
         write(8,*) Armtrax, Armtray, Armtraz 
         do im=1,nbox
            do imolty=1,nmolty
               write(8,*) rmtrax(imolty,im), rmtray(imolty,im)
     &              , rmtraz(imolty,im)
               write(8,*) rmrotx(imolty,im), rmroty(imolty,im)
     &              , rmrotz(imolty,im)
            enddo
         enddo
         do im=1, nbox
            write (8,*) (rmflcq(i,im),i=1,nmolty)
         enddo
         if ( lgibbs .or. lgrand .or. lnpt ) then
            write(8,*) (rmvol(ibox),ibox=1,nbox)
            do ibox = 1,nbox
               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  write(8,*) (rmhmat(ibox,i),i=1,9)
                  write(8,*) (hmat(ibox,i),i=1,9)
               else
                  write(8,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
               endif
            enddo
         endif
      endif
      write(8,*) nchain
      write(8,*) nmolty
      write(8,*) (nunit(i),i=1,nmolty)
      write(8,*) (moltyp(i),i=1,nchain)
      write(8,*) (nboxi(i),i=1,nchain)
      do i = 1, nmolty
         if ( lexpand(i) ) write(8,*) eetype(i)
      enddo
      do i = 1, nmolty
         if ( lexpand(i) ) write(8,*) rmexpc(i)
      enddo
      do  i = 1, nchain
         imolty = moltyp(i)
         do j = 1, nunit(imolty)
            write(8,*) rxu(i,j), ryu(i,j), rzu(i,j), qqu(i,j)
         enddo
      enddo
      close (8)
      return
      end
