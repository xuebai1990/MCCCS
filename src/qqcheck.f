      subroutine qqcheck(i,ibox,rxuu1,ryuu1,rzuu1)

c qqcheck
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

c     **************************************************************
c     ***  sets up the coulom() array for use with the group based *
c     ***  cutoff for charged molecules     Marcus Martin          *
c     **************************************************************

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'poten.inc'
      include 'system.inc'
      include 'qqlist.inc'

      integer i,ibox,jmolty,j

      double precision rcutsq,rxuu1,ryuu1,rzuu1,rxuij,ryuij,rzuij
     &     ,rijsq

      rcutsq = rcutchg(ibox)*rcutchg(ibox)

      if ( lpbc ) call setpbc(ibox)

      do j = 1,nchain
         lcoulom(j) = .false.
         jmolty = moltyp(j)

         if ( (nboxi(j) .eq. ibox) .and. (i .ne. j) ) then

c ---  check for the group based qq cutoff
            if ( lelect(jmolty) ) then
               rxuij = rxuu1 - rxu(j,1)
               ryuij = ryuu1 - ryu(j,1)
               rzuij = rzuu1 - rzu(j,1)

c --- minimum image the coulombic bead pair separations ***
               if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

               if ( (rijsq .lt. rcutsq) .or. lchgall) then
                  lcoulom(j) = .true.
               endif
            endif

         endif

      enddo

      return
      end
