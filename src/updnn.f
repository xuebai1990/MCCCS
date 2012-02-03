      subroutine updnn( i, ibox )

c updnn
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

      implicit logical(l)
      implicit character*40(c)
      implicit integer(f)
      implicit double precision(a,b,d,e,o-z)
c *** common blocks ***

      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'neigh.inc'
      include 'neigh2.inc' 

C -----------------------------------------------------------------------------
 
      rcnnsq = rcutnn(ibox)**2

      if ( lpbc ) call setpbc (ibox)

c -----------------------------------------------------------------
 
c --- set i-part of logical map to .false. ---

      do 9 j = 1, nchain
         lnn(i,j) = .false.
 9    continue

c -----------------------------------------------------------------
 
c --- loop over all beads ii of chain i 
      do 99 ii = 1, nunit(i)

         rxui = rxu(i,ii)
         ryui = ryu(i,ii)
         rzui = rzu(i,ii)

c --- loop over all chains j
         do 98 j = 1, nchain
               
            if ( i .eq. j ) go to 98

            if ( lnn(i,j) ) go to 98

c --- loop over all beads jj of chain j 
            do 97 jj = 1, nunit(j)

               rxuij = rxui - rxu(j,jj)
               ryuij = ryui - ryu(j,jj)
               rzuij = rzui - rzu(j,jj)

c *** minimum image the pair separations ***
               if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

               if ( rijsq .le. rcnnsq ) then
                  lnn(i,j) = .true.
                  go to 98
               endif
 
 97         continue
 98      continue
 99   continue
 
C -----------------------------------------------------------------------------
 
c *** set displacementvectors to zero ***
 
      do 1099 j = 1, 3
         disvec(1,i,j) = 0.0d0
         disvec(2,i,j) = 0.0d0
 1099 continue
 
c      write(2,*) '@@@ control updnn @@@'
c      write(2,*) 'lnn', i, '   ', (lnn(i,j),j=1,nchain)

      return
      end



