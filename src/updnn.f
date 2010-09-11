      subroutine updnn( i, ibox )

! updnn
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

      implicit logical(l)
      implicit character*40(c)
      implicit integer(f)
      implicit double precision(a,b,d,e,o-z)
! *** common blocks ***

      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'neigh.inc'
      include 'neigh2.inc' 

! -----------------------------------------------------------------------------
 
      rcnnsq = rcutnn(ibox)**2

      if ( lpbc ) call setpbc (ibox)

! -----------------------------------------------------------------
 
! --- set i-part of logical::map to .false. ---

      do 9 j = 1, nchain
         lnn(i,j) = .false.
 9    continue

! -----------------------------------------------------------------
 
! --- loop over all beads ii of chain i 
      do 99 ii = 1, nunit(i)

         rxui = rxu(i,ii)
         ryui = ryu(i,ii)
         rzui = rzu(i,ii)

! --- loop over all chains j
         do 98 j = 1, nchain
               
            if ( i .eq. j ) go to 98

            if ( lnn(i,j) ) go to 98

! --- loop over all beads jj of chain j 
            do 97 jj = 1, nunit(j)

               rxuij = rxui - rxu(j,jj)
               ryuij = ryui - ryu(j,jj)
               rzuij = rzui - rzu(j,jj)

! *** minimum image the pair separations ***
               if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

               if ( rijsq .le. rcnnsq ) then
                  lnn(i,j) = .true.
                  go to 98
               end if
 
 97         continue
 98      continue
 99   continue
 
! -----------------------------------------------------------------------------
 
! *** set displacementvectors to zero ***
 
      do 1099 j = 1, 3
         disvec(1,i,j) = 0.0d0
         disvec(2,i,j) = 0.0d0
 1099 continue
 
!      write(iou,*) '@@@ control updnn @@@'
!      write(iou,*) 'lnn', i, '   ', (lnn(i,j),j=1,nchain)

      return
      end



