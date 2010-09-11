      subroutine setnn ( ibox )

! setnn 
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

!    *******************************************************************
!    ** generates the complete near neighbour bitmap.                 **
!    **                                                               **
!    ** principal variables:                                          **
!    **                                                               **
!    ** logic   lnn               the bitmap                          **
!    ** dble    disvec            the displacementvector              **
!    **                                                               **
!    ** usage:                                                        **
!    **                                                               **
!    ** the subroutine generates the complete near neighbour bitmap   **
!    ** and sets all displacementvectors to zero at the beginning of  **
!    ** the run.                                                      **
!    *******************************************************************
 
      implicit logical(l)
      implicit character*40(c)
      implicit integer(f)
      implicit double precision(a,b,d,e,o-z)

      include 'control.inc'
      include 'coord.inc'
      include 'neigh2.inc'
      include 'neigh.inc'
      include 'system.inc' 


! -----------------------------------------------------------------------------
 
      rcnnsq = rcutnn(ibox)**2

      if ( lpbc ) call setpbc (ibox)

! --------------------------------------------------------------------
 
! --- set logical::map to .false. ---
      do 10 i = 1, nchain
         do 9 j = 1, nchain
            lnn(i,j) = .false.
 9       continue
 10   continue

! --------------------------------------------------------------------
 
! --- loop over all chains i 
      do 100 i = 1, nchain
 
         if ( nboxi(i) .ne. ibox ) go to 100

! --- loop over all beads ii of chain i 
         do 99 ii = 1, nunit(i)

            rxui = rxu(i,ii)
            ryui = ryu(i,ii)
            rzui = rzu(i,ii)

! --- loop over all chains j
            do 98 j = 1, nchain
               
               if ( ( i .eq. j ) .or. 
     &              ( nboxi(i) .ne. nboxi(j) ) ) go to 98

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
 
 97            continue
 98         continue
 99      continue
 100  continue
 
! -----------------------------------------------------------------------------
 
! *** set displacementvectors to zero ***
 
      do 1100 i = 1, nchain
 
         do 1099 j = 1, 3
            disvec(1,i,j) = 0.0d0
            disvec(2,i,j) = 0.0d0
!099     continue
!100  continue
 
!      write(iou,*) '@@@ control setnn @@@'
!      do 2000 i = 1, nchain
!         write(iou,*) 'lnn', i, '   ', (lnn(i,j),j=1,nchain)
!2000  continue
 
      write(iou,*)
      write(iou,*) 'neighbour table created by SETNN'

      return
      end

