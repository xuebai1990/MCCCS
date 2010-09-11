      subroutine update1(nblock,ipos,value,acmove,ibox,jbox)

! update
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

!
! *** this subroutine updates the block averages
!
      implicit none
      integer::nblock,ipos,ibox,jbox
      real(8)::acmove,dp,dn,value
      
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'neigh.inc'
      include 'neigh2.inc'
      include 'blkavg.inc'
      
      if (nblock.eq.1) then

! ---     first block

         
         dn = acmove
         dp = value
         if ( dn .lt. 0.5d0 ) then
            baver1(ipos,ibox,jbox,nblock) = 0.0d0
         else
            baver1(ipos,ibox,jbox,nblock) = dp / dn
         end if
         nccold1(ipos,ibox,jbox) = dn
         bccold1(ipos,ibox,jbox) = dp
         naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
         accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dp
      else      

! ---   other blocks

         dn = acmove - nccold1(ipos,ibox,jbox)
         dp = value - bccold1(ipos,ibox,jbox)
         if ( dn .lt. 0.5d0 ) then
            baver1(ipos,ibox,jbox,nblock) = 0.0d0
         else
            baver1(ipos,ibox,jbox,nblock) = dp / dn
         end if
         nccold1(ipos,ibox,jbox) = acmove
         bccold1(ipos,ibox,jbox) = value
         naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
         accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dp
      end if
      return
      end



















