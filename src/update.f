      subroutine update(nblock,ipos,ibox,value,acmove)

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
      integer::nblock,ipos,ibox
      real(8)::acmove,dp,dn,value
      
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'neigh.inc'
      include 'neigh2.inc'
      include 'blkavg.inc'
      
      if (nblock.eq.1) then

! ---     first block

!         write(iou,*) 'acmove', acmove       
  
         dn = acmove
         dp = value
         if ( dn .lt. 0.5d0 ) then
            baver(ipos,ibox,nblock) = 0.0d0
         else
            baver(ipos,ibox,nblock) = dp / dn
         end if
         nccold(ipos,ibox) = dn
         bccold(ipos,ibox) = dp
         naccu(ipos,ibox) = naccu(ipos,ibox) + dn
         accum(ipos,ibox) = accum(ipos,ibox) + dp
      else      

! ---   other blocks
!         write(iou,*) 'acmove', acmove
         dn = acmove - nccold(ipos,ibox)
         dp = value - bccold(ipos,ibox)
         if ( dn .lt. 0.5d0 ) then
            baver(ipos,ibox,nblock) = 0.0d0
         else
            baver(ipos,ibox,nblock) = dp / dn
         end if
         nccold(ipos,ibox) = acmove
         bccold(ipos,ibox) = value
         naccu(ipos,ibox) = naccu(ipos,ibox) + dn
         accum(ipos,ibox) = accum(ipos,ibox) + dp
      end if
      return
      end



















