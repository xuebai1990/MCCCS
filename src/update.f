      subroutine update(nblock,ipos,ibox,value,acmove)

c update
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

c
c *** this subroutine updates the block averages
c
      implicit none
      integer nblock,ipos,ibox
      double precision acmove,dp,dn,value
      
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'neigh.inc'
      include 'neigh2.inc'
      include 'blkavg.inc'
      
      if (nblock.eq.1) then

c ---     first block

!         write(2,*) 'acmove', acmove       
  
         dn = acmove
         dp = value
         if ( dn .lt. 0.5d0 ) then
            baver(ipos,ibox,nblock) = 0.0d0
         else
            baver(ipos,ibox,nblock) = dp / dn
         endif
         nccold(ipos,ibox) = dn
         bccold(ipos,ibox) = dp
         naccu(ipos,ibox) = naccu(ipos,ibox) + dn
         accum(ipos,ibox) = accum(ipos,ibox) + dp
      else      

c ---   other blocks
!         write(2,*) 'acmove', acmove
         dn = acmove - nccold(ipos,ibox)
         dp = value - bccold(ipos,ibox)
         if ( dn .lt. 0.5d0 ) then
            baver(ipos,ibox,nblock) = 0.0d0
         else
            baver(ipos,ibox,nblock) = dp / dn
         endif
         nccold(ipos,ibox) = acmove
         bccold(ipos,ibox) = value
         naccu(ipos,ibox) = naccu(ipos,ibox) + dn
         accum(ipos,ibox) = accum(ipos,ibox) + dp
      endif
      return
      end



















