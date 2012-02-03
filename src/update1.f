      subroutine update1(nblock,ipos,value,acmove,ibox,jbox)

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
      integer nblock,ipos,ibox,jbox
      double precision acmove,dp,dn,value
      
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'neigh.inc'
      include 'neigh2.inc'
      include 'blkavg.inc'
      
      if (nblock.eq.1) then

c ---     first block

         
         dn = acmove
         dp = value
         if ( dn .lt. 0.5d0 ) then
            baver1(ipos,ibox,jbox,nblock) = 0.0d0
         else
            baver1(ipos,ibox,jbox,nblock) = dp / dn
         endif
         nccold1(ipos,ibox,jbox) = dn
         bccold1(ipos,ibox,jbox) = dp
         naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
         accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dp
      else      

c ---   other blocks

         dn = acmove - nccold1(ipos,ibox,jbox)
         dp = value - bccold1(ipos,ibox,jbox)
         if ( dn .lt. 0.5d0 ) then
            baver1(ipos,ibox,jbox,nblock) = 0.0d0
         else
            baver1(ipos,ibox,jbox,nblock) = dp / dn
         endif
         nccold1(ipos,ibox,jbox) = acmove
         bccold1(ipos,ibox,jbox) = value
         naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
         accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dp
      endif
      return
      end



















