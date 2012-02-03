      subroutine mimage ( rxuij,ryuij,rzuij,ibox )

c mimage
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

      implicit none
      include 'peboco.inc'
      include 'control.inc'
      include 'system.inc'
      include 'cell.inc'

      integer ibox
      double precision rxuij,ryuij,rzuij,sx,sy,sz

c ----------------------------------------------------------------

      if (lsolid(ibox) .and. .not. lrect(ibox)) then
c$$$
c$$$c ***********************
c$$$c Sophia's (Bin's?) code
c$$$c *** non-rectangular box
c$$$         hsx = 0.5d0*hmat(1)
c$$$         hsy = 0.5d0*hmat(5)
c$$$         hsz = 0.5d0*hmat(9)
c$$$         sx = rxuij*hmati(1)+ryuij*hmati(4)+rzuij*hmati(7)
c$$$         sy = rxuij*hmati(2)+ryuij*hmati(5)+rzuij*hmati(8)
c$$$         sz = rxuij*hmati(3)+ryuij*hmati(6)+rzuij*hmati(9)
c$$$
c$$$c         print*,rxuij,ryuij,rzuij
c$$$
c$$$         if ( rzuij .gt. hsz ) then
c$$$            rzuij=rzuij-hmat(9)
c$$$            sz=sz-1d0
c$$$            if ( rzuij .gt. hsz ) then
c$$$               rzuij=rzuij-hmat(9)
c$$$               sz=sz-1d0
c$$$            endif
c$$$         elseif ( rzuij .lt. -hsz ) then
c$$$            rzuij=rzuij+hmat(9)
c$$$            sz=sz+1d0
c$$$            if ( rzuij .lt. -hsz ) then
c$$$               rzuij=rzuij+hmat(9)
c$$$               sz=sz+1d0
c$$$            endif
c$$$         endif
c$$$
c$$$         ryuij=sy*hmat(5)+sz*hmat(8)
c$$$         if ( ryuij .gt. hsy ) then
c$$$            ryuij=ryuij-hmat(5)
c$$$            sy=sy-1d0
c$$$            if ( ryuij .gt. hsy ) then
c$$$               ryuij=ryuij-hmat(5)
c$$$               sy=sy-1d0
c$$$            endif
c$$$         elseif ( ryuij .lt. -hsy ) then
c$$$            ryuij=ryuij+hmat(5)
c$$$            sy=sy+1d0
c$$$            if ( ryuij .lt. -hsy ) then
c$$$               ryuij=ryuij+hmat(5)
c$$$               sy=sy+1d0
c$$$            endif
c$$$         endif
c$$$
c$$$         rxuij=sx*hmat(1)+sy*hmat(4)+sz*hmat(7)
c$$$         if ( rxuij .gt. hsx ) then
c$$$            rxuij=rxuij-hmat(1)
c$$$            if ( rxuij .gt. hsx ) then
c$$$               rxuij=rxuij-hmat(1)
c$$$            endif
c$$$         elseif ( rxuij .lt. -hsx ) then
c$$$            rxuij=rxuij+hmat(1)
c$$$            if ( rxuij .lt. -hsx ) then
c$$$               rxuij=rxuij+hmat(1)
c$$$            endif
c$$$         endif
c$$$
c$$$         return


c ******************
c *** Collin's Code:
c *** non-rectangular box
         sx = rxuij*hmati(ibox,1)+ryuij*hmati(ibox,4)
     &        +rzuij*hmati(ibox,7)
         sy = rxuij*hmati(ibox,2)+ryuij*hmati(ibox,5)
     &        +rzuij*hmati(ibox,8)
         sz = rxuij*hmati(ibox,3)+ryuij*hmati(ibox,6)
     &        +rzuij*hmati(ibox,9)

         if ( sx .gt. 0.5d0 ) then
            sx = sx-1d0
         elseif ( sx .lt. -0.5d0 ) then
            sx = sx+1d0
         endif
         if ( sy .gt. 0.5d0 ) then
            sy = sy-1d0
         elseif ( sy .lt. -0.5d0 ) then
            sy = sy+1d0
         endif
         if ( sz .gt. 0.5d0 ) then
            sz = sz-1d0
         elseif ( sz .lt. -0.5d0 ) then
            sz = sz+1d0
         endif
c         sx = sx-nint(sx)
c         sy = sy-nint(sy)
c         sz = sz-nint(sz)
         rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
         ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)+sz*hmat(ibox,8)
         rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)+sz*hmat(ibox,9)

         return

      else

         if ( lpbcx ) then
            if ( lfold ) then
               if ( rxuij .gt. hbx ) then
                  rxuij=rxuij-bx
               else
                  if (rxuij.lt.-hbx) rxuij=rxuij+bx
               endif
            else
c            rxuij = rxuij - bx*anint(rxuij*bxi)
               rxuij = rxuij - bx*dint(rxuij*bxi+dsign(0.5d0,rxuij))
            endif
         endif
         
         if ( lpbcy ) then
            if ( lfold ) then
               if ( ryuij .gt. hby ) then
                  ryuij=ryuij-by
               else
                  if (ryuij.lt.-hby) ryuij=ryuij+by
               endif
            else
c            ryuij  = ryuij - by*anint(ryuij*byi)
               ryuij = ryuij - by*dint(ryuij*byi+dsign(0.5d0,ryuij))
            endif
         endif

         if ( lpbcz ) then
            if ( lfold ) then
               if (rzuij.gt.hbz) then
                  rzuij=rzuij-bz
               else
                  if (rzuij.lt.-hbz) rzuij=rzuij+bz
               endif
            else
c            rzuij  = rzuij - bz*anint(rzuij*bzi)
               rzuij = rzuij - bz*dint(rzuij*bzi+dsign(0.5d0,rzuij))
            endif
         endif

         return
         
      endif

      write(6,*) rxuij,ryuij,rzuij,ibox

      end
