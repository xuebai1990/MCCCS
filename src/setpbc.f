      subroutine setpbc (ibox)

c setpbc
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
      include 'control.inc'
      include 'system.inc'
      include 'peboco.inc'

      integer ibox

c ----------------------------------------------------------------

      if ( lpbcx ) then
         bx = boxlx(ibox)
         if ( lfold ) then
            hbx = 0.5d0 * bx
         else
            bxi = 1.0d0 / bx
         endif
      endif

      if ( lpbcy ) then
         by = boxly(ibox)
         if ( lfold ) then
            hby = 0.5d0 * by
         else
            byi = 1.0d0 / by
         endif
      endif

      if ( lpbcz ) then
         bz = boxlz(ibox)
         if ( lfold ) then
            hbz = 0.5d0 * bz
         else
            bzi = 1.0d0 / bz
         endif
      endif

      return
      end


