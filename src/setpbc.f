      subroutine setpbc (ibox)

! setpbc
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

      implicit none
      include 'control.inc'
      include 'system.inc'
      include 'peboco.inc'

      integer::ibox

! ----------------------------------------------------------------

      if ( lpbcx ) then
         bx = boxlx(ibox)
         if ( lfold ) then
            hbx = 0.5d0 * bx
         else
            bxi = 1.0d0 / bx
         end if
      end if

      if ( lpbcy ) then
         by = boxly(ibox)
         if ( lfold ) then
            hby = 0.5d0 * by
         else
            byi = 1.0d0 / by
         end if
      end if

      if ( lpbcz ) then
         bz = boxlz(ibox)
         if ( lfold ) then
            hbz = 0.5d0 * bz
         else
            bzi = 1.0d0 / bz
         end if
      end if

      return
      end


