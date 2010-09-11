       function ljpsur ( rijsq, ntij )
 
! ljpsur
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

!    *********************************************************
!    ** calculates energy for a polymeric surfactant bead.  **
!    *********************************************************
 
      implicit none

! *** common blocks ***
      include 'external.inc'

      real(8)::ljpsur, rijsq, sr, sr6
      integer::ntij

! --------------------------------------------------------------------
! AT PRESENT: all sigma = 1.0
!             epsilon   = 1.0
! --------------------------------------------------------------------

      sr = 1.0d0 / rijsq
      sr6 = sr**3

      if ( ntij .eq. 1 ) then
! *** nonpolar-nonpolar interaction ( LJ interaction )
         ljpsur = sr6*sr6 - sr6
      else
! *** polar-polar or polar-nonpolar interaction ( repulsive LJ interaction )
         if ( rijsq .le. 1.259921d0 ) then
            ljpsur = sr6*sr6 - sr6 + 0.25d0
         else
            ljpsur = 0.0d0
         end if
      end if

      return

! ----------------------------------------------------------------------------

      end
