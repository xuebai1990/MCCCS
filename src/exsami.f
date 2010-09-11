       function exsami ( z, ntj )

! exsami
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

!    **********************************************************
!    ** calculates the SAMI's external energy for a bead.    **
!    **********************************************************
 
      implicit none

! *** common blocks ***
      include 'external.inc'

      real(8)::exsami, z
      integer::ntj

! --------------------------------------------------------------------

      if ( ntj .eq. 1 ) then
! --- HEADgroup potential ---
         if ( z .le. alpha2 ) then
            exsami = 0.0d0
         else
            exsami = beta2 / ( 1.0d0 + ( (z/alpha2) - 1.0d0 )**tau2 )
         end if
      else
! --- TAILgroup potential ---
         if ( z .ge. alpha1 ) then
            exsami = 0.0d0
         else
            exsami = beta1 / ( 1.0d0 + ( 1.0d0 - (z/alpha1) )**tau1 )
         end if
      end if

      return

! ----------------------------------------------------------------------------

      end
