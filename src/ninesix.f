      function ninesix(rijsq,ntij)

! ninesix
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

!     ***************************************************
!     ***  calculates the energy of the 9-6 potential ***
!     ***  parameters defined in suijtab.f  JMS       ***
!     ***************************************************

      implicit none
      real(8)::rijsq,rij,ror,ninesix
      integer::ntij

      include 'control.inc'
      include 'nsix.inc'

      rij=dsqrt(rijsq)
      ror = rzero(ntij)/rij
      ninesix = 4.0d0*epsnx(ntij)*ror**6*(2.0d0*ror**3 - 3.0d0)
      if (lshift) then
         ninesix = ninesix - shiftnsix(ntij)
      end if

      return
      end

      
