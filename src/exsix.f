      function exsix(rijsq,ntij)

! exsix
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

!     **************************************************************
!     ***  calculates the energy using the exp-6 potential       ***
!     ***  parameters defined in suijtab.f  M.G. Martin          ***
!     **************************************************************
      implicit none
      real(8)::rijsq,rij,exsix
      integer::ntij

      include 'control.inc'
      include 'expsix.inc'

      rij=dsqrt(rijsq)
      exsix = aexsix(ntij)/(rijsq*rijsq*rijsq)
     &     + bexsix(ntij)*dexp(cexsix(ntij)*rij)
      if (lshift) exsix = exsix-sexsix(ntij)
      return
      end

      
