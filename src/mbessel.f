      function mbessel(z,nu)

! mbessel
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

        real(8)::z,nu
        real(8)::mbessel
        real(8)::pi
        parameter (pi = 3.14159265359d0)
        
!       mbessel = sqrt(pi/(2.0d0*z))*exp(-z)*
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))

! -- simple form     
        mbessel = sqrt(pi/(2.0d0*z))*exp(-z)
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))
        end function mbessel

