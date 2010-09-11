      subroutine bondlength( vibtype, requil, kvib, beta, length, vvib )

! bondlength
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
!     *** computes the bond length for a given vibration type    ***
!     *** vibtype is the vibration type                          ***
!     *** requil is the equilibrium bond length                  ***
!     *** kvib is the force constant for the bond length         ***
!     *** length is the bond length returned by this subroutine  ***
!     *** beta is 1/kT                                           ***
!     *** vvib is the vibration energy for this bond             ***
!     *** M.G. Martin  2-4-98                                    ***
!     **************************************************************

      implicit none

      include 'tabulated.inc'

      integer::vibtype
      real(8)::length, bond, random, bf, vvib, beta, kvib
     &     , requil, tabulated_vib

      vvib = 0.0d0

      if ( kvib .gt. 0.1d0 ) then
!     --- random bond length from Boltzmann distribution ---
 107     bond = (0.2d0*random()+0.9d0)
         
!     --- correct for jacobian by dividing by the max^3
!     --- 1.331  = (1.1)^3, the equilibrium bond length cancels out
         bf = bond*bond*bond/(1.331d0)
         if ( random() .ge. bf ) goto 107
         bond = bond * requil
         
!     --- correct for the bond energy
         vvib = kvib * (bond-requil )**2
         bf = dexp ( -(vvib * beta) )
         if ( random() .ge. bf ) goto 107
         length = bond

!  --- tabulated bond stretching potential
!  --- added 12/02/08 by KM
      elseif (L_vib_table) then
!     --- random bond length from Boltzmann distribution ---
!     ---  +/- 25% of equilibrium bond length
 108     bond = (0.5d0*random()+0.75d0)
         
!     --- correct for jacobian by dividing by the max^3
!     --- 1.953125  = (1.25)^3, the equilibrium bond length cancels out
         bf = bond*bond*bond/(1.953125d0)
         if ( random() .ge. bf ) goto 108
         
         bond = bond * requil
         call lininter_vib(bond, tabulated_vib, vibtype)
         vvib=tabulated_vib

         bf = dexp ( -(vvib * beta) )
         
         if ( random() .ge. bf ) goto 108
         length = bond

      else
!     --- fixed bond length
         length = requil
      end if
      
      return

      end
