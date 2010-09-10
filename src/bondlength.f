      subroutine bondlength( vibtype, requil, kvib, beta, length, vvib )

c bondlength
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

c     **************************************************************
c     *** computes the bond length for a given vibration type    ***
c     *** vibtype is the vibration type                          ***
c     *** requil is the equilibrium bond length                  ***
c     *** kvib is the force constant for the bond length         ***
c     *** length is the bond length returned by this subroutine  ***
c     *** beta is 1/kT                                           ***
c     *** vvib is the vibration energy for this bond             ***
c     *** M.G. Martin  2-4-98                                    ***
c     **************************************************************

      implicit none

      include 'tabulated.inc'

      integer::vibtype
      real(8)::length, bond, random, bf, vvib, beta, kvib
     &     , requil, tabulated_vib

      vvib = 0.0d0

      if ( kvib .gt. 0.1d0 ) then
c     --- random bond length from Boltzmann distribution ---
 107     bond = (0.2d0*random()+0.9d0)
         
c     --- correct for jacobian by dividing by the max^3
c     --- 1.331  = (1.1)^3, the equilibrium bond length cancels out
         bf = bond*bond*bond/(1.331d0)
         if ( random() .ge. bf ) goto 107
         bond = bond * requil
         
c     --- correct for the bond energy
         vvib = kvib * (bond-requil )**2
         bf = dexp ( -(vvib * beta) )
         if ( random() .ge. bf ) goto 107
         length = bond

c  --- tabulated bond stretching potential
c  --- added 12/02/08 by KM
      elseif (L_vib_table) then
c     --- random bond length from Boltzmann distribution ---
c     ---  +/- 25% of equilibrium bond length
 108     bond = (0.5d0*random()+0.75d0)
         
c     --- correct for jacobian by dividing by the max^3
c     --- 1.953125  = (1.25)^3, the equilibrium bond length cancels out
         bf = bond*bond*bond/(1.953125d0)
         if ( random() .ge. bf ) goto 108
         
         bond = bond * requil
         call lininter_vib(bond, tabulated_vib, vibtype)
         vvib=tabulated_vib

         bf = dexp ( -(vvib * beta) )
         
         if ( random() .ge. bf ) goto 108
         length = bond

      else
c     --- fixed bond length
         length = requil
      endif
      
      return

      end
