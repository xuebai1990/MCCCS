      function genlj (rijsq,sr2,epsilon2)

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

c  ********************************************************************
c  ***  calculates the energy of the Generalized Lennard Jones
c       potential                                               ***
c  ***      parameters defined  poten.inc           ***
c  *******************************************************************
c   Also you should make sure that the rcut you choose
c  is .gt. sigma*(2**(2/n0)) (because currently tail corrections
c  are evaluated with this assumption which is rarely incorrect.)
c  ********************************************************************

      implicit none
      real(8)::rijsq,rij,srij,sr2,epsilon2,genlj


      include 'control.inc'
      include 'poten.inc'

     
      rij=(dsqrt(rijsq))
      srij=dsqrt(sr2)

      if ( (rij) .le.(rij*srij)*2.0d0**(2.0d0/n0) ) then
         genlj = 4.0d0*epsilon2*(((srij)**n0)-((srij)**(n0/2.0d0)))
      else
         genlj =epsilon2*(((2.0d0**((4.0d0*n1/n0)))*((srij)**
     &  (2.0d0*n1)))-((2.0d0**((2.0d0*(n1/n0))+1.0d0))*((srij)**(n1))))
      end if

      
c In reduced units      
c      xij=(1.0d0/dsqrt(sr2))
c      
c      if ( (xij) .le. 2.0d0**(2.0d0/n0) ) then
c         genlj = 4.0d0 * (((1/xij)**n0)-((1/xij)**(n0/2.0d0)))
c      else
c         genlj =(( (2.0d0**((4.0d0*n1/n0)))*(1/xij)**(2.0d0*n1))-
c     & (2.0d0**((2.0d0*(n1/n0))+1.0d0)*((1/xij)**(n1))))
c      end if
c      vinter=vinter+e
       

!      not usually used in MC    
!      if (lshift) then
!         genlj = genlj - shiftgenlj()
!      endif
      
      return
      end

      
