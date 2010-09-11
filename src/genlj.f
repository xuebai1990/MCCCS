      function genlj (rijsq,sr2,epsilon2)

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

!  ********************************************************************
!  ***  calculates the energy of the Generalized Lennard Jones
!       potential                                               ***
!  ***      parameters defined  poten.inc           ***
!  *******************************************************************
!   Also you should make sure that the rcut you choose
!  is .gt. sigma*(2**(2/n0)) (because currently tail corrections
!  are evaluated with this assumption which is rarely incorrect.)
!  ********************************************************************

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

      
! In reduced units      
!      xij=(1.0d0/dsqrt(sr2))
!      
!      if ( (xij) .le. 2.0d0**(2.0d0/n0) ) then
!         genlj = 4.0d0 * (((1/xij)**n0)-((1/xij)**(n0/2.0d0)))
!      else
!         genlj =(( (2.0d0**((4.0d0*n1/n0)))*(1/xij)**(2.0d0*n1))-
!     & (2.0d0**((2.0d0*(n1/n0))+1.0d0)*((1/xij)**(n1))))
!      end if
!      vinter=vinter+e
       

!      not usually used in MC    
!      if (lshift) then
!         genlj = genlj - shiftgenlj()
!      end if
      
      return
      end

      
