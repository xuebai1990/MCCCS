      function garofalini(rijsq,ntij,qa,qb,aa,bb)

! garofalini
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
!     ***  calculates the energy using the garofalini (SiO2/H2O) ***
!     ***  exp-6 potential modified using articles in suijtab    ***
!     ***  parameters are defined in suijtab.f  KE ANDERSON      ***
!     **************************************************************
      implicit none

      include 'control.inc'
      include 'garofalini.inc'
      include 'poten.inc'

      real(8)::rijsq,rij,hterm,coul,erfunc,qa,qb
      real(8)::hterma
      integer::ntij,ntii,ntjj,aa,bb,i

!      write(6,*) 'input',rijsq,ntij,qa,qb
      rij = dsqrt(rijsq)
      hterm = 0.0d0
      coul = 0.0d0
      garofalini = 0.0d0


!      write(6,*) aa,bb,ntij,qa,qb,gbeta(ntij),galpha(ntij),grho(ntij)
!     H term
      do i=1,3
         hterma =  
     &        ga(ntij,i)/(1+dexp(gb(ntij,i)*(rij-gc(ntij,i))))
!         write(6,*) i,hterma,' (',ntij,')'
         hterm = hterm + hterma
      end do

      coul = qa*qb*erfunc(rij/gbeta(ntij))/rij
!      write(6,*) 'erfunc',coul*rij
      coul = coul * qqfact

      garofalini = galpha(ntij)*dexp(-rij/grho(ntij))
     &     + hterm 
     &     + coul

!      write(6,*) 'i,j,v2',aa,bb,'H term:',hterm,' Coul term:'
!     &     ,coul,' the rest:',garofalini-hterm-coul,
!     &     ' Total:',garofalini
!      write(6,*) '                     ',rij

!      write(6,*) 'leaving garofalini'

      return
      end

      
