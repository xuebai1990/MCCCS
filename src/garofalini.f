      function garofalini(rijsq,ntij,qa,qb,aa,bb)

c garofalini
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
c     ***  calculates the energy using the garofalini (SiO2/H2O) ***
c     ***  exp-6 potential modified using articles in suijtab    ***
c     ***  parameters are defined in suijtab.f  KE ANDERSON      ***
c     **************************************************************
      implicit none

      include 'control.inc'
      include 'garofalini.inc'
      include 'poten.inc'

      double precision rijsq,rij,hterm,coul,erfunc,qa,qb
      double precision hterma
      integer ntij,ntii,ntjj,aa,bb,i

c      write(6,*) 'input',rijsq,ntij,qa,qb
      rij = dsqrt(rijsq)
      hterm = 0.0d0
      coul = 0.0d0
      garofalini = 0.0d0


c      write(6,*) aa,bb,ntij,qa,qb,gbeta(ntij),galpha(ntij),grho(ntij)
c     H term
      do i=1,3
         hterma =  
     &        ga(ntij,i)/(1+dexp(gb(ntij,i)*(rij-gc(ntij,i))))
c         write(6,*) i,hterma,' (',ntij,')'
         hterm = hterm + hterma
      enddo

      coul = qa*qb*erfunc(rij/gbeta(ntij))/rij
c      write(6,*) 'erfunc',coul*rij
      coul = coul * qqfact

      garofalini = galpha(ntij)*dexp(-rij/grho(ntij))
     &     + hterm 
     &     + coul

c      write(6,*) 'i,j,v2',aa,bb,'H term:',hterm,' Coul term:'
c     &     ,coul,' the rest:',garofalini-hterm-coul,
c     &     ' Total:',garofalini
c      write(6,*) '                     ',rij

c      write(6,*) 'leaving garofalini'

      return
      end

      
