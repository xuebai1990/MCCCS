      function coruz(iunit,rho,ibox)

c coruz
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

c ***************************************************
c *** tail-corrections in energy with the zeolite ***
c ***************************************************

      implicit none
      include 'control.inc'
      include 'conver.inc'
      include 'system.inc'
      include 'zeopoten.inc'
      double precision coruz,eps,rci3,rho
      integer iunit,ibox
c --- note works only for alkanes!!!
      if (iunit.ne.1) then
        eps=zeps(1,4)+(iunit-2)*zeps(2,4)+zeps(3,4)
      else
        eps=zeps(1,4)
      endif
      rci3=zsig2(1,4)**(3.d0/2.d0)/rcut(ibox)**3 
      coruz=8.*onepi*eps*rho*(rci3*rci3*rci3/9.-rci3/3.)
      return
      end
