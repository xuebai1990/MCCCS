      function ninesix(rijsq,ntij)

c ninesix
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

c     ***************************************************
c     ***  calculates the energy of the 9-6 potential ***
c     ***  parameters defined in suijtab.f  JMS       ***
c     ***************************************************

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
      endif

      return
      end

      
