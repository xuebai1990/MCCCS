      function mmff(rijsq,ntij)

c mmff
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

c   ******************************************************************
c   ***  calculates the energy using the Buffered 14-7 potential   ***
c   ***  parameters defined in suijtab.f          Bin Chen         *** 
c   ******************************************************************

      implicit none
      double precision rijsq,mmff,rs2,rs1,sr1,sr7,rs7
      integer ntij

      include 'control.inc'
      include 'merck.inc'

        rs2 = rijsq / (sigisq(ntij))
        rs1 = dsqrt(rs2)
        rs7 = rs1*rs2*rs2*rs2
        sr1 = 1.07d0/(rs1+0.07d0)
        sr7 = sr1**7.0d0
        mmff = sr7*(1.12d0/(rs7+0.12d0) - 2.0d0)
     +           *epsimmff(ntij)

      if (lshift) mmff = mmff-smmff(ntij)
      return
      end
