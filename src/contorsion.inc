c CONTORSION.INC
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

      integer ntormax
      parameter (ntormax=700)
      double precision vtt0,vtt1,vtt2,vtt3,vtt4,vtt5,vtt6,vtt7,vtt8,vtt9
      common /contor/ vtt0(ntormax),vtt1(ntormax),vtt2(ntormax),
     +		      vtt3(ntormax),vtt4(ntormax),vtt5(ntormax),
     +		      vtt6(ntormax),vtt7(ntormax),vtt8(ntormax),
     +                vtt9(ntormax) 
      

