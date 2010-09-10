c ROSEN.INC
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

      integer::growfrom,growprev,growlist,grownum
      dimension growfrom(numax),growprev(numax),grownum(numax)
      dimension growlist(numax,numax)

      real(8)::weight,  weiold, voldt,voldbb,voldtg,voldext
     +     ,voldintra,voldinter,voldbvib,voldelect,voldewald
     +     ,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib
     +     ,vnewinter,vnewelect,vnewewald,vneworient
     +     ,vnewinterr,vnewextr,vnewelectr,voldinterr
     +     ,voldextr,voldelectr,voldorient
      common /rosena/ growfrom,growprev,growlist,grownum

      common /rosenb/ weight,  weiold, voldt,voldbb,voldtg,voldext
     +     ,voldintra,voldinter,voldbvib,voldelect,voldewald
     +     ,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib
     +     ,vnewinter,vnewelect,vnewewald,voldelectr
     +     ,vnewinterr,vnewextr,vnewelectr,voldinterr
     +     ,voldextr,vneworient,voldorient
    

