! ROSEN.INC
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

      integer::growfrom,growprev,growlist,grownum
      dimension growfrom(numax),growprev(numax),grownum(numax)
      dimension growlist(numax,numax)

      real(8)::weight,  weiold, voldt,voldbb,voldtg,voldext
     &     ,voldintra,voldinter,voldbvib,voldelect,voldewald
     &     ,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib
     &     ,vnewinter,vnewelect,vnewewald,vneworient
     &     ,vnewinterr,vnewextr,vnewelectr,voldinterr
     &     ,voldextr,voldelectr,voldorient
      common /rosena/ growfrom,growprev,growlist,grownum

      common /rosenb/ weight,  weiold, voldt,voldbb,voldtg,voldext
     &     ,voldintra,voldinter,voldbvib,voldelect,voldewald
     &     ,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib
     &     ,vnewinter,vnewelect,vnewewald,voldelectr
     &     ,vnewinterr,vnewextr,vnewelectr,voldinterr
     &     ,voldextr,vneworient,voldorient
    

