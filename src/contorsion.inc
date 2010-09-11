! CONTORSION.INC
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

      integer::ntormax
      parameter (ntormax=700)
      real(8)::vtt0,vtt1,vtt2,vtt3,vtt4,vtt5,vtt6,vtt7,vtt8,vtt9
      common /contor/ vtt0(ntormax),vtt1(ntormax),vtt2(ntormax),
     &		      vtt3(ntormax),vtt4(ntormax),vtt5(ntormax),
     &		      vtt6(ntormax),vtt7(ntormax),vtt8(ntormax),
     &                vtt9(ntormax) 
      

