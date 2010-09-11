! ENSEMBLE.INC
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

      real(8)::vbox,wbox,vinterb,vtailb,vintrab,vvibb,vbendb
     &   ,vtgb, vextb, velectb,vflucqb,v3garob
!      real(8)::velectb_intra, velectb_inter
      common /vwbox/  vbox(nbxmax),wbox(nbxmax),vinterb(nbxmax)
     &   ,vtailb(nbxmax),vintrab(nbxmax),vvibb(nbxmax)
     &   ,vbendb(nbxmax),vtgb(nbxmax),vextb(nbxmax)
     &   ,velectb(nbxmax),vflucqb(nbxmax),v3garob(nbxmax)
!      common /vabox/ velectb_intra(nbxmax),velectb_inter(nbxmax)
