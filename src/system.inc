! SYSTEM.INC
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

      logical::lsolid,lrect,ltwice,lideal
      integer::icell,nicell,iucell,solcount
      real(8)::boxlx, boxly, boxlz, rcut, rmin
     &   ,rcutnn, softcut,softlog,rcutin
     &   ,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc
      dimension icell(nmax),nicell(cmax),iucell(cmax,cmaxa)
     &   ,solcount(nbxmax,ntmax)
      dimension avsolinter(nbxmax,ntmax),avsolintra(nbxmax,ntmax)
     &   ,avsolbend(nbxmax,ntmax),avsoltor(nbxmax,ntmax)
     &   ,avsolelc(nbxmax,ntmax)

      common /boxpha/ lsolid(nbxmax),lrect(nbxmax),ltwice(nbxmax)
      common /boxdim/ boxlx(nbxmax),boxly(nbxmax),boxlz(nbxmax),
     &                 rcut(nbxmax),rcutnn(nbxmax),lideal(nbxmax)

      common /solav/ avsolinter,avsolintra,avsolbend,avsoltor,avsolelc
      common /input1/ rmin, softcut,softlog
     &  ,rcutin,solcount

      common /linkc/ icell,nicell,iucell
