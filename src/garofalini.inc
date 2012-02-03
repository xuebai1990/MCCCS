c GAROFALINI.INC
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
      integer pair_max
      parameter(pair_max=50000)

      double precision galpha(6),grho(6),gbeta(6),ga(6,3),
     &	 gb(6,3),gc(6,3),glambda(4),grij(4,2),ggamma(4,2),
     &   gtheta(4),grijsq(4,2)

      double precision dxij(pair_max),dyij(pair_max),
     &	 dzij(pair_max),dij(pair_max),dik(pair_max),dxik(pair_max)
     &   ,dyik(pair_max),dzik(pair_max)

      integer itr1(pair_max),itr2(pair_max),itr3(pair_max),
     &   ntr,tagged

      double precision v3garo,garofalini

      common /garo/ ga,gb,gc,ggamma,grij,galpha,grho,gbeta,
     &              glambda,gtheta,grijsq,tagged

      common /Ddist/ dij,dxij,dyij,dzij,dik,dxik,dyik,dzik
      common /triads/ itr1,itr2,itr3,ntr

      common /garoenrgy/ v3garo
