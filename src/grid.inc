c GRID.INC
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

      integer::maxtab,ngrx,ngry,ngrz,maxl,maxp
      real(8)::egrid,eps=1d-5
      real(8)::factx,facty,factz,dgrx,dgry,dgrz,xzz,yzz,zzz
c - parameters currently set low since zeolites not used
      parameter (maxtab=10,maxl=20,maxp=10)
      dimension egrid(0:maxtab)
     +    ,xzz(-maxp:maxl+maxp),yzz(-maxp:maxl+maxp),zzz(-10:maxl+maxp)
      common/gridd/factx,facty,factz,dgrx,dgry,dgrz,
     +             xzz,yzz,zzz,egrid
     +            ,ngrx,ngry,ngrz
