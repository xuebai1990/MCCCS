C FIX.INC
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

      integer::  iend,ipast,endnum,pastnum,fclose,fcount,iwbef
     +   ,ibef,befnum,wbefnum,nextnum,inext

      logical::lcrank

      real(8)::xx,yy,zz,xvec,yvec,zvec,distij,vtgtr,vtbend
     +   ,bsum_tor,vtvib

      dimension iend(numax),ipast(numax,numax)
     +   ,pastnum(numax),fclose(numax,numax)
     +   ,fcount(numax),iwbef(numax),ibef(numax,numax)
     +   ,befnum(numax),xx(numax),yy(numax),zz(numax)
     +   ,xvec(numax,numax),yvec(numax,numax)
     +   ,zvec(numax,numax),distij(numax,numax)
     +   ,nextnum(numax),inext(numax,numax),vtvib(nchmax)
     +   ,vtgtr(nchmax),vtbend(nchmax),bsum_tor(nchmax)

      common /fecba/ iend,ipast,endnum,pastnum,fclose,fcount
     +   ,befnum,wbefnum,iwbef,ibef,nextnum
     +   ,inext,lcrank
   
      common /fecbb/ xx,yy,zz,xvec,yvec,zvec,distij,vtgtr,vtbend
     +   ,bsum_tor,vtvib

