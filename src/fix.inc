! FIX.INC
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

      integer::  iend,ipast,endnum,pastnum,fclose,fcount,iwbef
     &   ,ibef,befnum,wbefnum,nextnum,inext

      logical::lcrank

      real(8)::xx,yy,zz,xvec,yvec,zvec,distij,vtgtr,vtbend
     &   ,bsum_tor,vtvib

      dimension iend(numax),ipast(numax,numax)
     &   ,pastnum(numax),fclose(numax,numax)
     &   ,fcount(numax),iwbef(numax),ibef(numax,numax)
     &   ,befnum(numax),xx(numax),yy(numax),zz(numax)
     &   ,xvec(numax,numax),yvec(numax,numax)
     &   ,zvec(numax,numax),distij(numax,numax)
     &   ,nextnum(numax),inext(numax,numax),vtvib(nchmax)
     &   ,vtgtr(nchmax),vtbend(nchmax),bsum_tor(nchmax)

      common /fecba/ iend,ipast,endnum,pastnum,fclose,fcount
     &   ,befnum,wbefnum,iwbef,ibef,nextnum
     &   ,inext,lcrank
   
      common /fecbb/ xx,yy,zz,xvec,yvec,zvec,distij,vtgtr,vtbend
     &   ,bsum_tor,vtvib

