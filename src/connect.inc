c CONNECT.INC
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

      logical linclu,lqinclu,lainclu,lexclu,wschvib,wschben,wschtor
      integer invib,itvib,ijvib,inben,itben,ijben2,ijben3
     +   ,intor,ittor,ijtor2,ijtor3,ijtor4
     +   ,wsched
      double precision brvib,brvibk,brvibmin,brvibmax,brben,brbenk

      common /connev/ brvib(nvmax),brvibk(nvmax),brvibmin(60)
     &   ,brvibmax(60)
     +   ,itvib(ntmax,numax,6),ijvib(ntmax,numax,6)
     +   ,invib(ntmax,numax)
      common /conneb/ brben(nvmax),brbenk(nvmax)
     +   ,itben(ntmax,numax,10),ijben2(ntmax,numax,10)
     +   ,ijben3(ntmax,numax,10)
     +   ,inben(ntmax,numax)
      common /connet/ ittor(ntmax,numax,10),ijtor2(ntmax,numax,10)
     +   ,ijtor3(ntmax,numax,10),ijtor4(ntmax,numax,10)
     +   ,intor(ntmax,numax)
      common /conneg/ wsched(numax),wschvib(numax,6)
     +    ,wschben(numax,10),wschtor(numax,10)
      common /includ/ linclu(ntmax,numax,numax)
     &    ,lqinclu(ntmax,numax,numax),lainclu(ntmax,numax,numax)
      common /exclud/ lexclu(ntmax,numax,ntmax,numax)




