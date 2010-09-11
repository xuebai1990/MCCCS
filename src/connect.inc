! CONNECT.INC
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

      logical::linclu,lqinclu,lainclu,lexclu,wschvib,wschben,wschtor
      integer::invib,itvib,ijvib,inben,itben,ijben2,ijben3
     &   ,intor,ittor,ijtor2,ijtor3,ijtor4
     &   ,wsched
      real(8)::brvib,brvibk,brvibmin,brvibmax,brben,brbenk
      
      common /connev/ brvib(nvmax),brvibk(nvmax),brvibmin(60)
     &   ,brvibmax(60)
     &   ,itvib(ntmax,numax,6),ijvib(ntmax,numax,6)
     &   ,invib(ntmax,numax)
      common /conneb/ brben(nvmax),brbenk(nvmax)
     &   ,itben(ntmax,numax,12),ijben2(ntmax,numax,12)
     &   ,ijben3(ntmax,numax,12)
     &   ,inben(ntmax,numax)
      common /connet/ ittor(ntmax,numax,12),ijtor2(ntmax,numax,12)
     &   ,ijtor3(ntmax,numax,12),ijtor4(ntmax,numax,12)
     &   ,intor(ntmax,numax)
      common /conneg/ wsched(numax),wschvib(numax,6)
     &    ,wschben(numax,12),wschtor(numax,12)
      common /includ/ linclu(ntmax,numax,numax)
     &    ,lqinclu(ntmax,numax,numax),lainclu(ntmax,numax,numax)
      common /exclud/ lexclu(ntmax,numax,ntmax,numax)




