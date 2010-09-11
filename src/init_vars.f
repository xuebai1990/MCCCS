      subroutine init_vars
! init_vars
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
 
!    *******************************************************************
!    *** initializes some variables that otherwise cause errors      ***
!    *** Matt McGrath, October 1, 2009                               ***
!    *******************************************************************
 
      implicit none

! *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'poten.inc'
      include 'conver.inc' 
      include 'external.inc'
      include 'zeolite.inc'
      include 'connect.inc'
      include 'ewaldsum.inc'
      include 'fepsi.inc'
      include 'inputdata.inc'
      include 'qqlist.inc'
      include 'clusterbias.inc'
      include 'neigh.inc'
      include 'cell.inc'
      include 'nsix.inc'
      include 'peboco.inc'     
      include 'ipswpar.inc'
      include 'eepar.inc'
      include 'cbmc.inc'
! KM 01/10 remove analysis
!      include 'gor.inc'
      include 'blkavg.inc'
! kea include for garofalini potential
      include 'garofalini.inc'
      include 'tabulated.inc'
! **********************************************************************
! local variables

      integer::ibend,ibox,i,imolty,iunit,iprop,jblock,itype


! **********************************************************************
      nmolty1=0
      leemove=.false.
      iratipsw=0
      acipsw=0.0d0
! KM remove analysis
!      nframe=0.0d0
      DO ibend=1,nvmax
         brben(ibend)=0.0d0
      end do
      DO itype=1,nntype
         lij(itype)=.FALSE.
      end do
      DO itype=1,nntype*nntype
         sig2ij(itype)=0.0d0
         epsij(itype)=0.0d0
      end do
      DO iprop=1,nprop
         DO ibox=1,nbxmax
            DO jblock=1,blockm
               baver(iprop,ibox,jblock)=0.0d0
            end do
         end do
      end do
      DO ibox=1,nbxmax
         rmvol(ibox)=0.0d0
         boxlx(ibox)=0.0d0
         boxly(ibox)=0.0d0
         boxlz(ibox)=0.0d0
         DO i=1,9
            hmat(ibox,i)=0.0d0
         end do
         lsolid(ibox)=.FALSE.
         lrect(ibox)=.FALSE.
      end do
      upnn=0.0d0

      DO imolty=1,ntmax
         llplace(imolty)=.FALSE.
         lwell(imolty)=.FALSE.
         nrig(imolty)=0
         lbias(imolty)=.false.
      end do

      DO iunit=1,numax
         lsave(iunit)=.false.
      end do

      return
      end 

