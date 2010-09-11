! BLKAVG.INC
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

      integer::nener,nprop,blockm
      real(8)::naccu,nccold,accum,bccold,aver,baver
     &     ,acsvol,acnvol,acshmat,acnhmat
! ** Neeraj adding for solubility parameter and heat of vaporization

      real(8)::naccu1,nccold1,accum1,bccold1,aver1,baver1
      integer::nprop1
      parameter (nener=18,nprop=4+nener+(4*ntmax)+3,blockm=100)
      parameter (nprop1=11)
      dimension naccu(nprop,nbxmax),nccold(nprop,nbxmax)
      dimension accum(nprop,nbxmax),bccold(nprop,nbxmax)
      dimension aver(nprop,nbxmax), baver(nprop,nbxmax,blockm)
      dimension acsvol(nbxmax),acnvol(nbxmax)
      dimension acshmat(nbxmax,9),acnhmat(nbxmax,9)
      
      dimension naccu1(nprop1,nbxmax,nbxmax), 
     &  nccold1(nprop1,nbxmax,nbxmax), aver1(nprop1,nbxmax,nbxmax)
      dimension accum1(nprop1,nbxmax,nbxmax),
     &  bccold1(nprop1,nbxmax,nbxmax), 
     &  baver1(nprop1,nbxmax,nbxmax,blockm)

      double precision, dimension (ntmax,nbxmax) :: acntray,acntrax,
     &     acntraz,acnrotx,acnroty,acnrotz,acstrax,acstray,acstraz,
     &     acsrotx,acsroty,acsrotz

      common /bokavg/ baver,naccu,nccold,accum,bccold,aver,acsvol
     &                ,acnvol,acshmat,acnhmat
      common /blkavg/ baver1,naccu1,nccold1,accum1,bccold1,aver1

      common /trarotavg/ acntrax,acntray,acntraz,acnrotx,acnroty,
     &     acnrotz,acstrax,acstray,acstraz,acsrotx,acsroty,
     &     acsrotz
