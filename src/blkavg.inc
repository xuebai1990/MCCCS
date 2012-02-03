c BLKAVG.INC
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

      integer nener,nprop,blockm
      double precision naccu,nccold,accum,bccold,aver,baver
     &                ,acsvol,acnvol,acshmat,acnhmat
      parameter (nener=18,nprop=4+nener+(4*ntmax),blockm=100)
      dimension naccu(nprop,nbxmax),nccold(nprop,nbxmax)
      dimension accum(nprop,nbxmax),bccold(nprop,nbxmax)
      dimension aver(nprop,nbxmax), baver(nprop,nbxmax,blockm)
      dimension acsvol(nbxmax),acnvol(nbxmax)
      dimension acshmat(nbxmax,9),acnhmat(nbxmax,9)

      common /bokavg/ baver,naccu,nccold,accum,bccold,aver,acsvol
     &                ,acnvol,acshmat,acnhmat

