      subroutine init_vars
c init_vars
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
 
c    *******************************************************************
c    *** initializes some variables that otherwise cause errors      ***
c    *** Matt McGrath, October 1, 2009                               ***
c    *******************************************************************
 
      implicit none

c *** common blocks ***
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
c KM remove analysis
c      include 'gor.inc'
      include 'blkavg.inc'
c kea include for garofalini potential
      include 'garofalini.inc'
      include 'tabulated.inc'
c **********************************************************************
c local variables

      integer ibend,ibox,i,imolty,iunit,iprop,jblock,itype


c **********************************************************************
      nmolty1=0
      leemove=.false.
      iratipsw=0
      acipsw=0.0d0
c KM remove analysis
c      nframe=0.0d0
      DO ibend=1,nvmax
         brben(ibend)=0.0d0
      ENDDO
      DO itype=1,nntype
         lij(itype)=.FALSE.
      ENDDO
      DO itype=1,nntype*nntype
         sig2ij(itype)=0.0d0
         epsij(itype)=0.0d0
      ENDDO
      DO iprop=1,nprop
         DO ibox=1,nbxmax
            DO jblock=1,blockm
               baver(iprop,ibox,jblock)=0.0d0
            ENDDO
         ENDDO
      ENDDO
      DO ibox=1,nbxmax
         rmvol(ibox)=0.0d0
         boxlx(ibox)=0.0d0
         boxly(ibox)=0.0d0
         boxlz(ibox)=0.0d0
         DO i=1,9
            hmat(ibox,i)=0.0d0
         ENDDO
         lsolid(ibox)=.FALSE.
         lrect(ibox)=.FALSE.
      ENDDO
      upnn=0.0d0

      DO imolty=1,ntmax
         llplace(imolty)=.FALSE.
         lwell(imolty)=.FALSE.
         nrig(imolty)=0
         lbias(imolty)=.false.
      ENDDO

      DO iunit=1,numax
         lsave(iunit)=.false.
      ENDDO

      return
      end 

