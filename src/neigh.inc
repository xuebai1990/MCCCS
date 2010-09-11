! NEIGH.INC
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

      logical::lnn(1,1),lneighbor
      integer::maxneigh
      parameter (maxneigh = 20)
      integer::neighbor(maxneigh,nmax),neigh_cnt(nmax),neigh_icnt,
     &        neighi(maxneigh),neighboro(maxneigh,nmax),neigh_o(nmax)
      real(8)::ndij(maxneigh,nmax),nxij(maxneigh,nmax),
     &        nyij(maxneigh,nmax),nxijo(maxneigh,nmax),
     &        nzij(maxneigh,nmax),ndiji(nmax),nxiji(maxneigh),
     &        nyiji(maxneigh),nziji(maxneigh),ndijo(maxneigh,nmax),
     &        nyijo(maxneigh,nmax),nzijo(maxneigh,nmax)
!      common /nn2/ lnn(1,1)
!	replace the 1s with nmax to use neighbor list 
      common /nn2/ lnn,lneighbor
      common /nn4/ neighbor
      common /nn3/ neigh_cnt,neigh_icnt,neighi
      common /nn5/ ndij,nxij,nyij,nzij
      common /nn6/ ndiji,nxiji,nyiji,nziji
      common /nn7/ ndijo,nxijo,nyijo,nzijo,neigh_o,neighboro
