c NEIGH.INC
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
c      common /nn2/ lnn(1,1)
c	replace the 1s with nmax to use neighbor list 
      common /nn2/ lnn,lneighbor
      common /nn4/ neighbor
      common /nn3/ neigh_cnt,neigh_icnt,neighi
      common /nn5/ ndij,nxij,nyij,nzij
      common /nn6/ ndiji,nxiji,nyiji,nziji
      common /nn7/ ndijo,nxijo,nyijo,nzijo,neigh_o,neighboro
