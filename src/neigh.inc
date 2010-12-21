! NEIGH.INC

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
