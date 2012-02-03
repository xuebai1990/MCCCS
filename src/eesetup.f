      subroutine eesetup(qelect)
c
c sets up EE. contains some EE stuff, see eemove.f for more details
c
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'eepar.inc'
      include 'inputdata.inc'
      include 'system.inc'
      include 'nsix.inc'
      include 'poten.inc'

      logical ovrlap,ldum
      integer i,m,j,ntii,ntij,ntjj,ntjjs,ii,jj,ntijs,imolty,irem
     &       ,isv,cnt
      double precision qelect,random,vmstate,volde,vnewe,dum,vintra
     &  ,vinter,vext,velect,vewald
      dimension qelect(nntype)

c --- initialize a few things

      leemove = .false.
      if ((pmexpc1.gt.1.0d-6).and.(.not.lexpee)) then
         write(2,*) 'pmexp nonzero but no lexpee?'
         stop
      elseif ((pmexpc1.lt.1.0d-6).and.lexpee) then
         write(2,*) 'pmexp zero but lexpee?'
         stop
      endif

c --- read necessary stuff

c --- moltyp (of fort.4) on which EE is performed
      read(44,*)
      read(44,*) imolty

c --- number of actual types of molecules (i.e. nmolty minus the
c --- types that identify intermediate states.
      read(44,*)
      read(44,*) nmolty1

c --- the number for final state (first one is 1)
      read(44,*)
      read(44,*) fmstate

      if (fmstate.lt.3) stop 'EE when no intermediate state'

c --- weight (psi) associated with each state
      read(44,*)
      read(44,*) (psi(i), i = 1, fmstate)

c --- the two states between which 'swap' is. ensure that sstate1 is
c --- sstate2-1
      read(44,*)
      read(44,*) sstate1, sstate2
      if (sstate1.ne.(sstate2-1)) stop 'choose sstates in order'

c --- once an ee move is performed, the prob that it will be
c --- ee_index_swap move (keep is quite low)
      read(44,*)
      read(44,*) eeratio

c --- read the starting mstate
      read(44,*)
      read(44,*) mstate

c --- check with temtyp
      cnt = 0
      do i = nmolty1, nmolty
         if (temtyp(i).gt.0) then
            if (temtyp(i).ne.1) stop 'ee must be on one molecule only'
            isv = i
            cnt = cnt+1
         endif
      enddo
      if (cnt.gt.1) stop 'only one state should be present in ee'
      if ((nmolty1+mstate-1).ne.isv) stop
     &   'initial mstate inconsistent with temtyp'
      if ((mstate.eq.1).or.(mstate.eq.fmstate)) lmstate = .true.

c --- setup rminee for each unit. for fully grown units (same as in
c --- the full molecules), rminee is to be set to rmin (read from
c --- fort.4). for the partially grown units scale rmin by equating
c --- the 12th power potential values of the partially grown beads
c --- with the 12th power of the equivalent full beads. this choice
c --- is more or less arbitrary - but is consistent.

      do i = 1, nmolty1-1
         do ii = 1, nunit(i)
            ntii = ntype(i,ii)
            do j = 1, nmolty1-1
               do jj = 1, nunit(j)
                  ntjj = ntype(j,jj)
                  if (lexpsix.or.lmmff) then
                     ntij = (ntii+ntjj)/2
                  elseif (lninesix) then
                     ntij = (ntii-1)*nxatom+ntjj
                  else
                     ntij = (ntii-1)*nntype+ntjj
                  endif
                  rminee(ntij) = rmin
c	write(2,*) i,ii,j,jj,rminee(ntij)
               enddo
            enddo
         enddo
      enddo
      do i = 1, nmolty1-1
         do ii = 1, nunit(i)
            ntii = ntype(i,ii)
            do j = nmolty1, nmolty
               do jj = 1, nunit(j)
                  ntjj = ntype(j,jj)
                  ntjjs = ntype(nmolty,jj)
                  if (lexpsix.or.lmmff) then
                     ntij = (ntii+ntjj)/2
                     ntijs = (ntii+ntjjs)/2
                  elseif (lninesix) then
                     ntij = (ntii-1)*nxatom+ntjj
                     ntijs = (ntii-1)*nxatom+ntjjs
                  else
                     ntij = (ntii-1)*nntype+ntjj
                     ntijs = (ntii-1)*nntype+ntjjs
                  endif
                  if (epsij(ntijs).ge.1.0d-6.and.sig2ij(ntijs).ge.
     &                1.0d-6) then
                     rminee(ntij) = (epsij(ntij)/epsij(ntijs))**
     &                (1.0d0/12.0d0)*dsqrt(sig2ij(ntij)/sig2ij(ntijs))*
     &                rmin
                  elseif ((dabs(qelect(ntii)*qelect(ntjj)))
     &                   .ge.1.0d-6) then
                     rminee(ntij) = rmin
                  else
                     rminee(ntij) = 0.0d0
                  endif
c	write(2,*) i,ii,j,jj,rminee(ntij)
c	write(2,*) 'nt', ntii,ntjj,ntij,ntjjs,ntijs
c	write(2,*) 'eps', sig2ij(ntij),sig2ij(ntijs),epsij(ntij),
c     &              epsij(ntijs),qelect(ntii),qelect(ntjj)
               enddo
            enddo
         enddo
      enddo
      do i = nmolty1, nmolty
         do ii = 1, nunit(i)
            ntii = ntype(i,ii)
            ntjjs = ntype(nmolty,ii)
            do j = 1, nmolty1-1
               do jj = 1, nunit(j)
                  ntjj = ntype(j,jj)
                  if (lexpsix.or.lmmff) then
                     ntij = (ntii+ntjj)/2
                     ntijs = (ntjjs+ntjj)/2
                  elseif (lninesix) then
                     ntij = (ntii-1)*nxatom+ntjj
                     ntijs = (ntjjs-1)*nxatom+ntjj
                  else
                     ntij = (ntii-1)*nntype+ntjj
                     ntijs = (ntjjs-1)*nntype+ntjj
                  endif
                  if (epsij(ntijs).ge.1.0d-6.and.sig2ij(ntijs).ge.
     &                1.0d-6) then
                     rminee(ntij) = (epsij(ntij)/epsij(ntijs))**
     &                (1.0d0/12.0d0)*dsqrt(sig2ij(ntij)/sig2ij(ntijs))*
     &                rmin
                  elseif ((dabs(qelect(ntii)*qelect(ntjj)))
     &                   .ge.1.0d-6) then
                     rminee(ntij) = rmin
                  else
                     rminee(ntij) = 0.0d0
                  endif
c	write(2,*) i,ii,j,jj,rminee(ntij)
c	write(2,*) 'nt', ntii,ntjj,ntij,ntjjs,ntijs
c	write(2,*) 'eps', sig2ij(ntij),sig2ij(ntijs),epsij(ntij),
c     &              epsij(ntijs)
               enddo
            enddo
         enddo
      enddo
      do i = nmolty1, nmolty
         do ii = 1, nunit(i)
            ntii = ntype(i,ii)
            do jj = 1, nunit(i)
               ntjj = ntype(i,jj)
               ntjjs = ntype(nmolty,jj)
               if (lexpsix.or.lmmff) then
                  ntij = (ntii+ntjj)/2
                  ntijs = (ntii+ntjjs)/2
               elseif (lninesix) then
                  ntij = (ntii-1)*nxatom+ntjj
                  ntijs = (ntii-1)*nxatom+ntjjs
               else
                  ntij = (ntii-1)*nntype+ntjj
                  ntijs = (ntii-1)*nntype+ntjjs
               endif
               if (epsij(ntijs).ge.1.0d-6.and.sig2ij(ntijs).ge.
     &             1.0d-6) then
                  rminee(ntij) = (epsij(ntij)/epsij(ntijs))**
     &             (1.0d0/12.0d0)*dsqrt(sig2ij(ntij)/sig2ij(ntijs))*
     &             rmin
               elseif ((dabs(qelect(ntii)*qelect(ntjj)))
     &                   .ge.1.0d-6) then
                  rminee(ntij) = rmin
               else
                  rminee(ntij) = 0.0d0
               endif
c	write(2,*) i,ii,i,jj,rminee(ntij)
            enddo
         enddo
      enddo

c	write(2,*) 'enumerate'
c	do i = 1, nmolty
c	do ii = 1, nunit(i)
c	ntii = ntype(i,ii)
c	do j = 1, nmolty
c	do jj = 1, nunit(j)
c	ntjj = ntype(j,jj)
c                  ntij = (ntii-1)*nntype+ntjj
c	write(2,999) i,ii,j,jj,rminee(ntij),epsij(ntij),sig2ij(ntij)
c	enddo
c	enddo
c	enddo
c	enddo
c	stop

c --- associate moltyp with mstate

      do m = 1, fmstate
         ee_moltyp(m) = nmolty1+m-1
      enddo
c      ee_moltyp(fmstate) = imolty
      do i = 1, nunit(imolty)
         do m = 1, fmstate
            ee_qqu(i,m) = qelect(ntype(nmolty1+m-1,i))
c	write(2,*) i,m,ee_qqu(i,m)
         enddo
c         ee_qqu(i,fmstate) = qelect(ntype(imolty,i))
c	write(2,*) i,1,ee_qqu(i,1)
c	write(2,*) i,fmstate,ee_qqu(i,fmstate)
      enddo

c --- associate a box with each state (convention: box 2 with states
c --- 1 to sstate1, and box 1 with sstate2 to fmstate)

      do m = 1, fmstate
         box_state(m) = 1
      enddo
c      do m = 1, sstate1
c         box_state(m) = 2
c      enddo
c      do m = sstate2, fmstate
c         box_state(m) = 1
c      enddo

c --- the underlying matix of the markov chain (nonsymmetric if one of
c --- the state is an end state)

      do m = 2, fmstate-1
         um_markov(m,m+1) = 0.5d0
         um_markov(m,m-1) = 0.5d0
      enddo
      um_markov(1,2) = 1.0d0
      um_markov(fmstate,fmstate-1) = 1.0d0

c --- pick a random chain at m=1 (i.e, in boxstate 1 ) to start off things
c --- if chain not present in boxstate 1, start with m = 6. for brute
c --- force method, there is always a unique tagged one (the tag doesn't
c --- change)

       eepointp = 1

c      if (dble(ncmt(box_state(1),imolty)).gt.0) then
c         eepointp = idint(dble(ncmt(box_state(1),imolty))*random())+1
c         mstate = 1
c      elseif (dble(ncmt(box_state(2),imolty)).gt.0) then
c         eepointp = idint(dble(ncmt(box_state(2),imolty))*random())+1
c         mstate = fmstate
c      else
c         write(2,*)'the type is in neither box, imolty:',imolty
c         stop
c      endif
c      lmstate = .true.
c	write(2,*) 'starting point', eepointp, mstate

c --- probability accumulators

c	write(2,*) 'prob check start'
      do m = 1, fmstate
         ee_prob(m) = 0
c	write(2,*) m, ee_prob(m)
      enddo
c	write(2,*) 'prob check end'

 999	format(4(i4,1x),3(1x,e17.8))
      return
      end
