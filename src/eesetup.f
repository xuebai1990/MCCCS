      subroutine eesetup
!
! sets up EE. contains some EE stuff, see eemove.f for more details
!
      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_runtime,only:err_exit
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'eepar.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'system.inc'
!$$$      include 'nsix.inc'
!$$$      include 'poten.inc'

      integer(KIND=normal_int)::i,m,j,ntii,ntij,ntjj,ntjjs,ii,jj,ntijs ,imolty,isv,cnt

! --- initialize a few things

      leemove = .false.
      if ((pmexpc1.gt.1.0d-6).and.(.not.lexpee)) then
         write(io_output,*) 'pmexp nonzero but no lexpee?'
         call err_exit('')
      else if ((pmexpc1.lt.1.0d-6).and.lexpee) then
         write(io_output,*) 'pmexp zero but lexpee?'
         call err_exit('')
      end if

! --- read necessary stuff

! --- moltyp (of fort.4) on which EE is performed
      read(44,*)
      read(44,*) imolty

! --- number of actual types of molecules (i.e. nmolty minus the
! --- types that identify intermediate states.
      read(44,*)
      read(44,*) nmolty1

! --- the number for final state (first one is 1)
      read(44,*)
      read(44,*) fmstate

      if (fmstate.lt.3) call err_exit('EE when no intermediate state')

! --- weight (psi) associated with each state
      read(44,*)
      read(44,*) (psi(i), i = 1, fmstate)

! --- the two states between which 'swap' is. ensure that sstate1 is
! --- sstate2-1
      read(44,*)
      read(44,*) sstate1, sstate2
      if (sstate1.ne.(sstate2-1)) call err_exit('choose sstates in order')

! --- once an ee move is performed, the prob that it will be
! --- ee_index_swap move (keep is quite low)
      read(44,*)
      read(44,*) eeratio

! --- read the starting mstate
      read(44,*)
      read(44,*) mstate

! --- check with temtyp
      cnt = 0
      do i = nmolty1, nmolty
         if (temtyp(i).gt.0) then
            if (temtyp(i).ne.1) call err_exit('ee must be on one molecule only')
            isv = i
            cnt = cnt+1
         end if
      end do
      if (cnt.gt.1) call err_exit('only one state should be present in ee')
      if ((nmolty1+mstate-1).ne.isv) call err_exit('initial mstate inconsistent with temtyp')
      if ((mstate.eq.1).or.(mstate.eq.fmstate)) lmstate = .true.

! --- setup rminee for each unit. for fully grown units (same as in
! --- the full molecules), rminee is to be set to rmin (read from
! --- fort.4). for the partially grown units scale rmin by equating
! --- the 12th power potential values of the partially grown beads
! --- with the 12th power of the equivalent full beads. this choice
! --- is more or less arbitrary - but is consistent.

      do i = 1, nmolty1-1
         do ii = 1, nunit(i)
            ntii = ntype(i,ii)
            do j = 1, nmolty1-1
               do jj = 1, nunit(j)
                  ntjj = ntype(j,jj)
                  if (lexpsix.or.lmmff) then
                     ntij = (ntii+ntjj)/2
                  else if (lninesix) then
                     ntij = (ntii-1)*nxatom+ntjj
                  else if (lgenlj) then
                     ntij = (ntii-1)*nntype+ntjj
                  else
                     ntij = (ntii-1)*nntype+ntjj
                  end if
                  rminee(ntij) = rmin
!	write(io_output,*) i,ii,j,jj,rminee(ntij)
               end do
            end do
         end do
      end do
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
                  else if (lninesix) then
                     ntij = (ntii-1)*nxatom+ntjj
                     ntijs = (ntii-1)*nxatom+ntjjs
                  else if (lgenlj) then
                     ntij = (ntii-1)*nntype+ntjj
                     ntijs = (ntii-1)*nntype+ntjjs
                  else
                     ntij = (ntii-1)*nntype+ntjj
                     ntijs = (ntii-1)*nntype+ntjjs
                  end if
                  if (epsij(ntijs).ge.1.0d-6.and.sig2ij(ntijs).ge. 1.0d-6) then
                     rminee(ntij) = (epsij(ntij)/epsij(ntijs))** (1.0d0/12.0d0)*dsqrt(sig2ij(ntij)/sig2ij(ntijs))* rmin
                  else if ((dabs(qelect(ntii)*qelect(ntjj))) .ge.1.0d-6) then
                     rminee(ntij) = rmin
                  else
                     rminee(ntij) = 0.0d0
                  end if
!	write(io_output,*) i,ii,j,jj,rminee(ntij)
!	write(io_output,*) 'nt', ntii,ntjj,ntij,ntjjs,ntijs
!	write(io_output,*) 'eps', sig2ij(ntij),sig2ij(ntijs),epsij(ntij),
!     &              epsij(ntijs),qelect(ntii),qelect(ntjj)
               end do
            end do
         end do
      end do
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
                  else if (lninesix) then
                     ntij = (ntii-1)*nxatom+ntjj
                     ntijs = (ntjjs-1)*nxatom+ntjj
                  else if (lgenlj) then
                     ntij = (ntii-1)*nntype+ntjj
                     ntijs = (ntii-1)*nntype+ntjjs
                  else
                     ntij = (ntii-1)*nntype+ntjj
                     ntijs = (ntjjs-1)*nntype+ntjj
                  end if
                  if (epsij(ntijs).ge.1.0d-6.and.sig2ij(ntijs).ge. 1.0d-6) then
                     rminee(ntij) = (epsij(ntij)/epsij(ntijs))** (1.0d0/12.0d0)*dsqrt(sig2ij(ntij)/sig2ij(ntijs))* rmin
                  else if ((dabs(qelect(ntii)*qelect(ntjj))) .ge.1.0d-6) then
                     rminee(ntij) = rmin
                  else
                     rminee(ntij) = 0.0d0
                  end if
!	write(io_output,*) i,ii,j,jj,rminee(ntij)
!	write(io_output,*) 'nt', ntii,ntjj,ntij,ntjjs,ntijs
!	write(io_output,*) 'eps', sig2ij(ntij),sig2ij(ntijs),epsij(ntij),
!     &              epsij(ntijs)
               end do
            end do
         end do
      end do
      do i = nmolty1, nmolty
         do ii = 1, nunit(i)
            ntii = ntype(i,ii)
            do jj = 1, nunit(i)
               ntjj = ntype(i,jj)
               ntjjs = ntype(nmolty,jj)
               if (lexpsix.or.lmmff) then
                  ntij = (ntii+ntjj)/2
                  ntijs = (ntii+ntjjs)/2
               else if (lninesix) then
                  ntij = (ntii-1)*nxatom+ntjj
                  ntijs = (ntii-1)*nxatom+ntjjs
               else if (lgenlj) then
                  ntij = (ntii-1)*nntype+ntjj
                  ntijs = (ntii-1)*nntype+ntjjs
               else
                  ntij = (ntii-1)*nntype+ntjj
                  ntijs = (ntii-1)*nntype+ntjjs
               end if
               if (epsij(ntijs).ge.1.0d-6.and.sig2ij(ntijs).ge. 1.0d-6) then
                  rminee(ntij) = (epsij(ntij)/epsij(ntijs))** (1.0d0/12.0d0)*dsqrt(sig2ij(ntij)/sig2ij(ntijs))* rmin
               else if ((dabs(qelect(ntii)*qelect(ntjj))) .ge.1.0d-6) then
                  rminee(ntij) = rmin
               else
                  rminee(ntij) = 0.0d0
               end if
!	write(io_output,*) i,ii,i,jj,rminee(ntij)
            end do
         end do
      end do

!	write(io_output,*) 'enumerate'
!	do i = 1, nmolty
!	do ii = 1, nunit(i)
!	ntii = ntype(i,ii)
!	do j = 1, nmolty
!	do jj = 1, nunit(j)
!	ntjj = ntype(j,jj)
!                  ntij = (ntii-1)*nntype+ntjj
!	write(io_output,'(4(i4,1x),3(1x,e17.8))') i,ii,j,jj,rminee(ntij),epsij(ntij),sig2ij(ntij)
!	end do
!	end do
!	end do
!	end do
!	call err_exit('')

! --- associate moltyp with mstate

      do m = 1, fmstate
         ee_moltyp(m) = nmolty1+m-1
      end do
!      ee_moltyp(fmstate) = imolty
      do i = 1, nunit(imolty)
         do m = 1, fmstate
            ee_qqu(i,m) = qelect(ntype(nmolty1+m-1,i))
!	write(io_output,*) i,m,ee_qqu(i,m)
         end do
!         ee_qqu(i,fmstate) = qelect(ntype(imolty,i))
!	write(io_output,*) i,1,ee_qqu(i,1)
!	write(io_output,*) i,fmstate,ee_qqu(i,fmstate)
      end do

! --- associate a box with each state (convention: box 2 with states
! --- 1 to sstate1, and box 1 with sstate2 to fmstate)

      do m = 1, fmstate
         box_state(m) = 1
      end do
!      do m = 1, sstate1
!         box_state(m) = 2
!      end do
!      do m = sstate2, fmstate
!         box_state(m) = 1
!      end do

! --- the underlying matix of the markov chain (nonsymmetric if one of
! --- the state is an end state)

      do m = 2, fmstate-1
         um_markov(m,m+1) = 0.5d0
         um_markov(m,m-1) = 0.5d0
      end do
      um_markov(1,2) = 1.0d0
      um_markov(fmstate,fmstate-1) = 1.0d0

! --- pick a random chain at m=1 (i.e, in boxstate 1 ) to start off things
! --- if chain not present in boxstate 1, start with m = 6. for brute
! --- force method, there is always a unique tagged one (the tag doesn't
! --- change)

       eepointp = 1

!      if (dble(ncmt(box_state(1),imolty)).gt.0) then
!         eepointp = idint(dble(ncmt(box_state(1),imolty))*random())+1
!         mstate = 1
!      else if (dble(ncmt(box_state(2),imolty)).gt.0) then
!         eepointp = idint(dble(ncmt(box_state(2),imolty))*random())+1
!         mstate = fmstate
!      else
!         write(io_output,*)'the type is in neither box, imolty:',imolty
!         call err_exit('')
!      end if
!      lmstate = .true.
!	write(io_output,*) 'starting point', eepointp, mstate

! --- probability accumulators

!	write(io_output,*) 'prob check start'
      do m = 1, fmstate
         ee_prob(m) = 0
!	write(io_output,*) m, ee_prob(m)
      end do
!	write(io_output,*) 'prob check end'

      return
      end
