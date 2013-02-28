MODULE moves_ee
  use util_random,only:random
  use util_runtime,only:err_exit
  use sim_system
  use energy_pairwise,only:energy
  use energy_kspace,only:ee_recip
  implicit none
  private
  save
  public::eesetup,eemove,ee_index_swap,expand,read_expand,init_ee,output_ee_stats

! EXPAND.INC
  real,allocatable,public::epsil(:,:,:),sigm(:,:,:)
  real,allocatable::qcharge(:,:,:),bnexpc(:,:),bsexpc(:,:),eta(:,:,:)
  integer,allocatable,public::numcoeff(:)

contains
!> \brief sets up EE. contains some EE stuff, see eemove.f for more details
  subroutine eesetup
    use energy_pairwise,only:type_2body

      integer::i,m,j,ntii,ntij,ntjj,ntjjs,ii,jj,ntijs ,imolty,isv,cnt

! initialize a few things

      leemove = .false.
      if ((pmexpc1.gt.1.0d-6).and.(.not.lexpee)) then
         write(io_output,*) 'pmexp nonzero but no lexpee?'
         call err_exit(__FILE__,__LINE__,'',myid+1)
      else if ((pmexpc1.lt.1.0d-6).and.lexpee) then
         write(io_output,*) 'pmexp zero but lexpee?'
         call err_exit(__FILE__,__LINE__,'',myid+1)
      end if

! read necessary stuff

! moltyp (of fort.4) on which EE is performed
      read(44,*)
      read(44,*) imolty

! number of actual types of molecules (i.e. nmolty minus the
! types that identify intermediate states.
      read(44,*)
      read(44,*) nmolty1

! the number for final state (first one is 1)
      read(44,*)
      read(44,*) fmstate

      if (fmstate.lt.3) call err_exit(__FILE__,__LINE__,'EE when no intermediate state',myid+1)

! weight (psi) associated with each state
      read(44,*)
      read(44,*) (psi(i), i = 1, fmstate)

! the two states between which 'swap' is. ensure that sstate1 is
! sstate2-1
      read(44,*)
      read(44,*) sstate1, sstate2
      if (sstate1.ne.(sstate2-1)) call err_exit(__FILE__,__LINE__,'choose sstates in order',myid+1)

! once an ee move is performed, the prob that it will be
! ee_index_swap move (keep is quite low)
      read(44,*)
      read(44,*) eeratio

! read the starting mstate
      read(44,*)
      read(44,*) mstate

! check with temtyp
      cnt = 0
      do i = nmolty1, nmolty
         if (temtyp(i).gt.0) then
            if (temtyp(i).ne.1) call err_exit(__FILE__,__LINE__,'ee must be on one molecule only',myid+1)
            isv = i
            cnt = cnt+1
         end if
      end do
      if (cnt.gt.1) call err_exit(__FILE__,__LINE__,'only one state should be present in ee',myid+1)
      if ((nmolty1+mstate-1).ne.isv) call err_exit(__FILE__,__LINE__,'initial mstate inconsistent with temtyp',myid+1)
      if ((mstate.eq.1).or.(mstate.eq.fmstate)) lmstate = .true.

! setup rminee for each unit. for fully grown units (same as in
! the full molecules), rminee is to be set to rmin (read from
! fort.4). for the partially grown units scale rmin by equating
! the 12th power potential values of the partially grown beads
! with the 12th power of the equivalent full beads. this choice
! is more or less arbitrary - but is consistent.

      do i = 1, nmolty1-1
         do ii = 1, nunit(i)
            ntii = ntype(i,ii)
            do j = 1, nmolty1-1
               do jj = 1, nunit(j)
                  ntjj = ntype(j,jj)
                  ntij=type_2body(ntii,ntjj)
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
                  ntij=type_2body(ntii,ntjj)
                  ntijs=type_2body(ntii,ntjjs)
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
                  ntij=type_2body(ntii,ntjj)
                  ntijs=type_2body(ntjjs,ntjj)
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
               ntij=type_2body(ntii,ntjj)
               ntijs=type_2body(ntii,ntjjs)
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
! ntij = (ntii-1)*nntype+ntjj
!	write(io_output,'(4(i4,1x),3(1x,e17.8))') i,ii,j,jj,rminee(ntij),epsij(ntij),sig2ij(ntij)
!	end do
!	end do
!	end do
!	end do
!	call err_exit(__FILE__,__LINE__,'',myid+1)

! associate moltyp with mstate

      do m = 1, fmstate
         ee_moltyp(m) = nmolty1+m-1
      end do
! ee_moltyp(fmstate) = imolty
      do i = 1, nunit(imolty)
         do m = 1, fmstate
            ee_qqu(i,m) = qelect(ntype(nmolty1+m-1,i))
!	write(io_output,*) i,m,ee_qqu(i,m)
         end do
! ee_qqu(i,fmstate) = qelect(ntype(imolty,i))
!	write(io_output,*) i,1,ee_qqu(i,1)
!	write(io_output,*) i,fmstate,ee_qqu(i,fmstate)
      end do

! associate a box with each state (convention: box 2 with states
! 1 to sstate1, and box 1 with sstate2 to fmstate)

      do m = 1, fmstate
         box_state(m) = 1
      end do
! do m = 1, sstate1
! box_state(m) = 2
! end do
! do m = sstate2, fmstate
! box_state(m) = 1
! end do

! the underlying matix of the markov chain (nonsymmetric if one of
! the state is an end state)

      do m = 2, fmstate-1
         um_markov(m,m+1) = 0.5d0
         um_markov(m,m-1) = 0.5d0
      end do
      um_markov(1,2) = 1.0d0
      um_markov(fmstate,fmstate-1) = 1.0d0

! pick a random chain at m=1 (i.e, in boxstate 1 ) to start off things
! if chain not present in boxstate 1, start with m = 6. for brute
! force method, there is always a unique tagged one (the tag doesn't
! change)

       eepointp = 1

! if (dble(ncmt(box_state(1),imolty)).gt.0) then
! eepointp = idint(dble(ncmt(box_state(1),imolty))*random())+1
! mstate = 1
! else if (dble(ncmt(box_state(2),imolty)).gt.0) then
! eepointp = idint(dble(ncmt(box_state(2),imolty))*random())+1
! mstate = fmstate
! else
! write(io_output,*)'the type is in neither box, imolty:',imolty
! call err_exit(__FILE__,__LINE__,'',myid+1)
! end if
! lmstate = .true.
!	write(io_output,*) 'starting point', eepointp, mstate

! probability accumulators

!	write(io_output,*) 'prob check start'
      do m = 1, fmstate
         ee_prob(m) = 0
!	write(io_output,*) m, ee_prob(m)
      end do
!	write(io_output,*) 'prob check end'

      return
  end subroutine eesetup

! peforms EE move. currently works when EE move is performed on
! only one type of molecule. put all states in fort.77 and add
! corresponding bead types in suijtab. can only do EE if the number
! of beads between states remain constant.
  subroutine eemove
      use transfer_swap,only:swap
      use sim_cell

      logical::ovrlap
      integer::i,ibox,iunit,imolty,imolty1,j,idummy(nmax)
      real::dum,vrecipn,vrecipo,vnew,vold,vintern ,vintero,vintran,vintrao,velectn,velecto,vewaldn,vewaldo,vextn ,vexto,deltv,deltvb,wdeltvb,vtailn,vtailo

! --------------------------------------------------------------------

! write(*,*) 'START EEMOVE'
! write(11,*) '1:',neigh_cnt(18)

      leemove = .true.
      leeacc = .false.

! choose nstate, given mstate

      if (mstate.eq.1) then
         nstate = mstate + 1
      else if (mstate.eq.fmstate) then
         nstate = mstate - 1
      else
         if (random().le.0.5d0) then
            nstate = mstate - 1
         else
            nstate = mstate + 1
         end if
      end if

!	write(io_output,*) 'typ', mstate, ee_moltyp(mstate),eepointp,
!     &   box_state(mstate), ncmt(box_state(mstate),ee_moltyp(mstate))

      if (ncmt(box_state(mstate),ee_moltyp(mstate)).eq.0) then
         write(io_output,*)'problem: mstate, but no molecule in mstate',mstate
         call err_exit(__FILE__,__LINE__,'',myid+1)
      end if

! type of move depending upon mstate and nstate. one type of move
! is 'swap', other is usual ee

      wee_ratio = dexp(psi(nstate)-psi(mstate))* um_markov(nstate,mstate)/um_markov(mstate,nstate)

!	write(io_output,*) 'mstate', mstate, ee_moltyp(mstate)
!	write(io_output,*) 'nstate', nstate, ee_moltyp(nstate)

      if ((mstate.eq.sstate1.and.nstate.eq.sstate2).or. (mstate.eq.sstate2.and.nstate.eq.sstate1)) then

         boxrem1 = box_state(mstate)
         boxins1 = box_state(nstate)

         eeirem = parbox(eepointp,boxrem1,ee_moltyp(mstate))
!	write(io_output,*) 'eeirem', eeirem, eepointp

         call swap

         if (.not.leeacc) goto 100

      else

! energy for the new state

         ibox = box_state(mstate)
         eeirem = parbox(eepointp,ibox,ee_moltyp(mstate))
!	write(io_output,*) 'eeirem1', eeirem, eepointp
         imolty = ee_moltyp(nstate)

         iunit = nunit(imolty)
         do i = 1, iunit
            rxuion(i,2) = rxu(eeirem,i)
            ryuion(i,2) = ryu(eeirem,i)
            rzuion(i,2) = rzu(eeirem,i)
            qquion(i,2) = ee_qqu(i,nstate)
         end do

         moltion(1) = imolty
         call energy(eeirem,imolty,vnew,vintran,vintern,vextn,velectn,vewaldn,2,ibox,1,iunit,.false.,ovrlap,.false.,dum,.false.,.false.,.false.)
         if (ovrlap) goto 100
         if (ltailc) then
! add tail corrections for the Lennard-Jones energy
            vtailn=ee_coru(ibox,imolty,2)
            vintern = vintern + vtailn
         end if

! energy for the old state

         imolty = ee_moltyp(mstate)
         do i = 1, iunit
            rxuion(i,1) = rxu(eeirem,i)
            ryuion(i,1) = ryu(eeirem,i)
            rzuion(i,1) = rzu(eeirem,i)
            qquion(i,1) = qqu(eeirem,i)
         end do
         moltion(2) = imolty
         call energy(eeirem,imolty,vold,vintrao,vintero,vexto,velecto,vewaldo,1,ibox,1,iunit,.false.,ovrlap,.false.,dum,.false.,.false.,.false.)
         if (ovrlap) then
            call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf eemove',myid+1)
         end if
         if (ltailc) then
! add tail corrections for the Lennard-Jones energy
            vtailo=ee_coru(ibox,imolty,1)
            vintero = vintero + vtailo
         end if

         if (lewald.and.(lelect(moltion(2)).or.lelect(moltion(1))))then
            call ee_recip(ibox,vrecipn,vrecipo,1)
            velectn = velectn+vrecipn+vewaldn
            velecto = velecto+vrecipo+vewaldo
            vnew = vnew + vrecipn
            vold = vold + vrecipo
         end if

! check for acceptance

         deltv = (vnew - vold)
         deltvb = beta*deltv
         wdeltvb = wee_ratio*dexp(-deltvb)

         if ((deltvb-dlog(wee_ratio)).le.0.0d0) then
! accept move
         else if (wdeltvb.gt.random()) then
! accept move
         else
! reject move
            goto 100
         end if

! update new

         imolty1 = ee_moltyp(nstate)
         parbox(ncmt(ibox,imolty1)+1,ibox,imolty1) = eeirem
!	write(io_output,*) 'update new', eepointp,ibox,imolty1,
!     &          parbox(eepointp,ibox,imolty1),
!     &          parbox(eepointp,ibox,imolty)

         parbox(eepointp,ibox,imolty) = parbox(ncmt(ibox,imolty),ibox,imolty)

!	write(io_output,*) 'update new1', imolty,parbox(eepointp,ibox,imolty),
!     &              ncmt(ibox,imolty)

         parbox(ncmt(ibox,imolty),ibox,imolty) = 0
         ncmt(ibox,imolty1) = ncmt(ibox,imolty1) + 1
         ncmt(ibox,imolty) = ncmt(ibox,imolty) - 1
         eepointp = ncmt(ibox,imolty1)

         vbox(ibox) = vbox(ibox) + deltv
         vinterb(ibox) = vinterb(ibox) + (vintern-vintero)
         vintrab(ibox) = vintrab(ibox) + (vintran-vintrao)
         vextb(ibox) = vextb(ibox) + (vextn-vexto)
         vtailb(ibox) = vtailb(ibox) + (vtailn-vtailo)
         velectb(ibox) = velectb(ibox) + (velectn-velecto)
!	write(io_output,*) vtailn,vtailo,vintern,vintero

! update reciprocal space term

         call ee_recip(ibox,dum,dum,2)

      end if

! update the present state

!	ovr = .true.

      temtyp(ee_moltyp(mstate)) = temtyp(ee_moltyp(mstate)) - 1
      temtyp(ee_moltyp(nstate)) = temtyp(ee_moltyp(nstate)) + 1
!	write(io_output,*) 'ttym', mstate,ee_moltyp(mstate),
!     &               temtyp(ee_moltyp(mstate))
!	write(io_output,*) 'ttyn', nstate,ee_moltyp(nstate),
!     &               temtyp(ee_moltyp(nstate))
      mstate = nstate

!	write(io_output,*) 'eemove eeirem', eeirem
      moltyp(eeirem) = ee_moltyp(nstate)
      do i = 1, nunit(ee_moltyp(nstate))
         qqu(eeirem,i) = ee_qqu(i,nstate)
!	write(io_output,*) 'charges', qqu(eeirem,i)
!	write(io_output,*) 'charges1'
      end do

      if (mstate.eq.1.or.mstate.eq.fmstate) then
         lmstate = .true.
      else
         lmstate = .false.
      end if

! update parall (since one molecule has changed its state/type)

      do  i = 1, nmolty
         idummy(i) = 0
      end do
      do j = 1, nchain
         i = moltyp(j)
         idummy(i) = idummy(i) + 1
         parall(i,idummy(i)) = j
      end do

!	write(io_output,*) 'accept', eepointp, mstate, nstate, eeirem, boxins1
!	write(io_output,*)'accept1',parbox(eepointp,boxins1,ee_moltyp(nstate))
!	write(io_output,*) '1line', leemove,lmstate,leeacc
!	write(io_output,*) '2line', fmstate,sstate1,sstate2,wee_ratio,eeratio
!	write(io_output,*) '3line',(ee_moltyp(i), i = 1, fmstate)
!	write(io_output,*) '4line',(box_state(i), i = 1, fmstate)
!	write(io_output,*) '5line',(psi(i), i = 1, fmstate)
!	write(io_output,*) '6line',(ee_qqu(1,i), i = 1, fmstate)
!	write(io_output,*) '7line',(um_markov(i,i+1),i = 1,fmstate-1)
!	write(io_output,*) '8line',(um_markov(i,i-1),i = 2,fmstate)

 100    continue
        leemove = .false.

      return
  end subroutine eemove

!
! swaps the tagged index for ee moves
!
  subroutine ee_index_swap

      integer::imolty,ibox,ibox1
      real::accr

! if mstate = 1, with equal probability change the tagged index
! to another one in the same box (m = 1, still), or with the other
! box (m = 6). note that the acceptance prob of m = 1 to m = 6
! move involves a permutation factor

      imolty = ee_moltyp(mstate)
      ibox = box_state(1)
      ibox1 = box_state(fmstate)
!	write(io_output,*) 'index swap old', mstate, eepointp, box_state(mstate),
!     &             ncmt(box_state(mstate),imolty)
      if (mstate.eq.1) then
         if (random().le.0.5d0) then
            eepointp = idint(dble(ncmt(ibox,imolty))*random())+1
         else
            accr=dble(ncmt(ibox1,imolty))/(dble(ncmt(ibox,imolty))+1.0) *dexp(psi(fmstate)-psi(1))
            if (random().le.accr) then
               mstate = fmstate
               eepointp = idint(dble(ncmt(ibox1,imolty))*random())+1
!	write(io_output,*) 'mstate,nstate', 1, 6, accr
            end if
         end if
      else if (mstate.eq.fmstate) then
         if (random().le.0.5d0) then
            eepointp = idint(dble(ncmt(ibox1,imolty))*random())+1
         else
            accr=dble(ncmt(ibox,imolty))/(dble(ncmt(ibox1,imolty))+1.0) *dexp(psi(1)-psi(fmstate))
            if (random().le.accr) then
               mstate = 1
               eepointp = idint(dble(ncmt(ibox,imolty))*random())+1
            end if
!	write(io_output,*) 'mstate,nstate', 6, 1, accr
         end if
      end if

!	write(io_output,*) ibox, ibox1,parbox(eepointp,ibox,imolty)
!     &    ,parbox(eepointp,ibox1,imolty)
!	write(io_output,*) 'index swap new', mstate, eepointp, box_state(mstate),
!     &             ncmt(box_state(mstate),imolty)

      return
  end subroutine ee_index_swap

!    **********************************************************************
! make a transition of a selected molecule from state i to state j **
! in an expanded-ensemble sampling.                                **
! written on Aug. 4/99 by Bin Chen.                                **
!    **********************************************************************
  subroutine expand
    use energy_kspace,only:recip
    use energy_pairwise,only:coru

      logical::ovrlap
      integer::i,ibox,iunit,flagon,itype,j,imolty,icbu ,ic,imt,jmt,itype2,disp
      real::dchain,vnew,vold ,vintran ,vintrao,deltv,deltvb,vintern,vintero,vextn,vexto ,velectn,velecto,vdum ,vrecipo,vrecipn,vexpta,vexptb,volume,rho

! write(io_output,*) 'start expand-ensemble move'
! select a chain at random ***
      dchain  = random()
      do icbu = 1,nmolty
         if ( dchain .lt. pmeemt(icbu) ) then
            imolty = icbu
            dchain = 2.0d0
         end if
      end do
      if ( .not. lexpand(imolty) )  call err_exit(__FILE__,__LINE__,'wrong type of molecule for the ES-move',myid+1)

      if (lgrand) then
! select a chain at random in box 1!
!         (in box 2 is an ideal gas!)
         ibox = 1
         if (nchbox(ibox).eq.0) then
            bnexpc(imolty,ibox) = bnexpc(imolty,ibox) + 1.0d0
            return
         end if
         i = dint( dble(ncmt(1,imolty))*random() ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(io_output,*) 'screwup'

      else

         dchain = dble(temtyp(imolty))
         i = int( dchain*random() ) + 1
         i = parall(imolty,i)
         ibox = nboxi(i)
      end if

      iunit = nunit(imolty)

! perform a move in the expanded coefficients

 10   disp = int( rmexpc(imolty)*(2.0d0*random()-1.0d0) )
      itype = mod(eetype(imolty)+disp+numcoeff(imolty) , numcoeff(imolty))
      if ( disp .eq. 0 ) goto 10
      if ( itype .eq. 0 ) itype = numcoeff(imolty)
      do ic = 1,2
         if ( ic .eq. 1 ) then
            itype2 = eetype(imolty)
         else
            itype2 = itype
         end if
         do j = 1,iunit
            rxuion(j,ic) = rxu(i,j)
            ryuion(j,ic) = ryu(i,j)
            rzuion(j,ic) = rzu(i,j)
            qquion(j,ic) = qcharge(imolty,j,itype2)
         end do
         moltion(ic) = imolty
      end do

! calculate the energy of i in the new configuration ***

      flagon = 2
      do j = 1,iunit
         epsilon_f(imolty,j) = epsil(imolty,j,itype)
         sigma_f(imolty,j) = sigm(imolty,j,itype)
      end do
      call energy(i,imolty, vnew,vintran, vintern,vextn,velectn ,vdum,flagon, ibox,1, iunit,.false.,ovrlap,.false. ,vdum,.false.,.false.,.false.)
      if (ovrlap) return

! Start of intermolecular tail correction for new

      if ( ltailc ) then

         volume = boxlx(ibox)*boxly(ibox)*boxlz(ibox)

         vexpta = 0.0d0

         do imt = 1, nmolty
            do jmt = 1, nmolty
               rho = dble(ncmt(ibox,jmt))/volume
               vexpta = vexpta + dble( ncmt(ibox,imt) ) * coru(imt,jmt,rho,ibox)
            end do
         end do
         vnew = vnew + vexpta
         vintern = vintern + vexpta
      end if

! calculate the energy of i in the old configuration ***
      flagon = 1
      do j = 1, iunit
         epsilon_f(imolty,j) = epsil(imolty,j,eetype(imolty))
         sigma_f(imolty,j) = sigm(imolty,j,eetype(imolty))
      end do
      call energy(i,imolty,vold,vintrao,vintero,vexto,velecto ,vdum,flagon,ibox,1, iunit,.false.,ovrlap,.false. ,vdum,.false.,.false.,.false.)


! Start of intermolecular tail correction for old

      if (ltailc ) then
         vexptb = 0.0d0
         do imt = 1, nmolty
            do jmt = 1, nmolty
               rho = ncmt(ibox,jmt) / volume
               vexptb = vexptb +  dble(ncmt(ibox,imt)) * coru(imt,jmt,rho,ibox)
            end do
         end do
         vold = vold + vexptb
         vintero = vintero + vexptb

      end if

      bnexpc(imolty,ibox) = bnexpc(imolty,ibox) + 1.0d0

      if (ovrlap) call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf of TRAXYZ',myid+1)

      if ( lewald ) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      end if

! check for acceptance ***

      deltv  = vnew - vold + eta(ibox,imolty,itype)  - eta(ibox,imolty,eetype(imolty))
      deltvb = beta * deltv

      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
! accept move
      else if ( dexp(-deltvb) .gt. random() ) then
! accept move
      else
! move rejected
         return
      end if

! write(io_output,*) 'expanded move accepted i',i,exp_cion(2)
      vbox(ibox)     = vbox(ibox) + vnew - vold
      vinterb(ibox)  = vinterb(ibox) + (vintern - vintero)
      vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
      vextb(ibox)    = vextb(ibox)   + (vextn   - vexto)
      velectb(ibox)   = velectb(ibox)  + (velectn - velecto)
      vtailb(ibox) = vtailb(ibox) + vexpta - vexptb

      ncmt2(ibox,imolty,itype) = ncmt2(ibox,imolty,itype) + 1
      ncmt2(ibox,imolty,eetype(imolty)) =  ncmt2(ibox,imolty,eetype(imolty)) - 1
      eetype(imolty) = itype
      do j = 1,iunit
         qqu(i,j) = qquion(j,2)
         epsilon_f(imolty,j) = epsil(imolty,j,itype)
         sigma_f(imolty,j) = sigm(imolty,j,itype)
      end do


      if (lewald) then
! update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      end if

      if (ldielect) then
          call dipole(ibox,1)
      end if

      bsexpc(imolty,ibox) = bsexpc(imolty,ibox) + 1.0d0

! write(io_output,*) 'end expand-ensemble move'

      return
  end subroutine expand

  subroutine read_expand
    use util_runtime,only:err_exit
    integer::imol,j,ii
    do imol=1,nmolty
       if ( lexpand(imol) ) then
          if ( temtyp(imol) .gt. 1 ) then
             call err_exit(__FILE__,__LINE__,'Only one molecule of this type is allowed!',myid+1)
          end if
          read(7,*)
          read(7,*) numcoeff(imol)
          do j = 1,numcoeff(imol)
             read(7,*)
             read(7,*) (epsil(imol,ii,j),ii=1,nunit(imol))
             write(io_output,*) 'itype:',j
             write(io_output,*) (epsil(imol,ii,j),ii=1,nunit(imol))
             read(7,*) (sigm(imol,ii,j),ii=1,nunit(imol))
             write(io_output,*) (sigm(imol,ii,j),ii=1,nunit(imol))
             read(7,*) (qcharge(imol,ii,j),ii=1,nunit(imol))
             write(io_output,*) (qcharge(imol,ii,j),ii=1,nunit(imol))
             read(7,*)
             read(7,*) (eta(ii,imol,j),ii=1,2)
             write(io_output,*) 'eta:',(eta(ii,imol,j),ii=1,2)
          end do
       end if
    end do
  end subroutine read_expand

  subroutine init_ee
    integer::jerr
    allocate(epsil(ntmax,numax,100),sigm(ntmax,numax,100),qcharge(ntmax,numax,100),bnexpc(ntmax,nbxmax),bsexpc(ntmax,nbxmax),eta(nbxmax,ntmax,20),numcoeff(ntmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'init_ee: allocation failed',jerr)
    end if
    bsexpc = 0.0d0
    bnexpc = 0.0d0

    call read_expand()
  end subroutine init_ee

  subroutine output_ee_stats(io_output)
    integer,intent(in)::io_output
    integer::i,j

    write(io_output,*)
    write(io_output,*)    '### Expanded Ensemble Move  ###'
    write(io_output,*)
    do i = 1, nmolty
       do j = 1,nbox
          if (lexpand(i) .and. bnexpc(i,j) .gt. 0.5) then
             write(io_output,*) 'molecule typ =',i,'  box =',j
             write(io_output,"(' attempts =',f8.1,'   accepted =',f8.1, ' accepted % =',f7.3)") bnexpc(i,j),bsexpc(i,j), bsexpc(i,j)/bnexpc(i,j)
          end if
       end do
    end do
  end subroutine output_ee_stats

!DEC$ ATTRIBUTES FORCEINLINE :: ee_coru
  function ee_coru(ibox,imolty,flagon) result(vtail)
    use sim_cell,only:cell_vol
    use energy_pairwise,only:coru
    real::vtail
    integer,intent(in)::ibox,imolty,flagon

    real::vol,rho
    integer::kmolty,jmolty

    vtail=0.0d0

    if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
       vol = cell_vol(ibox)
    else
       vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
    end if
    do kmolty = 1, nmolty
       do jmolty = 1, nmolty
! rho = ncmt(ibox,jmolty) /
!     &               ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
          if (flagon.eq.1) then
             rho = ncmt(ibox,jmolty) / vol
             vtail = vtail + ncmt(ibox,kmolty) * coru(kmolty,jmolty,rho,ibox)
! write(io_output,*) 'vtail',vtail
          else
             if (jmolty.eq.ee_moltyp(mstate)) then
                rho = (ncmt(ibox,jmolty)-1) / vol
             else if (jmolty.eq.imolty) then
                rho = (ncmt(ibox,jmolty)+1) / vol
             else
                rho = ncmt(ibox,jmolty) / vol
             end if
             if (kmolty.eq.imolty) then
                vtail = vtail + (ncmt(ibox,kmolty)+1) * coru(kmolty,jmolty,rho ,ibox)
             else if (kmolty.eq.ee_moltyp(mstate)) then
                vtail = vtail + (ncmt(ibox,kmolty)-1) * coru(kmolty,jmolty,rho ,ibox)
             else
                vtail = vtail + ncmt(ibox,kmolty) * coru(kmolty,jmolty,rho,ibox)
             end if
          end if
       end do
    end do

  end function ee_coru
end MODULE moves_ee
