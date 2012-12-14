      subroutine eemove
      
! peforms EE move. currently works when EE move is performed on
! only one type of molecule. put all states in fort.77 and add
! corresponding bead types in suijtab. can only do EE if the number
! of beads between states remain constant.

      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_runtime,only:err_exit
      use util_math
      use util_string
      use util_files
      use util_timings
      use transfer_swap,only:swap
      implicit none
      include 'common.inc'
      
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'conver.inc'
!$$$      include 'coord2.inc'
!$$$      include 'system.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'rosen.inc' 
!$$$      include 'boltzmann.inc'
!$$$      include 'external.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'poten.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'clusterbias.inc'
!$$$      include 'neigh.inc'
!$$$      include 'cell.inc'
!$$$      include 'eepar.inc'

      logical::ovrlap
      integer(KIND=normal_int)::i,ibox,iunit,imolty,imolty1,j ,idummy(nmax)
      real(KIND=double_precision)::dum,vrecipn,vrecipo,vnew,vold,vintern ,vintero,vintran,vintrao,velectn,velecto,vewaldn,vewaldo,vextn ,vexto,deltv,deltvb,random,wdeltvb,vtailn,vtailo
      
! --------------------------------------------------------------------

!      write(*,*) 'START EEMOVE'
!      write(11,*) '1:',neigh_cnt(18)

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
         call err_exit('')
      end if

! --- type of move depending upon mstate and nstate. one type of move
! --- is 'swap', other is usual ee

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

! --- energy for the new state

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
         call ee_energy(eeirem,imolty,vnew,vintran,vintern,vextn,velectn ,vewaldn,vtailn,2,ibox,1,iunit,.false.,ovrlap ,.false.,dum,.false.,.false.)
!         call energy(eeirem,imolty,vnew,vintran,vintern,vextn,velectn
!     &               ,vewaldn,2,ibox,1,iunit,.false.,ovrlap
!     &               ,.false.,dum,.false.,.false.)

         if (ovrlap) goto 100

! --- energy for the old state

         imolty = ee_moltyp(mstate)
         do i = 1, iunit
            rxuion(i,1) = rxu(eeirem,i)
            ryuion(i,1) = ryu(eeirem,i)
            rzuion(i,1) = rzu(eeirem,i)
            qquion(i,1) = qqu(eeirem,i)
         end do
         moltion(2) = imolty
         call ee_energy(eeirem,imolty,vold,vintrao,vintero,vexto,velecto ,vewaldo,vtailo,1,ibox,1,iunit,.false.,ovrlap ,.false.,dum,.false.,.false.)
!         call energy(eeirem,imolty,vold,vintrao,vintero,vexto,velecto
!     &               ,vewaldo,1,ibox,1,iunit,.false.,ovrlap
!     &               ,.false.,dum,.false.,.false.)
!	write(io_output,*) vnew,vold
         if (ovrlap) then
            write(io_output,*) 'disaster ovrlap in old conf eemove'
            call err_exit('')
         end if

         if (lewald.and.(lelect(moltion(2)).or.lelect(moltion(1))))then
            call ee_recip(ibox,vrecipn,vrecipo,1)
            velectn = velectn+vrecipn+vewaldn
            velecto = velecto+vrecipo+vewaldo
            vnew = vnew + vrecipn
            vold = vold + vrecipo
         end if

! --- check for acceptance

         deltv = (vnew - vold)
         deltvb = beta*deltv
         wdeltvb = wee_ratio*dexp(-deltvb)

         if ((deltvb-dlog(wee_ratio)).le.0.0d0) then
!           --- accept move
         else if (wdeltvb.gt.random()) then
!           --- accept move
         else
!           --- reject move
            goto 100
         end if

! --- update new

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

! --- update reciprocal space term

         call ee_recip(ibox,dum,dum,2)

      end if

! --- update the present state

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

! --- update parall (since one molecule has changed its state/type)

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
      end
