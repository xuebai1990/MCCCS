      subroutine eemove
      
c eemove

c peforms EE move. currently works when EE move is performed on
c only one type of molecule. put all states in fort.77 and add
c corresponding bead types in suijtab. can only do EE if the number
c of beads between states remain constant.

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

      implicit none
      
      include 'control.inc'
      include 'coord.inc'
      include 'conver.inc'
      include 'coord2.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'cbmc.inc'
      include 'rosen.inc' 
      include 'boltzmann.inc'
      include 'external.inc'
      include 'inputdata.inc'
      include 'ewaldsum.inc'
      include 'poten.inc'
      include 'fepsi.inc'
      include 'clusterbias.inc'
      include 'neigh.inc'
      include 'cell.inc'
      include 'eepar.inc'

      logical::ovrlap
      integer::i,ibox,iunit,imolty,imolty1,j,idummy(nmax)
      real(8)::dum,vrecipn,vrecipo,vnew,vold,vintern,vintero
     &                ,vintran,vintrao,velectn,velecto,vewaldn,vewaldo
     &                ,vextn,vexto,deltv,deltvb,random,wdeltvb,vtailn
     &                ,vtailo
      
C --------------------------------------------------------------------

c      write(*,*) 'START EEMOVE'
c      write(11,*) '1:',neigh_cnt(18)

      leemove = .true.
      leeacc = .false.

c choose nstate, given mstate

      if (mstate.eq.1) then
         nstate = mstate + 1
      elseif (mstate.eq.fmstate) then
         nstate = mstate - 1
      else
         if (random().le.0.5d0) then
            nstate = mstate - 1
         else
            nstate = mstate + 1
         endif
      endif

c	write(iou,*) 'typ', mstate, ee_moltyp(mstate),eepointp,
c     &   box_state(mstate), ncmt(box_state(mstate),ee_moltyp(mstate))

      if (ncmt(box_state(mstate),ee_moltyp(mstate)).eq.0) then
         write(iou,*)'problem: mstate, but no molecule in mstate',mstate
         call cleanup('')
      endif

c --- type of move depending upon mstate and nstate. one type of move
c --- is 'swap', other is usual ee

      wee_ratio = dexp(psi(nstate)-psi(mstate))*
     &      um_markov(nstate,mstate)/um_markov(mstate,nstate)

c	write(iou,*) 'mstate', mstate, ee_moltyp(mstate)
c	write(iou,*) 'nstate', nstate, ee_moltyp(nstate)

      if ((mstate.eq.sstate1.and.nstate.eq.sstate2).or.
     &    (mstate.eq.sstate2.and.nstate.eq.sstate1)) then

         boxrem1 = box_state(mstate)
         boxins1 = box_state(nstate)

         eeirem = parbox(eepointp,boxrem1,ee_moltyp(mstate))
c	write(iou,*) 'eeirem', eeirem, eepointp

         call swap(dum,dum,dum,dum,dum,dum,dum,dum,dum)

         if (.not.leeacc) goto 100

      else

c --- energy for the new state

         ibox = box_state(mstate)
         eeirem = parbox(eepointp,ibox,ee_moltyp(mstate))
c	write(iou,*) 'eeirem1', eeirem, eepointp
         imolty = ee_moltyp(nstate)

         iunit = nunit(imolty)
         do i = 1, iunit
            rxuion(i,2) = rxu(eeirem,i)
            ryuion(i,2) = ryu(eeirem,i)
            rzuion(i,2) = rzu(eeirem,i)
            qquion(i,2) = ee_qqu(i,nstate)
         enddo

         moltion(1) = imolty
         call ee_energy(eeirem,imolty,vnew,vintran,vintern,vextn,velectn
     &               ,vewaldn,vtailn,2,ibox,1,iunit,.false.,ovrlap
     &               ,.false.,dum,.false.,.false.)
c         call energy(eeirem,imolty,vnew,vintran,vintern,vextn,velectn
c     &               ,vewaldn,2,ibox,1,iunit,.false.,ovrlap
c     &               ,.false.,dum,.false.,.false.)

         if (ovrlap) goto 100

c --- energy for the old state

         imolty = ee_moltyp(mstate)
         do i = 1, iunit
            rxuion(i,1) = rxu(eeirem,i)
            ryuion(i,1) = ryu(eeirem,i)
            rzuion(i,1) = rzu(eeirem,i)
            qquion(i,1) = qqu(eeirem,i)
         enddo
         moltion(2) = imolty
         call ee_energy(eeirem,imolty,vold,vintrao,vintero,vexto,velecto
     &               ,vewaldo,vtailo,1,ibox,1,iunit,.false.,ovrlap
     &               ,.false.,dum,.false.,.false.)
c         call energy(eeirem,imolty,vold,vintrao,vintero,vexto,velecto
c     &               ,vewaldo,1,ibox,1,iunit,.false.,ovrlap
c     &               ,.false.,dum,.false.,.false.)
c	write(iou,*) vnew,vold
         if (ovrlap) then
            write(iou,*) 'disaster ovrlap in old conf eemove'
            call cleanup('')
         endif

         if (lewald.and.(lelect(moltion(2)).or.lelect(moltion(1))))then
            call ee_recip(ibox,vrecipn,vrecipo,1)
            velectn = velectn+vrecipn+vewaldn
            velecto = velecto+vrecipo+vewaldo
            vnew = vnew + vrecipn
            vold = vold + vrecipo
         endif

c --- check for acceptance

         deltv = (vnew - vold)
         deltvb = beta*deltv
         wdeltvb = wee_ratio*dexp(-deltvb)

         if ((deltvb-dlog(wee_ratio)).le.0.0d0) then
c           --- accept move
         elseif (wdeltvb.gt.random()) then
c           --- accept move
         else
c           --- reject move
            goto 100
         endif

c --- update new

         imolty1 = ee_moltyp(nstate)
         parbox(ncmt(ibox,imolty1)+1,ibox,imolty1) = eeirem
c	write(iou,*) 'update new', eepointp,ibox,imolty1,
c     &          parbox(eepointp,ibox,imolty1),
c     &          parbox(eepointp,ibox,imolty)

         parbox(eepointp,ibox,imolty) =
     &      parbox(ncmt(ibox,imolty),ibox,imolty)

c	write(iou,*) 'update new1', imolty,parbox(eepointp,ibox,imolty),
c     &              ncmt(ibox,imolty)

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
c	write(iou,*) vtailn,vtailo,vintern,vintero

c --- update reciprocal space term

         call ee_recip(ibox,dum,dum,2)

      endif

c --- update the present state

c	ovr = .true.

      temtyp(ee_moltyp(mstate)) = temtyp(ee_moltyp(mstate)) - 1
      temtyp(ee_moltyp(nstate)) = temtyp(ee_moltyp(nstate)) + 1
c	write(iou,*) 'ttym', mstate,ee_moltyp(mstate),
c     &               temtyp(ee_moltyp(mstate))
c	write(iou,*) 'ttyn', nstate,ee_moltyp(nstate),
c     &               temtyp(ee_moltyp(nstate))
      mstate = nstate

c	write(iou,*) 'eemove eeirem', eeirem
      moltyp(eeirem) = ee_moltyp(nstate)
      do i = 1, nunit(ee_moltyp(nstate))
         qqu(eeirem,i) = ee_qqu(i,nstate)
c	write(iou,*) 'charges', qqu(eeirem,i)
c	write(iou,*) 'charges1'
      enddo

      if (mstate.eq.1.or.mstate.eq.fmstate) then
         lmstate = .true.
      else
         lmstate = .false.
      endif

c --- update parall (since one molecule has changed its state/type)

      do  i = 1, nmolty
         idummy(i) = 0
      enddo
      do j = 1, nchain
         i = moltyp(j)
         idummy(i) = idummy(i) + 1
         parall(i,idummy(i)) = j
      enddo

c	write(iou,*) 'accept', eepointp, mstate, nstate, eeirem, boxins1
c	write(iou,*)'accept1',parbox(eepointp,boxins1,ee_moltyp(nstate))
c	write(iou,*) '1line', leemove,lmstate,leeacc
c	write(iou,*) '2line', fmstate,sstate1,sstate2,wee_ratio,eeratio
c	write(iou,*) '3line',(ee_moltyp(i), i = 1, fmstate)
c	write(iou,*) '4line',(box_state(i), i = 1, fmstate)
c	write(iou,*) '5line',(psi(i), i = 1, fmstate)
c	write(iou,*) '6line',(ee_qqu(1,i), i = 1, fmstate)
c	write(iou,*) '7line',(um_markov(i,i+1),i = 1,fmstate-1)
c	write(iou,*) '8line',(um_markov(i,i-1),i = 2,fmstate)

 100    continue
        leemove = .false.

      return
      end
