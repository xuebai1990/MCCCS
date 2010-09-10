      subroutine ee_index_swap
c
c swaps the tagged index for ee moves
c
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'eepar.inc'

      integer::imolty,ibox,ibox1
      real(8)::accr,random

c --- if mstate = 1, with equal probability change the tagged index
c --- to another one in the same box (m = 1, still), or with the other
c --- box (m = 6). note that the acceptance prob of m = 1 to m = 6
c --- move involves a permutation factor

      imolty = ee_moltyp(mstate)
      ibox = box_state(1)
      ibox1 = box_state(fmstate)
c	write(iou,*) 'index swap old', mstate, eepointp, box_state(mstate),
c     &             ncmt(box_state(mstate),imolty)
      if (mstate.eq.1) then
         if (random().le.0.5d0) then
            eepointp = idint(dble(ncmt(ibox,imolty))*random())+1
         else
            accr=dble(ncmt(ibox1,imolty))/(dble(ncmt(ibox,imolty))+1.0)
     &           *dexp(psi(fmstate)-psi(1))
            if (random().le.accr) then
               mstate = fmstate
               eepointp = idint(dble(ncmt(ibox1,imolty))*random())+1
c	write(iou,*) 'mstate,nstate', 1, 6, accr
            endif
         endif
      elseif (mstate.eq.fmstate) then
         if (random().le.0.5d0) then
            eepointp = idint(dble(ncmt(ibox1,imolty))*random())+1
         else
            accr=dble(ncmt(ibox,imolty))/(dble(ncmt(ibox1,imolty))+1.0)
     &           *dexp(psi(1)-psi(fmstate))
            if (random().le.accr) then
               mstate = 1
               eepointp = idint(dble(ncmt(ibox,imolty))*random())+1
            endif
c	write(iou,*) 'mstate,nstate', 6, 1, accr
         endif
      endif      

c	write(iou,*) ibox, ibox1,parbox(eepointp,ibox,imolty)
c     &    ,parbox(eepointp,ibox1,imolty)
c	write(iou,*) 'index swap new', mstate, eepointp, box_state(mstate),
c     &             ncmt(box_state(mstate),imolty)

      return
      end
