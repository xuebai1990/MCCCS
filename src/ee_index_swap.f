      subroutine ee_index_swap
!
! swaps the tagged index for ee moves
!
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'eepar.inc'

      integer(KIND=int)::imolty,ibox,ibox1
      real(KIND=double_precision)::accr,random

! --- if mstate = 1, with equal probability change the tagged index
! --- to another one in the same box (m = 1, still), or with the other
! --- box (m = 6). note that the acceptance prob of m = 1 to m = 6
! --- move involves a permutation factor

      imolty = ee_moltyp(mstate)
      ibox = box_state(1)
      ibox1 = box_state(fmstate)
!	write(iou,*) 'index swap old', mstate, eepointp, box_state(mstate),
!     &             ncmt(box_state(mstate),imolty)
      if (mstate.eq.1) then
         if (random().le.0.5d0) then
            eepointp = idint(dble(ncmt(ibox,imolty))*random())+1
         else
            accr=dble(ncmt(ibox1,imolty))/(dble(ncmt(ibox,imolty))+1.0)
     &           *dexp(psi(fmstate)-psi(1))
            if (random().le.accr) then
               mstate = fmstate
               eepointp = idint(dble(ncmt(ibox1,imolty))*random())+1
!	write(iou,*) 'mstate,nstate', 1, 6, accr
            end if
         end if
      elseif (mstate.eq.fmstate) then
         if (random().le.0.5d0) then
            eepointp = idint(dble(ncmt(ibox1,imolty))*random())+1
         else
            accr=dble(ncmt(ibox,imolty))/(dble(ncmt(ibox1,imolty))+1.0)
     &           *dexp(psi(1)-psi(fmstate))
            if (random().le.accr) then
               mstate = 1
               eepointp = idint(dble(ncmt(ibox,imolty))*random())+1
            end if
!	write(iou,*) 'mstate,nstate', 6, 1, accr
         end if
      end if      

!	write(iou,*) ibox, ibox1,parbox(eepointp,ibox,imolty)
!     &    ,parbox(eepointp,ibox1,imolty)
!	write(iou,*) 'index swap new', mstate, eepointp, box_state(mstate),
!     &             ncmt(box_state(mstate),imolty)

      return
      end
