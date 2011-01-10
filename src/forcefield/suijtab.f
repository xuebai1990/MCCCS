      subroutine suijtab( lmixlb,lmixjo)

! suijtab
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

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'external.inc'
      include 'poten.inc'
      include 'system.inc'
      include 'cell.inc'
      include 'expsix.inc'
      include 'merck.inc'
      include 'nsix.inc'
!kea
      include 'garofalini.inc'
      include 'conver.inc'

      logical::lmixlb, lmixjo

      integer::i, j, ij, ji,ibox,nmix,imix,jerr
      real(8)::rzeronx(nxatom),epsilonnx(nxatom)

      real(8)::rcheck, sr2, sr6, adum, bdum, rs1, rs7, sr7,
     &     pi,djay,sigmaTmp,epsilonTmp

! ----------------------------------------------------------------
      open(5,file='lj.poten')
      read(5,*)
      read(5,*) nntype
      read(5,*)
      allocate(sigi(nntype),epsi(nntype),mass(nntype),qelect(nntype)
     &     ,xiq(nntype),lqchg(nntype),lij(nntype),chname(nntype)
     &     ,chemid(nntype),jayself(nntype),extc12(nntype)
     &     ,extc3(nntype),extz0(nntype),q1(nntype),lpl(nntype)
     &     ,q1(nntype),sig2ij(nntype*nntype),epsij(nntype*nntype)
     &     ,ecut(nntype*nntype),jayq(nntype*nntype)
     &     ,aspecial(nntype*nntype),bspecial(nntype*nntype)
     &     ,lspecial(nntype*nntype)
     &     ,stat=jerr)
      if (jerr.ne.0) call cleanup('suijtab: allocation failed')

      sig2ij=0.
      epsij=0.
      lpl = .false.

      pi = 4.0d0*datan(1.0d0)
     
      if(lshift) then
         do ibox = 2,nbox
           if (dabs(rcut(1)-rcut(ibox)).gt.1.0d-10) then
               call cleanup('Keep rcut for each box same')
           end if
         end do
      end if

! KEA adding garofalini silica/water potential
!     see J. Phys. Chem. 94 5351 (1990)
      if (lgaro) then
         do i=1,6
            galpha(i) = 0.0d0
            grho(i) = 0.0d0
            gbeta(i) = 0.0d0
            ecut(i) = 0.0d0
            do j=1,3
               ga(i,j) = 0.0d0
               gb(i,j) = 0.0d0
               gc(i,j) = 0.0d0
            end do
         end do
         do i=1,4
            glambda(i) = 0.0d0
            do j=1,2
               grij(i,j) = 0.0d0
               grijsq(i,j) = 0.0d0
               ggamma(i,j) = 0.0d0
            end do
            gtheta(i) = 0.0d0
         end do
!     form:
!     U(2) = A exp(-rij/rho)+[zi zj erfc (rij/beta)]/rij + a/[1+exp(b/(rij-c))]
!     U(3) = h3(rij rik thetajik) + h3(rjk rji thetakji) + h3(rki rkj thetaijk)
!     h3(rij rik thetajik) = lambda exp[(gamma/(rij-rij*))+(gamma/(rik-rik*))]*
!             [cos(theta)-cos(theta*)]**2     for rij<rij* and rik<rik*
!                          = 0 otherwise

! conversion for coulomb potential (incl C/e, m-->A, J-->K, 4pi epsi naught)
       qqfact = 1.67125d5
! qqfact for bohr/hartree - from pot_KAng.f code
!         qqfact = 0.99865377d0

!     Parameters (galpha,grho,gbeta,ga,gb,gc; lambda,grij,ggamma,gtheta)
!     Si-Si
       galpha(1) = 13597175.7d0
       grho(1) = 0.29d0
       gbeta(1) = 2.29d0
       lqchg(1) = .true.
       qelect(1) = 4.0d0
       mass(1) = 28.09d0
       ecut(1) = garofalini(rcut(1)*rcut(1),1,qelect(1),qelect(1)
     &      ,1,1)
       chname(1) = ' Garo Si'
       chemid(1) = 'Si '

!     O-O
       galpha(2) = 5251972.5d0
       grho(2) = 0.29d0
       gbeta(2) = 2.34d0
       lqchg(2) = .true.
       qelect(2) = -2.0d0
       mass(2) = 16.00d0
       ecut(2) = garofalini(rcut(1)*rcut(1),2,qelect(2),qelect(2)
     &      ,2,2)
       chname(2) = ' Garo O'
       chemid(2) = 'O  '
       
!     H-H
       galpha(3) = 246299.4d0
       grho(3) = 0.35d0
       gbeta(3) = 2.1d0
       ga(3,1) = -38243.8d0
       gb(3,1) = 6.0d0 
       gc(3,1) = 1.51d0
       ga(3,2) = 2515.9d0
       gb(3,2) = 2.0d0 
       gc(3,2) = 2.42d0
       lqchg(3) = .true.
       qelect(3) = 1.0d0
       mass(3) = 1.0078d0
       ecut(3) = garofalini(rcut(1)*rcut(1),3,qelect(3),qelect(3)
     &      ,3,3)
       chname(3) = ' Garo H'
       chemid(3) = 'H  '
       
!     Si-O
       galpha(4) =  21457024.2d0
       grho(4) = 0.29d0
       gbeta(4) = 2.34d0
       ecut(4) = garofalini(rcut(1)*rcut(1),4,qelect(1),qelect(2)
     &      ,1,2)
       
!     Si-H
       galpha(5) = 499842.9d0
       grho(5) = 0.29d0
       gbeta(5) = 2.31d0
       ga(5,1) = -33715.5d0
       gb(5,1) = 6.0d0 
       gc(5,1) = 2.2d0
       ecut(5) = garofalini(rcut(1)*rcut(1),5,qelect(1),qelect(3)
     &      ,1,3)
       
!     0-H
       galpha(6) = 2886049.4d0
       grho(6) = 0.29d0
       gbeta(6) = 2.26d0
       ga(6,1) = -15096.7d0
       gb(6,1) = 15.0d0 
       gc(6,1) = 1.05d0
       ga(6,2) = 55353.6d0
       gb(6,2) = 3.2d0 
       gc(6,2) = 1.50d0
       ga(6,3) = -6038.7d0
       gb(6,3) = 5.0d0 
       gc(6,3) = 2.0d0
       ecut(6) = garofalini(rcut(1)*rcut(1),6,qelect(2),qelect(3)
     &      ,2,3)
       
!     Si-O-Si
       glambda(1) = 21732.3d0
       ggamma(1,1) = 2.0d0
       grij(1,1) = 2.6d0
       grijsq(1,1) = grij(1,1)*grij(1,1)
       gtheta(1) = dcos(109.5d0*degrad)
       
!     O-Si-O
       glambda(2) = 1376379.0d0
       ggamma(2,1) = 2.8d0
       grij(2,1) = 3.0d0
       grijsq(2,1) = grij(2,1)*grij(2,1)
       gtheta(2) = dcos(109.5d0*degrad)
       
!     H-O-H
       glambda(3) = 2535435.0d0
       ggamma(3,1) = 1.3d0
       grij(3,1) = 1.6d0
       grijsq(3,1) = grij(3,1)*grij(3,1)
       gtheta(3) = dcos(104.5d0*degrad)
       
!     Si-O-H
       glambda(4) = 362205.0d0
       ggamma(4,1) = 2.0d0
       ggamma(4,2) = 1.2d0
       grij(4,1) = grij(1,1)
       grij(4,2) = 1.5d0
       grijsq(4,1) = grij(4,1)*grij(4,1)
       grijsq(4,2) = grij(4,2)*grij(4,2)
       gtheta(4) = dcos(109.5d0*degrad)
       
       do i=1,6
          write(iou,*) 'garo ecut',i,ecut(i)
       end do
       return

      elseif(lexpsix) then
! --- Keep the rcut same for each box
         do ibox = 2,nbox
           if (dabs(rcut(1)-rcut(ibox)).gt.1.0d-10) then
              call cleanup('Keep rcut for each box same')
           end if
         end do

         do i = 1,natom
            aexsix(i) = 0.0d0
            bexsix(i) = 0.0d0
            cexsix(i) = 0.0d0
            qelect(i) = 0.0d0
            xiq(i) = 0.0d0
            lqchg(i) = .false.
            lij(i) = .true.     
         end do
! --- Explicit atom carbon Williams Exp-6 potential
! --- J.Chem.Phys 47 11 4680 (1967) paramter set IV
! --- note that the combination of C--H has to be the average of C 
!      and H the way it is set up in the energy subroutines
! --- natom is set in expsix.inc
!     U(r) = A*r^(-6) + B*exp[C*r]

! --- C---C nonbonded interaction
!         aexsix(1) = -2.858d5
!         bexsix(1) = 4.208d7
!         cexsix(1) = -3.60d0
!         aexsix(1) = -2.541d5
!         bexsix(1) = 3.115d7
!         cexsix(1) = -3.60d0
!         mass(1) = 12.011d0
!         aexsix(1) = -4.0d0
!         bexsix(1) = 1.0d0/dexp(-12.0d0)
!         cexsix(1) = -12.0d0/(2.0d0 ** (1.0d0/6.0d0))
         aexsix(1) = -0.590759e6
         bexsix(1) = 0.281431e9
         cexsix(1) = -0.396875e1
         mass(1) = 39.948d0

! --- C---H nonbonded interaction
!         aexsix(2) = -6.290d4
!         bexsix(2) = 4.411d6
!         cexsix(2) = -3.67d0
         aexsix(2) = -6.441d4
         bexsix(2) = 5.5355d6
         cexsix(2) = -3.67d0

! --- H---H nonbonded interaction
!         aexsix(3) = -1.374d4
!         bexsix(3) = 1.335d6
!         cexsix(3) = -3.74d0
         aexsix(3) = -1.6254d4
         bexsix(3) = 1.323d6
         cexsix(3) = -3.74d0
         mass(3) = 1.0079d0

         if (lshift) then
            do i=1,natom
               sexsix(i) = aexsix(i)*(rcut(1)**(-6))
     &              + bexsix(i)*Exp(cexsix(i)*rcut(1))
            end do
         else 
            do i=1,natom
               consp(i) = (2.0d0/3.0d0)*pi*(2.0d0*aexsix(i)/(rcut(1)
     &              *rcut(1)*rcut(1))+bexsix(i)*dexp(cexsix(i)*rcut(1))
     &              *(-6.0d0/(cexsix(i)*cexsix(i)*cexsix(i))+6.0d0
     &              *rcut(1)/(cexsix(i)*cexsix(i))-3.0d0*rcut(1)*
     &              rcut(1)/
     &              cexsix(i)+rcut(1)*rcut(1)*rcut(1)))
               consu(i) = 2.0d0*pi*(aexsix(i)/(3.0d0*rcut(1)*rcut(1)*
     &                     rcut(1))
     &              +(-rcut(1)*rcut(1)+2.0d0*rcut(1)/cexsix(i)-2.0d0/
     &                      (cexsix(i)*
     &               cexsix(i)))*bexsix(i)*dexp(cexsix(i)*rcut(1))/
     &               cexsix(i))
            end do
!            write(11,*) 'consp(i)',consp
!            write(11,*) 'consu(i)',consu
         end if

         write(iou,*) 
     &        ' i   aexsix       bexsix      cexsix     sexsix'
         do i = 1,natom
            write(iou,'(i3,2x,4e12.4)')i,aexsix(i),bexsix(i)
     &           ,cexsix(i),sexsix(i)
         end do

         return

      end if

      if (lmmff) then
! --- Keep the rcut same for each box
      do ibox = 2,nbox
         if (dabs(rcut(1)-rcut(ibox)).gt.1.0d-10) then
            call cleanup('Keep rcut for each box same')
         end if
      end do
! --- Merk Molecular Force Field (MMFF)
! --- J. Am. Chem. Soc. 1992 Vol.114 P7827-7843 (Thomas A. Halgren)
! --- natomtyp is set in mmff.inc
!     U(r) = epsi*(1.07/(rs+0.07))^7 * (1.12/(rs^7+0.12)-2)
!     rs = r / sigimmff

! --- C---C nonbonded interaction ***

         alphammff(1) = 1.050d0
         nmmff(1) = 2.490d0
         ammff(1) = 3.890d0
         gmmff(1) = 1.282d0
         sigimmff(1) = ammff(1)*dsqrt(dsqrt(alphammff(1)))
         epsimmff(1) = 45582.6d0*gmmff(1)*gmmff(1)*dsqrt(nmmff(1))
     &        /(ammff(1)**6.0d0)
!         sigimmff(1) = 1.126763255691509d0
!         epsimmff(1) = 0.9994354715851470d0
         sigisq(1) = sigimmff(1)*sigimmff(1)
         mass(1) = 12.011d0
!         mass(1) = 1.0d0
! --- H---H nonbonded interaction ***

         alphammff(3) = 0.250d0
         nmmff(3) = 0.800d0
         ammff(3) = 4.200d0
         gmmff(3) = 1.209d0
         sigimmff(3) = ammff(3)*dsqrt(dsqrt(alphammff(3)))
         epsimmff(3) = 45582.6d0*gmmff(3)*gmmff(3)*dsqrt(nmmff(3))
     &        /(ammff(3)**6.0d0)
         sigisq(3) = sigimmff(3)*sigimmff(3)
         mass(3) = 1.0078d0

! --- C---H nonbonded interaction by using cominbination rule ***
         sigimmff(2) = 0.5d0*(sigimmff(1)+sigimmff(3))*(1.0d0 + 
     &        0.2d0*(1.0d0-dexp(-12.d0*(((sigimmff(1)-sigimmff(3))/
     &        (sigimmff(3)+sigimmff(1)))**2.0d0))))
         epsimmff(2) = 91165.1d0*gmmff(1)*gmmff(3)*alphammff(1)*
     &        alphammff(3)/((dsqrt(alphammff(1)/nmmff(1))+
     &        dsqrt(alphammff(3)/nmmff(3)))*(sigimmff(2)**6.0d0))
         sigisq(2) = sigimmff(2)*sigimmff(2)


         if (lshift) then
            do i=1,natomtyp
               rs1 = rcut(1)/sigimmff(i)
               rs7 = rs1**7.0d0
               sr7 = (1.07d0/(rs1+0.07d0))**7.0d0               
               smmff(i) = epsimmff(i)*sr7*(1.12d0/(rs7+0.12)-2)
            end do
         else
            do i = 1,natomtyp
               smmff(i) = 0.0d0
            end do
            coru_cons(1) = -2.4837937263569310d-02
!            coru_cons(1) = -5.8244592746534724E-03
            coru_cons(2) = -1.7583010189381791d-02
            coru_cons(3) = -8.3770412792126582d-03
            corp_cons(1) = 0.1696349613545569d0
!            corp_cons(1) = 4.0098456560842058E-02
            corp_cons(2) = 0.1203650950025348d0
            corp_cons(3) = 5.7576802340310304d-02
         end if

         write(iou,*) ' i   epsimmff     sigimmff   smmff'
         do i = 1,natom
            write(iou,'(i3,2x,4e12.4)')i,epsimmff(i),sigimmff(i)
     &           ,smmff(i)
         end do

         return

      end if

      do i = 1, nntype
         sigi(i) = 0.0d0
         epsi(i) = 0.0d0
         mass(i) = 0.0d0
         qelect(i) = 0.0d0
         xiq(i) = 0.0d0
         lqchg(i) = .false.
         lij(i) = .true.
         chname(i) = '                  '
      end do

      if (lninesix) then
! * special potential for all-atom formic acid from llnl 4/6/04 jms 
      do ibox = 2,nbox
         if (dabs(rcut(1)-rcut(ibox)).gt.1.0d-10) then
            call cleanup('Keep rcut for each box same')
         end if
      end do


! *** sp2 carbon site H-[C](=O)-OH
         rzeronx(1) = 4.1834161d0
         epsilonnx(1) = 45.224d0
         qelect(1) = 0.44469d0
         lqchg(1) = .true.
         mass(1) = 12.011d0
         chname(1) = ' 9-6 formic acid C'

! *** hydroxyl oxygen site H-C(=O)-[O]H
         rzeronx(2) = 3.5694293136d0
         epsilonnx(2) = 47.151d0
         qelect(2) = -0.55296d0
         lqchg(2) = .true.
         mass(2) = 15.9996d0
         chname(2) = ' 9-6 formic acid O'

! *** carbonyl oxygen site H-C[=O]-OH
         rzeronx(3) = 3.0014635172d0
         epsilonnx(3) = 146.008d0
         qelect(3) = -0.43236d0
         lqchg(3) = .true.
         mass(3) = 15.9996d0
         chname(3) = ' 9-6 formic acid=O'

! *** hydrogen site [H]-C(=O)-OH
         rzeronx(4) = 0.8979696387d0
         epsilonnx(4) = 2.4054d0
         qelect(4) = 0.10732d0
         lqchg(4) = .true.
         mass(4) = 1.00794d0
         chname(4) = ' 9-6 formic acid H'

! *** acidic hydrogen site H-C(=O)-O[H] 
         rzeronx(5) = 1.115727276d0
         epsilonnx(5) = 12.027d0
         qelect(5) = 0.43331d0
         lqchg(5) = .true.
         mass(5) = 1.00794d0
         chname(5) = ' 9-6 formic acidOH'

! * calculate all site-site parameters via Lorentz-Berthelot rules
         do i = 1,nxatom
            do j = 1,nxatom
               ij = (i-1)*nxatom + j
               rzero(ij) = 0.5d0*(rzeronx(i) + rzeronx(j))
               epsnx(ij) = dsqrt(epsilonnx(i)*epsilonnx(j))
               if (lshift) then
                  shiftnsix(ij) = 4.0d0*epsnx(ij)*(rzero(ij)
     &                   /rcut(1))**6 *
     &               (2.0d0*(rzero(ij)/rcut(1))**3 - 3.0d0) 
               end if
            end do
         end do

         qqfact = 1.67d5

         return

      end if

! ================================================
! *** Begin Lennard-Jones Site Parameter Input ***
! ===================================================c
! Notes.  Throughout the listings the atoms or       c
! pseudoatoms that each entry is for is in brackets, c
! such as [CH3] for a united-atom methyl group.      c
! A double bond is "=" and a triple bond is "=-"     c
! PLEASE ADD NEW PARAMETERS TO END OF LISTING!!!     c
! ===================================================c

! * 1-5 correction term for unprotected hydrogen-oxygen interaction
! * for ether oxygens
      a15(1) = 4.0d7
! * for alcohol oxygens
      a15(2) = 7.5d7

! * OLD VALUES
! * this is 17^6
!      a15 = 24137569.0d0
! * this is 16^6
!      a15 = 16777216.0d0

      do j=1,nntype
         read(5,'(I,4F,2A)') i,sigi(i),epsi(i),qelect(i),mass(i),chemid(i),chname(i)
         if (qelect(i).ne.0) lqchg(i)=.true.
         if (sigi(i).eq.0 .and. epsi(i).eq.0) lij(i)=.false. 
      end do

! ===========================
! *** End Parameter Input ***
! ===========================

! --- Computation of un-like interactions
      if ( ljoe ) then
! --- STANDARD METHYL GROUP
         extc12(1) = 3.41d7
         extc3(1)  = 20800.0d0
         extz0(1)  = 0.86d0

! --- STANDARD METHYLENE GROUP
         extc12(5) = 2.80d7
         extc3(5)  = 17100.0d0
         extz0(5)  = 0.86d0

! --- Methane
         extc12(3) = 3.41d7
         extc3(3)  = 20800.0d0
         extz0(3)  = 0.86d0

! --- Martin's methyl (CH3)
         extc12(18) = 3.41d7
         extc3(18)  = 20800.0d0
         extz0(18)  = 0.86d0
      end if

! --- Assign jayq for pairs
      do i = 1,nntype
         do j = 1,nntype
            ij = (i-1)*nntype + j
            jayq(ij) = 0.0d0
         end do
      end do

! - CO2-FQ Carbon-Oxygen cross term (JCO)
      i = 131
      j = 132
      djay = (503.2d0)*(133.905d0)
      ij = (i-1)*nntype + j
      ji = (j-1)*nntype + i
      jayq(ij) = djay
      jayq(ji) = djay

! - CO2-FQ Oxygen-Oxygen cross term (JOO)
      i = 132
      j = 132
      djay = (503.2d0)*(1.09d0)
      ij = (i-1)*nntype + j
      jayq(ij) = djay

! --- SPC-FQ water Oxygen-Hydrogen cross term
      i = 109
      j = 110
      djay = (503.2d0)*(276.0d0)
      ij = (i-1)*nntype + j
      ji = (j-1)*nntype + i
      jayq(ij) = djay
      jayq(ji) = djay
      
! --- SPC-FQ water Hydrogen-Hydrogen cross term
      i = 110
      j = 110
      djay = (503.2d0)*(196.0d0)
      ij = (i-1)*nntype + j
      jayq(ij) = djay

! --- TIP4P water Charge-Hydrogen cross term
      i = 112
      j = 113
      djay = (503.2d0)*(286.4d0)
      ij = (i-1)*nntype + j
      ji = (j-1)*nntype + i
      jayq(ji) = djay
      jayq(ij) = djay

! --- TIP4P water Hydrogen-Hydrogen cross term
      i = 112
      j = 112
      djay = (503.2d0)*(203.6d0)
      ij = (i-1)*nntype + j
      jayq(ij) = djay

! --- Methane C-H cross term
      i = 28
      j = 29
      ij = (i-1)*nntype + j
      ji = (j-i)*nntype + i
      jayq(ji) = 114855.0d0
      jayq(ij) = 114855.0d0

! --- Methane H-H cross term
      i = 29
      j = 29
      ij = (i-1)*nntype + j
      jayq(ij) = 112537.0d0

! *** convert input data to program units ***
      if ( lsami ) then
      do ibox = 2,nbox
         if (dabs(rcut(1)-rcut(ibox)).gt.1.0d-10) then
            call cleanup('Keep rcut for each box same')
         end if
      end do
         call susami
         rcheck = 2.5d0 * 3.527d0
         if ( rcut(1) .ne. rcheck ) then
            write(iou,*) 'WARNING ### rcut set to 2.5sigma for SAMI'
            rcut(1) = rcheck
         end if
      else
! *** calculate square sigmas and epsilons for lj-energy subroutines ***

         do i = 1, nntype
            do j = 1, nntype
               ij = (i-1)*nntype + j
               if ( lspecial(ij) ) then
                  write(iou,*) 'ij,lspecial(ij)',ij,lspecial(ij)
                  adum = aspecial(ij)
                  bdum = bspecial(ij)
               else
                  adum = 1.0d0
                  bdum = 1.0d0
               end if
! --- Lorentz-Berthelot rules --- sig_ij = 0.5 [ sig_i + sig_j ]
               if ( lmixlb ) then
                  sig2ij(ij) =(adum* 0.5d0 * ( sigi(i) + sigi(j) ) )**2
                  if ( sigi(i) .eq. 0.0d0 .or. sigi(j) .eq. 0.0d0 )
     &                 sig2ij(ij) = 0.0d0
                  epsij(ij) = bdum*dsqrt( epsi(i) * epsi(j) )
               end if
! --- Jorgensen mixing rules --- sig_ij = [ sig_i * sig_j ]^(1/2)
               if(lmixjo) then
                  sig2ij(ij) = adum*adum*sigi(i) * sigi(j)
                  epsij(ij) = bdum*dsqrt( epsi(i) * epsi(j) )
               end if

	       if (lshift) then
                   sr2 = sig2ij(ij) / (rcut(1)*rcut(1))
                   sr6 = sr2 * sr2 * sr2
                   ecut(ij)= sr6*(sr6-1.0d0)*epsij(ij)
	       end if

            end do
         end do
      end if

      read(5,*)
      read(5,*) nmix 
      read(5,*)
      do imix=1,nmix
         read(5,*) i,j,sigmaTmp,epsilonTmp
         ij=(i-1)*nntype+j
         sig2ij(ij)=sigmaTmp*sigmaTmp
         sig2ij((j-1)*nntype+i)=sig2ij(ij)
         epsij(ij)=epsilonTmp
         epsij((j-1)*nntype+i)=epsilonTmp
      end do

! ---  Conversion factor for intermolecular coulomb interactions
      qqfact = 1.67d5

      return
      end





