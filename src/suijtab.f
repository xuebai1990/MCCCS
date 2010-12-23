      subroutine suijtab(lmixlb,lmixjo,ltab)

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'external.inc'
!$$$      include 'poten.inc'
!$$$      include 'system.inc'
!$$$      include 'cell.inc'
!$$$      include 'expsix.inc'
!$$$      include 'merck.inc'
!$$$      include 'nsix.inc'
!$$$!kea
!$$$      include 'garofalini.inc'
!$$$      include 'conver.inc'

      logical::lmixlb, lmixjo,ltab
      integer(KIND=int)::i, j, ij, ji,ibox,nntype5,nmix,imix
      real(KIND=double_precision)::rzeronx(nxatom),epsilonnx(nxatom)

      real(KIND=double_precision)::rcheck, sr2, sr6, adum, bdum, rs1, rs7, sr7,
     &     pi,djay,sigmaTmp,epsilonTmp

! ----------------------------------------------------------------
      do i = 1, nntype
         lpl(i) = .false.
      end do

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

! * ALKANES

! --- SKS-UA methyl group [CH3] standard (nature paper)
      sigi(1) = 3.93d0
      epsi(1) = 114.0d0
      mass(1) = 15.0347d0
      chname(1) = ' SKS-UA CH3 alkane'
      chemid(1)  = 'C  '        

! --- SKS-UA methylene group [CH2]
      sigi(2) = 3.93d0 
      epsi(2) = 47.0d0 
      mass(2) = 14.0268d0
      chname(2) = ' SKS-UA CH2 alkane'
      chemid(2)  = 'C  '


! --- TraPPE-UA Methane [CH4] sp3 (similar to OPLS) 
      sigi(3) = 3.73d0
      epsi(3) = 148.0d0
      mass(3) = 16.043d0
      chname(3) = ' Tr-UA CH4 alkane '
      chemid(3)  = 'C  '

! --- TraPPE-UA Methyl [CH3] sp3 (J. Phys. Chem.) (primary)
      sigi(4) = 3.75d0
      epsi(4) = 98.0d0
      mass(4) = 15.0347d0
      chname(4) = ' Tr-UA CH3 alkane '
      chemid(4)  = 'C  '

! --- TraPPE-UA Methylene [CH2] sp3 (secondary)
      sigi(5) = 3.95d0
      epsi(5) = 46.0d0
      mass(5) = 14.0268d0
      chname(5) = ' Tr-UA CH2 alkane '
      chemid(5)  = 'C  '

! --- TraPPE-UA Methine [CH] sp3 (ternary)
      sigi(6) = 4.68d0
      epsi(6) = 10.0d0
      mass(6) = 13.0191d0
      chname(6) = ' Tr-UA CH  alkane '
      chemid(6)  = 'C  '

! --- TraPPE-UA [C] (quaternary)  
      sigi(7) = 6.4d0 
      epsi(7) = 0.5d0
      mass(7) = 12.011d0
      chname(7) = ' Tr-UA C   alkane '
      chemid(7)  = 'C  '
!ccccccccccccccccccccccccccccccccccccccccccccccc
! --- coarse-grain end segment (CH3+CH2+CH2)
! --- sigi and epsi don't matter, just need to
! --- be read in so readdat doesn't get confused
      sigi(40) = 3.75d0
      epsi(40) = 98.0d0
      mass(40) = 43.0883d0
      chname(40) = ' Tr-UA CH3 alkane '
      chemid(40)  = 'C  '

! --- coarse-grain middle segement (CH2+CH2+CH2)
      sigi(50) = 3.95d0
      epsi(50) = 46.0d0
      mass(50) = 42.0804d0
      chname(50) = ' Tr-UA CH2 alkane '
      chemid(50)  = 'C  '
!ccccccccccccccccccccccccccccccccccccccccccccccc

! --- OPLS-UA ethane methyl [CH3] sp3
      sigi(8) = 3.775d0
      epsi(8) = 104.1d0
      mass(8) = 15.0347d0
      chname(8) = ' OPLSUA CH3 ethane'
      chemid(8)  = 'C  '

! --- OPLS-UA butane methyl [CH3] sp3
      sigi(9) = 3.905d0
      epsi(9) = 88.1d0
      mass(9) = 15.0347d0
      chname(9) = ' OPLSUA CH3 butane'
      chemid(9)  = 'C  '

! --- OPLS-UA methylene [CH2] sp3
      sigi(10) = 3.905d0
      epsi(10) = 59.4d0
      mass(10) = 14.0268d0
      chname(10) = ' OPLSUA CH2 alkane'
      chemid(10)  = 'C  '

! --- OPLS-UA Methine [CH] (ternary) CARBON GROUP Jorgensen 
      sigi(11) = 3.85d0
      epsi(11) = 32.0d0
      mass(11) = 13.0191d0
      chname(11) = ' OPLSUA CH alkane '
      chemid(11)  = 'C  '

! --- OPLS [C]  (quaternary) CARBON GROUP Jorgensen
      sigi(12) = 3.85d0
      epsi(12) = 25.0d0
      mass(12) = 12.0113d0
      chname(12) = ' OPLSUA C  alkane '
      chemid(12)  = 'C  '


! --- UA Methane [CH4] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(13) = 4.10d0
      epsi(13) = 140.0d0
      mass(13) = 16.043d0
      chname(13) = ' FreireUA CH4     '
      chemid(13)  = 'C  '

! --- UA Methyl [CH3] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(14) = 4.02d0
      epsi(14) = 96.0d0
      mass(14) = 15.0347d0
      chname(14) = ' FreireUA CH3 alkn'
      chemid(14)  = 'C  '

! --- UA Methylene [CH2] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(15) = 3.72d0
      epsi(15) = 57.0d0
      mass(15) = 14.0268d0
      chname(15) = ' FreireUA CH2 alkn'
      chemid(15)  = 'C  '

! --- UA Methine [CH] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(16) = 3.36d0
      epsi(16) = 36d0
      mass(16) = 13.019d0
      chname(16) = ' FreireUA CH alkn '
      chemid(16)  = 'C  '

! --- UA Quaternary [C] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(17) = 2.44d0
      epsi(17) = 9.0d0
      mass(17) = 12.011d0
      chname(17) = ' FreireUA C alkane'
      chemid(17)  = 'C  '


! --- Mol. Phys. UA Methyl [CH3] sp3
      sigi(18) = 3.77d0
      epsi(18) = 98.1d0
      mass(18) = 15.0347d0
      chname(18) = ' MPhysUA CH3 alkM '
      chemid(18)  = 'C  '

! --- Mol. Phys. UA METHYL-BRANCH Methyl [CH3] no tail correction
      sigi(19) = 3.93d0
      epsi(19) = 78.0d0
      mass(19) = 15.0347d0
      chname(19) = ' MPhysUA CH3 alkMB'
      chemid(19)  = 'C  '

! --- Mol. Phys. UA ETHYL-BRANCH Methyl [CH3] no tail correction
      sigi(20) = 3.93d0
      epsi(20) = 95.0d0
      mass(20) = 15.0347d0
      chname(20) = ' MPhysUA CH3 alkEB'
      chemid(20)  = 'C  '

! --- Mol. Phys. UA Methine [CH] sp3 (ternary)
      sigi(21) = 4.10d0
      epsi(21) = 12.0d0
      mass(21) = 13.0191d0
      chname(21) = ' MPhysUA CH alkane'
      chemid(21)  = 'C  '


! --- TraPPE-AA for alkane methane [C]-H4 carbon
      sigi(22) = 3.31d0
      epsi(22) = 0.01d0
      mass(22) = 12.011d0
      chname(22) = ' Tr-AA [C]H4 alkan'
      chemid(22)  = 'C  '

! --- TraPPE-AA for alkane methyl [C]-H3 carbon
      sigi(23) = 3.30d0
      epsi(23) = 4.0d0
      mass(23) = 12.011d0
      chname(23) = ' Tr-AA [C]H3 alkan'
      chemid(23)  = 'C  '

! --- TraPPE-AA for alkane methylene [C]-H2 carbon 
      sigi(24) = 3.65d0
      epsi(24) = 5.0d0
      mass(24) = 12.011d0
      chname(24) = ' Tr-AA [C]H2 alkan'
      chemid(24)  = 'C  '

! --- (TraPPE?)-AA for alkane methine [C]-H carbon
      sigi(25) = 4.0d0
      epsi(25) = 2.0d0
      mass(25) = 12.011d0
      chname(25) = ' Tr-AA [C]H alkane'
      chemid(25)  = 'C  '

! --- (TraPPE?)-AA for alkane quaternary [C] carbon
      sigi(26) = 4.35d0
      epsi(26) = 1.0d0
      mass(26) = 12.011d0
      chname(26) = ' Tr-AA C quat alkn'
      chemid(26)  = 'C  '

! --- TraPPE-AA for alkane carbon-hydrogen sigma bond H[-]C 
      sigi(27) = 3.31d0
      epsi(27) = 15.3d0
!      mass(27) = 0.0d0
      mass(27) = 1.0079d0
      chname(27) = ' Tr-AA H alkane   '
      chemid(27)  = 'H  '


! --- TraPPE-UA? Methane [CH4] sp3 charged with polarizability  
      sigi(28) = 3.73d0
      epsi(28) = 148.0d0
! is this correct?
      mass(28) = 16.043d0
      qelect(28) = -0.572d0
      lqchg(28) = .true.
      jayself(28) = 0.5d0*117403d0
      xiq(28) = 9449.3d0
      chname(28) = ' Tr C CH4 chg pol '
      chemid(28)  = 'C  '

! --- Methane hydrogen charged with polarizibility
      sigi(29) = 0.0d0
      epsi(29) = 0.0d0
      mass(29) = 1.0078d0
      qelect(29) = 0.143d0
      lqchg(29) = .true.
      jayself(29) = 0.5d0*177700d0
      xiq(29) = 0.0d0
      lij(29) = .false.
      chname(29) = ' Tr H CH4 chg pol '
      chemid(29)  = 'H  '


! --- OPLS-AA for alkane methylene [C]-H2 carbon 
      sigi(30) = 3.50d0
      epsi(30) = 33.2d0
      mass(30) = 12.011d0
      qelect(30) = -0.12d0
      lqchg(30) = .true.
      chname(30) = ' OPLSAA [C]H2 alkn'
      chemid(30)  = 'C  '

! --- OPLS AA for alkane hydrogen
      sigi(31) = 2.50d0
      epsi(31) = 15.1d0
      mass(31) = 1.0078d0
      qelect(31) = 0.06d0
      lqchg(31) = .true.
      chname(31) = ' OPLSAA H  alkane '
      chemid(31)  = 'H  '


! --- Tildesly explicit atom methyl [C]-H3 carbon
      sigi(32) = 3.367d0
      epsi(32) = 48.8d0
      mass(32) = 12.011d0
      chname(32) = ' TildAA [C]H3 alkn'
      chemid(32)  = 'C  '

! --- Tildesly explicit atom methylene [C]-H2 carbon
      sigi(33) = 3.367d0
      epsi(33) = 48.8d0
      mass(33) = 12.011d0
      chname(33) = ' TildAA [C]H2 alkn'
      chemid(33)  = 'C  '

! --- Tildesly explicit atom [H] hydrogen
      sigi(34) = 2.908d0
      epsi(34) = 6.84d0
      mass(34) = 1.0079d0
      chname(34) = ' TildAA H alkane  '
      chemid(34)  = 'H  '


! --- Lennard-Jonesium ethane [CH3-CH3]
      sigi(35) = 4.25d0
      epsi(35) = 236.0d0
      mass(35) = 30.070d0
      chname(35) = ' LJ ethane C2H6   '
      chemid(35)  = 'C  '

! --- Lennard-Jonesium heptane [CH3-(CH2)5-CH3]
      sigi(36) = 6.08d0
      epsi(36) = 418.0d0
      mass(36) = 100.203d0
      chname(36) = ' LJ heptane C7H16 '
      chemid(36)  = 'C  '

! --- SPECIAL LJ CHAIN fit to give phase diagram of octane J Phys Chem 98?
      sigi(37) = 2.91d0
      epsi(37) = 236.0d0
      mass(37) = 14.0268d0
      chname(37) = ' LJ CH2 octane fit'
      chemid(37)  = 'C  '

! --- Teja heptane at 366 K [CH3-(CH2)5-CH3]
      sigi(38) = 6.0471d0
      epsi(38) = 484.76d0
      mass(38) = 100.203d0
      chname(38) = ' Teja heptane 366K'
      chemid(38)  = 'C  '

! --- Teja heptane at 450 K [CH3-(CH2)5-CH3]
      sigi(39) = 6.0471d0
      epsi(39) = 456.82d0
      mass(39) = 100.203d0
      chname(39) = ' Teja heptane 450K'
      chemid(39)  = 'C  '

! * PERFLUOROALKANES

!$$$c --- perfluoromethane [CF4] 
!$$$      sigi(40) = 4.13d0
!$$$c      sigi(40) = 4.18d0 (fit for critical density)
!$$$      epsi(40) = 172.95d0
!$$$      mass(40) = 88.003d0
!$$$      chname(40) = ' UA CF4           '
!$$$      chemid(40)  = 'C  '
!$$$
!$$$c --- Bin's perfluoromethane [CF4] 
!$$$      sigi(41) = 4.15d0
!$$$      epsi(41) = 175.4d0
!$$$      mass(41) = 88.003d0
!$$$      chname(41) = ' Bin UA CF4       '
!$$$      chemid(41)  = 'C  '
!$$$
!$$$c --- TraPPE-UA (ilja email 4-14-99) [CF3] group
!$$$      sigi(42) = 4.36d0
!$$$      epsi(42) = 87.0d0
!$$$      mass(42) = 69.0065d0
!$$$      chname(42) = ' TrUA CF3 Ilja    '
!$$$      chemid(42)  = 'C  '
!$$$
!$$$c --- [CF3] group iterb      
!$$$      sigi(43) = 4.35d0
!$$$      epsi(43) = 87.0d0
!$$$      mass(43) = 69.006d0
!$$$      chname(43) = ' UA CF3 iterb     '
!$$$      chemid(43)  = 'C  '
!$$$
!$$$c --- TraPPE-UA (ilja email 4-14-99) [CF2] group
!$$$      sigi(44) = 4.73d0
!$$$      epsi(44) = 27.5d0
!$$$      mass(44) = 50.0081d0
!$$$      chname(44) = ' TrUA CF2 Ilja    '
!$$$      chemid(44)  = 'C  '
!$$$
!$$$
!$$$c --- Amber-AA for [C]F4 carbon (JCC 13(1992) P963)
!$$$      sigi(45) = 3.82d0/(2.0d0**(1.0d0/6.0d0))
!$$$      epsi(45) = 55.05d0
!$$$      mass(45) = 12.011d0
!$$$      qelect(45) = -0.756d0
!$$$      lqchg(45) = .true.
!$$$      chname(45) = ' AmberAA [C]F4    '
!$$$      chemid(45)  = 'C  '
!$$$
!$$$c --- Amber-AA for C[F]4 fluorine (JCC 13(1992) P963)
!$$$      sigi(46) = 3.50d0/(2.0d0**(1.0d0/6.0d0))
!$$$      epsi(46) = 30.70d0
!$$$      mass(46) = 18.9984d0
!$$$      qelect(46) = 0.189d0
!$$$      lqchg(46) = .true.
!$$$      chname(46) = ' AmberAA C[F]4    '
!$$$      chemid(46)  = 'F  '
!$$$
!$$$c --- AA for [C]F4 carbon (Surface Science 367(1996) P177)
!$$$      sigi(47) = 3.35d0
!$$$      epsi(47) = 32.73d0
!$$$      mass(47) = 12.011d0
!$$$c      qelect(47) = -0.808d0
!$$$      chname(47) = ' SurfSciAA [C]F4  '
!$$$      chemid(47)  = 'C  '
!$$$
!$$$c --- AA for C[F]4 fluorine (Surface Science 367(1996) P177)
!$$$c      sigi(48) = 2.95d0
!$$$c      epsi(48) = 37.0d0
!$$$      sigi(48) = 2.90d0
!$$$      epsi(48) = 34.3d0
!$$$      mass(48) = 18.9984d0
!$$$c      qelect(48) = 0.202d0
!$$$      chname(48) = ' SurfSciAA C[F]4  '
!$$$      chemid(48)  = 'F  '
!$$$
!$$$c --- AA for [C]F4 carbon (Nose and Klein J.Chem.Phys. 78(1983) 6928)
!$$$c      sigi(49) = 3.35d0
!$$$c      epsi(49) = 37.00d0
!$$$      sigi(49) = 3.35d0
!$$$      epsi(49) = 26.00d0
!$$$      mass(49) = 12.011d0
!$$$c      qelect(49) = -0.896d0
!$$$      chname(49) = ' NoseKleinAA [C]F4'
!$$$      chemid(49)  = 'C  '
!$$$
!$$$c --- AA for C[F]4 fluorine (Nose and Klein J.Chem.Phys. 78(1983) 6928)
!$$$      sigi(50) = 2.95d0
!$$$      epsi(50) = 38.50d0
!$$$      mass(50) = 18.9984d0
!$$$c      qelect(50) = 0.224d0
!$$$      chname(50) = ' NoseKleinAA C[F]4'
!$$$      chemid(50)  = 'F  '

! * ALKENES

! --- TraPPE-UA [CH2] sp2 alkene Try 2D 04-15-99 MGM
      sigi(51) = 3.675d0
      epsi(51) = 85.0d0
      mass(51) = 14.0269d0
      chname(51) = ' Tr-UA CH2 alkene '
      chemid(51)  = 'C  '

! --- TraPPE-UA [CH] sp2 alkene Try 3B 04-15-99 MGM
      sigi(52) = 3.73d0
      epsi(52) = 47.0d0
      mass(52) = 13.0191d0
      chname(52) = ' Tr-UA CH alkene  '
      chemid(52)  = 'C  '

! --- TraPPE-UA [C] sp2 Try A 04-15-99 MGM
      sigi(53) = 3.85d0
      epsi(53) = 20.0d0
      mass(53) = 12.011d0
      chname(53) = ' Tr-UA C alkene   '
      chemid(53)  = 'C  '


! --- OPLS-UA sp2 hybrid [CH2] group JACS 106, 6638-6646 (1984)
      sigi(54) = 3.85d0
      epsi(54) = 70.43d0
      mass(54) = 14.0269d0
      chname(54) = ' OPLSAA CH2 sp2   '
      chemid(54)  = 'C  '

! --- OPLS-UA sp2 hybrid [CH] group JACS 106, 6638-6646 (1984)
      sigi(55) = 3.800d0
      epsi(55) = 57.85d0
      mass(55) = 13.0191d0
      chname(55) = ' OPLSAA CH sp2    '
      chemid(55)  = 'C  '

! * AROMATICS

! --- TraPPE-UA [CH] benzene carbon
! *** maybe these three are for UA 9-site model?!?!?!?!?
      sigi(56) = 3.74d0
      epsi(56) = 48.0d0
      mass(56) = 13.0191d0
! * published CH(aro) for TraPPE-UA 6-site
!      sigi(56) = 3.695d0
!      epsi(56) = 50.5d0
!      mass(56) = 13.0191d0
      chname(56) = ' Tr-UA CH benzene6'
      chemid(56)  = 'C  '

! --- TraPPE-UA middle benzene site
      sigi(57) = 0.0d0
      epsi(57) = 0.0d0
      mass(57) = 0.0d0
      lqchg(57) = .true.
      qelect(57) = 2.42d0
      chname(57) = ' Tr-UA mid-q benz9'
      chemid(57)  = 'H  '
      
! --- TraPPE-UA pi electron benzene site
      sigi(58) = 0.0d0
      epsi(58) = 0.0d0
      mass(58) = 0.0d0
      lqchg(58) = .true.
      qelect(58) = -1.21d0
      chname(58) = ' Tr-UA pi-q benz9 '
      chemid(58)  = 'H  '
      
! --- TraPPE-UA [C] tertiary aromatic carbon for toluene
      sigi(59) = 3.88d0
      epsi(59) = 21.0d0
      mass(59) = 12.011d0
      chname(59) = ' Tr-UA C arom tolu'
      chemid(59)  = 'C  '

! --- TraPPE-UA [C] tertiary aromatic carbon for napthalene
      sigi(60) = 3.70d0
      epsi(60) = 30.0d0  
      mass(60) = 12.011d0
      chname(60) = ' Tr-UA C arom naph'
      chemid(60)  = 'C  '

! * ALCOHOLS      

! --- TraPPE-UA alkanol hydrogen [H]-O
      sigi(61) = 0.0d0
      epsi(61) = 0.0d0
      mass(61) = 1.0079d0
      qelect(61) = 0.435d0
      lij(61) = .false.
      lqchg(61) = .true.
      chname(61) = ' Tr-UA H alkanol  '
      chemid(61)  = 'H  '
      
! --- TraPPE-UA alkanol oxygen H-[O]-CHx
      sigi(62) = 3.02d0
      epsi(62) = 93.0d0
      mass(62) = 15.999d0
      qelect(62) = -0.700d0
      lqchg(62) = .true.
      chname(62) = ' Tr-UA O alkanol  '
      chemid(62)  = 'O  '

! --- TraPPE-UA methanol methyl [CH3]-OH 
      sigi(63) = sigi(4)
      epsi(63) = epsi(4)
      mass(63) = mass(4)
      qelect(63) = 0.265d0
      lqchg(63) = .true.
      chname(63) = ' Tr-UA CH3 alkanol'
      chemid(63)  = 'C  '

! --- TraPPE-UA alkanol methylene [CH2]-OH
      sigi(64) = sigi(5)
      epsi(64) = epsi(5)
      mass(64) = mass(5)
      qelect(64) = 0.265d0
      lqchg(64) = .true.
      chname(64) = ' Tr-UA CH2 alkanol'
      chemid(64)  = 'C  '

! --- TraPPE-UA alkanol methine [CH]-OH (from Bin 6-20-00)
      sigi(65) = 4.33d0
      epsi(65) = epsi(6)
      mass(65) = mass(6)
      qelect(65) = 0.265d0
      lqchg(65) = .true.
      chname(65) = ' Tr-UA CH alkanol '
      chemid(65)  = 'C  '

! --- TraPPE-UA alkanol quaternary carbon [C]-OH (from Bin 6-20-00)
      sigi(66) = 5.8d0
      epsi(66) = epsi(7)
      mass(66) = mass(7)
      qelect(66) = 0.265d0
      lqchg(66) = .true.
      chname(66) = ' Tr-UA C alkanol  '
      chemid(66)  = 'C  '


! --- OPLS-UA alkanol hydrogen [H]-O
      sigi(67) = 0.0d0
      epsi(67) = 0.0d0
      mass(67) = 1.0079d0
      qelect(67) = 0.435d0
      lij(67) = .false.
      lqchg(67) = .true.
      chname(67) = ' OPLSUA H alkanol '
      chemid(67)  = 'H  '

! --- OPLS-UA alkanol oxygen H-[O]-CHx 
      sigi(68) = 3.07d0
      epsi(68) = 85.578d0
      mass(68) = 15.999d0
      qelect(68) = -0.700d0
      lqchg(68) = .true.
      chname(68) = ' OPLSUA O alkanol '
      chemid(68)  = 'O  '

! --- OPLS-UA alkanol methyl [CH3]-OH 
      sigi(69) = sigi(8)
      epsi(69) = epsi(8)
      mass(69) = mass(8)
      qelect(69) = 0.265d0
      lqchg(69) = .true.
      chname(69) = ' OPLSUA CH3 alknol'
      chemid(69)  = 'C  '
    
! --- OPLS-UA alkanol methylene [CH2]-OH
      sigi(70) = sigi(10)
      epsi(70) = epsi(10)
      mass(70) = mass(10)
      qelect(70) = 0.265d0
      lqchg(70) = .true.
      chname(70) = ' OPLSUA CH2 alknol'
      chemid(70)  = 'C  '

! * ETHERS

! --- TraPPE-UA ether oxygen  CHx-[O]-CHx
      sigi(71) = 2.80d0
      epsi(71) = 55.0d0
      mass(71) = 16.00d0
      qelect(71) = -0.50d0
      lqchg(71) = .true.
      chname(71) = ' Tr-UA O ether    '
      chemid(71)  = 'O  '

! --- TraPPE-UA ether methyl [CH3]-O
      sigi(72) = sigi(4)
      epsi(72) = epsi(4)
      mass(72) = mass(4)
      qelect(72) = 0.25d0
      lqchg(72) = .true.
      chname(72) = ' Tr-UA CH3 ether  '
      chemid(72)  = 'C  '

! --- TraPPE-UA ether methylene [CH2]-O
      sigi(73) = sigi(5)
      epsi(73) = epsi(5)
      mass(73) = mass(5)
      qelect(73) = 0.25d0
      lqchg(73) = .true.
      chname(73) = ' Tr-UA CH2 ether  '
      chemid(73)  = 'C  '

! --- TraPPE-UA ether methine [CH]-O
      sigi(74) = sigi(65)
      epsi(74) = epsi(65)
      mass(74) = mass(65)
      qelect(74) = 0.25d0
      lqchg(74) = .true.
      chname(74) = ' Tr-UA CH ether   '
      chemid(74)  = 'C  '

! --- TraPPE-UA ether quaternary carbon [C]-O
      sigi(75) = sigi(66)
      epsi(75) = epsi(66)
      mass(75) = mass(66)
      qelect(75) = 0.25d0
      lqchg(75) = .true.
      chname(75) = ' Tr-UA C ether    '
      chemid(75)  = 'C  '

! --- TraPPE-UA Block copolymer ether oxygen next to carbonyl CH2-[O]-C=O
      sigi(76) = 2.80d0
      epsi(76) = 55.0d0
      mass(76) = 16.00d0
      qelect(76) = -0.25d0
      lqchg(76) = .true.
      chname(76) = ' Tr-UA O carbonate'
      chemid(76)  = 'O  '

! --- TraPPE-UA Block copolymer methylene O=C-O-[CH2]
      sigi(77) = sigi(5)
      epsi(77) = epsi(5)
      mass(77) = mass(5)
      qelect(77) = 0.30d0
      lqchg(77) = .true.
      chname(77) = ' Tr-UA CH2 carbnat'
      chemid(77)  = 'C  '


! --- OPLS-UA ether oxygen CHx-[O]-CHx (JCC 1990 vol 11, iss 8 958-971)
      sigi(78) = 3.00d0
      epsi(78) = 85.58d0
      mass(78) = 15.999d0
      qelect(78) = -0.50d0
      lqchg(78) = .true.
      chname(78) = ' OPLSUA O ether   '
      chemid(78)  = 'O  '

! --- OPLS-UA ether methyl [CH3]-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(79) = 3.80d0
      epsi(79) = 85.58d0
      mass(79) = 15.0347d0
      qelect(79) = 0.25d0
      lqchg(79) = .true.
      chname(79) = ' OPLSUA [CH3]-O   '
      chemid(79)  = 'C  '

! --- OPLS-UA ether methyl [CH3]-CH2-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(80) = 3.905d0
      epsi(80) = 88.06d0
      mass(80) = 15.0347d0
      chname(80) = ' OPLSUA [CH3]CH2-O'
      chemid(80)  = 'C  '

! --- OPLS-UA ether methylene [CH2]-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(81) = 3.80d0
      epsi(81) = 59.38d0
      mass(81) = 14.0268d0
      qelect(81) = 0.25d0
      lqchg(81) = .true.
      chname(81) = ' OPLSUA [CH2]-O   '
      chemid(81)  = 'C  '

! --- OPLS-UA THF methylene [CH2]-CH2-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(82) = 3.905d0
      epsi(82) = 59.38d0
      mass(82) = 14.0268d0
      chname(82) = ' OPLSUA [CH2]CH2-O'
      chemid(82)  = 'C  '

! * KETONES, ALDEHYDES AND ESTERS

!$$$c --- (TraPPE?)-UA ketone carbon [C]=O (jpotoff 12/13/99 + OPLS JPC v94 p1683 1990)
!$$$      sigi(83) = 3.82d0
!$$$      epsi(83) = 40.00d0
!$$$      mass(83) = 12.011d0
!$$$      qelect(83) = +0.424d0
!$$$      lqchg(83) = .true.
! --- TraPPE-UA ketone carbon [C]=O (TraPPE-6)
      sigi(83) = 3.82d0
      epsi(83) = 40.00d0
      mass(83) = 12.011d0
      qelect(83) = +0.424d0
      lqchg(83) = .true.
      chname(83) = ' Tr-UA [C]=O keton'
      chemid(83)  = 'C  '

!$$$c --- (TraPPE?)-UA ketone oxygen C=[O] (jpotoff 12/17/99 + OPLS JPC v94 p 1683 1990)
!$$$      sigi(84) = 3.04d0
!$$$      epsi(84) = 85.0d0
!$$$      mass(84) = 15.999d0
!$$$      qelect(84) = -0.424d0 
!$$$      lqchg(84) = .true. 
! --- TraPPE-UA ketone oxygen C=[O] (TraPPE-6)
      sigi(84) = 3.05d0
      epsi(84) = 79.0d0
      mass(84) = 15.999d0
      qelect(84) = -0.424d0 
      lqchg(84) = .true. 
      chname(84) = ' Tr-UA C=[O] keton'
      chemid(84)  = 'O  '

!$$$c --- (TraPPE?)-UA aldehyde carbon [CH]=O  (jpotoff 12/13/99 + OPLS JPC v94 p1683 1990)
!$$$      sigi(85) = 3.60d0
!$$$      epsi(85) = 58.00d0
!$$$      mass(85) = 13.011d0
!$$$      qelect(85) = +0.424d0
!$$$      lqchg(85) = .true.
! --- TraPPE-UA aldehyde carbon [CH]=O (TraPPE-6) 
      sigi(85) = 3.55d0
      epsi(85) = 65.00d0
      mass(85) = 13.019d0
      qelect(85) = +0.424d0
      lqchg(85) = .true.
      chname(85) = ' Tr-UA [CH]=O alde'
      chemid(85)  = 'C  '
      
! --- TraPPE-UA ester methylene group [CH2]-C=O 
      sigi(86) = sigi(5)
      epsi(86) = epsi(5)
      qelect(86) = 0.05d0
      mass(86) = mass(5)
      lqchg(86) = .true.
      chname(86) = ' Tr-UA [CH2]-C=O e'
      chemid(86)  = 'C  '

! --- TraPPE-UA ester methylene group C(=O)O-[CH2]
      sigi(87) = sigi(5)
      epsi(87) = epsi(5)
      qelect(87) = 0.25d0
      mass(87) = mass(5)
      lqchg(87) = .true.
      chname(87) = ' Tr-UA COO-[CH2] e'
      chemid(87)  = 'C  '

! --- TraPPE-UA ester oxygen C(=O)-[O]-CHx (uses TraPPE alcohol O)
      sigi(88) = sigi(62)
      epsi(88) = epsi(62)
      qelect(88) = -0.40d0
      mass(88) = mass(62)
      lqchg(88) = .true.
      chname(88) = ' Tr-UA CO[O]-CHx e'
      chemid(88)  = 'O  '

! --- TraPPE-UA ester oxygen in carbonyl C=[O] (uses TraPPE CO2 O)
      sigi(89) = 3.05d0
      epsi(89) = 79.0d0
      qelect(89) = -0.45d0
      mass(89) = 15.999d0
      lqchg(89) = .true.
      chname(89) = ' Tr-UA C=[O] ester'
      chemid(89)  = 'O  '

! --- TraPPE-UA ester carbon in carbonyl [C]=O 
      sigi(90) = 3.82d0
      epsi(90) = 40.0d0
      qelect(90) = 0.55d0
      mass(90) = 12.011d0
      lqchg(90) = .true.
      chname(90) = ' Tr-UA [C]=O ester'
      chemid(90)  = 'O  '

! * CARBOXYLIC ACIDS

! --- 91-94 are old parameters for TraPPE-UA carboxylic acids
! --- TraPPE-UA carboxylic acid hydrogen C(=O)-O-[H]
      sigi(91) = 0.0d0
      epsi(91) = 0.0d0
      mass(91) = 1.0079d0
      qelect(91) = 0.30d0
      lij(91) = .false.
      lqchg(91) = .true.
      chname(91) = 'oTr-UA COO[H] acid'
      chemid(91)  = 'H  '

! --- TraPPE iterB carbonyl oxygen  C[=O]-O-H
!      sigi(92) = 3.0d0
!      epsi(92) = 75.0d0
!      qelect(92) = -0.440d0
      sigi(92) = 3.04d0
      epsi(92) = 81.0d0
      qelect(92) = -0.424d0
      mass(92) = 15.999d0
      lqchg(92) = .true.
      chname(92) = 'oTr-UA C[O]OH acid'
      chemid(92)  = 'O  '

! --- TraPPE iterB carboxylic acid oxygen C(=O)-[O]-H
!      sigi(93) = 3.00d0
!      epsi(93) = 75.0d0
!      qelect(93) = -0.53d0
      sigi(93) = 3.02d0
      epsi(93) = 93.0d0
      qelect(93) = -0.30d0
      mass(93) = 15.999d0
      lqchg(93) = .true.
      chname(93) = 'oTr-UA CO[O]H acid'
      chemid(93)  = 'O  '

! --- TraPPE iterB carbonyl carbon  [C](=O)-O-H
!      sigi(94) = 4.0d0
!      epsi(94) = 42.0d0
!      qelect(94) = 0.52d0
      sigi(94) = 3.60d0
      epsi(94) = 65.0d0
      mass(94) = 12.011d0
      qelect(94) = 0.424d0
      lqchg(94) = .true.
      chname(94) = 'oTr-UA [C]OOH acid'
      chemid(94)  = 'C  '

! --- 95-98 new parameters for TraPPE-UA carboxylic acids 
! --- (jpotoff 12/17/99 + OPLS JPC v95 p 3315 1991)
! --- TraPPE-UA carboxylic acid carbonyl oxygen C=[O]-O-H 
      sigi(95) = 3.04d0
      epsi(95) = 81.0d0
      mass(95) = 15.999d0
      qelect(95) = -0.424d0 
      lqchg(95) = .true. 
      chname(95) = ' Tr-UA C[O]OH acid'
      chemid(95)  = 'O  '
      
! --- TraPPE-UA iterB carbonyl carbon  [C](=O)-O-H
      sigi(96) = 3.60d0
      epsi(96) = 65.00d0
      mass(96) = 12.011d0
      qelect(96) = 0.424d0
      lqchg(96) = .true.
      chname(96) = ' Tr-UA [C]OOH acid'
      chemid(96)  = 'C  '
      
! --- TraPPE-UA carboxylic acid hydrogen C(=O)-O-[H] (JPC v95 p. 3315, 1991)
      sigi(97) = 0.0d0
      epsi(97) = 0.0d0
      mass(97) = 1.0079d0
      qelect(97) = 0.30d0
      lij(97) = .false.
      lqchg(97) = .true.
      chname(97) = ' Tr-UA COO[H] acid'
      chemid(97)  = 'H  '

! --- TraPPE-UA carboxylic acid oxygen C(=O)-[O]-H
      sigi(98) = 3.02d0
      epsi(98) = 93.0d0
      mass(98) = 16.00d0
      qelect(98) = -0.30d0
      lqchg(98) = .true.
      chname(98) = ' Tr-UA CO[O]H acid'
      chemid(98)  = 'O  '

 
! --- OPLS-UA (1990) charged methyl [CH3]
      sigi(99) = 3.91d0
      epsi(99) = 80.6d0
      mass(99) = 15.0347d0
      qelect(99) = 0.080d0
      lqchg(99) = .true.
      chname(99) = ' OPLSUA CH3 acid? '
      chemid(99)  = 'C  '

! --- OPLS-UA (1990) carboxylic acid carbon [C](=O)-O-H
      sigi(100) = 3.75d0
      epsi(100) = 52.9d0
      mass(100) = 12.011d0
      qelect(100) = 0.55d0
      lqchg(100) = .true.
      chname(100) = ' OPLSUA [C]OOH acd'
      chemid(100)  = 'C  '

! --- OPLS-UA (1990) carboxylic acid oxygen C(=O)-[O]-H
      sigi(101) = 3.0d0
      epsi(101) = 85.6d0
      mass(101) = 15.999d0
      qelect(101) = -0.58d0
      lqchg(101) = .true.
      chname(101) = ' OPLSUA CO[O]H acd'
      chemid(101)  = 'O  '

! --- OPLS-UA (1990) carbonyl oxygen  C[=O]
      sigi(102) = 2.96d0
      epsi(102) = 105.7d0
      mass(102) = 15.999d0
      qelect(102) = -0.5d0
      lqchg(102) = .true.
      chname(102) = ' OPLSUA C[O]OH acd'
      chemid(102)  = 'O  '

! --- OPLS-UA (1990) and OPLS-AA (1995) carboxylic acid hydrogen C(=O)-O-[H]
      sigi(103) = 0.0d0
      epsi(103) = 0.0d0
      mass(103) = 1.0079d0
      qelect(103) = 0.45d0
      lij(103) = .false.
      lqchg(103) = .true.
      chname(103) = ' OPLSUA COO[H] acd'
      chemid(103)  = 'H  '

! --- OPLS-AA (1995) carboxylic acid carbonyl oxygen  C[=O]-O-H
      sigi(104) = 2.96d0
      epsi(104) = 105.8d0
      mass(104) = 15.999d0
      qelect(104) = -0.440d0
      lqchg(104) = .true.
      chname(104) = ' OPLSAA C[O]OH acd'
      chemid(104)  = 'O  '

! --- OPLS-AA (1995) carboxylic acid oxygen C(=O)-[O]-H
      sigi(105) = 3.00d0
      epsi(105) = 85.7d0
      mass(105) = 15.999d0
      qelect(105) = -0.53d0
      lqchg(105) = .true.
      chname(105) = ' OPLSAA CO[O]H acd'
      chemid(105)  = 'O  '

! --- OPLS-AA (1995) carboxylic acid carbon  [C](=O)-O-H
      sigi(106) = 3.75d0
      epsi(106) = 52.9d0
      mass(106) = 12.011d0
      qelect(106) = 0.52d0
      lqchg(106) = .true.
      chname(106) = ' OPLSAA [C]OOH acd'
      chemid(106)  = 'C  '

! * WATER

! --- SPC/E oxygen [O]   (simple point charge water oxygen)
      sigi(107) = 3.1655d0
      epsi(107) = 78.1958d0
      mass(107) = 16.000d0
      qelect(107) = -0.8476d0
      lqchg(107) = .true.
      chname(107) = ' SPC/E O water    '
      chemid(107)  = 'O  '

! --- SPC/E hydrogen [H] (simple point charge Enhanced water hydrogen)
      sigi(108) = 0.0d0
      epsi(108) = 0.0d0
      mass(108) = 1.0079d0
      qelect(108) = 0.4238d0      
      lqchg(108) = .true.
      chname(108) = ' SPC/E H water    '
      chemid(108)  = 'H  '

!$$$c --- TIP3P oxygen [O] 
!$$$      sigi(47) = 3.1506d0
!$$$      epsi(47) = 76.54d0
!$$$      mass(47) = 16.000d0
!$$$      qelect(47) = -0.834d0
!$$$      lqchg(47) = .true.
!$$$      chname(47) = ' TIP3P O water    '
!$$$      chemid(47)  = 'O  '
!$$$
!$$$c --- TIP3P hydrogen [H] 
!$$$      sigi(48) = 0.0d0
!$$$      epsi(48) = 0.0d0
!$$$      mass(48) = 1.0079d0
!$$$      qelect(48) = 0.417d0      
!$$$      lqchg(48) = .true.
!$$$      chname(48) = ' TIP3P H water    '
!$$$      chemid(48)  = 'H  '

! --- SPC-FQ oxygen [O]   S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(109) = 3.176
      epsi(109) = 148.0d0
      mass(109) = 15.999d0
      qelect(109) = -0.672123708
      lqchg(109) = .true.
      xiq(109) = 36899.0d0
      jayself(109) = (0.5d0)*(503.2d0)*(367.0d0)
      chname(109) = ' SPC-FQ O water   '
      chemid(109)  = '0  '

! --- SPC-FQ hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(110) = 0.0d0
      epsi(110) = 0.0d0
      mass(110) = 1.0079d0
      qelect(110) = 0.336061854
      lij(110) = .false.
      lqchg(110) = .true.
      xiq(110) = 0.0d0
      jayself(110) = (0.5d0)*(503.2d0)*(392.2d0)
      chname(110) = ' SPC-FQ H water   '
      chemid(110)  = 'H  '


! --- TIP4P-FQ Oxygen [O] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(111) = 3.159d0
      epsi(111) = 144.1d0
!      epsi(111) = 105.0d0
      mass(111) = 15.999d0
      chname(111) = ' TIP4P-FQ O water '
      chemid(111)  = 'O  '

! --- TIP4P-FQ Hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(112) = 0.0d0
      epsi(112) = 0.0d0
      mass(112) = 1.0079d0
      qelect(112) = 0.35d0
      lij(112) = .false.
      lqchg(112) = .true.
      xiq(112) = 0.0d0
      jayself(112) = (0.5d0)*(503.2d0)*(353.0d0)
      chname(112) = ' TIP4P-FQ H water '
      chemid(112)  = 'H  '

! --- TIP4P-FQ Charge [Q] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(113) = 0.0d0
      epsi(113) = 0.0d0
      mass(113) = 0.0d0
      qelect(113) = -0.70d0
      lij(113) = .false.
      lqchg(113) = .true.
      xiq(113) = 34464.0d0
      jayself(113) = (0.5d0)*(503.2d0)*(371.6d0)
      chname(113) = ' TIP4P-FQ M water '
      chemid(113)  = 'M  '


! --- TIP-4P water model --- [O] site
!      sigi(114) = 3.15365d0
      sigi(114) = 3.154d0
      epsi(114) = 78.0d0
      qelect(114) = 0.0d0
! the following value was listed under tip-4p water oxygen as well (type 152)
!      epsi(114) = 57.91d0
      mass(114) = 15.999d0
      chname(114) = ' TIP4P O water    '
      chemid(114)  = 'O  '

! --- TIP-4P water model --- [H] site
      sigi(115) = 0.0d0
      epsi(115) = 0.0d0
      mass(115) = 1.0079d0
      qelect(115) = 0.52d0
      lij(115) = .false.
      lqchg(115) = .true.
      chname(115) = ' TIP4P H water    '
      chemid(115)  = 'H  '

! --- TIP-4P water model --- [M] site
      sigi(116) = 0.0d0
      epsi(116) = 0.0d0
      mass(116) = 0.0d0
      qelect(116) = -1.04d0
      lij(116) = .false.
      lqchg(116) = .true.
      chname(116) = ' TIP4P M water    '
      chemid(116)  = 'M  '


! --- TIP5P oxygen [O]
      sigi(117) = 3.12d0
      epsi(117) = 80.512d0
      mass(117) = 15.999d0
      qelect(117) = 0.0d0
      lij(117) = .true.
      lqchg(117) = .false.
      chname(117) = ' TIP5P O water    '
      chemid(117)  = 'O  '

! --- TIP5P hydrogen [H]
      sigi(118) = 0.0d0
      epsi(118) = 0.0d0
      mass(118) = 1.0078d0
      qelect(118) = 0.241d0
      lij(118) = .false.
      lqchg(118) = .true.
      chname(118) = ' TIP5P H water    '
      chemid(118)  = 'H  '

! --- TIP5P lone-pair [L]
      sigi(119) = 0.0d0
      epsi(119) = 0.0d0
      mass(119) = 0.0d0
      qelect(119) = -0.241d0
      lij(119) = .false.
      lqchg(119) = .true.
      chname(119) = ' TIP5P L water    '
      chemid(119)  = 'L  '


! --- Fixed Charge Water oxygen [O] site
      sigi(120) = 3.34d0
      epsi(120) = 42.0d0
      mass(120) = 15.999d0
      qelect(120) = 6.0d0
      lqchg(120) = .true.
      chname(120) = ' FixedQ O water   '
      chemid(120)  = 'O  '

! --- Fixed Charge Water hydrogen [H] site
      sigi(121) = 0.0d0
      epsi(121) = 0.0d0
      mass(121) = 1.0079d0
      qelect(121) = 1.0d0
      lij(121) = .false.
      lqchg(121) = .true.
      chname(121) = ' FixedQ H water   '
      chemid(121)  = 'H  '

! --- Fixed Charge Water carbon-oxygen bond site (???hydrogen-oxygen bond???)
      sigi(122) = 2.2d0
      epsi(122) = 15.0d0
      mass(122) = 0.0d0
      qelect(122) = -2.16d0
      lqchg(122) = .true.
      chname(122) = ' FixedQ bond water'
      chemid(122)  = 'M  '

! --- Fixed Charge Water lone pair [L] site
      sigi(123) = 0.0d0
      epsi(123) = 0.0d0
      mass(123) = 0.0d0
      qelect(123) = -1.84d0
      lij(123) = .false.
      lqchg(123) = .true.
      chname(123) = ' FixedQ L water   '
      chemid(123)  = 'M  '

! * NOBLE GASES, CARBON MONOXIDE, CARBON DIOXIDE, NITROGEN, OXYGEN, HF

! --- TraPPE Helium (7-18-97 MGM)
      sigi(124) = 3.11d0
      epsi(124) = 4.0d0
!      sigi(124) = 3.065d0 used in JACS paper 1997
!      epsi(124) = 3.95d0
      mass(124) = 4.0026d0
!      sigi(124) = 2.556d0
!      epsi(124) = 10.2d0
!      mass(124) = 4.00d0
!      sigi(124) = 0.0d0
!      epsi(124) = 0.0d0
      chname(124) = ' TraPPE helium    '
      chemid(124)  = 'HE '

! --- TraPPE Argon (7-18-97 MGM)
      sigi(125) = 3.390d0
      epsi(125) = 116.0d0
      mass(125) = 39.948d0
      chname(125) = ' TraPPE argon     '
      chemid(125)  = 'Ar '

! --- Krypton  
      sigi(126) = 3.607d0
      epsi(126) = 161.0d0
      mass(126) = 83.80d0
      chname(126) = ' TraPPE? krypton  '
      chemid(126)  = 'Kr '


! --- carbon in carbon monoxide [C]=-O
      sigi(127) = 3.75d0
      epsi(127) = 52.9d0
      mass(127) = 12.011d0
      qelect(127) = -0.019d0
      lqchg(127) = .true.	
      chname(127) = ' carbon monoxide C'
      chemid(127)  = 'C  '

! --- oxygen in carbon monoxide C=-[O]
      sigi(128) = 2.96d0
      epsi(128) = 105.7d0
      mass(128) = 15.999d0
      qelect(128) = 0.019d0
      lqchg(128) = .true.
      chname(128) = ' carbon monoxide O'
      chemid(128)  = 'O  '


! --- Jeff's Amazing TraPPE CO2 model carbon [C]O2 (jpotoff 12/13/99)
      sigi(129) = 2.80d0
      epsi(129) = 27.0d0
      mass(129) = 12.011d0
      qelect(129) = 0.70d0
      lqchg(129) = .true.
      chname(129) = ' TraPPE C in CO2  '
      chemid(129)  = 'C  '

! --- Jeff's Amazing TraPPE CO2 model oxygen C[O]2 (jpotoff 12/13/99)
      sigi(130) = 3.05d0
      epsi(130) = 79.0d0
      mass(130) = 15.999d0
      qelect(130) = -0.350d0 
      lqchg(130) = .true. 
      chname(130) = ' TraPPE O in CO2  '
      chemid(130)  = 'O  '


! --- TraPPE carbon dioxide carbon in [C]O2-fq (jpotoff 2/15/00)
      sigi(131) = 2.80d0
      epsi(131) = 28.5d0
      mass(131) = 12.011d0
      qelect(131) = 0.6512d0
      lqchg(131) = .true.
!      xiq(131) = (503.2d0)*123.2d0
      xiq(131) = 0.0d0
      jayself(131) = (0.5d0)*(503.2d0)*(233.5d0)
      chname(131) = ' Tr-FQ C in CO2   '
      chemid(131)  = 'C  '

! --- TraPPE carbon dioxide oxygen in C[O]2-fq (jpotoff 2/15/00)
      sigi(132) = 3.06d0
      epsi(132) = 80.5d0
      mass(132) = 15.999d0
      qelect(132) = -0.3256d0
      lqchg(132) = .true.
!      xiq(132) = (503.2d0)*201.56d0
      xiq(132) = 39430.75d0
      jayself(132) = (0.5d0)*(503.2d0)*(308.17d0)
      chname(132) = ' Tr-FQ O in CO2   '
      chemid(132)  = 'O  '


! --- TraPPE nitrogen [N]2 (jpotoff 12/21/99)
      sigi(133) = 3.310d0
      epsi(133) = 36.00d0
      mass(133) = 14.00674d0
      qelect(133) = -0.50d0
      lqchg(133) = .true.
      chname(133) = ' TraPPE N in N2   '
      chemid(133)  = 'N  '

! --- TraPPE nitrogen COM charge cite for N2 (jpotoff 12/21/99)
      sigi(134) = 0.0d0
      epsi(134) = 0.0d0
      mass(134) = 0.0d0
      qelect(134) = 1.0d0
      lij(134) = .false.
      lqchg(134) = .true.
      chname(134) = ' TraPPE COM in N2 '
      chemid(134)  = 'M  '


! --- Tildesley nitrogen [N]2 
      sigi(135) = 3.31d0
      epsi(135) = 37.3d0
      mass(135) = 14.00674d0
      chname(135) = ' Tild. N in N2    '
      chemid(135)  = 'N  '


! --- TraPPE oxygen [O]2  Final parameter adjust 8-5-98
      sigi(136) = 3.07d0
      epsi(136) = 49.0d0
      mass(136) = 15.999d0
      chname(136) = ' TraPPE O in O2   '
      chemid(136)  = 'O  '


! --- OPLS hydrogen fluoride (HF) fluorine H-M-[F]
      sigi(137) = 2.984d0
      epsi(137) = 75.75d0
      mass(137) = 18.9984d0
      qelect(137) = 0.725d0
      lqchg(137) = .true.
      chname(137) = ' OPLS F in HMF    '
      chemid(137)  = 'F  '

! --- OPLS HF hydrogen [H]-M-F
      sigi(138) = 0.0d0
      epsi(138) = 0.0d0
      mass(138) = 1.0078d0
      qelect(138) = 0.725d0
      lij(138) = .false.
      lqchg(138) = .true.
      chname(138) = ' OPLS H in HMF    '
      chemid(138)  = 'H  '

! --- OPLS HF M site H-[M]-F
      sigi(139) = 0.0d0
      epsi(139) = 0.0d0
      mass(139) = 0.0d0
      qelect(139) = -1.45d0
      lij(139) = .false.
      lqchg(139) = .true.
      chname(139) = ' OPLS M in HMF    '
      chemid(139)  = 'M '

! * THIOLS, THIOETHERS

! --- TraPPE-UA dimethyl sulfide methyl group [CH3]-S-CH3
      sigi(140) = sigi(4)
      epsi(140) = epsi(4)
      mass(140) = mass(4)
      qelect(140) = 0.235d0
      lqchg(140) = .true.
      chemid(140)  = 'C  '
      chname(140) = ' Tr-UA CH3 thioeth'
      
! --- TraPPE-UA dimethyl sulfide sulfur CH3-[S]-CH3  
!     (1/25/00, based on JPC v90, p6379, 1986)
      sigi(141) = 3.52d0
      epsi(141) = 158.0d0
      mass(141) = 32.07d0
      qelect(141) = -0.47d0
      lqchg(141) = .true.
      chname(141) = ' Tr-UA S thioether'
      chemid(141)  = 'S  '

! --- TraPPE-UA methyl group [CH3]-S-H
      sigi(142) = sigi(4)
      epsi(142) = epsi(4)
      mass(142) = mass(4)
      qelect(142) = 0.18d0
      lqchg(142) = .true.
      chname(142) = ' Tr-UA CH3 thiol  '
      chemid(142)  = 'C  '
      
! --- TraPPE-UA sulfur CH3-[S]-H 
!     (1/25/00, based on JPC v90, p6379, 1986)
      sigi(143) = 3.62d0
      epsi(143) = 185.0d0
      mass(143) = 32.07d0
      qelect(143) = -0.45d0
      lqchg(143) = .true.
      chname(143) = ' Tr-UA S thiol    '
      chemid(143)  = 'S  '
      
! --- TraPPE-UA hydrogen CH3-S-[H]      
      sigi(144) = 0.0d0
      epsi(144) = 0.0d0
      mass(144) = 1.0079d0
      qelect(144) = 0.27d0
      lij(144) = .false.
      lqchg(144) = .true.
      chname(144) = ' Tr-UA H thiol    '
      chemid(144)  = 'H  '
 
! --- TraPPE-UA methylene group CH3-[CH2]-S-H
      sigi(145) = sigi(5)
      epsi(145) = epsi(5)
      mass(145) = mass(5)
      qelect(145) = 0.18d0
      lqchg(145) = .true.
      chname(145) = ' Tr-UA CH2 thiol  '
      chemid(145)  = 'C  '

! * AMINES 

! --- parameters for primary amines (2/28/00) based on JPC v94, p1683, 1990)
! --- TraPPE-UA methyl amine hydrogen CH3-N-[H]-H
      sigi(146) = 0.0d0
      epsi(146) = 0.0d0
      mass(146) = 1.0079d0
      qelect(146) = 0.275d0
      lij(146) = .false.
      lqchg(146) = .true.
      chname(146) = ' Tr-UA CH3-N[H]2  '
      chemid(146)  = 'H  '

! --- TraPPE-UA methyl amine nitrogen CH3-[N]-H2
      sigi(147) = 3.31d0
      epsi(147) = 165.0d0
      mass(147) = 14.00674d0
      qelect(147) = -0.65d0
      lqchg(147) = .true.
      chname(147) = ' Tr-UA CH3-[N]H2  '
      chemid(147)  = 'N  '

! --- TraPPE-UA methyl amine methyl [CH3]-N-H2
      sigi(148) = sigi(4)
      epsi(148) = epsi(4)
      mass(148) = mass(4)
      qelect(148) = 0.10d0
      lqchg(148) = .true.
      chname(148) = ' Tr-UA [CH3]-NH2  '
      chemid(148)  = 'C  '

! --- TraPPE-UA dimethyl amine nitrogen CH3-[N]-CH3-H 
      sigi(149) = 3.31d0
      epsi(149) = 115.0d0
      mass(149) = 14.00674d0
      qelect(149) = -0.75d0
      lqchg(149) = .true.
      chname(149) = ' Tr-UA (CH3)2[N]H '
      chemid(149)  = 'N  '

! --- TraPPE-UA trimethyl amine nitrogen CH3-[N]-CH3-CH3
      sigi(150) = 3.31d0
      epsi(150) = 115.0d0
      mass(150) = 14.00674d0
      qelect(150) = -0.60d0
      lqchg(150) = .true.           
      chname(150) = ' Tr-UA (CH3)3[N]  '
      chemid(150)  = 'N  '

! * NITRILES

! --- TraPPE-UA nitrile nitrogen C=-[N]
      sigi(151) = 2.95d0
      epsi(151) = 60.0d0
      mass(151) = 14.007d0
      qelect(151) = -0.398d0
      lqchg(151) = .true.
      chname(151) = ' Tr-UA N nitrile  '
      chemid(151)  = 'N  '

! --- TraPPE-UA nitrile carbon [C]=-N
      sigi(152) = 3.55d0
      epsi(152) = 60.0d0
      mass(152) = 12.011d0
      qelect(152) = 0.129d0
      lqchg(152) = .true.
      chname(152) = ' Tr-UA C nitrile  '
      chemid(152)  = 'C  '

! --- TraPPE hydrogen cyanide hydrogen [H]-C=-N
! TRIAL VALUES
      sigi(153) = 0.0d0
      epsi(153) = 0.0d0
      mass(153) = 1.0079d0
      qelect(153) = 0.269d0
      lij(153) = .false.
      lqchg(153) = .true.
      chname(153) = ' Tr-UA H in HCN   '
      chemid(153)  = 'H  '

! --- TraPPE-UA acetonitrile methyl [CH3]-C=-N
      sigi(154) = sigi(4)
      epsi(154) = epsi(4)
      mass(154) = mass(4)
      qelect(154) = 0.269d0
      lqchg(154) = .true.
      chname(154) = ' Tr-UA CH3 nitrile'
      chemid(154)  = 'C  '

! --- TraPPE-UA alkyl nitrile methylene R-[CH2]-C=-N
      sigi(155) = sigi(5)
      epsi(155) = epsi(5)
      mass(155) = mass(5)
      qelect(155) = 0.269d0
      lqchg(155) = .true.
      chname(155) = ' Tr-UA CH2 nitrile'
      chemid(155)  = 'C  '


! --- OPLS-UA nitrile nitrogen C=-[N]
      sigi(156) = 3.20d0
      epsi(156) = 85.51d0
      mass(156) = 14.007d0
      qelect(156) = -0.430d0
      lqchg(156) = .true.
      chname(156) = ' OPLSUA N nitrile '
      chemid(156)  = 'N  '

! --- OPLS-UA nitrile carbon R-[C]=-N
      sigi(157) = 3.65d0
      epsi(157) = 75.53d0
      mass(157) = 12.011d0
      qelect(157) = 0.280d0
      lqchg(157) = .true.
      chname(157) = ' OPLSUA C nitrile '
      chemid(157)  = 'C  '

! --- OPLS-UA acetonitrile methyl [CH3]-C=-N
      sigi(158) = 3.775d0
      epsi(158) = 104.16d0
      mass(158) = 15.035d0
      qelect(158) = 0.15d0
      lqchg(158) = .true.
      chname(158) = ' OPLSUA CH3 nitril'
      chemid(158)  = 'C  '


! --- McDonald UA nitrile nitrogen R-C=-[N]
      sigi(159) = 3.3d0
      epsi(159) = 50.0d0
      mass(159) = 14.007d0
      qelect(159) = -0.398d0
      lqchg(159) = .true.
      chname(159) = ' McDUA N nitrile  '
      chemid(159)  = 'N  '

! --- McDonald UA nitrile carbon R-[C]=-N
      sigi(160) = 3.4d0
      epsi(160) = 50.0d0
      mass(160) = 12.011d0
      qelect(160) = 0.129d0
      lqchg(160) = .true.
      chname(160) = ' McDUA C nitrile  '
      chemid(160)  = 'C  '

! --- McDonald UA acetonitrile methyl [CH3]-C=-N
      sigi(161) = 3.6d0
      epsi(161) = 191.0d0
      mass(161) = 15.035d0
      qelect(161) = 0.269d0
      lqchg(161) = .true.
      chname(161) = ' McDUA CH3 nitrile'
      chemid(161)  = 'C  '

! * CHARMM

! --- Charmm C2 (methylene group carbon)
      sigi(162) = 3.8754d0
      epsi(162) = 19.6257d0
      mass(162) = 0.003d0
      chname(162) = ' CHARMM C2 ???    '     
      chemid(162)  =  'C  '

! --- Charmm H (hydrogen)
      sigi(163) = 2.4500d0
      epsi(163) = 19.1225d0
      mass(163) = 1.0078d0
      chname(163) = ' CHARMM H  ???    '
      chemid(163)  = 'H  '

! --- Charmm O (bound with 2 single bonds)
      sigi(164) = 2.8598d0
      epsi(164) = 114.7348d0
      mass(164) = 16.00d0
      chemid(164)  = 'O  '
      chname(164) = ' CHARMM O sp3 ??? '

! --- Charmm P
      sigi(165) = 3.7418d0
      epsi(165) = 100.6446d0
      mass(165) = 0.003d0
      chname(165) = ' CHARMM P ???     '
      chemid(165)  = 'P  '

! --- Charmm O' (bound with a double bond)
      sigi(166) = 2.8598d0
      epsi(166) = 114.7348d0
      mass(166) = 16.00d0
      chname(166) = ' CHARMM P ???     '
      chemid(166)  = 'O  '

! --- Charmm N3 (tertiary ammonia)
      sigi(167) = 3.5012d0
      epsi(167) = 84.0382d0
      mass(167) = 0.003d0
      chname(167) = ' CHARMM P ???     '
      chemid(167)  = 'N  '

! --- Charmm C3 (methyl group carbon)
      sigi(168) = 3.8754d0
      epsi(168) = 19.6257d0
      mass(168) = 0.003d0
      chname(168) = ' CHARMM C3 ???    '
      chemid(168)  = 'C  '

! --- Charmm C1 (ternary carbon)
      sigi(169) = 3.8754d0
      epsi(169) = 19.6257d0
      mass(169) = 0.003d0
      chname(169) = ' CHARMM C1 ???    '
      chemid(169)  = 'C  '

! --- Charmm C' (carboxylic head group carbon)
      sigi(170) = 3.6170d0
      epsi(170) = 74.4770d0
      mass(170) = 0.003d0
      chname(170) = ' CHARMM C ???    '
      chemid(170)  = 'C  '


! * ALL-ATOM NITRILES

! --- TraPPE-AA nitrile nitrogen C=-[N]
      sigi(171) = 2.95d0
      epsi(171) = 60.0d0
      mass(171) = 14.007d0
      qelect(171) = -0.398d0
      lqchg(171) = .true.
      chname(171) = ' Tr-AA N nitrile  '
      chemid(171)  = 'N  '

! --- TraPPE-AA nitrile carbon [C]=-N
      sigi(172) = 3.55d0
      epsi(172) = 60.0d0
      mass(172) = 12.011d0
      qelect(172) = 0.129d0
      lqchg(172) = .true.
      chname(172) = ' Tr-AA C nitrile  '
      chemid(172)  = 'C  '

! --- TraPPE-AA acetonitrile methyl carbon H3[C]-C=-N
      sigi(173) = 3.3d0
      epsi(173) = 4.0d0
      mass(173) = 12.011d0
      qelect(173) = 0.269d0
      lqchg(173) = .true.
      chname(173) = ' Tr-AA H3[C]-C=-N '
      chemid(173)  = 'C  '

! --- TraPPE-AA acetonitrile methyl hydrogen C[H3]-C=-N
      sigi(174) = sigi(27)
      epsi(174) = epsi(27)
      mass(174) = mass(27)
      chname(174) = ' Tr-AA [H]3C-C=-N '
      chemid(174)  = 'H  '


!     * SILICA

!----[Si]-O-Si
       sigi(177) = 0.0d0
       epsi(177) = 0.0d0
       mass(177) = 28.0d0      
       qelect(177) =1.216d0
       lqchg(177) = .true.
       lij(177) = .false.
       chname(177) = ' [Si]-O-Si '
      chemid(177)  = 'Si '

!----Si-[O]-Si
       sigi(178) = 3.35d0
       epsi(178) = 70.0d0
       mass(178) = 16.0d0      
       qelect(178) =-0.608d0
       lqchg(178) = .true.
       lij(178) = .true.
       chname(178) = ' Si-[O]-Si '
      chemid(178)  = 'O  '

!----O-[Si]-CH2
       sigi(179) = 6.4d0
       epsi(179) = 0.5d0
       mass(179) = 28.0d0      
       qelect(179) = 0.304d0
       lqchg(179) = .true.
       lij(179) = .true.
       chname(179) = ' O-[Si]-CH2 '
      chemid(179)  = 'Si '

!---- [CH3]-Si-O 
       sigi(180) = sigi(4)
       epsi(180) = epsi(4)
       mass(180) = mass(4)     
       qelect(180) =0.0d0
       lqchg(180) = .false.
       lij(180) = .true.
       chname(180) = '[CH3]-Si-O'
      chemid(180)  = 'C '
 
!----[CH2]-Si-O
       sigi(181) = sigi(5)
       epsi(181) = epsi(5)
       mass(181) = mass(5)     
       qelect(181) =0.0d0
       lqchg(181) = .false.
       lij(181) = .true.
       chname(181) = ' [CH2]-Si-O '
      chemid(181)  = 'C  '

!-- [Si] in SiO2 substrate
       sigi(182) = 0.0d0
       epsi(182) = 0.0d0
       mass(182) = 28.0d0
       qelect(182) =1.216d0
       lqchg(182) = .true.
       lij(182) = .false.
       chname(182) = ' [Si] in SiO2 substrate '
      chemid(182)  = 'Si '

!---  silanol oxygen H-[O]-Si
       sigi(183) = 3.35d0
       epsi(183) = 70.0d0
       mass(183) = 16.00d0
       qelect(183) = -0.739d0
       lqchg(183) = .true.
       lij(183) = .true.      
       chname(183) = ' H-[O]-Si '
      chemid(183)  = 'O  '
       
!---  fullerene [C]
       sigi(184) = 3.469d0
       epsi(184) = 33.247d0
       mass(184) =  12.011d0
      chemid(184)  = 'C  '
       

!   12 site benzene model with hydrogen at the normal position

      sigi(185) = 3.60d0
      epsi(185) = 30.7d0
      mass(185) = 12.011d0
      qelect(185) = -0.095d0
      lqchg(185) = .true.
      lij(185) = .true.
      chname(185) = 'C Trappe AA benzene  '
      chemid(185)  = 'C  '

!   benzene 12 site model

      sigi(186) = 2.36d0
      epsi(186) = 25.44d0
      mass(186) = 1.0079d0
      qelect(186) = 0.095d0
      lqchg(186) = .true.
      lij(186) = .true.
      chname(186) = 'H Trappe AA benzene  '
      chemid(186)  = 'H  '

! -- MFI silicalite-1 oxygen
      sigi(190) = 3.0d0
      epsi(190) = 93.53d0
      mass(190) = 15.999d0
      chname(190) = ' silicalite-1 O '
      chemid(190) = 'O  '

! ---- added 7/12/06 for nitrotoluene

! -- TraPPE-UA [C] alpha aro carbon for nitro
      sigi(196) = 4.50d0
      epsi(196) = 15.0d0
      mass(196) = 12.011d0
      lqchg(196) = .true.
      qelect(196) = 0.14d0
      chname(196) = ' Tr-UA C aro nitro '
      chemid(196) = 'c  '

! -- TraPPE-UA [N] nitro
      sigi(197) = 3.31d0
      epsi(197) = 40.0d0
      mass(197) = 14.007d0
      lqchg(197) = .true.
      qelect(197) = 0.82d0
      chname(197) = 'Tr-UA N nitro '
      chemid(197) = 'N  '

! -- TraPPE-UA [O] nitro
      sigi(198) = 2.90d0
      epsi(198) = 80.0d0
      mass(198) = 15.999d0
      lqchg(198) = .true.
      qelect(198) = -0.48d0
      chname(198) = 'Tr-UA O nitro '
      chemid(198) = 'O  '

! --- TraPPE-UA [CH] benzene9 carbon also #56
      sigi(199) = 3.74d0
      epsi(199) = 48.0d0
      mass(199) = 13.0191d0
! * published CH(aro) for TraPPE-UA 6-site
!      sigi(56) = 3.695d0
!      epsi(56) = 50.5d0
!      mass(56) = 13.0191d0
      chname(199) = ' Tr-UA CH benzene9'
      chemid(199)  = 'C  '


! --- JLR 12-1-09 parameters for gradually growing in benzene
       sigi(203) = 1.6d0
       epsi(203) = 20.0d0
       mass(203) = 13.091d0
       lqchg(203) = .false.
       lij(203) = .true.
       chname(203) = 'stage 1 benzene'

       sigi(204) = 2.1d0
       epsi(204) = 30.0d0
       mass(204) = 13.091d0
       lqchg(204) = .false.
       lij(204) = .true.
       chname(204) = 'stage 2 benzene'

       sigi(205) = 2.5d0
       epsi(205) = 37.0d0
       mass(205) = 13.091d0
       lqchg(205) = .false.
       lij(205) = .true.
       chname(205) = 'stage 3 benzene'

       sigi(206) = 2.9d0
       epsi(206) = 42.0d0
       mass(206) = 13.091d0
       lqchg(206) = .false.
       lij(206) = .true.
       chname(206) = 'stage 4 benzene'

       sigi(207) = 3.3d0
       epsi(207) = 46.0d0
       mass(207) = 13.091d0
       lqchg(207) = .false.
       lij(207) = .true.
       chname(207) = 'stage 5 benzene'
! --- END JLR 12-1-09 ---

! - parameters for acrylates
! -- some are already listed; listed twice for convenience during fitting

! --- methyl group attached to ether oxygen (TraPPE 6) #72
      sigi(210) = 3.75d0
      epsi(210) = 98.0d0
      mass(210) = 15.0347d0
      qelect(210) = 0.25d0
      lij(210) = .true.
      lqchg(210) = .true.
      chname(210) = 'Tr-UA ether CH3'
      chemid(210) = 'C  '

! --- ether oxygen #71
      sigi(211) = 2.80d0
      epsi(211) = 55.0d0
      mass(211) = 15.999d0
      qelect(211) = -0.25d0
      lij(211) = .true.
      lqchg(211) = .true.
      chname(211) = 'Tr-UA ether O'
      chemid(211) = 'O  '

! --  ketone with CM4 charge
      sigi(212) = 3.82d0
      epsi(212) = 40.0d0
      mass(212) = 12.011d0
      qelect(212) = 0.4d0
      lij(212) = .true.
      lqchg(212) = .true.
      chname(212) = 'carbonyl C'
      chemid(212) = 'C  '

! --- C=O oxygen CM4 charge
      sigi(213) = 3.05d0
      epsi(213) = 79.0d0
      mass(213) = 15.999d0
      qelect(213) = -0.4d0
      lij(213) = .true.
      lqchg(213) = .true.
      chname(213) = 'C=O oxygen'
      chemid(213) = 'O  '

! -- TraPPE-UA sp2 butadiene
      sigi(214) = 3.71d0
      epsi(214) = 52.0d0
      mass(214) = 13.0191d0
      qelect(214) = 0.0d0
      lij(214) = .true.
      lqchg(214) = .true.
      chname(214) = 'Tr-UA sp2 CH w/charge'
      chemid(214) = 'C  '
      
! --  TraPPE-UA sp2 CH2 #51
      sigi(215) = 3.675d0
      epsi(215) = 85.0d0
      mass(215) = 14.0268d0
      qelect(215) = 0.0d0
      lij(215) = .true.
      lqchg(215) = .true.
      chname(215) = 'Tr-UA sp2 CH2'
      chemid(215) = 'C  '

! --  TraPPE-UA methyl CH3 #4
      sigi(216) = 3.75d0
      epsi(216) = 98.0d0
      mass(216) = 15.0347d0
      lij(216) = .true.
      chname(216) = 'Tr-UA CH3'
      chemid(216) = 'C  '

! --  TraPPE-UA C(sp2)
      sigi(217) = 3.85d0
      epsi(217) = 22.0d0
      mass(217) = 12.011d0
      lij(217) = .true.
      chname(217) = 'Tr-UA sp2 C'
      chemid(217) = 'C  '

! --  TraPPE-UA CH2-(ether O) with different charge
      sigi(218) = 3.95d0
      epsi(218) = 46.0d0
      mass(218) = 14.0268d0
      qelect(218) = 0.25d0
      lij(218) = .true.
      lqchg(218) = .true.
      chname(218) = 'Tr-UA ether CH2'
      chemid(218) = 'C  '

! -- parameters for primary amines

! -- CH3-N-[H]-H
      sigi(220) = 0.0d0
      epsi(220) = 0.0d0
      mass(220) = 1.0079d0
!      qelect(220) = 0.385d0

!     --- first degree
      qelect(220) = 0.356d0
      lqchg(220) = .true.
      lij(220) = .false.
      chname(220) = 'TraPPE-AA 2o H Amine'
      chemid(220)  = 'H   '

! -- CH3-[N]-H2

!     *** 2nd degree *******
!      sigi(221) = 3.52d0
!      epsi(221) = 58.0d0
!      qelect(221) = -0.745d0
!     ***********************


!     *** 3rd degree ****
!      sigi(221) = 3.78d0
!      epsi(221) = 12.0d0
!      qelect(221) = -0.54d0

!     ********************************
                                         
!     ********************************

!     *** 1st degree *****
      sigi(221) = 3.34d0
      epsi(221) = 111.0d0
      qelect(221) = -0.892d0
      mass(221) = 14.00674d0
      lqchg(221) = .true.
      lij(221) = .true.
      chname(221) = 'TraPPE-AA H Amine'
      chemid(221)  = 'N   '
        

! -- [C(methylene)]-N-H2
      sigi(222) = sigi(24)
      epsi(222) = epsi(24)
      mass(222) = mass(24)
      qelect(222) = 0.18d0
      lqchg(222) = .true.
      chname(222) = 'TraPPE-AA C Amine'
      chemid(222)  = 'C   '


! -- N-[C]-[H2]
      sigi(223) = sigi(24)
      epsi(223) = epsi(24)
      mass(223) =  mass(24)
      qelect(223) = 0.18d0
      lqchg(223) = .true.
      chname(223) = 'TraPPE-AA C Amine'
      chemid(223)  = 'C   '

!-- CH3-[N]-CH3-H, dimethylamine
       sigi(224) = 3.31d0
       epsi(224) = 115.0d0
       mass(224) = 14.00674d0
       qelect(224) = -0.75d0
       lqchg(224) = .true.
       chname(224) = 'TraPPE-AA N Amine'
       chemid(224)  = 'N   '

 

! -- starting for the carboxylic acid (Jeff's 2004 paper)

! --- TraPPE-UA carboxylic acid carbonyl oxygen C=[O]-O-H
      sigi(230) = 3.05d0
      epsi(230) = 79.0d0
      mass(230) = 15.999d0
      qelect(230) = -0.45d0
      lqchg(230) = .true.
      chname(230) = ' Tr-UA C[O]OH acid'
      chemid(230)  = 'O  '

! --- TraPPE-UA iterB carbonyl carbon  [C](=O)-O-H
      sigi(231) = 3.90d0
      epsi(231) = 41.00d0
      mass(231) = 12.011d0
      qelect(231) = 0.42d0
      lqchg(231) = .true.
      chname(231) = ' Tr-UA [C]OOH acid'
      chemid(231)  = 'C  '

! --- TraPPE-UA carboxylic acid hydrogen C(=O)-O-[H] (JPC v95 p. 3315, 1991)
      sigi(232) = 0.0d0
      epsi(232) = 0.0d0
      mass(232) = 1.0079d0
      qelect(232) = 0.37d0
      lij(232) = .false.
      lqchg(232) = .true.
      chname(232) = ' Tr-UA COO[H] acid'
      chemid(232)  = 'H  '

! --- TraPPE-UA carboxylic acid oxygen C(=O)-[O]-H
      sigi(233) = 3.02d0
      epsi(233) = 93.0d0
      mass(233) = 16.00d0
      qelect(233) = -0.46d0
      lqchg(233) = .true.
      chname(233) = ' Tr-UA CO[O]H acid'
      chemid(233)  = 'O  '

! --- TraPPE-UA carboxylic acid oxygen [CH3]-C(=O)-O-H
      sigi(234) = sigi(4)
      epsi(234) = epsi(4)
      mass(234) = mass(4)
      qelect(234) = 0.12d0
      lqchg(234) = .true.
      chname(234) = ' Tr-UA CO[O]H acid'
      chemid(234)  = 'C  '

! --- TraPPE-UA carboxylic acid oxygen [CH2-C(=O)-O-H
      sigi(235) = sigi(5)
      epsi(235) = epsi(5)
      mass(235) = mass(5)
      qelect(235) = 0.12d0
      lqchg(235) = .true.
      chname(235) = ' Tr-UA CO[O]H acid'
      chemid(235)  = 'C  '


! -- starting for fluoropropane

! --  [C]F3 Terminal Methyl

      sigi(250) = 3.55d0
      epsi(250) = 35.0d0
      mass(250) = 12.011d0
      qelect(250) = 0.245d0
      lqchg(250) = .true.
      lij(250) = .true.
      chname(250) = 'C methyl'
      chemid(250) ='C  '
 
! -- [C]F2 Methylene
      sigi(251) = 3.55d0
      epsi(251) = 35.0d0
      mass(251) = 12.011d0
      qelect(251) = 0.152d0
      lqchg(251) = .true.
      lij(251) = .true.
      chname(251) = 'C methylene'
      chemid(251) ='C  '
! -- C[F]3
      sigi(252) = 2.95d0
      epsi(252) = 25.0d0
      mass(252) = 18.9984d0
      qelect(252) = -0.078d0
      lqchg(252) = .true.
      lij(252) = .true.
      chname(252) = 'F in CF3'
      chemid(252) ='F  '

! -- C[F]2
      sigi(253) = 2.95d0
      epsi(253) = 25.0d0
      mass(253) = 18.9984d0
      qelect(253) = -0.087d0
      lqchg(253) = .true.
      lij(253) = .true.
      chname(253) = 'F CF2'
      chemid(253) ='F  '

! -- [H]-CF (bonded methylene type carbon) for HFC227
      sigi(254) = 2.36d0
      epsi(254) = 20.40d0
      mass(254) = 1.0079d0
      qelect(254) = 0.084d0
      lqchg(254) = .true.
      lij(254) = .true.
      chname(254) = 'H CHF'
      chemid(254) ='H  '

! --  [C]F3 Terminal Methyl
      sigi(255) = 3.65d0
      epsi(255) = 27.50d0
      mass(255) = 12.011d0
      qelect(255) = 0.256d0
      lqchg(255) = .true.
      lij(255) = .true.
      chname(255) = 'C methyl'
      chemid(255) ='C  '

! -- [C]F2 Methylene
      sigi(256) = 3.70d0
      epsi(256) = 28.0d0
      mass(256) = 12.011d0
      qelect(256) = 0.068d0
      lqchg(256) = .true.
      lij(256) = .true.
      chname(256) = 'C methylene'
      chemid(256) ='C  '

! -- C[F]3
      sigi(257) = 2.92d0
      epsi(257) = 32.50d0
      mass(257) = 18.9984d0
      qelect(257) = -0.09d0
      lqchg(257) = .true.
      lij(257) = .true.
      chname(257) = 'F in CF3'
      chemid(257) ='F  '

! -- C[F]H
      sigi(258) = 2.92d0
      epsi(258) = 32.50d0
      mass(258) = 18.9984d0
      qelect(258) = -0.124d0
      lqchg(258) = .true.
      lij(258) = .true.
      chname(258) = 'F CF2'
      chemid(258) ='F  '


! -- Starting all atom alkane. Starting Ethane and then Ethanol

! -- [C]H3 Methyl carbon
      sigi(275) = 3.55d0
      epsi(275) = 35.0d0
      mass(275) = 12.011d0
      qelect(275) = -0.24d0
      lqchg(275) = .true.
      lij(275) = .true.
      chname(275) = 'C methyl AA'
      chemid(275) ='C  '

! -- [H]-CH2 Hydrogen AA
      sigi(276) = 2.55d0
      epsi(276) = 17.50d0
      mass(276) = 1.0079d0
      qelect(276) = 0.08d0
      lqchg(276) = .true.
      lij(276) = .true.
      chname(276) = 'H in CH3'
      chemid(276) ='H  '


! -- [C]H3 in ethanol
      sigi(285) = 3.55d0
      epsi(285) = 35.0d0
      mass(285) = 12.011d0
      qelect(285) = -0.231d0
      lqchg(285) = .true.
      lij(285) = .true.
      chname(285) = 'C in CH3'
      chemid(285) ='C  '

! -- [C]H2 in ethanol
      sigi(286) = 3.55d0
      epsi(286) = 35.50d0
      mass(286) = 12.011d0
      qelect(286) = 0.026d0
      lqchg(286) = .true.
      lij(286) = .true.
      chname(286) = 'C in CH2'
      chemid(286) ='C  '

!  -- C[H]H2 in ethanol
      sigi(287) = 2.55d0
      epsi(287) = 15.50d0
      mass(287) = 1.0079d0
      qelect(287) = 0.086d0
      lqchg(287) = .true.
      lij(287) = .true.
      chname(287) = 'H in CH3'
      chemid(287) ='H  '

!  -- C[H]H in ethanol
      sigi(288) = 2.55d0
      epsi(288) = 15.50d0
      mass(288) = 1.0079d0
      qelect(288) = 0.055d0
      lqchg(288) = .true.
      lij(288) = .true.
      chname(288) = 'H in CH2'
      chemid(288) ='H  ' 

!  -- [O] in ethanol
      sigi(289) = 2.9d0
      epsi(289) = 80.50d0
      mass(289) = 15.9998d0
      qelect(289) = -0.478d0
      lqchg(289) = .true.
      lij(289) = .true.
      chname(289) = 'O in C2H5OH'
      chemid(289) ='O  '

! -- [H]-O in ethanol
      sigi(290) = 0.5d0
      epsi(290) = 12.00d0
      mass(290) = 1.0079d0
      qelect(290) = 0.315d0
      lqchg(290) = .true.
      lij(290) = .true.
      chname(290) = 'H in OH'
      chemid(290) ='H  '

!   EH m-nitrotoluene 4/2/09 KM
! -- #299-312 charge model 1
! --  simply combine H and C charges for CH3
! -- #313-319 charge model 2
! --  adjust CH3 and ring carbon charges
! --  to recover the CM4 dipole moment

! -- TraPPE-EH [O1] nitro m
! -- from TraPPE 10 
      sigi(299) = 2.70d0
      epsi(299) = 42.0d0
      mass(299) = 15.999d0
      qelect(299) = -0.202d0
      lqchg(299) = .true.
      lij(299) = .true.
      chname(299) = 'O in NO2'
      chemid(299 ) ='O  '

! -- TraPPE-EH [O2] nitro m
! -- from TraPPE 10 
      sigi(300) = 2.70d0
      epsi(300) = 42.0d0
      mass(300) = 15.999d0
      qelect(300) = -0.196d0
      lqchg(300) = .true.
      lij(300) = .true.
      chname(300) = 'O in NO2'
      chemid(300) ='O  '

! -- TraPPE-EH [N] nitro m
! -- from TraPPE 10 
      sigi(301) = 2.90d0
      epsi(301) = 30.0d0
      mass(301) = 14.007d0
      qelect(301) = 0.038d0
      lqchg(301) = .true.
      lij(301) = .true.
      chname(301) = 'N in NO2'
      chemid(301) ='N  '

! -- TraPPE-EH [C] aro alpha nitro m
! -- from TraPPE 9
! -- nothing special about alpha to nitro?
      sigi(302) = 3.60d0
      epsi(302) = 30.70d0
      mass(302) = 12.011d0
      qelect(302) = 0.096d0
      lqchg(302) = .true.
      lij(302) = .true.
      chname(302) = 'C aro'
      chemid(302) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(303) = sigi(302)
      epsi(303) = epsi(302)
      mass(303) = mass(302)
      qelect(303) = -0.053d0
      lqchg(303) = .true.
      lij(303) = .true.
      chname(303) = 'C aro'
      chemid(303) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(304) = sigi(302)
      epsi(304) = epsi(302)
      mass(304) = mass(302)
      qelect(304) = -0.059d0
      lqchg(304) = .true.
      lij(304) = .true.
      chname(304) = 'C aro'
      chemid(304) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(305) = sigi(302)
      epsi(305) = epsi(302)
      mass(305) = mass(302)
      qelect(305) = -0.048d0
      lqchg(305) = .true.
      lij(305) = .true.
      chname(305) = 'C aro'
      chemid(305) ='C  '

! -- TraPPE-EH [C] aro alpha CH3 m
! -- from TraPPE 9
      sigi(306) = sigi(302)
      epsi(306) = epsi(302)
      mass(306) = mass(302)
      qelect(306) = -0.015d0
      lqchg(306) = .true.
      lij(306) = .true.
      chname(306) = 'C aro'
      chemid(306) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(307) = sigi(302)
      epsi(307) = epsi(302)
      mass(307) = mass(302)
      qelect(307) = -0.057d0
      lqchg(307) = .true.
      lij(307) = .true.
      chname(307) = 'C aro'
      chemid(307) ='C  '

! -- TraPPE-EH [H] aro m
! -- from TraPPE 9
      sigi(308) = 2.36d0
      epsi(308) = 25.45d0
      mass(308) = 1.008d0
      qelect(308) = 0.086d0
      lqchg(308) = .true.
      lij(308) = .true.
      chname(308) = 'H aro'
      chemid(308) ='H  '

! -- TraPPE-EH [H] aro m
! -- from TraPPE 9
      sigi(309) = 2.36d0
      epsi(309) = 25.45d0
      mass(309) = 1.008d0
      qelect(309) = 0.10d0
      lqchg(309) = .true.
      lij(309) = .true.
      chname(309) = 'H aro'
      chemid(309) ='H  '

! -- TraPPE-EH [H] aro m
! -- from TraPPE 9
      sigi(310) = 2.36d0
      epsi(310) = 25.45d0
      mass(310) = 1.008d0
      qelect(310) = 0.106d0
      lqchg(310) = .true.
      lij(310) = .true.
      chname(310) = 'H aro'
      chemid(310) ='H  '

! -- TraPPE-EH [H] aro m
! -- from TraPPE 9
      sigi(311) = 2.36d0
      epsi(311) = 25.45d0
      mass(311) = 1.008d0
      qelect(311) = 0.10d0
      lqchg(311) = .true.
      lij(311) = .true.
      chname(311) = 'H aro'
      chemid(311) ='H  '

! -- TraPPE-UA CH3 for toluene m
      sigi(312) = sigi(4)
      epsi(312) = epsi(4)
      mass(312) = mass(4)
      qelect(312) = 0.104d0
      lqchg(312) = .true.
      lij(312) = .true.
      chname(312) = 'CH3 '
      chemid(312) ='C  '


! -- charge model 2

! -- TraPPE-EH [C] aro alpha nitro m
! -- from TraPPE 9
! -- nothing special about alpha to nitro?
      sigi(313) = 3.60d0
      epsi(313) = 30.70d0
      mass(313) = 12.011d0
      qelect(313) = 0.089d0
      lqchg(313) = .true.
      lij(313) = .true.
      chname(313) = 'C aro'
      chemid(313) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(314) = sigi(302)
      epsi(314) = epsi(302)
      mass(314) = mass(302)
      qelect(314) = -0.060d0
      lqchg(314) = .true.
      lij(314) = .true.
      chname(314) = 'C aro'
      chemid(314) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(315) = sigi(302)
      epsi(315) = epsi(302)
      mass(315) = mass(302)
      qelect(315) = -0.066d0
      lqchg(315) = .true.
      lij(315) = .true.
      chname(315) = 'C aro'
      chemid(315) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(316) = sigi(302)
      epsi(316) = epsi(302)
      mass(316) = mass(302)
      qelect(316) = -0.054d0
      lqchg(316) = .true.
      lij(316) = .true.
      chname(316) = 'C aro'
      chemid(316) ='C  '

! -- TraPPE-EH [C] aro alpha CH3 m
! -- from TraPPE 9
      sigi(317) = sigi(302)
      epsi(317) = epsi(302)
      mass(317) = mass(302)
      qelect(317) = -0.022d0
      lqchg(317) = .true.
      lij(317) = .true.
      chname(317) = 'C aro'
      chemid(317) ='C  '

! -- TraPPE-EH [C] aro m
! -- from TraPPE 9
      sigi(318) = sigi(302)
      epsi(318) = epsi(302)
      mass(318) = mass(302)
      qelect(318) = -0.064d0
      lqchg(318) = .true.
      lij(318) = .true.
      chname(318) = 'C aro'
      chemid(318) ='C  '

! -- TraPPE-UA CH3 for toluene m
      sigi(340) = sigi(4)
      epsi(340) = epsi(4)
      mass(340) = mass(4)
      qelect(340) = 0.145d0
      lqchg(340) = .true.
      lij(340) = .true.
      chname(340) = 'CH3 '
      chemid(340) ='C  '

!   EH o-nitrotoluene 4/2/09 KM
! -- #319-332 charge model 1
! --  simply combine H and C charges for CH3
! -- #333-339 charge model 2
! --  adjust CH3 and ring carbon charges
! --  to recover the CM4 dipole moment

! -- TraPPE-EH [O1] nitro o
! -- from TraPPE 10 
      sigi(319) = 2.70d0
      epsi(319) = 42.0d0
      mass(319) = 15.999d0
      qelect(319) = -0.173d0
      lqchg(319) = .true.
      lij(319) = .true.
      chname(319) = 'O in NO2'
      chemid(319) ='O  '

! -- TraPPE-EH [O2] nitro o
! -- from TraPPE 10 
      sigi(320) = 2.70d0
      epsi(320) = 42.0d0
      mass(320) = 15.999d0
      qelect(320) = -0.192d0
      lqchg(320) = .true.
      lij(320) = .true.
      chname(320) = 'O in NO2'
      chemid(320) ='O  '

! -- TraPPE-EH [N] nitro o
! -- from TraPPE 10 
      sigi(321) = 2.90d0
      epsi(321) = 30.0d0
      mass(321) = 14.007d0
      qelect(321) = 0.029d0
      lqchg(321) = .true.
      lij(321) = .true.
      chname(321) = 'N in NO2'
      chemid(321) ='N  '

! -- TraPPE-EH [C] aro alpha nitro o
! -- from TraPPE 9
! -- nothing special about alpha to nitro?
      sigi(322) = 3.60d0
      epsi(322) = 30.70d0
      mass(322) = 12.011d0
      qelect(322) = 0.059d0
      lqchg(322) = .true.
      lij(322) = .true.
      chname(322) = 'C aro'
      chemid(322) ='C  '

! -- TraPPE-EH [C] aro o
! -- from TraPPE 9
      sigi(323) = sigi(302)
      epsi(323) = epsi(302)
      mass(323) = mass(302)
      qelect(323) = -0.051d0
      lqchg(323) = .true.
      lij(323) = .true.
      chname(323) = 'C aro'
      chemid(323) ='C  '

! -- TraPPE-EH [C] aro o
! -- from TraPPE 9
      sigi(324) = sigi(302)
      epsi(324) = epsi(302)
      mass(324) = mass(302)
      qelect(324) = -0.071d0
      lqchg(324) = .true.
      lij(324) = .true.
      chname(324) = 'C aro'
      chemid(324) ='C  '

! -- TraPPE-EH [C] aro o
! -- from TraPPE 9
      sigi(325) = sigi(302)
      epsi(325) = epsi(302)
      mass(325) = mass(302)
      qelect(325) = -0.026d0
      lqchg(325) = .true.
      lij(325) = .true.
      chname(325) = 'C aro'
      chemid(325) ='C  '

! -- TraPPE-EH [C] aro CH3 o
! -- from TraPPE 9
      sigi(326) = sigi(302)
      epsi(326) = epsi(302)
      mass(326) = mass(302)
      qelect(326) = -0.074d0
      lqchg(326) = .true.
      lij(326) = .true.
      chname(326) = 'C aro'
      chemid(326) ='C  '

! -- TraPPE-EH [C] aro alpha CH3 o
! -- from TraPPE 9
      sigi(327) = sigi(302)
      epsi(327) = epsi(302)
      mass(327) = mass(302)
      qelect(327) = 0.003d0
      lqchg(327) = .true.
      lij(327) = .true.
      chname(327) = 'C aro'
      chemid(327) ='C  '

! -- TraPPE-EH [H] aro o
! -- from TraPPE 9
      sigi(328) = 2.36d0
      epsi(328) = 25.45d0
      mass(328) = 1.008d0
      qelect(328) = 0.095d0
      lqchg(328) = .true.
      lij(328) = .true.
      chname(328) = 'H aro'
      chemid(328) ='H  '

! -- TraPPE-EH [H] aro o
! -- from TraPPE 9
      sigi(329) = 2.36d0
      epsi(329) = 25.45d0
      mass(329) = 1.008d0
      qelect(329) = 0.109d0
      lqchg(329) = .true.
      lij(329) = .true.
      chname(329) = 'H aro'
      chemid(329) ='H  '

! -- TraPPE-EH [H] aro o
! -- from TraPPE 9
      sigi(330) = 2.36d0
      epsi(330) = 25.45d0
      mass(330) = 1.008d0
      qelect(330) = 0.105d0
      lqchg(330) = .true.
      lij(330) = .true.
      chname(330) = 'H aro'
      chemid(330) ='H  '

! -- TraPPE-EH [H] aro o
! -- from TraPPE 9
      sigi(331) = 2.36d0
      epsi(331) = 25.45d0
      mass(331) = 1.008d0
      qelect(331) = 0.099d0
      lqchg(331) = .true.
      lij(331) = .true.
      chname(331) = 'H aro'
      chemid(331) ='H  '

! -- TraPPE-UA CH3 for toluene o
      sigi(332) = sigi(4)
      epsi(332) = epsi(4)
      mass(332) = mass(4)
      qelect(332) = 0.088d0
      lqchg(332) = .true.
      lij(332) = .true.
      chname(332) = 'CH3 '
      chemid(332) ='C  '

! -- charge model 2 

! -- TraPPE-EH [C] aro alpha nitro o
! -- from TraPPE 9
! -- nothing special about alpha to nitro?
      sigi(333) = 3.60d0
      epsi(333) = 30.70d0
      mass(333) = 12.011d0
      qelect(333) = 0.0535d0
      lqchg(333) = .true.
      lij(333) = .true.
      chname(333) = 'C aro'
      chemid(333) ='C  '

! -- TraPPE-EH [C] aro o
! -- from TraPPE 9
      sigi(334) = sigi(302)
      epsi(334) = epsi(302)
      mass(334) = mass(302)
      qelect(334) = -0.0565d0
      lqchg(334) = .true.
      lij(334) = .true.
      chname(334) = 'C aro'
      chemid(334) ='C  '

! -- TraPPE-EH [C] aro o
! -- from TraPPE 9
      sigi(335) = sigi(302)
      epsi(335) = epsi(302)
      mass(335) = mass(302)
      qelect(335) = -0.0765d0
      lqchg(335) = .true.
      lij(335) = .true.
      chname(335) = 'C aro'
      chemid(335) ='C  '

! -- TraPPE-EH [C] aro o
! -- from TraPPE 9
      sigi(336) = sigi(302)
      epsi(336) = epsi(302)
      mass(336) = mass(302)
      qelect(336) = -0.0315d0
      lqchg(336) = .true.
      lij(336) = .true.
      chname(336) = 'C aro'
      chemid(336) ='C  '

! -- TraPPE-EH [C] aro CH3 o
! -- from TraPPE 9
      sigi(337) = sigi(302)
      epsi(337) = epsi(302)
      mass(337) = mass(302)
      qelect(337) = -0.0795d0
      lqchg(337) = .true.
      lij(337) = .true.
      chname(337) = 'C aro'
      chemid(337) ='C  '

! -- TraPPE-EH [C] aro alpha CH3 o
! -- from TraPPE 9
      sigi(338) = sigi(302)
      epsi(338) = epsi(302)
      mass(338) = mass(302)
      qelect(338) = -0.0025d0
      lqchg(338) = .true.
      lij(338) = .true.
      chname(338) = 'C aro'
      chemid(338) ='C  '

! -- TraPPE-UA CH3 for toluene o
      sigi(339) = sigi(4)
      epsi(339) = epsi(4)
      mass(339) = mass(4)
      qelect(339) = 0.121d0
      lqchg(339) = .true.
      lij(339) = .true.
      chname(339) = 'CH3 '
      chemid(339) ='C  '












! -- Chlorobenzene parameters all atom

      sigi(341) = 3.60
      epsi(341) = 30.7
      mass(341) = 12.011d0
      qelect(341) = -0.11d0
      lqchg(341) = .true.
      lij(341) = .true.
      chname(341) = 'C1 Chlorobenzene AA'
      chemid(341)  = 'C  '

      sigi(342) = 3.600d0
      epsi(342) = 30.7d0
      mass(342) =  12.011d0
      qelect(342) = -0.07d0
      lqchg(342) = .true.
      lij(342) = .true.
      chname(342) = 'C2 Chlorobenzene AA'
      chemid(342)  = 'C  '

      sigi(343) = 3.600d0
      epsi(343) = 30.70d0
      mass(343) = 12.011d0
      qelect(343) = -0.11d0
      lqchg(343) = .true.
      lij(343) = .true.
      chname(343) = 'C3 Chlorobenzene AA'
      chemid(343)  = 'C  '

      sigi(344) = 3.600d0
      epsi(344) = 30.70d0
      mass(344) = 12.011d0
      qelect(344) = -0.055d0
      lqchg(344) = .true.
      lij(344) = .true.
      chname(344) = 'C4 Chlorobenzene AA'
      chemid(344)  = 'C  '

      sigi(345) = 3.600d0
      epsi(345) = 30.7d0
      mass(345) = 12.011d0
      qelect(345) = 0.06d0
      lqchg(345) = .true.
      lij(345) = .true.
      chname(345) = 'C5 Chlorobenzene AA'
      chemid(345)  = 'C  '

      sigi(346) = 3.60d0
      epsi(346) = 30.70d0
      mass(346) = 12.011d0
      qelect(346) = -0.055d0
      lqchg(346) = .true.
      lij(346) = .true.
      chname(346) = 'C6 Chlorobenzene AA'
      chemid(346)  = 'C  '

      sigi(347) = 2.36d0
      epsi(347) = 25.44d0
      mass(347) = 1.0079d0
      qelect(347) = 0.1d0
      lqchg(347) = .true.
      lij(347) = .true.
      chname(347) = 'H7 Chlorobenzene AA'
      chemid(347)  = 'H  '

      sigi(348) = 2.36d0
      epsi(348) = 25.44d0
      mass(348) = 1.0079d0
      qelect(348) = 0.09d0
      lqchg(348) = .true.
      lij(348) = .true.
      chname(349) = 'H8 Chlorobenzene AA'
      chemid(349)  = 'H  '

      sigi(350) = 2.360d0
      epsi(350) = 25.44d0
      mass(350) = 1.0079d0
      qelect(350) = 0.1d0
      lqchg(350) = .true.
      lij(350) = .true.
      chname(351) = 'H9 Chlorobenzene AA'
      chemid(351)  = 'H  '

      sigi(352) = 2.36d0
      epsi(352) = 25.44d0
      mass(352) = 1.0079d0
      qelect(353) = 0.09d0
      lqchg(353) = .true.
      lij(353) = .true.
      chname(353) = 'H10 Chlorobenzene AA'
      chemid(353)  = 'H  '

      sigi(354) = 2.36d0
      epsi(354) = 25.44d0
      mass(354) = 1.0079d0
      qelect(354) = 0.09d0
      lqchg(354) = .true.
      lij(354) = .true.
      chname(354) = 'H11 Chlorobenzene AA'
      chemid(354)  = 'H  '

      sigi(355) = 3.5d0
      epsi(355) = 158.0d0
      mass(355) = 35.4527d0
      qelect(355) = -0.13d0
      lqchg(355) = .true.
      lij(355) = .true.
      chname(355) = 'Cl12 Chlorobenzene AA'
      chemid(355)  = 'Cl '

   
! -- O-dichlorobenznene


      sigi(356) = 3.60d0
      epsi(356) = 30.7d0
      mass(356) = 12.011d0
      qelect(356) = -0.11d0
      lqchg(356) = .true.
      lij(356) = .true.
      chname(356) = 'C1 OChlorobenzene AA'
      chemid(356)  = 'C  '

      sigi(357) = 3.600d0
      epsi(357) = 30.7d0
      mass(357) =  12.011d0
      qelect(357) = -0.08d0
      lqchg(357) = .true.
      lij(357) = .true.
      chname(357) = 'C2 OChlorobenzene AA'
      chemid(357)  = 'C  '

      sigi(358) = 3.600d0
      epsi(358) = 30.70d0
      mass(358) = 12.011d0
      qelect(358) = 0.08d0
      lqchg(358) = .true.
      lij(358) = .true.
      chname(358) = 'C3 OChlorobenzene AA'
      chemid(358)  = 'C  '

      sigi(359) = 3.600d0
      epsi(359) = 30.70d0
      mass(359) = 12.011d0
      qelect(359) = 0.08d0
      lqchg(359) = .true.
      lij(359) = .true.
      chname(359) = 'C4 OChlorobenzene AA'
      chemid(359)  = 'C  '

      sigi(360) = 3.600d0
      epsi(360) = 30.7d0
      mass(360) = 12.011d0
      qelect(360) = -0.08d0
      lqchg(360) = .true.
      lij(360) = .true.
      chname(360) = 'C5 OChlorobenzene AA'
      chemid(360)  = 'C  '

      sigi(361) = 3.60d0
      epsi(361) = 30.70d0
      mass(361) = 12.011d0
      qelect(361) = -0.11d0
      lqchg(361) = .true.
      lij(361) = .true.
      chname(361) = 'C6 OChlorobenzene AA'
      chemid(361)  = 'C  '

      sigi(362) = 2.36d0
      epsi(362) = 25.44d0
      mass(362) = 1.0079d0
      qelect(362) = 0.11d0
      lqchg(362) = .true.
      lij(362) = .true.
      chname(362) = 'H7 OChlorobenzene AA'
      chemid(362)  = 'H  '

      sigi(363) = 2.36d0
      epsi(363) = 25.44d0
      mass(363) = 1.0079d0
      qelect(363) = 0.10d0
      lqchg(363) = .true.
      lij(363) = .true.
      chname(363) = 'H8 OChlorobenzene AA'
      chemid(363)  = 'H  '

      sigi(364) = 2.360d0
      epsi(364) = 25.44d0
      mass(364) = 1.0079d0
      qelect(364) = 0.1d0
      lqchg(364) = .true.
      lij(364) = .true.
      chname(364) = 'H9 OChlorobenzene AA'
      chemid(364)  = 'H  '

      sigi(365) = 2.36d0
      epsi(365) = 25.44d0
      mass(365) = 1.0079d0
      qelect(365) = 0.11d0
      lqchg(365) = .true.
      lij(365) = .true.
      chname(365) = 'H10 OChlorobenzene AA'
      chemid(365)  = 'H  '

      sigi(366) = 3.5d0
      epsi(366) = 158.0d0
      mass(366) = 35.4527d0
      qelect(366) = -0.1d0
      lqchg(366) = .true.
      lij(366) = .true.
      chname(366) = 'Cl11 OChlorobenzene AA'
      chemid(366)  = 'Cl '

      sigi(367) = 3.5d0
      epsi(367) = 158.0d0
      mass(367) = 35.4527d0
      qelect(367) = -0.1d0
      lqchg(367) = .true.
      lij(367) = .true.
      chname(367) = 'Cl12 OChlorobenzene AA'
      chemid(367)  = 'Cl '

! --Adding for ethylene and propylene carbonate
      
! -- [CH2]-O-

      sigi(370) = sigi(5) 
      epsi(370) = epsi(5)
      mass(370) = mass(5)
      qelect(370) = 0.25d0
      lqchg(370) = .true.
      lij(370) = .true.
      chname(370) = 'Ch2 ether'
      chemid(370) ='C  '

! -- CH2-[O]-

      sigi(371) = 2.85d0
      epsi(371) = 55.0d0
      mass(371) = 15.9998d0
      qelect(371) = -0.50d0
      lqchg(371) = .true.
      lij(371) = .true.
      chname(371) = 'O ether'
      chemid(371) ='O  '

! -- CH2-O-[C]=O

      sigi(372) = 3.1d0
      epsi(372) = 35.0d0
      mass(372) = 12.011d0
      qelect(372) = 1.150d0
      lqchg(372) = .true.
      lij(372) = .true.
      chname(372) = 'C carbonate'
      chemid(372) ='C  '

! -- CH2-O-C=[O]

      sigi(373) = 3.04d0
      epsi(373) = 85.0d0
      mass(373) = 15.9998d0
      qelect(373) = -0.65d0
      lqchg(373) = .true.
      lij(373) = .true.
      chname(373) = 'O ketone'
      chemid(373) ='O  '

! TATB JCP 2004 120 7059
! -- [C]-NO2
      sigi(380) = 3.60d0
      epsi(380) = 30.7d0
      mass(380) = 12.011d0
      qelect(380) = -0.242d0
      lqchg(380) = .true.
      lij(380) = .true.
      chname(380) = 'C in TATB'
      chemid(380) ='C  '

! -- [C]-NH2
      sigi(381) = 3.60d0
      epsi(381) = 30.70d0
      mass(381) = 12.011d0
      qelect(381) = 0.408d0
      lqchg(381) = .true.
      lij(381) = .true.
      chname(381) = 'C in TATB'
      chemid(381) ='C  '

! -- [N]-O2 Nitro group
      sigi(382) = 2.90d0
      epsi(382) = 30.0d0
      mass(382) = 14.00747d0
      qelect(382) = 0.008d0
      lqchg(382) = .true.
      lij(382) = .true.
      chname(382) = 'N in TATB'
      chemid(382) ='N  '

! -- [N]-H2 Nitro group
      sigi(383) = 3.25d0
      epsi(383) = 160.0d0
      mass(383) = 14.00747d0
      qelect(383) = -0.738d0
      lqchg(383) = .true.
      lij(383) = .true.
      chname(383) = 'N in TATB'
      chemid(383) ='N  '

!c -- [O]- Nitro group
      sigi(384) = 2.70d0
      epsi(384) = 42.0d0
      mass(384) = 15.9998d0
      qelect(384) = -0.104d0
      lqchg(384) = .true.
      lij(384) = .true.
      chname(384) = 'O in TATB'
      chemid(384) ='O  '

! -- [H]- Nitro group
      sigi(385) = 0.50d0
      epsi(385) = 12.0d0
      mass(385) = 1.0079d0
      qelect(385) = 0.386d0
      lqchg(385) = .true.
      lij(385) = .true.
      chname(385) = 'H in TATB'
      chemid(385) ='H  '


!  LEFTOVER PIECES

! --- Monica's alcohol methyl Not bonded to O (-CH2-) - CH3
!     (van Leeuwen JPC 99, 1831 (1995))
!      sigi() = 3.93d0
!      epsi() = 110.0d0
!      mass() = 15.0347d0
!
!- Dummy methylene for Sciece Paper MGM 12-17-97
!      sigi() = 3.95d0
!      epsi() = 55.0d0
!      mass() = 14.0268d0
!
! --- AA for alkane methyl (CH3) carbon ---
! --- if 0.45 as C-H bond length
!      sigi() = 3.44d0
!      epsi() = 13.0d0
!      qelect() = -0.572d0
!      lqchg() = .true.
!
!c --- AA for alkane carbon ---
!      sigi() = 3.65d0
!      epsi() = 5.0d0
!      mass() = 12.011d0
!      qelect() = 0.265d0
!      lqchg() = .true.
!      lij() = .true.

! ===========================
! *** End Parameter Input ***
! ===========================

! --- Computation of un-like interactions
      if ( ljoe ) then
! --- STANDARD METHYL GROUP
         extc12(4) = 3.41d7
         extc3(4)  = 20800.0d0
         extz0(4)  = 0.86d0

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

      if (ltab) then
         open(5,file='lj.poten')
         read(5,*)
         read(5,*) nntype5
         read(5,*)
         do j=1,nntype5
            read(5,*) i,sigi(i),epsi(i),qelect(i)
!         read(5,'(I,4F,2A)') i,sigi(i),epsi(i),qelect(i),mass(i)
!     &        ,chemid(i),chname(i)
            if (qelect(i).ne.0) then
               lqchg(i)=.true.
            else
               lqchg(i)=.false.
            end if
            if (sigi(i).eq.0 .and. epsi(i).eq.0) then
               lij(i)=.false. 
            else
               lij(i)=.true.
            end if
         end do
      end if

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
! ---  Conversion factor for intermolecular coulomb interactions
      qqfact = 1.67d5

      if (ltab) then
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
            if (lshift) then
               sr2 = sig2ij(ij) / (rcut(1)*rcut(1))
               sr6 = sr2 * sr2 * sr2
               ecut(ij)= sr6*(sr6-1.0d0)*epsij(ij)
               ecut((j-1)*nntype+i)=ecut(ij)
            end if
         end do
         close(5)
      end if

      return
      end
