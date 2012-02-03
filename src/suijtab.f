      subroutine suijtab( lmixlb,lmixjo,qelect )

c suijtab
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
      include 'external.inc'
      include 'poten.inc'
      include 'system.inc'
      include 'cell.inc'
      include 'expsix.inc'
      include 'merck.inc'
      include 'nsix.inc'

      logical lmixlb, lmixjo

      integer i, j, ij, ji
      double precision qelect(nntype)
      double precision rzeronx(nxatom),epsilonnx(nxatom)

      double precision rcheck, sr2, sr6, adum, bdum, rs1, rs7, sr7,
     +     pi,djay

c ----------------------------------------------------------------
      do i = 1, nntype
         lpl(i) = .false.
      enddo

      pi = 4.0d0*datan(1.0d0)
      if(lexpsix) then
         do i = 1,natom
            aexsix(i) = 0.0d0
            bexsix(i) = 0.0d0
            cexsix(i) = 0.0d0
            qelect(i) = 0.0d0
            xiq(i) = 0.0d0
            lqchg(i) = .false.
            lij(i) = .true.     
         enddo
c --- Explicit atom carbon Williams Exp-6 potential
c --- J.Chem.Phys 47 11 4680 (1967) paramter set IV
c --- note that the combination of C--H has to be the average of C 
c      and H the way it is set up in the energy subroutines
c --- natom is set in expsix.inc
c     U(r) = A*r^(-6) + B*exp[C*r]

c --- C---C nonbonded interaction
c         aexsix(1) = -2.858d5
c         bexsix(1) = 4.208d7
c         cexsix(1) = -3.60d0
c         aexsix(1) = -2.541d5
c         bexsix(1) = 3.115d7
c         cexsix(1) = -3.60d0
c         mass(1) = 12.011d0
c         aexsix(1) = -4.0d0
c         bexsix(1) = 1.0d0/dexp(-12.0d0)
c         cexsix(1) = -12.0d0/(2.0d0 ** (1.0d0/6.0d0))
         aexsix(1) = -0.590759e6
         bexsix(1) = 0.281431e9
         cexsix(1) = -0.396875e1
         mass(1) = 39.948d0

c --- C---H nonbonded interaction
c         aexsix(2) = -6.290d4
c         bexsix(2) = 4.411d6
c         cexsix(2) = -3.67d0
         aexsix(2) = -6.441d4
         bexsix(2) = 5.5355d6
         cexsix(2) = -3.67d0

c --- H---H nonbonded interaction
c         aexsix(3) = -1.374d4
c         bexsix(3) = 1.335d6
c         cexsix(3) = -3.74d0
         aexsix(3) = -1.6254d4
         bexsix(3) = 1.323d6
         cexsix(3) = -3.74d0
         mass(3) = 1.0079d0

         if (lshift) then
            do i=1,natom
               sexsix(i) = aexsix(i)*(rcut**(-6))
     &              + bexsix(i)*Exp(cexsix(i)*rcut)
            enddo
         else 
            do i=1,natom
               consp(i) = (2.0d0/3.0d0)*pi*(2.0d0*aexsix(i)/(rcut
     &              *rcut*rcut)+bexsix(i)*dexp(cexsix(i)*rcut)
     &              *(-6.0d0/(cexsix(i)*cexsix(i)*cexsix(i))+6.0d0
     &              *rcut/(cexsix(i)*cexsix(i))-3.0d0*rcut*rcut/
     &              cexsix(i)+rcut*rcut*rcut))
               consu(i) = 2.0d0*pi*(aexsix(i)/(3.0d0*rcut*rcut*rcut)
     &              +(-rcut*rcut+2.0d0*rcut/cexsix(i)-2.0d0/(cexsix(i)*
     &               cexsix(i)))*bexsix(i)*dexp(cexsix(i)*rcut)/
     &               cexsix(i))
            enddo
c            write(11,*) 'consp(i)',consp
c            write(11,*) 'consu(i)',consu
         endif

         write(6,*) 
     &        ' i   aexsix       bexsix      cexsix     sexsix'
         do i = 1,natom
            write(6,'(i3,2x,4e12.4)')i,aexsix(i),bexsix(i)
     &           ,cexsix(i),sexsix(i)
         enddo

         return

      endif

      if (lmmff) then

c --- Merk Molecular Force Field (MMFF)
c --- J. Am. Chem. Soc. 1992 Vol.114 P7827-7843 (Thomas A. Halgren)
c --- natomtyp is set in mmff.inc
c     U(r) = epsi*(1.07/(rs+0.07))^7 * (1.12/(rs^7+0.12)-2)
c     rs = r / sigimmff

c --- C---C nonbonded interaction ***

         alphammff(1) = 1.050d0
         nmmff(1) = 2.490d0
         ammff(1) = 3.890d0
         gmmff(1) = 1.282d0
         sigimmff(1) = ammff(1)*dsqrt(dsqrt(alphammff(1)))
         epsimmff(1) = 45582.6d0*gmmff(1)*gmmff(1)*dsqrt(nmmff(1))
     +        /(ammff(1)**6.0d0)
c         sigimmff(1) = 1.126763255691509d0
c         epsimmff(1) = 0.9994354715851470d0
         sigisq(1) = sigimmff(1)*sigimmff(1)
         mass(1) = 12.011d0
c         mass(1) = 1.0d0
c --- H---H nonbonded interaction ***

         alphammff(3) = 0.250d0
         nmmff(3) = 0.800d0
         ammff(3) = 4.200d0
         gmmff(3) = 1.209d0
         sigimmff(3) = ammff(3)*dsqrt(dsqrt(alphammff(3)))
         epsimmff(3) = 45582.6d0*gmmff(3)*gmmff(3)*dsqrt(nmmff(3))
     +        /(ammff(3)**6.0d0)
         sigisq(3) = sigimmff(3)*sigimmff(3)
         mass(3) = 1.0078d0

c --- C---H nonbonded interaction by using cominbination rule ***
         sigimmff(2) = 0.5d0*(sigimmff(1)+sigimmff(3))*(1.0d0 + 
     +        0.2d0*(1.0d0-dexp(-12.d0*(((sigimmff(1)-sigimmff(3))/
     +        (sigimmff(3)+sigimmff(1)))**2.0d0))))
         epsimmff(2) = 91165.1d0*gmmff(1)*gmmff(3)*alphammff(1)*
     +        alphammff(3)/((dsqrt(alphammff(1)/nmmff(1))+
     +        dsqrt(alphammff(3)/nmmff(3)))*(sigimmff(2)**6.0d0))
         sigisq(2) = sigimmff(2)*sigimmff(2)


         if (lshift) then
            do i=1,natomtyp
               rs1 = rcut/sigimmff(i)
               rs7 = rs1**7.0d0
               sr7 = (1.07d0/(rs1+0.07d0))**7.0d0               
               smmff(i) = epsimmff(i)*sr7*(1.12d0/(rs7+0.12)-2)
            enddo
         else
            do i = 1,natomtyp
               smmff(i) = 0.0d0
            enddo
            coru_cons(1) = -2.4837937263569310d-02
c            coru_cons(1) = -5.8244592746534724E-03
            coru_cons(2) = -1.7583010189381791d-02
            coru_cons(3) = -8.3770412792126582d-03
            corp_cons(1) = 0.1696349613545569d0
c            corp_cons(1) = 4.0098456560842058E-02
            corp_cons(2) = 0.1203650950025348d0
            corp_cons(3) = 5.7576802340310304d-02
         endif

         write(6,*) ' i   epsimmff     sigimmff   smmff'
         do i = 1,natom
            write(6,'(i3,2x,4e12.4)')i,epsimmff(i),sigimmff(i)
     &           ,smmff(i)
         enddo

         return

      endif

      do i = 1, nntype
         sigi(i) = 0.0d0
         epsi(i) = 0.0d0
         mass(i) = 0.0d0
         qelect(i) = 0.0d0
         xiq(i) = 0.0d0
         lqchg(i) = .false.
         lij(i) = .true.
         chname(i) = '                  '
      enddo

      if (lninesix) then
c * special potential for all-atom formic acid from llnl 4/6/04 jms 

c *** sp2 carbon site H-[C](=O)-OH
         rzeronx(1) = 4.1834161d0
         epsilonnx(1) = 45.224d0
         qelect(1) = 0.44469d0
         lqchg(1) = .true.
         mass(1) = 12.011d0
         chname(1) = ' 9-6 formic acid C'

c *** hydroxyl oxygen site H-C(=O)-[O]H
         rzeronx(2) = 3.5694293136d0
         epsilonnx(2) = 47.151d0
         qelect(2) = -0.55296d0
         lqchg(2) = .true.
         mass(2) = 15.9996d0
         chname(2) = ' 9-6 formic acid O'

c *** carbonyl oxygen site H-C[=O]-OH
         rzeronx(3) = 3.0014635172d0
         epsilonnx(3) = 146.008d0
         qelect(3) = -0.43236d0
         lqchg(3) = .true.
         mass(3) = 15.9996d0
         chname(3) = ' 9-6 formic acid=O'

c *** hydrogen site [H]-C(=O)-OH
         rzeronx(4) = 0.8979696387d0
         epsilonnx(4) = 2.4054d0
         qelect(4) = 0.10732d0
         lqchg(4) = .true.
         mass(4) = 1.00794d0
         chname(4) = ' 9-6 formic acid H'

c *** acidic hydrogen site H-C(=O)-O[H] 
         rzeronx(5) = 1.115727276d0
         epsilonnx(5) = 12.027d0
         qelect(5) = 0.43331d0
         lqchg(5) = .true.
         mass(5) = 1.00794d0
         chname(5) = ' 9-6 formic acidOH'

c * calculate all site-site parameters via Lorentz-Berthelot rules
         do i = 1,nxatom
            do j = 1,nxatom
               ij = (i-1)*nxatom + j
               rzero(ij) = 0.5d0*(rzeronx(i) + rzeronx(j))
               epsnx(ij) = dsqrt(epsilonnx(i)*epsilonnx(j))
               if (lshift) then
                  shiftnsix(ij) = 4.0d0*epsnx(ij)*(rzero(ij)/rcut)**6 *
     &               (2.0d0*(rzero(ij)/rcut)**3 - 3.0d0) 
               endif
            enddo
         enddo

         qqfact = 1.67d5

         return

      endif

c ================================================
c *** Begin Lennard-Jones Site Parameter Input ***
c ===================================================c
c Notes.  Throughout the listings the atoms or       c
c pseudoatoms that each entry is for is in brackets, c
c such as [CH3] for a united-atom methyl group.      c
c A double bond is "=" and a triple bond is "=-"     c
c PLEASE ADD NEW PARAMETERS TO END OF LISTING!!!     c
c ===================================================c

c * 1-5 correction term for unprotected hydrogen-oxygen interaction
c * for ether oxygens
      a15(1) = 4.0d7
c * for alcohol oxygens
      a15(2) = 7.5d7

c * OLD VALUES
c * this is 17^6
c      a15 = 24137569.0d0
c * this is 16^6
c      a15 = 16777216.0d0

c * ALKANES

c --- SKS-UA methyl group [CH3] standard (nature paper)
      sigi(1) = 3.93d0
      epsi(1) = 114.0d0
      mass(1) = 15.0347d0
      chname(1) = ' SKS-UA CH3 alkane'

c --- SKS-UA methylene group [CH2]
      sigi(2) = 3.93d0 
      epsi(2) = 47.0d0 
      mass(2) = 14.0268d0
      chname(2) = ' SKS-UA CH2 alkane'


c --- TraPPE-UA Methane [CH4] sp3 (similar to OPLS) 
      sigi(3) = 3.73d0
      epsi(3) = 148.0d0
      mass(3) = 16.043d0
      chname(3) = ' Tr-UA CH4 alkane '

c --- TraPPE-UA Methyl [CH3] sp3 (J. Phys. Chem.) (primary)
      sigi(4) = 3.75d0
      epsi(4) = 98.0d0
      mass(4) = 15.0347d0
      chname(4) = ' Tr-UA CH3 alkane '

c --- TraPPE-UA Methylene [CH2] sp3 (secondary)
      sigi(5) = 3.95d0
      epsi(5) = 46.0d0
      mass(5) = 14.0268d0
      chname(5) = ' Tr-UA CH2 alkane '

c --- TraPPE-UA Methine [CH] sp3 (ternary)
      sigi(6) = 4.68d0
      epsi(6) = 10.0d0
      mass(6) = 13.0191d0
      chname(6) = ' Tr-UA CH  alkane '

c --- TraPPE-UA [C] (quaternary)  
      sigi(7) = 6.4d0 
      epsi(7) = 0.5d0
      mass(7) = 12.011d0
      chname(7) = ' Tr-UA C   alkane '


c --- OPLS-UA ethane methyl [CH3] sp3
      sigi(8) = 3.775d0
      epsi(8) = 104.1d0
      mass(8) = 15.0347d0
      chname(8) = ' OPLSUA CH3 ethane'

c --- OPLS-UA butane methyl [CH3] sp3
      sigi(9) = 3.905d0
      epsi(9) = 88.1d0
      mass(9) = 15.0347d0
      chname(9) = ' OPLSUA CH3 butane'

c --- OPLS-UA methylene [CH2] sp3
      sigi(10) = 3.905d0
      epsi(10) = 59.4d0
      mass(10) = 14.0268d0
      chname(10) = ' OPLSUA CH2 alkane'

c --- OPLS-UA Methine [CH] (ternary) CARBON GROUP Jorgensen 
      sigi(11) = 3.85d0
      epsi(11) = 32.0d0
      mass(11) = 13.0191d0
      chname(11) = ' OPLSUA CH alkane '

c --- OPLS [C]  (quaternary) CARBON GROUP Jorgensen
      sigi(12) = 3.85d0
      epsi(12) = 25.0d0
      mass(12) = 12.0113d0
      chname(12) = ' OPLSUA C  alkane '


c --- UA Methane [CH4] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(13) = 4.10d0
      epsi(13) = 140.0d0
      mass(13) = 16.043d0
      chname(13) = ' FreireUA CH4     '

c --- UA Methyl [CH3] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(14) = 4.02d0
      epsi(14) = 96.0d0
      mass(14) = 15.0347d0
      chname(14) = ' FreireUA CH3 alkn'

c --- UA Methylene [CH2] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(15) = 3.72d0
      epsi(15) = 57.0d0
      mass(15) = 14.0268d0
      chname(15) = ' FreireUA CH2 alkn'

c --- UA Methine [CH] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(16) = 3.36d0
      epsi(16) = 36d0
      mass(16) = 13.019d0
      chname(16) = ' FreireUA CH alkn '

c --- UA Quaternary [C] Freire Mol. Phys. 91, (2), 189-201 (1997)
      sigi(17) = 2.44d0
      epsi(17) = 9.0d0
      mass(17) = 12.011d0
      chname(17) = ' FreireUA C alkane'


c --- Mol. Phys. UA Methyl [CH3] sp3
      sigi(18) = 3.77d0
      epsi(18) = 98.1d0
      mass(18) = 15.0347d0
      chname(18) = ' MPhysUA CH3 alkM '

c --- Mol. Phys. UA METHYL-BRANCH Methyl [CH3] no tail correction
      sigi(19) = 3.93d0
      epsi(19) = 78.0d0
      mass(19) = 15.0347d0
      chname(19) = ' MPhysUA CH3 alkMB'

c --- Mol. Phys. UA ETHYL-BRANCH Methyl [CH3] no tail correction
      sigi(20) = 3.93d0
      epsi(20) = 95.0d0
      mass(20) = 15.0347d0
      chname(20) = ' MPhysUA CH3 alkEB'

c --- Mol. Phys. UA Methine [CH] sp3 (ternary)
      sigi(21) = 4.10d0
      epsi(21) = 12.0d0
      mass(21) = 13.0191d0
      chname(21) = ' MPhysUA CH alkane'


c --- TraPPE-AA for alkane methane [C]-H4 carbon
      sigi(22) = 3.31d0
      epsi(22) = 0.01d0
      mass(22) = 12.011d0
      chname(22) = ' Tr-AA [C]H4 alkan'

c --- TraPPE-AA for alkane methyl [C]-H3 carbon
      sigi(23) = 3.30d0
      epsi(23) = 4.0d0
      mass(23) = 12.011d0
      chname(23) = ' Tr-AA [C]H3 alkan'

c --- TraPPE-AA for alkane methylene [C]-H2 carbon 
      sigi(24) = 3.65d0
      epsi(24) = 5.0d0
      mass(24) = 12.011d0
      chname(24) = ' Tr-AA [C]H2 alkan'

c --- (TraPPE?)-AA for alkane methine [C]-H carbon
      sigi(25) = 4.0d0
      epsi(25) = 2.0d0
      mass(25) = 12.011d0
      chname(25) = ' Tr-AA [C]H alkane'

c --- (TraPPE?)-AA for alkane quaternary [C] carbon
      sigi(26) = 4.35d0
      epsi(26) = 1.0d0
      mass(26) = 12.011d0
      chname(26) = ' Tr-AA C quat alkn'

c --- TraPPE-AA for alkane carbon-hydrogen sigma bond H[-]C 
      sigi(27) = 3.31d0
      epsi(27) = 15.3d0
c      mass(27) = 0.0d0
      mass(27) = 1.0079d0
      chname(27) = ' Tr-AA H alkane   '


c --- TraPPE-UA? Methane [CH4] sp3 charged with polarizability  
      sigi(28) = 3.73d0
      epsi(28) = 148.0d0
c is this correct?
      mass(28) = 16.043d0
      qelect(28) = -0.572d0
      lqchg(28) = .true.
      jayself(28) = 0.5d0*117403d0
      xiq(28) = 9449.3d0
      chname(28) = ' Tr C CH4 chg pol '

c --- Methane hydrogen charged with polarizibility
      sigi(29) = 0.0d0
      epsi(29) = 0.0d0
      mass(29) = 1.0078d0
      qelect(29) = 0.143d0
      lqchg(29) = .true.
      jayself(29) = 0.5d0*177700d0
      xiq(29) = 0.0d0
      lij(29) = .false.
      chname(29) = ' Tr H CH4 chg pol '


c --- OPLS-AA for alkane methylene [C]-H2 carbon 
      sigi(30) = 3.50d0
      epsi(30) = 33.2d0
      mass(30) = 12.011d0
      qelect(30) = -0.12d0
      lqchg(30) = .true.
      chname(30) = ' OPLSAA [C]H2 alkn'

c --- OPLS AA for alkane hydrogen
      sigi(31) = 2.50d0
      epsi(31) = 15.1d0
      mass(31) = 1.0078d0
      qelect(31) = 0.06d0
      lqchg(31) = .true.
      chname(31) = ' OPLSAA H  alkane '


c --- Tildesly explicit atom methyl [C]-H3 carbon
      sigi(32) = 3.367d0
      epsi(32) = 48.8d0
      mass(32) = 12.011d0
      chname(32) = ' TildAA [C]H3 alkn'

c --- Tildesly explicit atom methylene [C]-H2 carbon
      sigi(33) = 3.367d0
      epsi(33) = 48.8d0
      mass(33) = 12.011d0
      chname(33) = ' TildAA [C]H2 alkn'

c --- Tildesly explicit atom [H] hydrogen
      sigi(34) = 2.908d0
      epsi(34) = 6.84d0
      mass(34) = 1.0079d0
      chname(34) = ' TildAA H alkane  '


c --- Lennard-Jonesium ethane [CH3-CH3]
      sigi(35) = 4.25d0
      epsi(35) = 236.0d0
      mass(35) = 30.070d0
      chname(35) = ' LJ ethane C2H6   '

c --- Lennard-Jonesium heptane [CH3-(CH2)5-CH3]
      sigi(36) = 6.08d0
      epsi(36) = 418.0d0
      mass(36) = 100.203d0
      chname(36) = ' LJ heptane C7H16 '

c --- SPECIAL LJ CHAIN fit to give phase diagram of octane J Phys Chem 98?
      sigi(37) = 2.91d0
      epsi(37) = 236.0d0
      mass(37) = 14.0268d0
      chname(37) = ' LJ CH2 octane fit'

c --- Teja heptane at 366 K [CH3-(CH2)5-CH3]
      sigi(38) = 6.0471d0
      epsi(38) = 484.76d0
      mass(38) = 100.203d0
      chname(38) = ' Teja heptane 366K'

c --- Teja heptane at 450 K [CH3-(CH2)5-CH3]
      sigi(39) = 6.0471d0
      epsi(39) = 456.82d0
      mass(39) = 100.203d0
      chname(39) = ' Teja heptane 450K'

c * PERFLUOROALKANES

c --- perfluoromethane [CF4] 
      sigi(40) = 4.13d0
c      sigi(40) = 4.18d0 (fit for critical density)
      epsi(40) = 172.95d0
      mass(40) = 88.003d0
      chname(40) = ' UA CF4           '

c --- Bin's perfluoromethane [CF4] 
      sigi(41) = 4.15d0
      epsi(41) = 175.4d0
      mass(41) = 88.003d0
      chname(41) = ' Bin UA CF4       '

c --- TraPPE-UA (ilja email 4-14-99) [CF3] group
      sigi(42) = 4.36d0
      epsi(42) = 87.0d0
      mass(42) = 69.0065d0
      chname(42) = ' TrUA CF3 Ilja    '

c --- [CF3] group iterb      
      sigi(43) = 4.35d0
      epsi(43) = 87.0d0
      mass(43) = 69.006d0
      chname(43) = ' UA CF3 iterb     '

c --- TraPPE-UA (ilja email 4-14-99) [CF2] group
      sigi(44) = 4.73d0
      epsi(44) = 27.5d0
      mass(44) = 50.0081d0
      chname(44) = ' TrUA CF2 Ilja    '


c --- Amber-AA for [C]F4 carbon (JCC 13(1992) P963)
      sigi(45) = 3.82d0/(2.0d0**(1.0d0/6.0d0))
      epsi(45) = 55.05d0
      mass(45) = 12.011d0
      qelect(45) = -0.756d0
      lqchg(45) = .true.
      chname(45) = ' AmberAA [C]F4    '

c --- Amber-AA for C[F]4 fluorine (JCC 13(1992) P963)
      sigi(46) = 3.50d0/(2.0d0**(1.0d0/6.0d0))
      epsi(46) = 30.70d0
      mass(46) = 18.9984d0
      qelect(46) = 0.189d0
      lqchg(46) = .true.
      chname(46) = ' AmberAA C[F]4    '

c --- AA for [C]F4 carbon (Surface Science 367(1996) P177)
      sigi(47) = 3.35d0
      epsi(47) = 32.73d0
      mass(47) = 12.011d0
c      qelect(47) = -0.808d0
      chname(47) = ' SurfSciAA [C]F4  '

c --- AA for C[F]4 fluorine (Surface Science 367(1996) P177)
c      sigi(48) = 2.95d0
c      epsi(48) = 37.0d0
      sigi(48) = 2.90d0
      epsi(48) = 34.3d0
      mass(48) = 18.9984d0
c      qelect(48) = 0.202d0
      chname(48) = ' SurfSciAA C[F]4  '

c --- AA for [C]F4 carbon (Nose and Klein J.Chem.Phys. 78(1983) 6928)
c      sigi(49) = 3.35d0
c      epsi(49) = 37.00d0
      sigi(49) = 3.35d0
      epsi(49) = 26.00d0
      mass(49) = 12.011d0
c      qelect(49) = -0.896d0
      chname(49) = ' NoseKleinAA [C]F4'

c --- AA for C[F]4 fluorine (Nose and Klein J.Chem.Phys. 78(1983) 6928)
      sigi(50) = 2.95d0
      epsi(50) = 38.50d0
      mass(50) = 18.9984d0
c      qelect(50) = 0.224d0
      chname(50) = ' NoseKleinAA C[F]4'

c * ALKENES

c --- TraPPE-UA [CH2] sp2 alkene Try 2D 04-15-99 MGM
      sigi(51) = 3.675d0
      epsi(51) = 85.0d0
      mass(51) = 14.0269d0
      chname(51) = ' Tr-UA CH2 alkene '

c --- TraPPE-UA [CH] sp2 alkene Try 3B 04-15-99 MGM
      sigi(52) = 3.73d0
      epsi(52) = 47.0d0
      mass(52) = 13.0191d0
      chname(52) = ' Tr-UA CH alkene  '

c --- TraPPE-UA [C] sp2 Try A 04-15-99 MGM
      sigi(53) = 3.85d0
      epsi(53) = 20.0d0
      mass(53) = 12.011d0
      chname(53) = ' Tr-UA C alkene   '


c --- OPLS-UA sp2 hybrid [CH2] group JACS 106, 6638-6646 (1984)
      sigi(54) = 3.85d0
      epsi(54) = 70.43d0
      mass(54) = 14.0269d0
      chname(54) = ' OPLSAA CH2 sp2   '

c --- OPLS-UA sp2 hybrid [CH] group JACS 106, 6638-6646 (1984)
      sigi(55) = 3.800d0
      epsi(55) = 57.85d0
      mass(55) = 13.0191d0
      chname(55) = ' OPLSAA CH sp2    '

c * AROMATICS

c --- TraPPE-UA [CH] benzene carbon
c *** maybe these three are for UA 9-site model?!?!?!?!?
c      sigi(56) = 3.74d0
c      epsi(56) = 48.0d0
c      mass(56) = 13.0191d0
c * published CH(aro) for TraPPE-UA 6-site
      sigi(56) = 3.695d0
      epsi(56) = 50.5d0
      mass(56) = 13.0191d0
      chname(56) = ' Tr-UA CH benzene6'

c --- TraPPE-UA middle benzene site
      sigi(57) = 0.0d0
      epsi(57) = 0.0d0
      mass(57) = 0.0d0
      lqchg(57) = .true.
      qelect(57) = 2.42d0
      chname(57) = ' Tr-UA mid-q benz9'
      
c --- TraPPE-UA pi electron benzene site
      sigi(58) = 0.0d0
      epsi(58) = 0.0d0
      mass(58) = 0.0d0
      lqchg(58) = .true.
      qelect(58) = -1.21d0
      chname(58) = ' Tr-UA pi-q benz9 '
      
c --- TraPPE-UA [C] tertiary aromatic carbon for toluene
      sigi(59) = 3.88d0
      epsi(59) = 21.0d0
      mass(59) = 12.011d0
      chname(59) = ' Tr-UA C arom tolu'

c --- TraPPE-UA [C] tertiary aromatic carbon for napthalene
      sigi(60) = 3.70d0
      epsi(60) = 30.0d0  
      mass(60) = 12.011d0
      chname(60) = ' Tr-UA C arom naph'

c * ALCOHOLS      

c --- TraPPE-UA alkanol hydrogen [H]-O
      sigi(61) = 0.0d0
      epsi(61) = 0.0d0
      mass(61) = 1.0079d0
      qelect(61) = 0.435d0
      lij(61) = .false.
      lqchg(61) = .true.
      chname(61) = ' Tr-UA H alkanol  '
      
c --- TraPPE-UA alkanol oxygen H-[O]-CHx
      sigi(62) = 3.02d0
      epsi(62) = 93.0d0
      mass(62) = 16.00d0
      qelect(62) = -0.700d0
      lqchg(62) = .true.
      chname(62) = ' Tr-UA O alkanol  '

c --- TraPPE-UA methanol methyl [CH3]-OH 
      sigi(63) = sigi(4)
      epsi(63) = epsi(4)
      mass(63) = mass(4)
      qelect(63) = 0.265d0
      lqchg(63) = .true.
      chname(63) = ' Tr-UA CH3 alkanol'

c --- TraPPE-UA alkanol methylene [CH2]-OH
      sigi(64) = sigi(5)
      epsi(64) = epsi(5)
      mass(64) = mass(5)
      qelect(64) = 0.265d0
      lqchg(64) = .true.
      chname(64) = ' Tr-UA CH2 alkanol'

c --- TraPPE-UA alkanol methine [CH]-OH (from Bin 6-20-00)
      sigi(65) = 4.33d0
      epsi(65) = epsi(6)
      mass(65) = mass(6)
      qelect(65) = 0.265d0
      lqchg(65) = .true.
      chname(65) = ' Tr-UA CH alkanol '

c --- TraPPE-UA alkanol quaternary carbon [C]-OH (from Bin 6-20-00)
      sigi(66) = 5.8d0
      epsi(66) = epsi(7)
      mass(66) = mass(7)
      qelect(66) = 0.265d0
      lqchg(66) = .true.
      chname(66) = ' Tr-UA C alkanol  '


c --- OPLS-UA alkanol hydrogen [H]-O
      sigi(67) = 0.0d0
      epsi(67) = 0.0d0
      mass(67) = 1.0079d0
      qelect(67) = 0.435d0
      lij(67) = .false.
      lqchg(67) = .true.
      chname(67) = ' OPLSUA H alkanol '

c --- OPLS-UA alkanol oxygen H-[O]-CHx 
      sigi(68) = 3.07d0
      epsi(68) = 85.578d0
      mass(68) = 15.999d0
      qelect(68) = -0.700d0
      lqchg(68) = .true.
      chname(68) = ' OPLSUA O alkanol '

c --- OPLS-UA alkanol methyl [CH3]-OH 
      sigi(69) = sigi(8)
      epsi(69) = epsi(8)
      mass(69) = mass(8)
      qelect(69) = 0.265d0
      lqchg(69) = .true.
      chname(69) = ' OPLSUA CH3 alknol'
    
c --- OPLS-UA alkanol methylene [CH2]-OH
      sigi(70) = sigi(10)
      epsi(70) = epsi(10)
      mass(70) = mass(10)
      qelect(70) = 0.265d0
      lqchg(70) = .true.
      chname(70) = ' OPLSUA CH2 alknol'

c * ETHERS

c --- TraPPE-UA ether oxygen  CHx-[O]-CHx
      sigi(71) = 2.80d0
      epsi(71) = 55.0d0
      mass(71) = 16.00d0
      qelect(71) = -0.50d0
      lqchg(71) = .true.
      chname(71) = ' Tr-UA O ether    '

c --- TraPPE-UA ether methyl [CH3]-O
      sigi(72) = sigi(4)
      epsi(72) = epsi(4)
      mass(72) = mass(4)
      qelect(72) = 0.25d0
      lqchg(72) = .true.
      chname(72) = ' Tr-UA CH3 ether  '

c --- TraPPE-UA ether methylene [CH2]-O
      sigi(73) = sigi(5)
      epsi(73) = epsi(5)
      mass(73) = mass(5)
      qelect(73) = 0.25d0
      lqchg(73) = .true.
      chname(73) = ' Tr-UA CH2 ether  '

c --- TraPPE-UA ether methine [CH]-O
      sigi(74) = sigi(65)
      epsi(74) = epsi(65)
      mass(74) = mass(65)
      qelect(74) = 0.25d0
      lqchg(74) = .true.
      chname(74) = ' Tr-UA CH ether   '

c --- TraPPE-UA ether quaternary carbon [C]-O
      sigi(75) = sigi(66)
      epsi(75) = epsi(66)
      mass(75) = mass(66)
      qelect(75) = 0.25d0
      lqchg(75) = .true.
      chname(75) = ' Tr-UA C ether    '

c --- TraPPE-UA Block copolymer ether oxygen next to carbonyl CH2-[O]-C=O
      sigi(76) = 2.80d0
      epsi(76) = 55.0d0
      mass(76) = 16.00d0
      qelect(76) = -0.25d0
      lqchg(76) = .true.
      chname(76) = ' Tr-UA O carbonate'

c --- TraPPE-UA Block copolymer methylene O=C-O-[CH2]
      sigi(77) = sigi(5)
      epsi(77) = epsi(5)
      mass(77) = mass(5)
      qelect(77) = 0.30d0
      lqchg(77) = .true.
      chname(77) = ' Tr-UA CH2 carbnat'


c --- OPLS-UA ether oxygen CHx-[O]-CHx (JCC 1990 vol 11, iss 8 958-971)
      sigi(78) = 3.00d0
      epsi(78) = 85.58d0
      mass(78) = 15.999d0
      qelect(78) = -0.50d0
      lqchg(78) = .true.
      chname(78) = ' OPLSUA O ether   '

c --- OPLS-UA ether methyl [CH3]-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(79) = 3.80d0
      epsi(79) = 85.58d0
      mass(79) = 15.0347d0
      qelect(79) = 0.25d0
      lqchg(79) = .true.
      chname(79) = ' OPLSUA [CH3]-O   '

c --- OPLS-UA ether methyl [CH3]-CH2-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(80) = 3.905d0
      epsi(80) = 88.06d0
      mass(80) = 15.0347d0
      chname(80) = ' OPLSUA [CH3]CH2-O'

c --- OPLS-UA ether methylene [CH2]-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(81) = 3.80d0
      epsi(81) = 59.38d0
      mass(81) = 14.0268d0
      qelect(81) = 0.25d0
      lqchg(81) = .true.
      chname(81) = ' OPLSUA [CH2]-O   '

c --- OPLS-UA THF methylene [CH2]-CH2-O (JCC 1990 vol 11, iss 8 958-971)
      sigi(82) = 3.905d0
      epsi(82) = 59.38d0
      mass(82) = 14.0268d0
      chname(82) = ' OPLSUA [CH2]CH2-O'

c * KETONES, ALDEHYDES AND ESTERS

c$$$c --- (TraPPE?)-UA ketone carbon [C]=O (jpotoff 12/13/99 + OPLS JPC v94 p1683 1990)
c$$$      sigi(83) = 3.82d0
c$$$      epsi(83) = 40.00d0
c$$$      mass(83) = 12.011d0
c$$$      qelect(83) = +0.424d0
c$$$      lqchg(83) = .true.
c --- TraPPE-UA ketone carbon [C]=O (TraPPE-6)
      sigi(83) = 3.82d0
      epsi(83) = 40.00d0
      mass(83) = 12.011d0
      qelect(83) = +0.424d0
      lqchg(83) = .true.
      chname(83) = ' Tr-UA [C]=O keton'

c$$$c --- (TraPPE?)-UA ketone oxygen C=[O] (jpotoff 12/17/99 + OPLS JPC v94 p 1683 1990)
c$$$      sigi(84) = 3.04d0
c$$$      epsi(84) = 85.0d0
c$$$      mass(84) = 15.999d0
c$$$      qelect(84) = -0.424d0 
c$$$      lqchg(84) = .true. 
c --- TraPPE-UA ketone oxygen C=[O] (TraPPE-6)
      sigi(84) = 3.05d0
      epsi(84) = 79.0d0
      mass(84) = 15.999d0
      qelect(84) = -0.424d0 
      lqchg(84) = .true. 
      chname(84) = ' Tr-UA C=[O] keton'

c$$$c --- (TraPPE?)-UA aldehyde carbon [CH]=O  (jpotoff 12/13/99 + OPLS JPC v94 p1683 1990)
c$$$      sigi(85) = 3.60d0
c$$$      epsi(85) = 58.00d0
c$$$      mass(85) = 13.011d0
c$$$      qelect(85) = +0.424d0
c$$$      lqchg(85) = .true.
c --- TraPPE-UA aldehyde carbon [CH]=O (TraPPE-6) 
      sigi(85) = 3.55d0
      epsi(85) = 65.00d0
      mass(85) = 13.019d0
      qelect(85) = +0.424d0
      lqchg(85) = .true.
      chname(85) = ' Tr-UA [CH]=O alde'
      
c --- TraPPE-UA ester methylene group [CH2]-C=O 
      sigi(86) = sigi(5)
      epsi(86) = epsi(5)
      qelect(86) = 0.05d0
      mass(86) = mass(5)
      lqchg(86) = .true.
      chname(86) = ' Tr-UA [CH2]-C=O e'

c --- TraPPE-UA ester methylene group C(=O)O-[CH2]
      sigi(87) = sigi(5)
      epsi(87) = epsi(5)
      qelect(87) = 0.25d0
      mass(87) = mass(5)
      lqchg(87) = .true.
      chname(87) = ' Tr-UA COO-[CH2] e'

c --- TraPPE-UA ester oxygen C(=O)-[O]-CHx (uses TraPPE alcohol O)
      sigi(88) = sigi(62)
      epsi(88) = epsi(62)
      qelect(88) = -0.40d0
      mass(88) = mass(62)
      lqchg(88) = .true.
      chname(88) = ' Tr-UA CO[O]-CHx e'

c --- TraPPE-UA ester oxygen in carbonyl C=[O] (uses TraPPE CO2 O)
      sigi(89) = 3.05d0
      epsi(89) = 79.0d0
      qelect(89) = -0.45d0
      mass(89) = 15.999d0
      lqchg(89) = .true.
      chname(89) = ' Tr-UA C=[O] ester'

c --- TraPPE-UA ester carbon in carbonyl [C]=O 
      sigi(90) = 3.82d0
      epsi(90) = 40.0d0
      qelect(90) = 0.55d0
      mass(90) = 12.011d0
      lqchg(90) = .true.
      chname(90) = ' Tr-UA [C]=O ester'

c * CARBOXYLIC ACIDS

c --- 91-94 are old parameters for TraPPE-UA carboxylic acids
c --- TraPPE-UA carboxylic acid hydrogen C(=O)-O-[H]
      sigi(91) = 0.0d0
      epsi(91) = 0.0d0
      mass(91) = 1.0079d0
      qelect(91) = 0.30d0
      lij(91) = .false.
      lqchg(91) = .true.
      chname(91) = 'oTr-UA COO[H] acid'

c --- TraPPE iterB carbonyl oxygen  C[=O]-O-H
c      sigi(92) = 3.0d0
c      epsi(92) = 75.0d0
c      qelect(92) = -0.440d0
      sigi(92) = 3.04d0
      epsi(92) = 81.0d0
      qelect(92) = -0.424d0
      mass(92) = 15.999d0
      lqchg(92) = .true.
      chname(92) = 'oTr-UA C[O]OH acid'

c --- TraPPE iterB carboxylic acid oxygen C(=O)-[O]-H
c      sigi(93) = 3.00d0
c      epsi(93) = 75.0d0
c      qelect(93) = -0.53d0
      sigi(93) = 3.02d0
      epsi(93) = 93.0d0
      qelect(93) = -0.30d0
      mass(93) = 15.999d0
      lqchg(93) = .true.
      chname(93) = 'oTr-UA CO[O]H acid'

c --- TraPPE iterB carbonyl carbon  [C](=O)-O-H
c      sigi(94) = 4.0d0
c      epsi(94) = 42.0d0
c      qelect(94) = 0.52d0
      sigi(94) = 3.60d0
      epsi(94) = 65.0d0
      mass(94) = 12.011d0
      qelect(94) = 0.424d0
      lqchg(94) = .true.
      chname(94) = 'oTr-UA [C]OOH acid'

c --- 95-98 new parameters for TraPPE-UA carboxylic acids 
c --- (jpotoff 12/17/99 + OPLS JPC v95 p 3315 1991)
c --- TraPPE-UA carboxylic acid carbonyl oxygen C=[O]-O-H 
      sigi(95) = 3.04d0
      epsi(95) = 81.0d0
      mass(95) = 15.999d0
      qelect(95) = -0.424d0 
      lqchg(95) = .true. 
      chname(95) = ' Tr-UA C[O]OH acid'
      
c --- TraPPE-UA iterB carbonyl carbon  [C](=O)-O-H
      sigi(96) = 3.60d0
      epsi(96) = 65.00d0
      mass(96) = 12.011d0
      qelect(96) = 0.424d0
      lqchg(96) = .true.
      chname(96) = ' Tr-UA [C]OOH acid'
      
c --- TraPPE-UA carboxylic acid hydrogen C(=O)-O-[H] (JPC v95 p. 3315, 1991)
      sigi(97) = 0.0d0
      epsi(97) = 0.0d0
      mass(97) = 1.0079d0
      qelect(97) = 0.30d0
      lij(97) = .false.
      lqchg(97) = .true.
      chname(97) = ' Tr-UA COO[H] acid'

c --- TraPPE-UA carboxylic acid oxygen C(=O)-[O]-H
      sigi(98) = 3.02d0
      epsi(98) = 93.0d0
      mass(98) = 16.00d0
      qelect(98) = -0.30d0
      lqchg(98) = .true.
      chname(98) = ' Tr-UA CO[O]H acid'

 
c --- OPLS-UA (1990) charged methyl [CH3]
      sigi(99) = 3.91d0
      epsi(99) = 80.6d0
      mass(99) = 15.0347d0
      qelect(99) = 0.080d0
      lqchg(99) = .true.
      chname(99) = ' OPLSUA CH3 acid? '

c --- OPLS-UA (1990) carboxylic acid carbon [C](=O)-O-H
      sigi(100) = 3.75d0
      epsi(100) = 52.9d0
      mass(100) = 12.011d0
      qelect(100) = 0.55d0
      lqchg(100) = .true.
      chname(100) = ' OPLSUA [C]OOH acd'

c --- OPLS-UA (1990) carboxylic acid oxygen C(=O)-[O]-H
      sigi(101) = 3.0d0
      epsi(101) = 85.6d0
      mass(101) = 15.999d0
      qelect(101) = -0.58d0
      lqchg(101) = .true.
      chname(101) = ' OPLSUA CO[O]H acd'

c --- OPLS-UA (1990) carbonyl oxygen  C[=O]
      sigi(102) = 2.96d0
      epsi(102) = 105.7d0
      mass(102) = 15.999d0
      qelect(102) = -0.5d0
      lqchg(102) = .true.
      chname(102) = ' OPLSUA C[O]OH acd'

c --- OPLS-UA (1990) and OPLS-AA (1995) carboxylic acid hydrogen C(=O)-O-[H]
      sigi(103) = 0.0d0
      epsi(103) = 0.0d0
      mass(103) = 1.0079d0
      qelect(103) = 0.45d0
      lij(103) = .false.
      lqchg(103) = .true.
      chname(103) = ' OPLSUA COO[H] acd'

c --- OPLS-AA (1995) carboxylic acid carbonyl oxygen  C[=O]-O-H
      sigi(104) = 2.96d0
      epsi(104) = 105.8d0
      mass(104) = 15.999d0
      qelect(104) = -0.440d0
      lqchg(104) = .true.
      chname(104) = ' OPLSAA C[O]OH acd'

c --- OPLS-AA (1995) carboxylic acid oxygen C(=O)-[O]-H
      sigi(105) = 3.00d0
      epsi(105) = 85.7d0
      mass(105) = 15.999d0
      qelect(105) = -0.53d0
      lqchg(105) = .true.
      chname(105) = ' OPLSAA CO[O]H acd'

c --- OPLS-AA (1995) carboxylic acid carbon  [C](=O)-O-H
      sigi(106) = 3.75d0
      epsi(106) = 52.9d0
      mass(106) = 12.011d0
      qelect(106) = 0.51d0
      lqchg(106) = .true.
      chname(106) = ' OPLSAA [C]OOH acd'

c * WATER

c --- SPC/E oxygen [O]   (simple point charge water oxygen)
      sigi(107) = 3.1655d0
      epsi(107) = 78.1958d0
      mass(107) = 16.000d0
      qelect(107) = -0.8476d0
      lqchg(107) = .true.
      chname(107) = ' SPC/E O water    '

c --- SPC/E hydrogen [H] (simple point charge Enhanced water hydrogen)
      sigi(108) = 0.0d0
      epsi(108) = 0.0d0
      mass(108) = 1.0079d0
      qelect(108) = 0.4238d0      
      lqchg(108) = .true.
      chname(108) = ' SPC/E H water    '


c --- SPC-FQ oxygen [O]   S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(109) = 3.176
      epsi(109) = 148.0d0
      mass(109) = 15.999d0
      qelect(109) = -0.672123708
      lqchg(109) = .true.
      xiq(109) = 36899.0d0
      jayself(109) = (0.5d0)*(503.2d0)*(367.0d0)
      chname(109) = ' SPC-FQ O water   '

c --- SPC-FQ hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(110) = 0.0d0
      epsi(110) = 0.0d0
      mass(110) = 1.0079d0
      qelect(110) = 0.336061854
      lij(110) = .false.
      lqchg(110) = .true.
      xiq(110) = 0.0d0
      jayself(110) = (0.5d0)*(503.2d0)*(392.2d0)
      chname(110) = ' SPC-FQ H water   '


c --- TIP4P-FQ Oxygen [O] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(111) = 3.159d0
      epsi(111) = 144.1d0
c      epsi(111) = 105.0d0
      mass(111) = 15.999d0
      chname(111) = ' TIP4P-FQ O water '

c --- TIP4P-FQ Hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(112) = 0.0d0
      epsi(112) = 0.0d0
      mass(112) = 1.0079d0
      qelect(112) = 0.35d0
      lij(112) = .false.
      lqchg(112) = .true.
      xiq(112) = 0.0d0
      jayself(112) = (0.5d0)*(503.2d0)*(353.0d0)
      chname(112) = ' TIP4P-FQ H water '

c --- TIP4P-FQ Charge [Q] S.W. Rick et al JCP 101 (7), 1 1994 6141
      sigi(113) = 0.0d0
      epsi(113) = 0.0d0
      mass(113) = 0.0d0
      qelect(113) = -0.70d0
      lij(113) = .false.
      lqchg(113) = .true.
      xiq(113) = 34464.0d0
      jayself(113) = (0.5d0)*(503.2d0)*(371.6d0)
      chname(113) = ' TIP4P-FQ M water '


c --- TIP-4P water model --- [O] site
c      sigi(114) = 3.15365d0
      sigi(114) = 3.154d0
      epsi(114) = 78.0d0
c the following value was listed under tip-4p water oxygen as well (type 152)
c      epsi(114) = 57.91d0
      mass(114) = 15.999d0
      chname(114) = ' TIP4P O water    '

c --- TIP-4P water model --- [H] site
      sigi(115) = 0.0d0
      epsi(115) = 0.0d0
      mass(115) = 1.0079d0
      qelect(115) = 0.52d0
      lij(115) = .false.
      lqchg(115) = .true.
      chname(115) = ' TIP4P H water    '

c --- TIP-4P water model --- [M] site
      sigi(116) = 0.0d0
      epsi(116) = 0.0d0
      mass(116) = 0.0d0
      qelect(116) = -1.04d0
      lij(116) = .false.
      lqchg(116) = .true.
      chname(116) = ' TIP4P M water    '


c --- TIP5P oxygen [O]
      sigi(117) = 3.12d0
      epsi(117) = 80.512d0
      mass(117) = 15.999d0
      qelect(117) = 0.0d0
      lij(117) = .true.
      lqchg(117) = .false.
      chname(117) = ' TIP5P O water    '

c --- TIP5P hydrogen [H]
      sigi(118) = 0.0d0
      epsi(118) = 0.0d0
      mass(118) = 1.0078d0
      qelect(118) = 0.241d0
      lij(118) = .false.
      lqchg(118) = .true.
      chname(118) = ' TIP5P H water    '

c --- TIP5P lone-pair [L]
      sigi(119) = 0.0d0
      epsi(119) = 0.0d0
      mass(119) = 0.0d0
      qelect(119) = -0.241d0
      lij(119) = .false.
      lqchg(119) = .true.
      chname(119) = ' TIP5P L water    '


c --- Fixed Charge Water oxygen [O] site
      sigi(120) = 3.34d0
      epsi(120) = 42.0d0
      mass(120) = 15.999d0
      qelect(120) = 6.0d0
      lqchg(120) = .true.
      chname(120) = ' FixedQ O water   '

c --- Fixed Charge Water hydrogen [H] site
      sigi(121) = 0.0d0
      epsi(121) = 0.0d0
      mass(121) = 1.0079d0
      qelect(121) = 1.0d0
      lij(121) = .false.
      lqchg(121) = .true.
      chname(121) = ' FixedQ H water   '

c --- Fixed Charge Water carbon-oxygen bond site (???hydrogen-oxygen bond???)
      sigi(122) = 2.2d0
      epsi(122) = 15.0d0
      mass(122) = 0.0d0
      qelect(122) = -2.16d0
      lqchg(122) = .true.
      chname(122) = ' FixedQ bond water'

c --- Fixed Charge Water lone pair [L] site
      sigi(123) = 0.0d0
      epsi(123) = 0.0d0
      mass(123) = 0.0d0
      qelect(123) = -1.84d0
      lij(123) = .false.
      lqchg(123) = .true.
      chname(123) = ' FixedQ L water   '

c * NOBLE GASES, CARBON MONOXIDE, CARBON DIOXIDE, NITROGEN, OXYGEN, HF

c --- TraPPE Helium (7-18-97 MGM)
      sigi(124) = 3.11d0
      epsi(124) = 4.0d0
c      sigi(124) = 3.065d0 used in JACS paper 1997
c      epsi(124) = 3.95d0
      mass(124) = 4.0026d0
c      sigi(124) = 2.556d0
c      epsi(124) = 10.2d0
c      mass(124) = 4.00d0
      chname(124) = ' TraPPE helium    '

c --- TraPPE Argon (7-18-97 MGM)
      sigi(125) = 3.390d0
      epsi(125) = 116.0d0
      mass(125) = 39.948d0
      chname(125) = ' TraPPE argon     '

c --- Krypton  
      sigi(126) = 3.607d0
      epsi(126) = 161.0d0
      mass(126) = 83.80d0
      chname(126) = ' TraPPE? krypton  '


c --- carbon in carbon monoxide [C]=-O
      sigi(127) = 3.75d0
      epsi(127) = 52.9d0
      mass(127) = 12.011d0
      qelect(127) = -0.019d0
      lqchg(127) = .true.
      chname(127) = ' carbon monoxide C'

c --- oxygen in carbon monoxide C=-[O]
      sigi(128) = 2.96d0
      epsi(128) = 105.7d0
      mass(128) = 15.999d0
      qelect(128) = 0.019d0
      lqchg(128) = .true.
      chname(128) = ' carbon monoxide O'


c --- Jeff's Amazing TraPPE CO2 model carbon [C]O2 (jpotoff 12/13/99)
      sigi(129) = 2.80d0
      epsi(129) = 27.0d0
      mass(129) = 12.011d0
      qelect(129) = 0.70d0
      lqchg(129) = .true.
      chname(129) = ' TraPPE C in CO2  '

c --- Jeff's Amazing TraPPE CO2 model oxygen C[O]2 (jpotoff 12/13/99)
      sigi(130) = 3.05d0
      epsi(130) = 79.0d0
      mass(130) = 15.999d0
      qelect(130) = -0.350d0 
      lqchg(130) = .true. 
      chname(130) = ' TraPPE O in CO2  '


c --- TraPPE carbon dioxide carbon in [C]O2-fq (jpotoff 2/15/00)
      sigi(131) = 2.80d0
      epsi(131) = 28.5d0
      mass(131) = 12.011d0
      qelect(131) = 0.6512d0
      lqchg(131) = .true.
c      xiq(131) = (503.2d0)*123.2d0
      xiq(131) = 0.0d0
      jayself(131) = (0.5d0)*(503.2d0)*(233.5d0)
      chname(131) = ' Tr-FQ C in CO2   '

c --- TraPPE carbon dioxide oxygen in C[O]2-fq (jpotoff 2/15/00)
      sigi(132) = 3.06d0
      epsi(132) = 80.5d0
      mass(132) = 15.999d0
      qelect(132) = -0.3256d0
      lqchg(132) = .true.
c      xiq(132) = (503.2d0)*201.56d0
      xiq(132) = 39430.75d0
      jayself(132) = (0.5d0)*(503.2d0)*(308.17d0)
      chname(132) = ' Tr-FQ O in CO2   '


c --- TraPPE nitrogen [N]2 (jpotoff 12/21/99)
      sigi(133) = 3.310d0
      epsi(133) = 36.00d0
      mass(133) = 14.00674d0
      qelect(133) = -0.50d0
      lqchg(133) = .true.
      chname(133) = ' TraPPE N in N2   '

c --- TraPPE nitrogen COM charge cite for N2 (jpotoff 12/21/99)
      sigi(134) = 0.0d0
      epsi(134) = 0.0d0
      mass(134) = 0.0d0
      qelect(134) = 1.0d0
      lij(134) = .false.
      lqchg(134) = .true.
      chname(134) = ' TraPPE COM in N2 '


c --- Tildesley nitrogen [N]2 
      sigi(135) = 3.31d0
      epsi(135) = 37.3d0
      mass(135) = 14.00674d0
      chname(135) = ' Tild. N in N2    '


c --- TraPPE oxygen [O]2  Final parameter adjust 8-5-98
      sigi(136) = 3.07d0
      epsi(136) = 49.0d0
      mass(136) = 15.999d0
      chname(136) = ' TraPPE O in O2   '


c --- OPLS hydrogen fluoride (HF) fluorine H-M-[F]
      sigi(137) = 2.984d0
      epsi(137) = 75.75d0
      mass(137) = 18.9984d0
      qelect(137) = 0.725d0
      lqchg(137) = .true.
      chname(137) = ' OPLS F in HMF    '

c --- OPLS HF hydrogen [H]-M-F
      sigi(138) = 0.0d0
      epsi(138) = 0.0d0
      mass(138) = 1.0078d0
      qelect(138) = 0.725d0
      lij(138) = .false.
      lqchg(138) = .true.
      chname(138) = ' OPLS H in HMF    '

c --- OPLS HF M site H-[M]-F
      sigi(139) = 0.0d0
      epsi(139) = 0.0d0
      mass(139) = 0.0d0
      qelect(139) = -1.45d0
      lij(139) = .false.
      lqchg(139) = .true.
      chname(139) = ' OPLS M in HMF    '

c * THIOLS, THIOETHERS

c --- TraPPE-UA dimethyl sulfide methyl group [CH3]-S-CH3
      sigi(140) = sigi(4)
      epsi(140) = epsi(4)
      mass(140) = mass(4)
      qelect(140) = 0.235d0
      lqchg(140) = .true.
      chname(140) = ' Tr-UA CH3 thioeth'
      
c --- TraPPE-UA dimethyl sulfide sulfur CH3-[S]-CH3  
c     (1/25/00, based on JPC v90, p6379, 1986)
      sigi(141) = 3.52d0
      epsi(141) = 158.0d0
      mass(141) = 32.07d0
      qelect(141) = -0.47d0
      lqchg(141) = .true.
      chname(141) = ' Tr-UA S thioether'

c --- TraPPE-UA methyl group [CH3]-S-H
      sigi(142) = sigi(4)
      epsi(142) = epsi(4)
      mass(142) = mass(4)
      qelect(142) = 0.18d0
      lqchg(142) = .true.
      chname(142) = ' Tr-UA CH3 thiol  '
      
c --- TraPPE-UA sulfur CH3-[S]-H 
c     (1/25/00, based on JPC v90, p6379, 1986)
      sigi(143) = 3.62d0
      epsi(143) = 185.0d0
      mass(143) = 32.07d0
      qelect(143) = -0.45d0
      lqchg(143) = .true.
      chname(143) = ' Tr-UA S thiol    '
      
c --- TraPPE-UA hydrogen CH3-S-[H]      
      sigi(144) = 0.0d0
      epsi(144) = 0.0d0
      mass(144) = 1.0079d0
      qelect(144) = 0.27d0
      lij(144) = .false.
      lqchg(144) = .true.
      chname(144) = ' Tr-UA H thiol    '
 
c --- TraPPE-UA methylene group CH3-[CH2]-S-H
      sigi(145) = sigi(5)
      epsi(145) = epsi(5)
      mass(145) = mass(5)
      qelect(145) = 0.18d0
      lqchg(145) = .true.
      chname(145) = ' Tr-UA CH2 thiol  '

c * AMINES 

c --- parameters for primary amines (2/28/00) based on JPC v94, p1683, 1990)
c --- TraPPE-UA methyl amine hydrogen CH3-N-[H]-H
      sigi(146) = 0.0d0
      epsi(146) = 0.0d0
      mass(146) = 1.0079d0
      qelect(146) = 0.275d0
      lij(146) = .false.
      lqchg(146) = .true.
      chname(146) = ' Tr-UA CH3-N[H]2  '

c --- TraPPE-UA methyl amine nitrogen CH3-[N]-H2
      sigi(147) = 3.31d0
      epsi(147) = 165.0d0
      mass(147) = 14.00674d0
      qelect(147) = -0.65d0
      lqchg(147) = .true.
      chname(147) = ' Tr-UA CH3-[N]H2  '

c --- TraPPE-UA methyl amine methyl [CH3]-N-H2
      sigi(148) = sigi(4)
      epsi(148) = epsi(4)
      mass(148) = mass(4)
      qelect(148) = 0.10d0
      lqchg(148) = .true.
      chname(148) = ' Tr-UA [CH3]-NH2  '

c --- TraPPE-UA dimethyl amine nitrogen CH3-[N]-CH3-H 
      sigi(149) = 3.31d0
      epsi(149) = 115.0d0
      mass(149) = 14.00674d0
      qelect(149) = -0.75d0
      lqchg(149) = .true.
      chname(149) = ' Tr-UA (CH3)2[N]H '

c --- TraPPE-UA trimethyl amine nitrogen CH3-[N]-CH3-CH3
      sigi(150) = 3.31d0
      epsi(150) = 115.0d0
      mass(150) = 14.00674d0
      qelect(150) = -0.60d0
      lqchg(150) = .true.           
      chname(150) = ' Tr-UA (CH3)3[N]  '

c * NITRILES

c --- TraPPE-UA nitrile nitrogen C=-[N]
      sigi(151) = 2.95d0
      epsi(151) = 60.0d0
      mass(151) = 14.007d0
      qelect(151) = -0.398d0
      lqchg(151) = .true.
      chname(151) = ' Tr-UA N nitrile  '

c --- TraPPE-UA nitrile carbon [C]=-N
      sigi(152) = 3.55d0
      epsi(152) = 60.0d0
      mass(152) = 12.011d0
      qelect(152) = 0.129d0
      lqchg(152) = .true.
      chname(152) = ' Tr-UA C nitrile  '

c --- TraPPE hydrogen cyanide hydrogen [H]-C=-N
c TRIAL VALUES
      sigi(153) = 0.0d0
      epsi(153) = 0.0d0
      mass(153) = 1.0079d0
      qelect(153) = 0.269d0
      lij(153) = .false.
      lqchg(153) = .true.
      chname(153) = ' Tr-UA H in HCN   '

c --- TraPPE-UA acetonitrile methyl [CH3]-C=-N
      sigi(154) = sigi(4)
      epsi(154) = epsi(4)
      mass(154) = mass(4)
      qelect(154) = 0.269d0
      lqchg(154) = .true.
      chname(154) = ' Tr-UA CH3 nitrile'

c --- TraPPE-UA alkyl nitrile methylene R-[CH2]-C=-N
      sigi(155) = sigi(5)
      epsi(155) = epsi(5)
      mass(155) = mass(5)
      qelect(155) = 0.150d0
      lqchg(155) = .true.
      chname(155) = ' Tr-UA CH2 nitrile'


c --- OPLS-UA nitrile nitrogen C=-[N]
      sigi(156) = 3.20d0
      epsi(156) = 85.51d0
      mass(156) = 14.007d0
      qelect(156) = -0.430d0
      lqchg(156) = .true.
      chname(156) = ' OPLSUA N nitrile '

c --- OPLS-UA nitrile carbon R-[C]=-N
      sigi(157) = 3.65d0
      epsi(157) = 75.53d0
      mass(157) = 12.011d0
      qelect(157) = 0.280d0
      lqchg(157) = .true.
      chname(157) = ' OPLSUA C nitrile '

c --- OPLS-UA acetonitrile methyl [CH3]-C=-N
      sigi(158) = 3.775d0
      epsi(158) = 104.16d0
      mass(158) = 15.035d0
      qelect(158) = 0.15d0
      lqchg(158) = .true.
      chname(158) = ' OPLSUA CH3 nitril'


c --- McDonald UA nitrile nitrogen R-C=-[N]
      sigi(159) = 3.3d0
      epsi(159) = 50.0d0
      mass(159) = 14.007d0
      qelect(159) = -0.398d0
      lqchg(159) = .true.
      chname(159) = ' McDUA N nitrile  '

c --- McDonald UA nitrile carbon R-[C]=-N
      sigi(160) = 3.4d0
      epsi(160) = 50.0d0
      mass(160) = 12.011d0
      qelect(160) = 0.129d0
      lqchg(160) = .true.
      chname(160) = ' McDUA C nitrile  '

c --- McDonald UA acetonitrile methyl [CH3]-C=-N
      sigi(161) = 3.6d0
      epsi(161) = 191.0d0
      mass(161) = 15.035d0
      qelect(161) = 0.269d0
      lqchg(161) = .true.
      chname(161) = ' McDUA CH3 nitrile'

c * CHARMM

c --- Charmm C2 (methylene group carbon)
      sigi(162) = 3.8754d0
      epsi(162) = 19.6257d0
      mass(162) = 0.003d0
      chname(162) = ' CHARMM C2 ???    '     

c --- Charmm H (hydrogen)
      sigi(163) = 2.4500d0
      epsi(163) = 19.1225d0
      mass(163) = 1.0078d0
      chname(163) = ' CHARMM H  ???    '

c --- Charmm O (bound with 2 single bonds)
      sigi(164) = 2.8598d0
      epsi(164) = 114.7348d0
      mass(164) = 16.00d0
      chname(164) = ' CHARMM O sp3 ??? '

c --- Charmm P
      sigi(165) = 3.7418d0
      epsi(165) = 100.6446d0
      mass(165) = 0.003d0
      chname(165) = ' CHARMM P ???     '

c --- Charmm O' (bound with a double bond)
      sigi(166) = 2.8598d0
      epsi(166) = 114.7348d0
      mass(166) = 16.00d0
      chname(166) = ' CHARMM P ???     '

c --- Charmm N3 (tertiary ammonia)
      sigi(167) = 3.5012d0
      epsi(167) = 84.0382d0
      mass(167) = 0.003d0
      chname(167) = ' CHARMM P ???     '

c --- Charmm C3 (methyl group carbon)
      sigi(168) = 3.8754d0
      epsi(168) = 19.6257d0
      mass(168) = 0.003d0
      chname(168) = ' CHARMM C3 ???    '

c --- Charmm C1 (ternary carbon)
      sigi(169) = 3.8754d0
      epsi(169) = 19.6257d0
      mass(169) = 0.003d0
      chname(169) = ' CHARMM C1 ???    '

c --- Charmm C' (carboxylic head group carbon)
      sigi(170) = 3.6170d0
      epsi(170) = 74.4770d0
      mass(170) = 0.003d0
      chname(170) = ' CHARMM C\' ???    '


c * ALL-ATOM NITRILES

c --- TraPPE-AA nitrile nitrogen C=-[N]
      sigi(171) = 2.95d0
      epsi(171) = 60.0d0
      mass(171) = 14.007d0
      qelect(171) = -0.398d0
      lqchg(171) = .true.
      chname(171) = ' Tr-AA N nitrile  '

c --- TraPPE-AA nitrile carbon [C]=-N
      sigi(172) = 3.55d0
      epsi(172) = 60.0d0
      mass(172) = 12.011d0
      qelect(172) = 0.129d0
      lqchg(172) = .true.
      chname(172) = ' Tr-AA C nitrile  '

c --- TraPPE-AA acetonitrile methyl carbon H3[C]-C=-N
      sigi(173) = 3.3d0
      epsi(173) = 4.0d0
      mass(173) = 12.011d0
      qelect(173) = 0.269d0
      lqchg(173) = .true.
      chname(173) = ' Tr-AA H3[C]-C=-N '

c --- TraPPE-AA acetonitrile methyl hydrogen C[H3]-C=-N
      sigi(174) = sigi(27)
      epsi(174) = epsi(27)
      mass(174) = mass(27)
      chname(174) = ' Tr-AA [H]3C-C=-N '

c * LEFTOVER PIECES

C --- Monica's alcohol methyl Not bonded to O (-CH2-) - CH3
C     (van Leeuwen JPC 99, 1831 (1995))
c      sigi() = 3.93d0
c      epsi() = 110.0d0
c      mass() = 15.0347d0
c
C- Dummy methylene for Sciece Paper MGM 12-17-97
c      sigi() = 3.95d0
c      epsi() = 55.0d0
c      mass() = 14.0268d0
c
c --- AA for alkane methyl (CH3) carbon ---
c --- if 0.45 as C-H bond length
c      sigi() = 3.44d0
c      epsi() = 13.0d0
c      qelect() = -0.572d0
c      lqchg() = .true.
c
cc --- AA for alkane carbon ---
c      sigi() = 3.65d0
c      epsi() = 5.0d0
c      mass() = 12.011d0
c      qelect() = 0.265d0
c      lqchg() = .true.
c      lij() = .true.

c ===========================
c *** End Parameter Input ***
c ===========================

c --- Computation of un-like interactions
      if ( ljoe ) then
C --- STANDARD METHYL GROUP
         extc12(1) = 3.41d7
         extc3(1)  = 20800.0d0
         extz0(1)  = 0.86d0

C --- STANDARD METHYLENE GROUP
         extc12(5) = 2.80d7
         extc3(5)  = 17100.0d0
         extz0(5)  = 0.86d0

C --- Methane
         extc12(3) = 3.41d7
         extc3(3)  = 20800.0d0
         extz0(3)  = 0.86d0

C --- Martin's methyl (CH3)
         extc12(18) = 3.41d7
         extc3(18)  = 20800.0d0
         extz0(18)  = 0.86d0
      endif

c --- Assign jayq for pairs
      do i = 1,nntype
         do j = 1,nntype
            ij = (i-1)*nntype + j
            jayq(ij) = 0.0d0
         enddo
      enddo

c - CO2-FQ Carbon-Oxygen cross term (JCO)
      i = 131
      j = 132
      djay = (503.2d0)*(133.905d0)
      ij = (i-1)*nntype + j
      ji = (j-1)*nntype + i
      jayq(ij) = djay
      jayq(ji) = djay

c - CO2-FQ Oxygen-Oxygen cross term (JOO)
      i = 132
      j = 132
      djay = (503.2d0)*(1.09d0)
      ij = (i-1)*nntype + j
      jayq(ij) = djay

c --- SPC-FQ water Oxygen-Hydrogen cross term
      i = 109
      j = 110
      djay = (503.2d0)*(276.0d0)
      ij = (i-1)*nntype + j
      ji = (j-1)*nntype + i
      jayq(ij) = djay
      jayq(ji) = djay
      
c --- SPC-FQ water Hydrogen-Hydrogen cross term
      i = 110
      j = 110
      djay = (503.2d0)*(196.0d0)
      ij = (i-1)*nntype + j
      jayq(ij) = djay

c --- TIP4P water Charge-Hydrogen cross term
      i = 112
      j = 113
      djay = (503.2d0)*(286.4d0)
      ij = (i-1)*nntype + j
      ji = (j-1)*nntype + i
      jayq(ji) = djay
      jayq(ij) = djay

c --- TIP4P water Hydrogen-Hydrogen cross term
      i = 112
      j = 112
      djay = (503.2d0)*(203.6d0)
      ij = (i-1)*nntype + j
      jayq(ij) = djay

c --- Methane C-H cross term
      i = 28
      j = 29
      ij = (i-1)*nntype + j
      ji = (j-i)*nntype + i
      jayq(ji) = 114855.0d0
      jayq(ij) = 114855.0d0

c --- Methane H-H cross term
      i = 29
      j = 29
      ij = (i-1)*nntype + j
      jayq(ij) = 112537.0d0

c *** convert input data to program units ***
      if ( lsami ) then
         call susami
         rcheck = 2.5d0 * 3.527d0
         if ( rcut .ne. rcheck ) then
            write(6,*) 'WARNING ### rcut set to 2.5sigma for SAMI'
            rcut = rcheck
         endif
      else
c *** calculate square sigmas and epsilons for lj-energy subroutines ***
         do i = 1, nntype
            do j = 1, nntype
               ij = (i-1)*nntype + j
               if ( lspecial(ij) ) then
                  write(6,*) 'ij,lspecial(ij)',ij,lspecial(ij)
                  adum = aspecial(ij)
                  bdum = bspecial(ij)
               else
                  adum = 1.0d0
                  bdum = 1.0d0
               endif
c --- Lorentz-Berthelot rules --- sig_ij = 0.5 [ sig_i + sig_j ]
               if ( lmixlb ) then
                  sig2ij(ij) =(adum* 0.5d0 * ( sigi(i) + sigi(j) ) )**2
                  if ( sigi(i) .eq. 0.0d0 .or. sigi(j) .eq. 0.0d0 )
     &                 sig2ij(ij) = 0.0d0
                  epsij(ij) = bdum*dsqrt( epsi(i) * epsi(j) )
               endif
c --- Jorgensen mixing rules --- sig_ij = [ sig_i * sig_j ]^(1/2)
               if(lmixjo) then
                  sig2ij(ij) = adum*adum*sigi(i) * sigi(j)
                  epsij(ij) = bdum*dsqrt( epsi(i) * epsi(j) )
               endif

	       if (lshift) then
                   sr2 = sig2ij(ij) / (rcut*rcut)
                   sr6 = sr2 * sr2 * sr2
                   ecut(ij)= sr6*(sr6-1.0d0)*epsij(ij)
	       endif

            enddo
         enddo
      endif
c ---  Conversion factor for intermolecular coulomb interactions
      qqfact = 1.67d5

      return
      end





