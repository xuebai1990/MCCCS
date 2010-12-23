      subroutine suvibe

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
!$$$      include 'connect.inc'
!$$$      include 'conver.inc'
!$$$      include 'contorsion.inc'

      integer(KIND=normal_int)::i

! ----------------------------------------------------------------

! - TraPPE-UA- fixed bond length for C-C bonds -
      brvib(1) = 1.54d0
      brvibk(1) = 0.0d0

! - bond stretching for C-C bonds -
!      brvib(2) = 1.54d0
!      brvibk(2) = 226450.3d0
! - bond stretching for C-C bonds - CHARMM
      brvib(2) = 1.54d0
      brvibk(2) = 237018.0d0

! - fixed bond length for C - C double bonds (not aromatic)
! - from table 1.6 of McMurry Organic Chemistry 3rd edition
      brvib(3) = 1.33d0
      brvibk(3) = 0.0d0

! - fixed bond length for sp2 C - C sp3 single bonds (not aromatic)
! - from Jorgensen JACS 106 6638-6646 (1984)
      brvib(4) = 1.50d0
      brvibk(4) = 0.0d0

! - fixed bond length for H-C bond Williams 
      brvib(5) = 1.04d0
      brvibk(5) = 0.0d0

! - fixed bond length for C-C bond Williams, also Freire UA C-C, OPLS-UA 
      brvib(6) = 1.53d0
      brvibk(6) = 0.0d0

! - fixed bond length for C-CH3 in neopentane Freire Mol Phys 1997
      brvib(7) = 1.55d0
      brvibk(7) = 0.0d0

! - fixed bond length for C-CH2 in neohexane and CHCH 23 dimethyl butane
! -  Freire Mol Phys 1997
      brvib(8) = 1.56d0
      brvibk(8) = 0.0d0

! - bond stretching for C-COOH bonds - CHARMM
      brvib(10) = 1.52d0
      brvibk(10) = 0.0d0
!      brvibk(10) = 239627.0d0

! - bond stretching for C=O bonds in COOH - CHARMM
      brvib(11) = 1.23d0
      brvibk(11) = 0.0d0
!      brvibk(11) = 619288.2d0

! - bond stretching for C(CO)-OH bonds - CHARMM
      brvib(12) = 1.37d0
      brvibk(12) = 402578.3d0

!     AMBER94 Cornell et. at.  JACS 117, 5179-5197 (1995)

! - fixed bond length for C==O bonds in COOH - AMBER94
      brvib(13) = 1.229d0
      brvibk(13) = 0.0d0

! - fixed bond length for C--O bonds in COOH - AMBER 94
      brvib(14) = 1.364d0
      brvibk(14) = 0.0d0

! - fixed bond length for sp2C-sp3C bonds - AMBER 94
      brvib(15) = 1.522d0
      brvibk(15) = 0.0d0

! - fixed bond length for O-H bonds in COOH and COH - AMBER 94
      brvib(16) = 0.96d0
      brvibk(16) = 0.0d0

! - fixed bond length for C-O bonds in alkanol 
      brvib(17) = 1.41d0
      brvibk(17) = 0.0d0

! - fixed bond length for H-O bonds in water 
      brvib(18) = 0.9572d0
      brvibk(18) = 0.0d0

! - OPLS AA bond length for c-c in alkanes
!      brvib(19) = 1.529d0
!      brvibk(19) = 0.0d0
      brvib(19) = 1.535d0
      brvibk(19) = 0.0d0

! - MMFF bond length for c-c in alkanes
!      brvib(19) = 1.508d0
!      brvibk(19) = 0.0d0

! - OPLS AA bond length for c-h in alkanes
!      brvib(20) = 1.09d0
!      brvibk(20) = 0.0d0
      brvib(20) = 0.55d0
      brvibk(20) = 0.0d0

! - MMFF bond length for c-h in alkanes
!      brvib(20) = 1.093
!      brvibk(20) = 0.0d0

! - OPLS AA bond length for h-h from CRC 72nd Ed.
!      brvib(21) = 0.74611d0
      brvib(21) = 1.535d0/2.0d0
      brvibk(21) = 0.0d0

! --- AA for CF4 bond length C-F (Surface Science 367(1996) P177) ---
      brvib(22) = 1.37d0
      brvibk(22) = 0.0d0

! --- AA for CF4 bond length C-F(Nose and Klein J.Chem.Phys. 78(1983) 6928) ---
      brvib(23) = 1.323d0
      brvibk(23) = 0.0d0

! --- Amber AA for CF4 bond length C-F (JCC 13(1992) P963) ---
      brvib(24) = 1.38d0
      brvibk(24) = 0.0d0

! --- SPC-FQ JCP 101, (7) 1 1994 6141
      brvib(25) = 1.0d0
      brvibk(25) = 0.0d0

! --- TIP4P OH bond length
      brvib(26) = 0.9572d0
      brvibk(26) = 0.0d0

! --- TIP4P OM length
      brvib(27) = 0.15d0
      brvibk(27) = 0.0d0

! --- Fixed bond length for O-O in dioxygen
!     J Chem Phys 98 (12) 9895--9904 Muller-Plathe et al
      brvib(28) = 1.21d0
      brvibk(28) = 0.0d0

! - TraPPE fixed bond length for O-H bonds in COH - from OPLS
      brvib(29) = 0.945d0
      brvibk(29) = 0.0d0

! - TraPPE fixed bond length for C-O bonds in alkanol - from OPLS
      brvib(30) = 1.43d0
      brvibk(30) = 0.0d0

! --- Fixed bond length for N-N in dinitrogen
!     53rd ed CRC Handbook page F-180, also used for C-H distance
      brvib(31) = 1.10d0
      brvibk(31) = 0.0d0

! --- TraPPE C-O bond length for CO2
       brvib(32) = 1.160d0
       brvibk(32) = 0.0d0

! --- SPECIAL LJ CHAIN J Phys Chem fit to give octane phase diagram
      brvib(33) = 3.2664d0
      brvibk(33) = 0.0d0

! --- N - charge site bond length for N2 w/quadrupole
      brvib(34) = 0.55d0
      brvibk(34) = 0.0d0

! -- CH3-C in ketones, also carboxylic acids
      brvib(35) = 1.52d0
      brvibk(35) = 0.0d0

! --- C=O bond length in carboxylic acid (OPLS)
      brvib(36) = 1.214d0
      brvibk(36) = 0.0d0
      
! -- C-O in carboxylic acids
      brvib(37) = 1.364d0
      brvibk(37) = 0.0d0

! -- O-H in carboxylic acids
      brvib(38) = 0.970d0
      brvibk(38) = 0.0d0

! -- CH3-S in thiols
      brvib(39) = 1.82d0
      brvibk(39) = 0.0d0

! -- S-H in thiols
      brvib(40) = 1.34d0
      brvibk(40) = 0.0d0

! -- C-N in amines (TraPPE-7)
      brvib(41) = 1.448d0
      brvibk(41) = 0.0d0

! -- N-H in amines (TraPPE-7)
      brvib(42) = 1.01d0
      brvibk(42) = 0.0d0

! -- C-N tertiary amines
      brvib(43) = 1.47d0
      brvibk(43) = 0.0d0

! -- TraPPE-UA C-O in ethers from OPLS
      brvib(44) = 1.41d0
      brvibk(44) = 0.0d0

! --- H-F bond length in HF
      brvib(45) = 0.917d0
      brvibk(45) = 0.0d0

! --- M site bond length in HF
      brvib(46) = 0.166d0
      brvibk(46) = 0.0d0

! --- C(arom)-C(arom) fixed bond length
      brvib(47) = 1.4d0
      brvibk(47) = 0.0d0

! --- TraPPE-UA CH3-[C=-N] triple bond length (from OPLS)
      brvib(48) = 1.157d0
      brvibk(48) = 0.0d0

! --- OPLS [CH3-C]=-N bond length -not- for TraPPE-UA
      brvib(49) = 1.458d0
      brvibk(49) = 0.0d0

! --- CH3-CH2 bond length in OPLS-UA Ether model
      brvib(50) = 1.516d0
      brvibk(50) = 0.0d0

! --- dinitrogen bond length (Tildesley?) From Nose + Klein Mol Phys 1983, 50, 1055.
      brvib(51) = 1.098d0
      brvibk(51) = 0.0d0

! * formic acid model from llnl 4/6/04 jms
! --- C-O bond length in 9-6 formic acid model
      brvib(52) = 1.376d0
      brvibk(52) = 237518.0d0

! --- C=O bond length in 9-6 formic acid model
      brvib(53) = 1.188d0
      brvibk(53) = 558190.5d0

! --- C-H bond length in 9-6 formic acid model
      brvib(54) = 1.075d0
      brvibk(54) = 227599.8d0

! --- O-H bond length in 9-6 formic acid model
      brvib(55) = 1.038d0
      brvibk(55) = 317466.8d0

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! -- added 7/12/06 for nitro group and rigid aromatic ring

! -- C-N bond length for nitro
      brvib(56) = 1.49d0
      brvibk(56) = 0.0d0

! -- N-O bond length for nitro
      brvib(57) = 1.225d0
      brvibk(57) = 0.0d0

! -- center-pi site, 9 site benzene
      brvib(58) = 0.785d0
      brvibk(58) = 0.0d0

! -- C4-C5 oNT
      brvib(59) = 1.396d0
      brvibk(59) = 0.0d0

! -- C4-C9 oNT
      brvib(60) = 1.4075d0
      brvibk(60) = 0.0d0

! -- C4-C5 mNT
      brvib(61) = 1.3908d0
      brvibk(61) = 0.0d0

! -- C4-C9 mNT
      brvib(62) = 1.3927d0
      brvibk(62) = 0.0d0

! -- N-O1 oNT
      brvib(63) = 1.2271d0
      brvibk(63) = 0.0d0

! -- N-O2 oNT
      brvib(64) = 1.2272d0
      brvibk(64) = 0.0d0

! -- N-O1 mNT
      brvib(65) = 1.2265d0
      brvibk(65) = 0.0d0

! -- N-O2 mNT
      brvib(66) = 1.2263d0
      brvibk(66) = 0.0d0

! -- C4-N oNT
      brvib(67) = 1.4692d0
      brvibk(67) = 0.0d0

! -- C4-N mNT
      brvib(68) = 1.4699d0
      brvibk(68) = 0.0d0

! -- bond lengths for Neimark DMMP JPCA v108, 1435 (2004)

!  -- P=O
      brvib(70) = 1.458d0
      brvibk(70) = 0.0d0

! -- P-CH3
      brvib(71) = 1.795d0
      brvibk(71) = 0.0d0

! -- P-O(CH3)
      brvib(72) = 1.586d0
      brvibk(72) = 0.0d0

! -- O-CH3
      brvib(73) = 1.418d0
      brvibk(73) = 0.0d0

! -- end parameters for Neimark DMMP


! -- Starting floro-alkanes (Neeraj)

! --  C-F All atom floroalkanes
      brvib(80)  = 1.33d0
      brvibk(80) = 0.0d0

! --  C-C All atom Floroalkanes

      brvib(81) = 1.54d0
      brvibk(81) = 0.0d0

! --  C-H all atom Fluoroalkanes
      brvib(82) = 1.08d0
      brvibk(82)= 0.0d0

! -- C-Cl bond for Chloroflouroalkanes
      brvib(83) = 1.754
      brvibk(83) = 0.0d0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! - parameters for acrylates
! - some values already listed elsewhere - will fix that later

! - bond length for CHx-O(ether) bonds
      brvib(100) = 1.41d0
      brvibk(100) = 0.0d0

! - bond length for O(ether)-C(carbonyl) bonds
      brvib(101) = 1.344d0
      brvibk(101) = 0.0d0

! - bond length for C=O for carboylate esters
      brvib(102) = 1.20d0
      brvibk(102) = 0.0d0

! - bond length for C-CH bonds
      brvib(103) = 1.52d0
      brvibk(103) = 0.0d0

! - bond length for C=C bonds
      brvib(104) = 1.33d0
      brvibk(104) = 0.0d0
      
! - bond length for CHx-CHy bonds (#1)
      brvib(105) = 1.54d0
      brvibk(105) = 0.0d0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! - test parameters for coarse-grain model
      brvib(110) = 3.25d0
      brvibk(110) = 0.0d0

      brben(110) = 150.0d0
      brbenk(110) = 1.0d-12
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!--- JLR 11-11-09
!--- adding vibrational parameters for RPLC simulation
! --- Si-O in silica
      brvib(111) = 1.61d0
      brvibk(111) = 0.0d0
! --- Si-CHx in alkylsilane
      brvib(112) = 1.90d0
      brvibk(112) = 0.0d0
!--- END JLR 11-11-09

! * unused 
!
! * methanol O-H from H2 parameter set from Monica's dissertation *
!      brvib() = 1.0285d0
!      brvibk() = 0.0d0
!
! * methanol C-O from H2 parameter set from Monica's dissertation *
!      brvib() = 1.4175d0
!      brvibk() = 0.0d0

! -- C--C TATB
       brvib(120) = 1.442d0
       brvibk(120) = 0.5d0*1400.0d0*503.25
! -- C--NO2 TATB
       brvib(121) = 1.419d0
       brvibk(121) = 0.5d0*1400.0d0*503.25
! -- C--NH2 TATB
       brvib(122) = 1.314d0
       brvibk(122) = 0.5d0*1400.0d0*503.25
! -- N--O TATB
       brvib(123) = 1.243d0
       brvibk(123) = 0.5d0*1400.0d0*503.25
! -- N--H TATB
       brvib(124) = 1.0d0
       brvibk(124) = 0.5d0*700.0d0*503.25


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! - TraPPE-UA bond angle for alkane segment centered at methylene (CH2 sp3)-
      brben(1) = 114.0d0
!     write(iou,*) '***** brben', brben(1) * raddeg
      brbenk(1) = 31250.0d0

! - TraPPE-UA bond angle for alkane segment centered at ternary (CH sp3)-
! - also used by Freire for most alkane bond angles
      brben(2) = 112.0d0
      brbenk(2) = 31250.0d0

! - TraPPE-UA bond angle for alkane segment centered at quaternary (C sp3)-
! - basically the same as used by Freire for neopentane and 
!                                        2,3-dimethylbutane
      brben(3) = 109.47d0
      brbenk(3) = 31250.0d0

! - bond angle for C=C-C -  segment centered at (CH sp2)
! - taken from Amber JACS (1995) 117, 5179-5197
      brben(4) = 119.70d0
      brbenk(4) = 35210.0d0

! - bond angle for cis C=C-C NO LONGER USED 4-13-99
      brben(5) = 127.4d0
      brbenk(5) = 31250.0d0

! - bond angle for H-C[methyl]-C[methylene] Ryckaert 
      brben(6) = 112.49d0
      brbenk(6) = 0.0d0

! - bond angle for C-C-C Van der Ploeg
      brben(7) = 112.0d0
      brbenk(7) = 31278.0d0

! - bond angle for H-C[methylene]-H Ryckaert
      brben(8) = 106.0d0
      brbenk(8) = 0.0d0 

! - bond angle for C-C=O (in carboxylic headgroup) - taken from Charmm
!      brben(10) = 120.0d0
! --- OPLS value
      brben(10) = 125.0d0
      brbenk(10) = 68438.0d0

! - bond angle for C-C-O (in carboxylic headgroup) - taken from Charmm
      brben(11) = 110.0d0
      brbenk(11) = 12359.0d0

! - bond angle for C2-P-O (in DPPC) - taken from Charmm
      brben(12) = 120.5d0
      brbenk(12) = 0.0d0

! - bond angle for C2-C1-H (in DPPC) - taken from Charmm
      brben(13) = 110.0d0
      brbenk(13) = 0.0d0

! - bond angle for H-C2-H (in DPPC) - taken from Charmm
      brben(14) = 106.4d0
      brbenk(14) = 0.0d0

! - bond angle for H-C2-O (in DPPC) - taken from Charmm
      brben(15) = 109.5d0
      brbenk(15) = 0.0d0

! - bond angle for H-C2-C1 (in DPPC) - taken from Charmm
      brben(16) = 110.0d0
      brbenk(16) = 0.0d0

! - bond angle for O-C2-C1 (in DPPC) - taken from Charmm
      brben(17) = 109.5d0
      brbenk(17) = 0.0d0

! - bond angle for O-P-O' (in DPPC) - taken from Charmm
      brben(18) = 108.2d0
      brbenk(18) = 0.0d0

! - bond angle for O-P-O (in DPPC) - taken from Charmm
      brben(19) = 102.6d0
      brbenk(19) = 0.0d0

! - bond angle for O'-P-O' (in DPPC) - taken from Charmm
      brben(20) = 119.9d0
      brbenk(20) = 0.0d0

! - bond angle for O-C2-C2 (in DPPC) - taken from Charmm
      brben(21) = 109.5d0
      brbenk(21) = 0.0d0

! - bond angle for C2-C2-H (in DPPC) - taken from Charmm
      brben(22) = 110.0d0
      brbenk(22) = 0.0d0

! - bond angle for C2-C2-N3 (in DPPC) - taken from Charmm
      brben(23) = 109.5d0
      brbenk(23) = 0.0d0

! - bond angle for C2-N3-C3 (in DPPC) - taken from Charmm
!   same parameters as OPLS-UA ether C-C-O and C-O-C
      brben(24) = 112.0d0
      brbenk(24) = 0.0d0

! - bond angle for H-C2-N3 (in DPPC) - taken from Charmm
      brben(25) = 109.5d0
      brbenk(25) = 0.0d0

! - bond angle for N3-C3-H (in DPPC) - taken from Charmm
      brben(26) = 109.5d0
      brbenk(26) = 0.0d0

! - bond angle for C3-N3-C3 (in DPPC) - taken from Charmm
      brben(27) = 112.0d0
      brbenk(27) = 0.0d0

! - bond angle for H-C3-H (in DPPC) - taken from Charmm
      brben(28) = 106.4d0
      brbenk(28) = 0.0d0

! - bond angle for O-C'-O' (in DPPC) - taken from Charmm
      brben(29) = 123.0d0
      brbenk(29) = 0.0d0

! - bond angle for C1-O-C' (in DPPC) - taken from Charmm
      brben(30) = 109.5d0
      brbenk(30) = 0.0d0

!   bond angles from the All-atom AMBER force field
!   Cornell et al JACS 117 19 5179-5197 1995
!   and from OPLS-AA JACS 118 45 11225-11236 (1996)

! - bond angle for C-C=O H3C-COOH  AMBER 
      brben(31) = 120.4d0
      brbenk(31) = 4.03d4

! - bond angle for C--O--H (carboxyllic acid) AMBER
      brben(32) = 113.0d0
      brbenk(32) = 1.76d4

! - bond angle for H3C--C--OH (carboxyllic acid) OPLS-AA 1995
      brben(33) = 108.0d0
      brbenk(33) = 3.53d4
 
! - bond angle for O==C--OH (carboxyllic acid) OPLS-AA 1995
      brben(34) = 121.0d0
      brbenk(34) = 4.03d4

! - TraPPE bond angle for C-O-H in alkanol - from AMBER (OPLS) flexible
      brben(35) = 108.5d0
      brbenk(35) = 27720.0d0

! - TraPPE bond angle for C-C-O in alkanol - from AMBER
      brben(36) = 109.47d0
      brbenk(36) = 25200.0d0

! - WATER - SPC-FQ JCP 101 (7) 1 1994
! --   SPC/E J. Phys. Chem. 91 6269-6271 (1987)
      brben(37) = 109.47d0
      brbenk(37) = 17640.0d0

! - OPLS AA for c-c-c in alkanes
      brben(38) = 112.7d0
      brbenk(38) = 29382.3d0

! - MMFF for c-c-c in alkanes
!      brben(38) = 109.608d0
!      brbenk(38) = 30818.6d0

! - OPLS AA for c-c-h in alkanes
      brben(39) = 110.7d0
!      brbenk(39) = 0.0d0
      brbenk(39) = 18883.2d0

! - MMFF for c-c-h in alkanes
!      brben(39) = 110.549d0
!      brbenk(39) = 0.0d0
!     brbenk(39) = 23032.5d0

! - OPLS AA for h-c-h in alkanes
      brben(40) = 107.8d0
!    brbenk(40) = 0.0d0
! - Methane h-c-h
!      brben(40) = 109.4712206344907d0
      brbenk(40) = 16617.2d0
! - MMFF for h-c-h in alkanes
!      brben(40) = 108.836d0
!      brbenk(40) = 0.0d0
!      brbenk(40) = 18686.7d0


! - AA F-C-F & F-C-C in Methyl group
      brben(41) = 109.4712206344907d0
      brbenk(41) = 0.0d0

! - AA C-C-C in C3F8 (JPC 95 1991 P3136)
      brben(42) = 115.9d0
      brbenk(42) = 31250.0d0

! - AA F-C-F in Methylene group (JPC 95 1991 P3136)
      brben(43) = 107.0
      brbenk(43) = 0.0d0

! - TIP4P HOH angle
      brben(44) = 104.52d0
      brbenk(44) = 0.0d0

! - ??? OPLS C-C-O angle
      brben(45) = 108.0d0
      brbenk(45) = 25200.0d0

! - bond angle for C-O-H in alkanol 
! - Monica's Alcohol Fixed bond angle
      brben(46) = 108.5d0
      brbenk(46) = 0.0d0
! *** methanol C-O-H angle from H2 param set Monica's dissertation ***
!      brben(46) = 108.63d0
!      brbenk(46) = 0.0d0

! - bond angle for O-C-C in alkanol 
! - Monica's Alcohol Fixed bond angle
      brben(47) = 108.0d0
      brbenk(47) = 0.0d0

! - bond angle for C-C-C in alkanol 
! - Monica's Alcohol Fixed bond angle
      brben(48) = 112.0d0
      brbenk(48) = 0.0d0

! - bond angle for O-C-O in CO2 and R-C=-N (nitriles)
      brben(49) = 180.0d0
      brbenk(49) = 0.0d0

! - bond angle for (O=)C-O-C for OPLS
      brben(50) = 115.0d0
      brbenk(50) = 31250d0

! - bond angle for O=C-O for ester
      brben(51) = 125d0
      brbenk(51) = 31250d0
 
! - bond angle for C-C=O in carboxylic acids
      brben(52) = 126.0d0
!      brbenk(52) = 31250.0d0
      brbenk(52) = brbenk(31)
      
! - bond angle for C-C-O in carboxylic acids
      brben(53) = 111.0d0
!      brbenk(53) = 31250.0d0
      brbenk(53) = brbenk(33)

! - bond angle for O=C-O in carboxylic acids
      brben(54) = 123.0d0
!      brbenk(54) = 31250.0d0
      brbenk(54) = brbenk(34)

! - bond angle for C-O-H in carboxylic acids
      brben(55) = 107.0d0
!      brbenk(55) = 0.0d0
      brbenk(55) = brbenk(32)

! --- TraPPE Ether C-O-C angle from OPLS
!   - constant from AMBER ether (off of website)
      brben(56) = 112.0d0
! wrong value!     brbenk(56) = 15102.0d0
      brbenk(56) = 30200.0d0

! --- TraPPE Ether C-C-O angle from OPLS
!   - constant from AMBER '94 website
      brben(57) = 112.0d0
! wrong value!     brbenk(57) = 12581.2d0
      brbenk(57) = 25150.0d0

! -- C-N-H in amines (TraPPE-7)
      brben(58) = 112.9d0
      brbenk(58) = 31250.0d0

! -- H-N-H in amines (TraPPE-7)
      brben(59) = 106.4d0
      brbenk(59) = 21955.0d0
 
! -- C-N-C in amines (TraPPE-7)
      brben(60) = 109.5d0
      brbenk(60) = 25178.0d0

! -- Ch3-C=O bond angle -> in Ketones
      brben(61) = 121.4d0
      brbenk(61) = 31250.0d0
!	brbenk(61) = 0.0d0

! -- Ch3-C-CH3 in ketones
      brben(62) = 117.2d0
      brbenk(62) = 31250.0d0 
!        brbenk(62) = 0.0d0

! -- C-S-H in thiols
      brben(63) = 96.0d0
      brbenk(63) = 31250.0d0 

! -- C-S-C in sulfides
      brben(64) = 99.0d0
      brbenk(64) = 31250.0d0

! -- C-O-C in hydrofuran
!        brben(65) = 111.0d0
      brben(65) = 110.0d0
      brbenk(65) = 0.0d0

! * formic acid model from llnl 4/6/04 jms
! --- O-C=O bend in 9-6 formic acid model
      brben(66) = 125.6d0
      brbenk(66) = 76992.0d0

! --- O-C-H bend in 9-6 formic acid model
      brben(67) = 112.3d0
      brbenk(67) = 23148.0d0

! --- O=C-H bend in 9-6 formic acid model
      brben(68) = 122.1d0
      brbenk(68) = 23148.0d0

! --- C-O-H bend in 9-6 formic acid model
      brben(69) = 103.0d0
      brbenk(69) = 26670.0d0
      
! --  C-N-O bond angle for nitro (UA)
      brben(70) = 111.5d0
      brbenk(70) = 40284.0d0

! -- O-N-O for nitro (UA + EH/TATB)
      brben(71) = 125.0d0
      brbenk(71) = 40284.0d0

! -- CH(aro)-CH(aro)-[CH(aro) or CHx]
      brben(72) = 120.0d0
      brbenk(72) = 0.0d0

! --  C-N-O bond angle for nitro (EH/TATB)
      brben(73) = 117.7d0
      brbenk(73) = 40284.0d0

! -- Parameters for DMMP (bending constants are not correct)
! --  O=P-O     
      brben(78) = 114.2d0
      brbenk(78) = 31250.0d0
! -- O=P-CH3
      brben(79) = 119.25d0
      brbenk(79) = 31250.0d0
! -- O-P-O
      brben(80) = 106.5d0
      brbenk(80) = 31250.0d0
! -- CH3-P-O
      brben(81)= 100.5d0
      brbenk(81) = 31250.0d0
! -- O-P-O in DMMP
      brben(82) = 159.0d0
      brbenk(82) = 31250.0d0	
!-  C(aro)-C(aro)-C(aro)
      brben(83) = 120.0d0
      brbenk(83) = 0.0d0

! -- BEGIN parameters for Neimark DMMP
! -- O=P-CH3
      brben(84) = 116.3d0
      brbenk(84) = 40293.0d0 

! -- O=P-O
      brben(85) = 116.5d0
      brbenk(85) = 50397.0d0

! -- CH3-P-O
      brben(86) = 104.3d0
      brbenk(86) = 20447.0d0

! -- P-O-CH3
      brben(87) = 121.0
      brbenk(87) = 40293.0d0

! -- END parameters for Neimark DMMP


! -- Starting All atom floro-alkanes (Neeraj)

! -- F-C-F bend angle for floroalkanes (Bending const. OPLS)

      brben(90)  = 109.47d0
      brbenk(90) = 38731.0d0

! -- F-C-C (OPLS)

      brben(91)  = 109.5d0
      brbenk(91) = 25150.0d0

! -- H-C-F (OPLS)
 
      brben(92)  = 107.0d0
      brbenk(92) = 20120.0d0

! -- C-C-C
      brben(93) = 115.90d0
      brbenk(93) = 31250.0d0

! -- C-C-H
      brben(94)  = 108.30d0
      brbenk(94) = 20120.0d0

! -- Cl-C-Cl Cf2Cl2
      brben(95)  = 111.70d0
      brbenk(95) = 39234.0d0

! -- Cl-C-F CF2Cl2
      brben(96) = 107.80
      brbenk(96) = 37725.0d0 

! - parameters for acrylates

! - bond angle for CHx-O-C
      brben(100) = 115.0d0
      brbenk(100) = 31250.0d0
   
! - bond angle for O-C=O
!c      brben(101) = 125.0d0
!c      brbenk(101) = 62500.0d0
      brben(101) = 123.0d0
      brbenk(101) = 20150.0d0

! - bond angle for O-CHx-CHy
      brben(102) = 111.0d0
      brbenk(102) = 17650.0d0

! - bond angle for O=C-CHx (carboxylic acids)
      brben(103) = 126.0d0
      brbenk(103) = 20150.0d0

! - bond angle for C-CH=CH2, C-C=CH2, and CHx-C(sp2)-CHx
      brben(104) = 119.7d0
      brbenk(104) = 35210.0d0


! -- TraPPE-7 Bending parameters Collin's part added by Neeraj

! -- C-N-0 in nitro
      brben(105) = 111.5d0
      brbenk(105) = 40284.0d0

! -- C-C-N in nitro
      brben(106) = 111.1d0
      brbenk(106) = 31724.0d0

! -- H-C-N in nitro
      brben(107) = 105.0d0
      brbenk(107) = 17624.0d0

! -- O-N-O in nitro
      brben(108) = 125.0d0
      brbenk(108) = 40284.0d0

! --- H - C - O in alcohols and ethers H-C-N and C-N-H for amines
      brben(109) = 109.5d0
      brbenk(109) = 17624.39d0

! --- H-N-H for amine
      brben(110) = 106.4d0
      brbenk(110) = 21955.0d0

! --- C - N - C for amine and O - C - H for alkanol
      brben(111) = 109.5d0
      brbenk(111) = 25178.0d0

! --- amine C-C-N
      brben(112) = 109.47d0
      brbenk(112) = 28300.0d0

! --- OPLS C - C - Cl
      brben(113) =  111.7d0
      brbenk(113) = 39276.9d0

! --- OPLS H - C - Cl
!      brben(114) = 114.20d0
!      brbenk(114) = 35248.5d0

! -- OPLS H-C-O in alcohols
      brben(115) = 109.5d0
      brbenk(115) = 17605d0 

! --Begin TATB
! -- C--C--N TATB
      brben(150) = 120.0d0
      brbenk(150) = 0.5d0*100.0d0*503.25

! -- C--C(NO2)--C, N--C--H, C--N--O TATB
      brben(151) = 122.0d0
      brbenk(151) = 0.5d0*100.0d0*503.25

! -- C--(NH2)--C, O--N--O
      brben(152) = 118.0d0
      brbenk(152) = 0.5d0*100.0d0*503.25

! -- H--N--H TATB
      brben(153) = 123.0d0
      brbenk(153) = 0.5d0*100.0d0*503.25
! -- END TATTB

! --- JLR 11-11-09
! --- Silica Si-O-Si
      brben(154) = 147.0d0
      brbenk(154) = 20.0d0*1000.0d0/1.9872
! --- END JLR 11-11-09

      do i = 1, nvmax
         brben(i) = brben(i) * 8.0d0 * datan(1.0d0) / 360.0d0
      end do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! - TraPPE-UA alkane torsional parameters for linear segment -
!     - OPLS torsion Jorgensen, Madura, Swenson JACS 1984 106, 813
      vtt0(1) = 0.0d0
      vtt1(1) = 355.03d0
      vtt2(1) = -68.19d0
      vtt3(1) = 791.32d0

! - TraPPE-UA alkane torsional parameters for segment containing a ternary -
! - carbon CH midsegment split OPLS torsion 
! - Siepmann et al Mol Phys 1997, 90, 687-693
!      vtt0(2) =  0.0d0
      vtt0(2) = -251.06d0
      vtt1(2) =  428.73d0
      vtt2(2) = -111.85d0
      vtt3(2) =  441.27d0

! - TraPPE-UA alkane torsional parameters for segment containing a -
! - quaternary carbon C as mid-segment -
! - Mundy et al, Faraday Disc 104, 17-36 (1996)
      vtt0(3) =    0.00d0
      vtt1(3) =    0.00d0
      vtt2(3) =    0.00d0
      vtt3(3) =  461.29d0

! - torsional parameters for X--(CH2 sp3)--(CH sp2)--Y
! - taken from jorgensen JACS 106 22 6638-6646 (1984)
      vtt0(4) =    686.1d0
      vtt1(4) =     86.4d0
      vtt2(4) =   -109.8d0
      vtt3(4) =   -282.2d0

! - TraPPE torsional parameters for alkanol (H-)-O-C-(-C) from OPLS
      vtt0(5) =    0.00d0
      vtt1(5) =  209.82d0
      vtt2(5) =  -29.17d0
      vtt3(5) =  187.93d0

! - TraPPE torsional parameters for alkanol (HO)-C-C-(-C) from OPLS
      vtt0(6) =    0.00d0
      vtt1(6) =  176.62d0
      vtt2(6) =  -53.34d0
      vtt3(6) =  769.93d0

! - TraPPE torsional parameters for alkanol (H-)-O-CH-(CH3)_2 from OPLS
      vtt0(7) =  215.89d0
      vtt1(7) =  197.33d0
      vtt2(7) =   31.46d0
      vtt3(7) =  -173.92d0

! - torsional parameters for TRANS alkene X-(CH sp2)-(CH sp2)-Y
! - used in harmonic potential E = vtt0 * (theta - vtt1)**2
! - Fitted from pcmodel trans 2-butene 4-14-99 MGM
      vtt0(11) =  13400.0d0
      vtt1(11) =  0.0d0

! - torsional parameters for CIS alkene X-(CH sp2)-(CH sp2)-Y
! - used in harmonic potential E = vtt0 * (theta - vtt1)**2
! - Fitted from pcmodel CIS 2-butene 4-14-99 MGM
      vtt0(12) =  12400.0d0
      vtt1(12) =  onepi

! *** torsional parameters for diethyl ether C-C-O-C from OPLS
! *** V = v1(1+cos()) - v2(1-cos(2*)) + v3(1+cos(3*))
! *** this is torsional type 25.

! *** type 27 is O-C-C-O

! - torsional parameters for acrylates

! - CHx-O-C=O
      vtt0(34) = 1820.74d0
      vtt1(34) = -417.41d0
      vtt2(34) = -1373.14d0
      vtt3(34) = -30.19d0
      vtt4(34) = 0.0d0

! - CH2=CH-C-O
      vtt0(35) = 823.03d0
      vtt1(35) = 47.91d0
      vtt2(35) = -773.13d0
      vtt3(35) = 1.99d0
      vtt4(35) = 0.0d0

! - CHx-O-C-CHy(sp2)
      vtt0(36) = 1820.74d0
      vtt1(36) = 417.41d0
      vtt2(36) = -1373.14d0
      vtt3(36) = 30.19d0
      vtt4(36) = 0.0d0

! - CH2=C-CH=O
      vtt0(37) = 823.03d0
      vtt1(37) = -47.91d0
      vtt2(37) = -773.13d0
      vtt3(37) = -1.99 d0
      vtt4(37) = 0.0d0

! - O=C-C(sp2)-CHx(sp3)
      vtt0(38) = 195.185d0
      vtt1(38) = -149.30d0
      vtt2(38) = 164.38d0
      vtt3(38) = 24.12d0
      vtt4(38) = 28.45d0

! - O-C(carbonyl)-C(sp2)-CHx(sp3)
      vtt0(39) =195.185d0
      vtt1(39) = 149.30d0
      vtt2(39) = 164.38d0
      vtt3(39) = -24.12d0
      vtt4(39) = 28.45d0

! - CHx-CHy-O-C(carbonyl)
      vtt0(40) = 2029.99d0
      vtt1(40) = -751.83d0
      vtt2(40) = -538.95d0
      vtt3(40) = -22.10d0
      vtt4(40) = -51.27d0

! - CH2=CH-CH=CH2 potential for butadiene
      vtt0(41) = 2034.577d0
      vtt1(41) = 531.571d0
      vtt2(41) = -1239.35d0
      vtt3(41) = 460.038d0
      vtt4(41) = 196.382d0

! - O=CH-CH-CH3 torsion for 2-methyl propanal
      vtt0(42) = 1063.29d0
      vtt1(42) = -736.9d0
      vtt2(42) = 57.84d0
      vtt3(42) = -293.23d0
      vtt4(42) = 0.0d0

! - O-CH2-CH2-CHx torsion for ether/acrylate w/long side chain (new functional form)
      vtt0(43) = 893.21d0
      vtt1(43) = 176.62d0
      vtt2(43) = 53.34d0
      vtt3(43) = 769.93d0
      vtt4(43) = 0.0d0

! - H-O-CH2-CH2 torsion from TraPPE 5 alcohols (new functional form)
      vtt0(44) = 368.58d0
      vtt1(44) = 209.8d0
      vtt2(44) = 29.17d0
      vtt3(44) = 187.93d0
      vtt4(44) = 0.0d0

! - O-CH2-CH2-O from TraPPE 6 glycols for HEA (new functional form)
      vtt0(45) = 1258.09d0
      vtt1(45) = 0.0d0
      vtt2(45) = 251.62d0
      vtt3(45) = 1006.47d0
      vtt4(45) = 0.0d0

! - CH2=CH-C(sp2)-CH3 potential for isoprene
      vtt0(46) = 1861.286d0
      vtt1(46) = -349.966d0
      vtt2(46) = -1048.70d0
      vtt3(46) = -580.535d0
      vtt4(46) = 117.915d0


! -- Torsional parameters for Neimark DMMP.  Six parameter
! -- torsional function
! -- CH3-P-O-CH3
        vtt0(48) = 33.80d0
        vtt1(48) = 317.0d0
        vtt2(48) = 38.0d0 
        vtt3(48) = -29.35d0
        vtt4(48) = 37.0d0
        vtt5(48) = -3.0d0
! -- O=P-O-CH3
        vtt0(49) = 0.0d0
        vtt1(49) = 0.0d0
        vtt2(49) = 50.5d0 
        vtt3(49) = 0.0d0
        vtt4(49) = 0.0d0
        vtt5(49) = 0.0d0
! -- O-P-O-CH3
        vtt0(50) = 0.0d0
        vtt1(50) = 480.0d0
        vtt2(50) = 252.6d0 
        vtt3(50) = 0.0d0
        vtt4(50) = 0.0d0
        vtt5(50) = 0.0d0

! *** Starting  fluorocarbons: Fit by Neeraj 6/24/2006 MP2/6-311+G**.
! *** All atom forcefield

! -- F-C-C-F
        vtt0(51) = 2543.43d0
        vtt1(51) = 1.25603d0
        vtt2(51) = -8.58713d0
        vtt3(51) = -1261.59

! -- F-C-C-C
        vtt0(52)  = 1985.58d0
        vtt1(52)  = -0.18585d0
        vtt2(52)  = 4.07924d0
        vtt3(52)  = -992.516d0

! -- CF-CF-CF-CF Fitted for perfluoropentane
        vtt0(53)   = 1124.71
        vtt1(53)   = 849.824
        vtt2(53)   = 331.375
        vtt3(53)   = 908.94
        vtt4(53)   = 434.207
        vtt5(53)   = 201.725
        vtt6(53)   = 127.14
        vtt7(53)   = -163.547
        vtt8(53)   = 46.9091
        vtt9(53)   = 60.7336

! -- H-C-C-F
        vtt0(54)  = 819.016
        vtt1(54)  = -6.14744
        vtt2(54)  = -43.9376
        vtt3(54)  = 895.398
        vtt4(54)  = 42.3686
        vtt5(54)  = -16.1043
        vtt6(54)  = 87.2734


! -- CF-CF-CF-CF for flourobutane
        vtt0(55)   = 1231.36d0 
        vtt1(55)   = 972.012d0
        vtt2(55)   = 365.577
        vtt3(55)   = 981.158
        vtt4(55)   = 364.393
        vtt5(55)   = 226.897
        vtt6(55)   = 121.712
        vtt7(55)   = -123.205
        vtt8(55)   = 28.2703
        vtt9(55)   = 44.6317
       
! -- Available for other Fluororcarbon. They have been shifted down 101 102 etc...
! -- Ethane H-C-C-H vtorso=a0+a1*(1-cosx)+a2*(1+cos2x)+a3*(1-cos3x)
        vtt0(56) = 1521.29d0
        vtt1(56) = -0.135221d0
        vtt2(56) = -0.545298d0
        vtt3(56) = -765.161d0

! -- Ethanol H-O-C-C vtorso=a0+a1*(1-cosx)+a2*(1+cos2x)+a3*(1-cos3x)+a4*(1+cos4x)
        vtt0(57) = 639.492d0
        vtt1(57) = -101.095d0
        vtt2(57) = 10.2389d0
        vtt3(57) = -321.075d0
        vtt4(57) = 89.8948d0

! -- Ethanol H-O-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x+a7*cos7x
        vtt0(58) = 262.743 
        vtt1(58) = -72.2022
        vtt2(58) = 25.3956
        vtt3(58) = 261.653
        vtt4(58) = -38.3658 
        vtt5(58) =  42.3685
        vtt6(58) = 7.93367
        vtt7(58) = 15.1805

!-- Ethanol O-C-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x
        vtt0(59) = 853.463 
        vtt1(59) = 11.4499
        vtt2(59) = -12.8932
        vtt3(59) = 887.455
        vtt4(59) = 12.9193
        vtt5(59) = -10.5521
        vtt6(59) = 35.1449

! - type 44 is under type 39


! -Hydrofluoroethers F-C-O-C vtorso=a0+a1*cosx+a2*cos(2*x)=a3*cos(3*x)
       vtt0(60) = 804.608 
       vtt1(60) = -6.3210
       vtt2(60) = 9.1809
       vtt3(60) = 785.878
 
! -Hydrofluoroethers H-C-O-C vtorso=a0+a1*cosx+a2*cos(2*x)+a3*cos(3*x)
       vtt0(61) = 327.282d0
       vtt1(61) = 5.29603d0
       vtt2(61) = 9.29972d0
       vtt3(61) = 324.084d0
   
! - Hydrofluoroethers F-C-C-O vtorso = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x) + a5*c
! os(5*x)
        vtt0(62) = 1738.42d0 
        vtt1(62) = -462.352d0
        vtt2(62) =  9.39616d0
        vtt3(62) =  1086.9d0
        vtt4(62) = 238.459d0
        vtt5(62) =  40.9771d0
      
! - Hydrofluorethers C-O-C-C vtorso = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x) + a4*cos(4*x) + a5
! *cos(5*x) + a6*cos(6*x) + a7*cos(7*x
        vtt0(63) = 1207.59d0
        vtt1(63) = 1146.14d0
        vtt2(63) = 90.5438d0
        vtt3(63) = 252.856d0
        vtt4(63) = 306.492d0
        vtt5(63) = 101.542d0
        vtt6(63) = 14.9379d0
        vtt7(63) = -104.586d0
      
! - Hydrofluorethers H-C-C-O vtorso = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x) + a4*cos(4*x) + a5
! *cos(5*x) + a6*cos(6*x) + a7*cos(7*x)
        vtt0(64) = 1434.36d0
        vtt1(64) = -56.6214d0
        vtt2(64) = 76.6241d0
        vtt3(64) = 767.995d0
        vtt4(64) = 360.543d0
        vtt5(64) = 307.889d0
        vtt6(64) = 555.559d0
        vtt7(64) = 354.076d0

! *** Ketones Jeff's Email

! --- Ch3-C(=O)-CH2-CH2
        vtt0(65) = -17.26d0
        vtt1(65) = 752.6d0
        vtt2(65) = 14.89d0
        vtt3(65) = 282.1d0
 
! --- O=CH2-CH2-CHx: Fit by me (Jeff) to HF/3-21G ab initio data

        vtt0(66) = 2035.5876d0
        vtt1(66) = -736.8992d0
        vtt2(66) = 57.8440d0
        vtt3(66) = -293.229d0


! -- From TraPPE-7 Added by Neeraj for amines and nitro
! *** 60 through 70 fit by the BEST TORSIONAL FITTING PROGRAM EVER!

! - torsional parameters for nitro H-C-C-N fit to OPLS
      vtt0(70) = 165.24d0
      vtt1(70) = -219.263d0
      vtt2(70) = 63.667d0
      vtt3(70) = 4.98368d0
      vtt4(70) = 8.18974d0
      vtt5(70) = -2.63063d0
      vtt6(70) = 0.78009d0

! - torsional parameters for nitro C-C-N-O fit to OPLS
      vtt0(71) = 69.1666d0
      vtt1(71) = -41.3563d0
      vtt2(71) = -14.5474d0
      vtt3(71) = -19.1091d0
      vtt4(71) = 8.02837d0
      vtt5(71) = -2.91134d0
      vtt6(71) = 0.954035d0

! - torsional parameters for nitro H-C-N-O fit to OPLS
      vtt0(72) = 75.4217d0
      vtt1(72) = -40.797d0
      vtt2(72) = 80.445d0
      vtt3(72) = -41.0586d0
      vtt4(72) = 16.0442d0
      vtt5(72) = -5.44169d0
      vtt6(72) = 1.67867d0

! - torsional parameters for ether H-C-O-H fit to OPLS
      vtt0(73) = 192.557d0
      vtt1(73) = -88.3325d0
      vtt2(73) = 10.2361d0
      vtt3(73) = -114.617d0
      vtt4(73) = 0.177971d0
      vtt5(73) = -0.0247285d0
      vtt6(73) = 0.00349949d0

! - torsional parameters for alkanol and ether H-C-C-O OPLS
      vtt0(74) = 215.758d0
      vtt1(74) = 94.6829d0
      vtt2(74) = 40.9651d0
      vtt3(74) = -144.295d0
      vtt4(74) = 10.7712d0
      vtt5(74) = -3.6513d0
      vtt6(74) = 1.11172d0

! - torsional parameters for alkanol H-C-O-C fit to OPLS
      vtt0(75) = 351.912d0
      vtt1(75) = -289.934d0
      vtt2(75) = 195.209d0
      vtt3(75) = -284.436d0
      vtt4(75) = 37.1459d0
      vtt5(75) = -13.173d0
      vtt6(75) = 4.28713d0

! - torisonal parameters for amine H-C-N-H
      vtt0(76) = 198.768d0
      vtt1(76) = -109.123d0
      vtt2(76) = 12.4603d0
      vtt3(76) = -102.29d0
      vtt4(76) = 0.210352d0
      vtt5(76) = -0.0287978d0
      vtt6(76) = 0.00401613

! - torisonal parameters for amine H-C-N-C
      vtt0(77) = 173.871d0
      vtt1(77) = -36.9908d0
      vtt2(77) = 4.69016d0
      vtt3(77) = -141.655d0
      vtt4(7) = 0.0976006d0
      vtt5(77) = -0.0148358d0
      vtt6(77) = 0.0022968d0

! - torisonal parameters for amine H-N-C-C
      vtt0(78) = 189.877d0
      vtt1(78) =47.8376d0
      vtt2(78) =104.991d0
      vtt3(78) =-105.243d0
      vtt4(78) = 0.0d0
      vtt5(78) = 0.0d0
      vtt6(78) = 0.0d0

! - torisonal parameters for amine C-N-C-C
      vtt0(79) = 1466.12d0
      vtt1(79) = -2188.07d0
      vtt2(79) = 1380.77d0
      vtt3(79) = -889.694d0
      vtt4(79) = 329.24d0
      vtt5(79) = -136.897d0
      vtt6(79) = 52.6532d0

!      vtt0(79) = 864.411d0
!      vtt1(79) = -1029.11d0
!      vtt2(79) = 718.434d0
!      vtt3(79) = -43.7331d0
!      vtt4(79) = 8.24626d0
!      vtt5(79) = -1.59901d0
!      vtt6(79) = 0.315767d0

! - torisonal parameters for amine H-C-C-N same as C-C-C-N
      vtt0(80) = 438.11d0
      vtt1(80) = 480.681d0
      vtt2(80) = 150.364d0
      vtt3(80) = -115.192d0
      vtt4(80) = -0.566972d0
      vtt5(80) = 0.0847927d0
      vtt6(80) = -0.0129149d0

! -- Starting methyl, dimethyl, diethyl acetamide torsions
! -- Mp2/6-311+g**//HF/6-311+g**

! --------------Beging All atom alkane potentials---------------

! -- All atom alkane. This torsion is fit to the C5H12.

! -- torsion H-C3-C2-C2 (Linear Alkanes) vtorso=a0 + a3*cos(3x)       
        vtt0(100) = 750.517d0+22.0d0
        vtt1(100) = 0.0d0
        vtt2(100) = 0.0d0
        vtt3(100) = 772.345d0
        vtt4(100) = 0.0d0
        vtt5(100) = 0.0d0 
        vtt6(100) = 0.0d0
        vtt7(100) = 0.0d0
        vtt8(100) = 0.0d0
        vtt9(100) = 0.0d0

! -- Ethane H-C3-C3-H vtorso=a0+a1*(cosx)+a2*(cos2x)+a3*(cos3x)
! --- Using the same for H-C3-C2-H
        vtt0(101) = 755.453+10.0
        vtt1(101) = 0.135223
        vtt2(101) = -0.545309
        vtt3(101) =  765.161d0
        vtt4(101) = 0.0d0
        vtt5(101) = 0.0d0
        vtt6(101) = 0.0d0
        vtt7(101) = 0.0d0
        vtt8(101) = 0.0d0
        vtt9(101) = 0.0d0

! -- torsion H-C3-C0-C3 (NeoPentane) V = a0+a1*cosx+a2*cos(2x)+a3*cos(3*x)
        vtt0(102) = 929.219d0 + 7.8d0
        vtt1(102) = -0.00692145d0
        vtt2(102) = 0.314788d0
        vtt3(102) = 961.663d0
        vtt4(102) = 0.0d0
        vtt5(102) = 0.0d0
        vtt6(102) = 0.0d0
        vtt7(102) = 0.0d0
        vtt8(102) = 0.0d0
        vtt9(102) = 0.0d0

! -- torsion for C2-C2-C2-C3 MP2/6-311+G** for C5H12
        vtt0(103) = 1124.71
        vtt1(103)  = 849.824
        vtt2(103)  = 331.375
        vtt3(103)  = 908.94
        vtt4(103)  = 434.207
        vtt5(103)  = 201.725
        vtt6(103)  = 127.14
        vtt7(103)  = -163.547
        vtt8(103)  = 46.9091
        vtt9(103)  = 60.7336

! Type 104   C*-C1-C1-C*


! Type 105   C*-C1-C0-C*


! Type 106   C*-C0-C0-C*



! --HF/6-311+g** C3-C0-C2-C1 fitted for 224trimethylhexane
        vtt0(107) = 1107.68+100.0
        vtt1(107) = 36.2517
        vtt2(107) = 121.795
        vtt3(107) = 1149.47
        vtt4(107) = -97.2671
        vtt5(107) = 168.586
        vtt6(107) = 0.0d0
        vtt7(107) = 0.0d0
        vtt8(107) = 0.0d0
        vtt9(107) = 0.0d0
! --HF/6-311+g** C3-C2-C1-C2 fitted for 224trimethylhexane
! y = a0 + a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)
        vtt0(108) = 1837.67
        vtt1(108) = 190.996
        vtt2(108) = 51.5605
        vtt3(108) = 1539.46
        vtt4(108) = -133.516
        vtt5(108) = 0.0d0
        vtt6(108) = 0.0d0
        vtt7(108) = 0.0d0
        vtt8(108) = 0.0d0
        vtt9(108) = 0.0d0

! --MP2/6-311+G** C3-C1-C2-C2 2Methylhexane
! y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)
        vtt0(109) = 788.503
        vtt1(109) = 410.738
        vtt2(109) = 283.868
        vtt3(109) = 874.922
        vtt4(109) = 155.877
        vtt5(109) = 31.9489
        vtt6(109) = 92.541
        vtt7(109) = 0.0d0
        vtt8(109) = 0.0d0
        vtt9(109) = 0.0d0

!--MP2/6-311+G** C2-C2-C1-C2 4Methyl Hexane
! y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)
        vtt0(110) = 1326.64+5.0
        vtt1(110) = 129.728
        vtt2(110) = 154.814
        vtt3(110) = 1239.51
        vtt4(110) = -114.036
        vtt5(110) = 0.0d0
        vtt6(110) = 0.0d0
        vtt7(110) = 0.0d0
        vtt8(110) = 0.0d0
        vtt9(110) = 0.0d0
 
! -- MP2/6-311+G** H-C2-C2-H
! -- y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)
        vtt0(111) = 952.756+17
        vtt1(111) = -337.392
        vtt2(111) = -508.6
        vtt3(111) = 1031.73
        vtt4(111) = 37.6621
        vtt5(111) = -116.642
        vtt6(111) = 105.495
        vtt7(111) = 0.0d0
        vtt8(111) = 0.0d0
        vtt9(111) = 0.0d0
 
! MP2/6-311+G** H-C2-C2-C2
! y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)
        vtt0(112) = 759.436+3
        vtt1(112) = -76.9545
        vtt2(112) = 2.98074
        vtt3(112) = 721.005
        vtt4(112) = 74.762
        vtt5(112) = -20.6381
        vtt6(112) = 52.1718
        vtt7(112) = 0.0d0
        vtt8(112) = 0.0d0
        vtt9(112) = 0.0d0

! Type 113   H-C2-C1-H
        vtt0(113) = 1315.38
        vtt1(113) = 118.146
        vtt2(113) = -210.975
        vtt3(113) = 1217.28
        vtt4(113) = 61.743
        vtt5(113) = -26.8647
        vtt6(113) = 180.745
        vtt7(113) = 0.0d0
        vtt8(113) = 0.0d0
        vtt9(113) = 0.0d0  


! Type 114   H-C3-C1-H
        vtt0(114) = 799.89+55
        vtt1(114) = 3.09152
        vtt2(114) = 2.50577
        vtt3(114) = 853.972
        vtt4(114) = 0.0d0
        vtt5(114) = 0.0d0
        vtt6(114) = 0.0d0
        vtt7(114) = 0.0d0
        vtt8(114) = 0.0d0
        vtt9(114) = 0.0d0

! Type 115   H-C1-C1-H 


! Type 116   H-C3-C1-C*
        vtt0(116) = 773.277
        vtt1(116) = 1.9856
        vtt2(116) = 138.805
        vtt3(116) = 775.89
        vtt4(116) = -97.855
        vtt5(116) = 2.75463
        vtt6(116) = 0.0d0
        vtt7(116) = 0.0d0
        vtt8(116) = 0.0d0
        vtt9(116) = 0.0d0
 

! Type 117   H-C2-C1-C*
        vtt0(117) = 1079.02+30
        vtt1(117) = -449.079
        vtt2(117) = -450.254
        vtt3(117) = 1043.84
        vtt4(117) = 101.841
        vtt5(117) = -137.144
        vtt6(117) = 119.671
        vtt7(117) = 0.0d0
        vtt8(117) = 0.0d0
        vtt9(117) = 0.0d0


! Type 118   H-C2-C0-C*
        vtt0(118) = 1112.86
        vtt1(118) = -14.327
        vtt2(118) = 123.98
        vtt3(118) = 1275.8
        vtt4(118) = -88.994
        vtt5(118) = 51.9279
        vtt6(118) = 168.154
        vtt7(118) = 0.0d0
        vtt8(118) = 0.0d0
        vtt9(118) = 0.0d0


! Type 119   H-C1-C2-C*  (Hc1c2c3)
        vtt0(119) = 1115.19
        vtt1(119) = -513.756
        vtt2(119) = -231.342
        vtt3(119) = 1093.66
        vtt4(119) = -98.7663
        vtt5(119) = -50.8918
        vtt6(119) = 117.262
        vtt7(119) = 0.0d0
        vtt8(119) = 0.0d0
        vtt9(119) = 0.0d0 

! Type 120   H-C1-C1-C*


! Type 121   H-C1-C0-C*







! ----------End All atom alkane potentials-----------------------------

! ---------Begin Amide potential------------------------------------

! -- Adding for Formamide
! -- MP2/6-311+G** H-C(=O)-N-H
! y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)
        vtt0(131) = 1914.53
        vtt1(131) = -361.234
        vtt2(131) = -2574.16
        vtt3(131) = 540.275
        vtt4(131) = 743.347
        vtt5(131) = -228.483
        vtt6(131) = 0.0d0 
        vtt7(131) = 0.0d0
        vtt8(131) = 0.0d0
        vtt9(131) = 0.0d0
! MP2/6-311+G** O-C(=O)-N-H
! y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7*cos(7*x)+a8*cos(8*x)
!omputed values:
        vtt0(132) = 2014.47+75
        vtt1(132) = 184.106
        vtt2(132) = -2728.3
        vtt3(132) = -301.602
        vtt4(132) = 811.395
        vtt5(132) = 162.795
        vtt6(132) = -157.987
        vtt7(132) = -90.937
        vtt8(132) = 125.165
        vtt9(132) = 0.0d0

! 1) H-C(=O)-N-C

! Fitting with formula: y =
! a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7*
! cos(7*x)+a8*cos(8*x)
        vtt0(133) = 2667.38+50
        vtt1(133) = -566.191
        vtt2(133) = -3339.08
        vtt3(133) = 713.659
        vtt4(133) = 766.582
        vtt5(133) = -205.182
        vtt6(133) = -188.234
        vtt7(133) = 114.007
        vtt8(133) = 114.365
        vtt9(133) = 0.0d0 

!  H-C-N-C(=O)

! y=a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7*
! cos(7*x)+a8*cos(8*x)+a9*cos(9*x)
        vtt0(134) = 253.627
        vtt1(134) = -5.62679
        vtt2(134) = -2.32438
        vtt3(134) = -369.616
        vtt4(134) = 4.03741
        vtt5(134) = 8.73542
        vtt6(134) = 154.936
        vtt7(134) = -0.441499
        vtt8(134) = -8.09161
        vtt9(134) = -62.3274

! H-C-N-C
        vtt0(135) = 532.285+40
        vtt1(135) = 11.956
        vtt2(135) = -7.03127
        vtt3(135) = 550.882
        vtt4(135) = 0.0d0
        vtt5(135) = 0.0
        vtt6(135) = 0.0
        vtt7(135) = 0.0d0
        vtt8(135) = 0.0d0
        vtt9(135) = 0.0d0
 
! C-N-C(=O)-O

! y= a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7* cos(7*x)+a8*cos(8*x)+a9*cos(9*x)
        vtt0(136) = 2835.23+20
        vtt1(136) = 411.734
        vtt2(136) = -3616.28
        vtt3(136) = -552.388
        vtt4(136) = 904.816
        vtt5(136) = 191.827
        vtt6(136) = -284.345
        vtt7(136) = -99.6365
        vtt8(136) = 232.172
        vtt9(136) = 68.6528
! ---------------------------------------------------------------------
! Starting Alkanol
! -- Ethanol H-O-C-C vtorso=a0+a1*(1-cosx)+a2*(1+cos2x)+a3*(1-cos3x)+a4*(1+cos4x)
        vtt0(144) = 639.492d0
        vtt1(144) = -101.095d0
        vtt2(144) = 10.2389d0
        vtt3(144) = -321.075d0
        vtt4(144) = 89.8948d0

! -- Ethanol H-O-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x+a7*cos7x
        vtt0(145) = 262.743
        vtt1(145) = -72.2022
        vtt2(145) = 25.3956
        vtt3(145) = 261.653
        vtt4(145) = -38.3658
        vtt5(145) =  42.3685
        vtt6(145) = 7.93367
        vtt7(145) = 15.1805

!-- Ethanol O-C-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x
        vtt0(146) = 853.463
        vtt1(146) = 11.4499
        vtt2(146) = -12.8932
        vtt3(146) = 887.455
        vtt4(146) = 12.9193
        vtt5(146) = -10.5521
        vtt6(146) = 35.1449

! Starting TATB part

!    H--N--C--C
       vtt0(200) = 0.5d0*17.0d0*503.25d0

!    O--N--C--C
       vtt0(201) = 0.5d0*5.6d0*503.25d0

!    C--C--C--C, C--C--C--N, N--C--C--N
       vtt0(202) = 0.5d0*25.0d0*503.25d0

      return
      end





