      subroutine suvibe

c suvibe
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
      include 'connect.inc'
      include 'conver.inc'
      include 'contorsion.inc'

      integer i

c ----------------------------------------------------------------

c - TraPPE-UA- fixed bond length for C-C bonds -
      brvib(1) = 1.54d0
      brvibk(1) = 0.0d0

c - bond stretching for C-C bonds -
c      brvib(2) = 1.54d0
c      brvibk(2) = 226450.3d0
c - bond stretching for C-C bonds - CHARMM
      brvib(2) = 1.54d0
      brvibk(2) = 237018.0d0

c - fixed bond length for C - C double bonds (not aromatic)
c - from table 1.6 of McMurry Organic Chemistry 3rd edition
      brvib(3) = 1.33d0
      brvibk(3) = 0.0d0

c - fixed bond length for sp2 C - C sp3 single bonds (not aromatic)
c - from Jorgensen JACS 106 6638-6646 (1984)
      brvib(4) = 1.50d0
      brvibk(4) = 0.0d0

c - fixed bond length for H-C bond Williams 
      brvib(5) = 1.04d0
      brvibk(5) = 0.0d0

c - fixed bond length for C-C bond Williams, also Freire UA C-C, OPLS-UA 
      brvib(6) = 1.53d0
      brvibk(6) = 0.0d0

c - fixed bond length for C-CH3 in neopentane Freire Mol Phys 1997
      brvib(7) = 1.55d0
      brvibk(7) = 0.0d0

c - fixed bond length for C-CH2 in neohexane and CHCH 23 dimethyl butane
c -  Freire Mol Phys 1997
      brvib(8) = 1.56d0
      brvibk(8) = 0.0d0

c - bond stretching for C-COOH bonds - CHARMM
      brvib(10) = 1.52d0
      brvibk(10) = 0.0d0
c      brvibk(10) = 239627.0d0

c - bond stretching for C=O bonds in COOH - CHARMM
      brvib(11) = 1.23d0
      brvibk(11) = 0.0d0
c      brvibk(11) = 619288.2d0

c - bond stretching for C(CO)-OH bonds - CHARMM
      brvib(12) = 1.37d0
      brvibk(12) = 402578.3d0

c     AMBER94 Cornell et. at.  JACS 117, 5179-5197 (1995)

c - fixed bond length for C==O bonds in COOH - AMBER94
      brvib(13) = 1.229d0
      brvibk(13) = 0.0d0

c - fixed bond length for C--O bonds in COOH - AMBER 94
      brvib(14) = 1.364d0
      brvibk(14) = 0.0d0

c - fixed bond length for sp2C-sp3C bonds - AMBER 94
      brvib(15) = 1.522d0
      brvibk(15) = 0.0d0

c - fixed bond length for O-H bonds in COOH and COH - AMBER 94
      brvib(16) = 0.96d0
      brvibk(16) = 0.0d0

c - fixed bond length for C-O bonds in alkanol 
      brvib(17) = 1.41d0
      brvibk(17) = 0.0d0

c - fixed bond length for H-O bonds in water 
      brvib(18) = 0.9572d0
      brvibk(18) = 0.0d0

c - OPLS AA bond length for c-c in alkanes
c      brvib(19) = 1.529d0
c      brvibk(19) = 0.0d0
      brvib(19) = 1.535d0
      brvibk(19) = 0.0d0

c - MMFF bond length for c-c in alkanes
c      brvib(19) = 1.508d0
c      brvibk(19) = 0.0d0

c - OPLS AA bond length for c-h in alkanes
c      brvib(20) = 1.09d0
c      brvibk(20) = 0.0d0
      brvib(20) = 0.55d0
      brvibk(20) = 0.0d0

c - MMFF bond length for c-h in alkanes
c      brvib(20) = 1.093
c      brvibk(20) = 0.0d0

c - OPLS AA bond length for h-h from CRC 72nd Ed.
c      brvib(21) = 0.74611d0
      brvib(21) = 1.535d0/2.0d0
      brvibk(21) = 0.0d0

c --- AA for CF4 bond length C-F (Surface Science 367(1996) P177) ---
      brvib(22) = 1.37d0
      brvibk(22) = 0.0d0

c --- AA for CF4 bond length C-F(Nose and Klein J.Chem.Phys. 78(1983) 6928) ---
      brvib(23) = 1.323d0
      brvibk(23) = 0.0d0

c --- Amber AA for CF4 bond length C-F (JCC 13(1992) P963) ---
      brvib(24) = 1.38d0
      brvibk(24) = 0.0d0

c --- SPC-FQ JCP 101, (7) 1 1994 6141
      brvib(25) = 1.0d0
      brvibk(25) = 0.0d0

c --- TIP4P OH bond length
      brvib(26) = 0.9572d0
      brvibk(26) = 0.0d0

c --- TIP4P OM length
      brvib(27) = 0.15d0
      brvibk(27) = 0.0d0

c --- Fixed bond length for O-O in dioxygen
c     J Chem Phys 98 (12) 9895--9904 Muller-Plathe et al
      brvib(28) = 1.21d0
      brvibk(28) = 0.0d0

c - TraPPE fixed bond length for O-H bonds in COH - from OPLS
      brvib(29) = 0.945d0
      brvibk(29) = 0.0d0

c - TraPPE fixed bond length for C-O bonds in alkanol - from OPLS
      brvib(30) = 1.43d0
      brvibk(30) = 0.0d0

c --- Fixed bond length for N-N in dinitrogen
c     53rd ed CRC Handbook page F-180, also used for C-H distance
      brvib(31) = 1.10d0
      brvibk(31) = 0.0d0

c --- TraPPE C-O bond length for CO2
       brvib(32) = 1.160d0
       brvibk(32) = 0.0d0

c --- SPECIAL LJ CHAIN J Phys Chem fit to give octane phase diagram
      brvib(33) = 3.2664d0
      brvibk(33) = 0.0d0

c --- N - charge site bond length for N2 w/quadrupole
      brvib(34) = 0.55d0
      brvibk(34) = 0.0d0

c -- CH3-C in ketones, also carboxylic acids
      brvib(35) = 1.52d0
      brvibk(35) = 0.0d0

c --- C=O bond length in carboxylic acid (OPLS)
      brvib(36) = 1.214d0
      brvibk(36) = 0.0d0
      
c -- C-O in carboxylic acids
      brvib(37) = 1.364d0
      brvibk(37) = 0.0d0

c -- O-H in carboxylic acids
      brvib(38) = 0.970d0
      brvibk(38) = 0.0d0

c -- CH3-S in thiols
      brvib(39) = 1.82d0
      brvibk(39) = 0.0d0

c -- S-H in thiols
      brvib(40) = 1.34d0
      brvibk(40) = 0.0d0

c -- C-N in amines (TraPPE-7)
      brvib(41) = 1.448d0
      brvibk(41) = 0.0d0

c -- N-H in amines (TraPPE-7)
      brvib(42) = 1.01d0
      brvibk(42) = 0.0d0

c -- C-N tertiary amines
      brvib(43) = 1.47d0
      brvibk(43) = 0.0d0

c -- TraPPE-UA C-O in ethers from OPLS
      brvib(44) = 1.41d0
      brvibk(44) = 0.0d0

c --- H-F bond length in HF
      brvib(45) = 0.917d0
      brvibk(45) = 0.0d0

c --- M site bond length in HF
      brvib(46) = 0.166d0
      brvibk(46) = 0.0d0

c --- C(arom)-C(arom) fixed bond length
      brvib(47) = 1.4d0
      brvibk(47) = 0.0d0

c --- TraPPE-UA CH3-[C=-N] triple bond length (from OPLS)
      brvib(48) = 1.157d0
      brvibk(48) = 0.0d0

c --- OPLS [CH3-C]=-N bond length -not- for TraPPE-UA
      brvib(49) = 1.458d0
      brvibk(49) = 0.0d0

c --- CH3-CH2 bond length in OPLS-UA Ether model
      brvib(50) = 1.516d0
      brvibk(50) = 0.0d0

c --- dinitrogen bond length (Tildesley?) From Nose + Klein Mol Phys 1983, 50, 1055.
      brvib(51) = 1.098d0
      brvibk(51) = 0.0d0

c * formic acid model from llnl 4/6/04 jms
c --- C-O bond length in 9-6 formic acid model
      brvib(52) = 1.376d0
      brvibk(52) = 237518.0d0

c --- C=O bond length in 9-6 formic acid model
      brvib(53) = 1.188d0
      brvibk(53) = 558190.5d0

c --- C-H bond length in 9-6 formic acid model
      brvib(54) = 1.075d0
      brvibk(54) = 227599.8d0

c --- O-H bond length in 9-6 formic acid model
      brvib(55) = 1.038d0
      brvibk(55) = 317466.8d0

c -- bond lengths for Neimark DMMP JPCA v108, 1435 (2004)

c  -- P=O
      brvib(70) = 1.458d0
      brvibk(70) = 0.0d0

c -- P-CH3
      brvib(71) = 1.795d0
      brvibk(71) = 0.0d0

c -- P-O(CH3)
      brvib(72) = 1.586d0
      brvibk(72) = 0.0d0

c -- O-CH3
      brvib(73) = 1.418d0
      brvibk(73) = 0.0d0
c -- end parameters for Neimark DMMP


C -- Starting floro-alkanes (Neeraj)

c --  C-F All atom floroalkanes
      brvib(80)  = 1.33d0
      brvibk(80) = 0.0d0

c --  C-C All atom Floroalkanes

      brvib(81) = 1.54d0
      brvibk(81) = 0.0d0

c --  C-H all atom Fluoroalkanes
      brvib(82) = 1.08d0
      brvibk(82)= 0.0d0

c -- C-Cl bond for Chloroflouroalkanes
      brvib(83) = 1.754
      brvibk(83) = 0.0d0


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c - test parameters for coarse-grain model
      brvib(110) = 3.25d0
      brvibk(110) = 0.0d0

      brben(110) = 150.0d0
      brbenk(110) = 1.0d-12
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c * unused 
c
c * methanol O-H from H2 parameter set from Monica's dissertation *
c      brvib() = 1.0285d0
c      brvibk() = 0.0d0
c
c * methanol C-O from H2 parameter set from Monica's dissertation *
c      brvib() = 1.4175d0
c      brvibk() = 0.0d0

c -- C--C TATB
       brvib(120) = 1.442d0
       brvibk(120) = 0.5d0*1400.0d0*503.25
c -- C--NO2 TATB
       brvib(121) = 1.419d0
       brvibk(121) = 0.5d0*1400.0d0*503.25
c -- C--NH2 TATB
       brvib(122) = 1.314d0
       brvibk(122) = 0.5d0*1400.0d0*503.25
c -- N--O TATB
       brvib(123) = 1.243d0
       brvibk(123) = 0.5d0*1400.0d0*503.25
c -- N--H TATB
       brvib(124) = 1.0d0
       brvibk(124) = 0.5d0*700.0d0*503.25


c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c - TraPPE-UA bond angle for alkane segment centered at methylene (CH2 sp3)-
      brben(1) = 114.0d0
c     write(iou,*) '***** brben', brben(1) * raddeg
      brbenk(1) = 31250.0d0

c - TraPPE-UA bond angle for alkane segment centered at ternary (CH sp3)-
c - also used by Freire for most alkane bond angles
      brben(2) = 112.0d0
      brbenk(2) = 31250.0d0

c - TraPPE-UA bond angle for alkane segment centered at quaternary (C sp3)-
c - basically the same as used by Freire for neopentane and 
c                                        2,3-dimethylbutane
      brben(3) = 109.47d0
      brbenk(3) = 31250.0d0

c - bond angle for C=C-C -  segment centered at (CH sp2)
c - taken from Amber JACS (1995) 117, 5179-5197
      brben(4) = 119.70d0
      brbenk(4) = 35210.0d0

c - bond angle for cis C=C-C NO LONGER USED 4-13-99
      brben(5) = 127.4d0
      brbenk(5) = 31250.0d0

c - bond angle for H-C[methyl]-C[methylene] Ryckaert 
      brben(6) = 112.49d0
      brbenk(6) = 0.0d0

c - bond angle for C-C-C Van der Ploeg
      brben(7) = 112.0d0
      brbenk(7) = 31278.0d0

c - bond angle for H-C[methylene]-H Ryckaert
      brben(8) = 106.0d0
      brbenk(8) = 0.0d0 

c - bond angle for C-C=O (in carboxylic headgroup) - taken from Charmm
c      brben(10) = 120.0d0
c --- OPLS value
      brben(10) = 125.0d0
      brbenk(10) = 68438.0d0

c - bond angle for C-C-O (in carboxylic headgroup) - taken from Charmm
      brben(11) = 110.0d0
      brbenk(11) = 12359.0d0

c - bond angle for C2-P-O (in DPPC) - taken from Charmm
      brben(12) = 120.5d0
      brbenk(12) = 0.0d0

c - bond angle for C2-C1-H (in DPPC) - taken from Charmm
      brben(13) = 110.0d0
      brbenk(13) = 0.0d0

c - bond angle for H-C2-H (in DPPC) - taken from Charmm
      brben(14) = 106.4d0
      brbenk(14) = 0.0d0

c - bond angle for H-C2-O (in DPPC) - taken from Charmm
      brben(15) = 109.5d0
      brbenk(15) = 0.0d0

c - bond angle for H-C2-C1 (in DPPC) - taken from Charmm
      brben(16) = 110.0d0
      brbenk(16) = 0.0d0

c - bond angle for O-C2-C1 (in DPPC) - taken from Charmm
      brben(17) = 109.5d0
      brbenk(17) = 0.0d0

c - bond angle for O-P-O' (in DPPC) - taken from Charmm
      brben(18) = 108.2d0
      brbenk(18) = 0.0d0

c - bond angle for O-P-O (in DPPC) - taken from Charmm
      brben(19) = 102.6d0
      brbenk(19) = 0.0d0

c - bond angle for O'-P-O' (in DPPC) - taken from Charmm
      brben(20) = 119.9d0
      brbenk(20) = 0.0d0

c - bond angle for O-C2-C2 (in DPPC) - taken from Charmm
      brben(21) = 109.5d0
      brbenk(21) = 0.0d0

c - bond angle for C2-C2-H (in DPPC) - taken from Charmm
      brben(22) = 110.0d0
      brbenk(22) = 0.0d0

c - bond angle for C2-C2-N3 (in DPPC) - taken from Charmm
      brben(23) = 109.5d0
      brbenk(23) = 0.0d0

c - bond angle for C2-N3-C3 (in DPPC) - taken from Charmm
c   same parameters as OPLS-UA ether C-C-O and C-O-C
      brben(24) = 112.0d0
      brbenk(24) = 0.0d0

c - bond angle for H-C2-N3 (in DPPC) - taken from Charmm
      brben(25) = 109.5d0
      brbenk(25) = 0.0d0

c - bond angle for N3-C3-H (in DPPC) - taken from Charmm
      brben(26) = 109.5d0
      brbenk(26) = 0.0d0

c - bond angle for C3-N3-C3 (in DPPC) - taken from Charmm
      brben(27) = 112.0d0
      brbenk(27) = 0.0d0

c - bond angle for H-C3-H (in DPPC) - taken from Charmm
      brben(28) = 106.4d0
      brbenk(28) = 0.0d0

c - bond angle for O-C'-O' (in DPPC) - taken from Charmm
      brben(29) = 123.0d0
      brbenk(29) = 0.0d0

c - bond angle for C1-O-C' (in DPPC) - taken from Charmm
      brben(30) = 109.5d0
      brbenk(30) = 0.0d0

c   bond angles from the All-atom AMBER force field
c   Cornell et al JACS 117 19 5179-5197 1995
c   and from OPLS-AA JACS 118 45 11225-11236 (1996)

c - bond angle for C-C=O H3C-COOH  AMBER 
      brben(31) = 120.4d0
      brbenk(31) = 4.03d4

c - bond angle for C--O--H (carboxyllic acid) AMBER
      brben(32) = 113.0d0
      brbenk(32) = 1.76d4

c - bond angle for H3C--C--OH (carboxyllic acid) OPLS-AA 1995
      brben(33) = 108.0d0
      brbenk(33) = 3.53d4
 
c - bond angle for O==C--OH (carboxyllic acid) OPLS-AA 1995
      brben(34) = 121.0d0
      brbenk(34) = 4.03d4

c - TraPPE bond angle for C-O-H in alkanol - from AMBER (OPLS) flexible
      brben(35) = 108.5d0
      brbenk(35) = 27720.0d0

c - TraPPE bond angle for C-C-O in alkanol - from AMBER
      brben(36) = 109.47d0
      brbenk(36) = 25200.0d0

c - WATER - SPC-FQ JCP 101 (7) 1 1994
c --   SPC/E J. Phys. Chem. 91 6269-6271 (1987)
      brben(37) = 109.47d0
      brbenk(37) = 17640.0d0

c - OPLS AA for c-c-c in alkanes
      brben(38) = 112.7d0
      brbenk(38) = 29382.3d0

c - MMFF for c-c-c in alkanes
c      brben(38) = 109.608d0
c      brbenk(38) = 30818.6d0

c - OPLS AA for c-c-h in alkanes
      brben(39) = 110.7d0
c      brbenk(39) = 0.0d0
      brbenk(39) = 18883.2d0

c - MMFF for c-c-h in alkanes
c      brben(39) = 110.549d0
c      brbenk(39) = 0.0d0
c     brbenk(39) = 23032.5d0

c - OPLS AA for h-c-h in alkanes
      brben(40) = 107.8d0
c    brbenk(40) = 0.0d0
c - Methane h-c-h
c      brben(40) = 109.4712206344907d0
      brbenk(40) = 16617.2d0
c - MMFF for h-c-h in alkanes
c      brben(40) = 108.836d0
c      brbenk(40) = 0.0d0
c      brbenk(40) = 18686.7d0


c - AA F-C-F & F-C-C in Methyl group
      brben(41) = 109.4712206344907d0
      brbenk(41) = 0.0d0

c - AA C-C-C in C3F8 (JPC 95 1991 P3136)
      brben(42) = 115.9d0
      brbenk(42) = 31250.0d0

c - AA F-C-F in Methylene group (JPC 95 1991 P3136)
      brben(43) = 107.0
      brbenk(43) = 0.0d0

c - TIP4P HOH angle
      brben(44) = 104.52d0
      brbenk(44) = 0.0d0

c - ??? OPLS C-C-O angle
      brben(45) = 108.0d0
      brbenk(45) = 25200.0d0

c - bond angle for C-O-H in alkanol 
c - Monica's Alcohol Fixed bond angle
      brben(46) = 108.5d0
      brbenk(46) = 0.0d0
c *** methanol C-O-H angle from H2 param set Monica's dissertation ***
c      brben(46) = 108.63d0
c      brbenk(46) = 0.0d0

c - bond angle for O-C-C in alkanol 
c - Monica's Alcohol Fixed bond angle
      brben(47) = 108.0d0
      brbenk(47) = 0.0d0

c - bond angle for C-C-C in alkanol 
c - Monica's Alcohol Fixed bond angle
      brben(48) = 112.0d0
      brbenk(48) = 0.0d0

c - bond angle for O-C-O in CO2 and R-C=-N (nitriles)
      brben(49) = 180.0d0
      brbenk(49) = 0.0d0

c - bond angle for (O=)C-O-C for OPLS
      brben(50) = 115.0d0
      brbenk(50) = 31250d0

c - bond angle for O=C-O for ester
      brben(51) = 125d0
      brbenk(51) = 31250d0
 
c - bond angle for C-C=O in carboxylic acids
      brben(52) = 126.0d0
c      brbenk(52) = 31250.0d0
      brbenk(52) = brbenk(31)
      
c - bond angle for C-C-O in carboxylic acids
      brben(53) = 111.0d0
c      brbenk(53) = 31250.0d0
      brbenk(53) = brbenk(33)

c - bond angle for O=C-O in carboxylic acids
      brben(54) = 123.0d0
c      brbenk(54) = 31250.0d0
      brbenk(54) = brbenk(34)

c - bond angle for C-O-H in carboxylic acids
      brben(55) = 107.0d0
c      brbenk(55) = 0.0d0
      brbenk(55) = brbenk(32)

c --- TraPPE Ether C-O-C angle from OPLS
c   - constant from AMBER ether (off of website)
      brben(56) = 112.0d0
c wrong value!     brbenk(56) = 15102.0d0
      brbenk(56) = 30200.0d0

c --- TraPPE Ether C-C-O angle from OPLS
c   - constant from AMBER '94 website
      brben(57) = 112.0d0
c wrong value!     brbenk(57) = 12581.2d0
      brbenk(57) = 25150.0d0

c -- C-N-H in amines (TraPPE-7)
      brben(58) = 112.9d0
      brbenk(58) = 31250.0d0

c -- H-N-H in amines (TraPPE-7)
      brben(59) = 106.4d0
      brbenk(59) = 21955.0d0
 
c -- C-N-C in amines (TraPPE-7)
      brben(60) = 109.5d0
      brbenk(60) = 25178.0d0

c -- Ch3-C=O bond angle -> in Ketones
      brben(61) = 121.4d0
      brbenk(61) = 31250.0d0
c	brbenk(61) = 0.0d0

c -- Ch3-C-CH3 in ketones
      brben(62) = 117.2d0
      brbenk(62) = 31250.0d0 
c        brbenk(62) = 0.0d0

c -- C-S-H in thiols
      brben(63) = 96.0d0
      brbenk(63) = 31250.0d0 

c -- C-S-C in sulfides
      brben(64) = 99.0d0
      brbenk(64) = 31250.0d0

c -- C-O-C in hydrofuran
c        brben(65) = 111.0d0
      brben(65) = 110.0d0
      brbenk(65) = 0.0d0

c * formic acid model from llnl 4/6/04 jms
c --- O-C=O bend in 9-6 formic acid model
      brben(66) = 125.6d0
      brbenk(66) = 76992.0d0

c --- O-C-H bend in 9-6 formic acid model
      brben(67) = 112.3d0
      brbenk(67) = 23148.0d0

c --- O=C-H bend in 9-6 formic acid model
      brben(68) = 122.1d0
      brbenk(68) = 23148.0d0

c --- C-O-H bend in 9-6 formic acid model
      brben(69) = 103.0d0
      brbenk(69) = 26670.0d0
      
c -- Parameters for DMMP (bending constants are not correct)
c --  O=P-O     
      brben(78) = 114.2d0
      brbenk(78) = 31250.0d0
c -- O=P-CH3
      brben(79) = 119.25d0
      brbenk(79) = 31250.0d0
c -- O-P-O
      brben(80) = 106.5d0
      brbenk(80) = 31250.0d0
c -- CH3-P-O
      brben(81)= 100.5d0
      brbenk(81) = 31250.0d0
c -- O-P-O in DMMP
      brben(82) = 159.0d0
      brbenk(82) = 31250.0d0	
c-  C(aro)-C(aro)-C(aro)
      brben(83) = 120.0d0
      brbenk(83) = 0.0d0

c -- BEGIN parameters for Neimark DMMP
c -- O=P-CH3
      brben(84) = 116.3d0
      brbenk(84) = 40293.0d0 

c -- O=P-O
      brben(85) = 116.5d0
      brbenk(85) = 50397.0d0

c -- CH3-P-O
      brben(86) = 104.3d0
      brbenk(86) = 20447.0d0

c -- P-O-CH3
      brben(87) = 121.0
      brbenk(87) = 40293.0d0

c -- END parameters for Neimark DMMP


C -- Starting All atom floro-alkanes (Neeraj)

c -- F-C-F bend angle for floroalkanes (Bending const. OPLS)

      brben(90)  = 109.47d0
      brbenk(90) = 38731.0d0

c -- F-C-C (OPLS)

      brben(91)  = 109.5d0
      brbenk(91) = 25150.0d0

c -- H-C-F (OPLS)
 
      brben(92)  = 107.0d0
      brbenk(92) = 20120.0d0

c -- C-C-C
      brben(93) = 115.90d0
      brbenk(93) = 31250.0d0

c -- C-C-H
      brben(94)  = 108.30d0
      brbenk(94) = 20120.0d0

c -- Cl-C-Cl Cf2Cl2
      brben(95)  = 111.70d0
      brbenk(95) = 39234.0d0

c -- Cl-C-F CF2Cl2
      brben(96) = 107.80
      brbenk(96) = 37725.0d0 



c -- TraPPE-7 Bending parameters Collin's part added by Neeraj

c -- C-N-0 in nitro
      brben(101) = 111.5d0
      brbenk(101) = 40284.0d0

c -- C-C-N in nitro
      brben(102) = 111.1d0
      brbenk(102) = 31724.0d0

c -- H-C-N in nitro
      brben(103) = 105.0d0
      brbenk(103) = 17624.0d0

c -- O-N-O in nitro
      brben(104) = 125.0d0
      brbenk(104) = 40284.0d0

c --- H - C - O in alcohols and ethers H-C-N and C-N-H for amines
      brben(105) = 109.5d0
      brbenk(105) = 17624.39d0

c --- H-N-H for amine
      brben(106) = 106.4d0
      brbenk(106) = 21955.0d0

c --- C - N - C for amine and O - C - H for alkanol
      brben(107) = 109.5d0
      brbenk(107) = 25178.0d0

c --- amine C-C-N
      brben(108) = 109.47d0
      brbenk(108) = 28300.0d0

c --- OPLS C - C - Cl
      brben(109) =  111.7d0
      brbenk(109) = 39276.9d0

c --- OPLS H - C - Cl
c      brben(110) = 114.20d0
c      brbenk(110) = 35248.5d0

c -- OPLS H-C-O in alcohols
      brben(111) = 109.5d0
      brbenk(111) = 17605d0 

c --Begin TATB
c -- C--C--N TATB
      brben(150) = 120.0d0
      brbenk(150) = 0.5d0*100.0d0*503.25

c -- C--C(NO2)--C, N--C--H, C--N--O TATB
      brben(151) = 122.0d0
      brbenk(151) = 0.5d0*100.0d0*503.25

c -- C--(NH2)--C, O--N--O
      brben(152) = 118.0d0
      brbenk(152) = 0.5d0*100.0d0*503.25

c -- H--N--H TATB
      brben(153) = 123.0d0
      brbenk(153) = 0.5d0*100.0d0*503.25
c -- END TATTB


      do i = 1, nvmax
         brben(i) = brben(i) * 8.0d0 * datan(1.0d0) / 360.0d0
      enddo
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c - TraPPE-UA alkane torsional parameters for linear segment -
c     - OPLS torsion Jorgensen, Madura, Swenson JACS 1984 106, 813
      vtt0(1) = 0.0d0
      vtt1(1) = 355.03d0
      vtt2(1) = -68.19d0
      vtt3(1) = 791.32d0

c - TraPPE-UA alkane torsional parameters for segment containing a ternary -
c - carbon CH midsegment split OPLS torsion 
c - Siepmann et al Mol Phys 1997, 90, 687-693
c      vtt0(2) =  0.0d0
      vtt0(2) = -251.06d0
      vtt1(2) =  428.73d0
      vtt2(2) = -111.85d0
      vtt3(2) =  441.27d0

c - TraPPE-UA alkane torsional parameters for segment containing a -
c - quaternary carbon C as mid-segment -
c - Mundy et al, Faraday Disc 104, 17-36 (1996)
      vtt0(3) =    0.00d0
      vtt1(3) =    0.00d0
      vtt2(3) =    0.00d0
      vtt3(3) =  461.29d0

c - torsional parameters for X--(CH2 sp3)--(CH sp2)--Y
c - taken from jorgensen JACS 106 22 6638-6646 (1984)
      vtt0(4) =    686.1d0
      vtt1(4) =     86.4d0
      vtt2(4) =   -109.8d0
      vtt3(4) =   -282.2d0

c - TraPPE torsional parameters for alkanol (H-)-O-C-(-C) from OPLS
      vtt0(5) =    0.00d0
      vtt1(5) =  209.82d0
      vtt2(5) =  -29.17d0
      vtt3(5) =  187.93d0

c - TraPPE torsional parameters for alkanol (HO)-C-C-(-C) from OPLS
      vtt0(6) =    0.00d0
      vtt1(6) =  176.62d0
      vtt2(6) =  -53.34d0
      vtt3(6) =  769.93d0

c - TraPPE torsional parameters for alkanol (H-)-O-CH-(CH3)_2 from OPLS
      vtt0(7) =  215.89d0
      vtt1(7) =  197.33d0
      vtt2(7) =   31.46d0
      vtt3(7) =  -173.92d0

c - torsional parameters for TRANS alkene X-(CH sp2)-(CH sp2)-Y
c - used in harmonic potential E = vtt0 * (theta - vtt1)**2
c - Fitted from pcmodel trans 2-butene 4-14-99 MGM
      vtt0(11) =  13400.0d0
      vtt1(11) =  0.0d0

c - torsional parameters for CIS alkene X-(CH sp2)-(CH sp2)-Y
c - used in harmonic potential E = vtt0 * (theta - vtt1)**2
c - Fitted from pcmodel CIS 2-butene 4-14-99 MGM
      vtt0(12) =  12400.0d0
      vtt1(12) =  onepi

c *** torsional parameters for diethyl ether C-C-O-C from OPLS
c *** V = v1(1+cos()) - v2(1-cos(2*)) + v3(1+cos(3*))
c *** this is torsional type 25.

c *** type 27 is O-C-C-O

c -- Torsional parameters for Neimark DMMP.  Six parameter
c -- torsional function
c -- CH3-P-O-CH3
        vtt0(33) = 33.80d0
        vtt1(33) = 317.0d0
        vtt2(33) = 38.0d0 
        vtt3(33) = -29.35d0
        vtt4(33) = 37.0d0
        vtt5(33) = -3.0d0
c -- O=P-O-CH3
        vtt0(34) = 0.0d0
        vtt1(34) = 0.0d0
        vtt2(34) = 50.5d0 
        vtt3(34) = 0.0d0
        vtt4(34) = 0.0d0
        vtt5(34) = 0.0d0
c -- O-P-O-CH3
        vtt0(35) = 0.0d0
        vtt1(35) = 480.0d0
        vtt2(35) = 252.6d0 
        vtt3(35) = 0.0d0
        vtt4(35) = 0.0d0
        vtt5(35) = 0.0d0

c *** Starting  fluorocarbons: Fit by Neeraj 6/24/2006 MP2/6-311+G**.
c *** All atom forcefield

c -- F-C-C-F
        vtt0(36) = 2543.43d0
        vtt1(36) = 1.25603d0
        vtt2(36) = -8.58713d0
        vtt3(36) = -1261.59

c -- F-C-C-C
        vtt0(37)  = 1985.58d0
        vtt1(37)  = -0.18585d0
        vtt2(37)  = 4.07924d0
        vtt3(37)  = -992.516d0

c -- CF-CF-CF-CF Fitted for perfluoropentane
        vtt0(38)   = 1124.71
        vtt1(38)   = 849.824
        vtt2(38)   = 331.375
        vtt3(38)   = 908.94
        vtt4(38)   = 434.207
        vtt5(38)   = 201.725
        vtt6(38)   = 127.14
        vtt7(38)   = -163.547
        vtt8(38)   = 46.9091
        vtt9(38)   = 60.7336

c -- H-C-C-F
        vtt0(39)  = 819.016
        vtt1(39)  = -6.14744
        vtt2(39)  = -43.9376
        vtt3(39)  = 895.398
        vtt4(39)  = 42.3686
        vtt5(39)  = -16.1043
        vtt6(39)  = 87.2734


c -- CF-CF-CF-CF for flourobutane
        vtt0(44)   = 1231.36d0 
        vtt1(44)   = 972.012d0
        vtt2(44)   = 365.577
        vtt3(44)   = 981.158
        vtt4(44)   = 364.393
        vtt5(44)   = 226.897
        vtt6(44)   = 121.712
        vtt7(44)   = -123.205
        vtt8(44)   = 28.2703
        vtt9(44)   = 44.6317
       
C -- Available for other Fluororcarbon. They have been shifted down 101 102 etc...
c -- Ethane H-C-C-H vtorso=a0+a1*(1-cosx)+a2*(1+cos2x)+a3*(1-cos3x)
        vtt0(40) = 1521.29d0
        vtt1(40) = -0.135221d0
        vtt2(40) = -0.545298d0
        vtt3(40) = -765.161d0

c -- Ethanol H-O-C-C vtorso=a0+a1*(1-cosx)+a2*(1+cos2x)+a3*(1-cos3x)+a4*(1+cos4x)
        vtt0(41) = 639.492d0
        vtt1(41) = -101.095d0
        vtt2(41) = 10.2389d0
        vtt3(41) = -321.075d0
        vtt4(41) = 89.8948d0

c -- Ethanol H-O-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x+a7*cos7x
        vtt0(42) = 262.743 
        vtt1(42) = -72.2022
        vtt2(42) = 25.3956
        vtt3(42) = 261.653
        vtt4(42) = -38.3658 
        vtt5(42) =  42.3685
        vtt6(42) = 7.93367
        vtt7(42) = 15.1805

c-- Ethanol O-C-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x
        vtt0(43) = 853.463 
        vtt1(43) = 11.4499
        vtt2(43) = -12.8932
        vtt3(43) = 887.455
        vtt4(43) = 12.9193
        vtt5(43) = -10.5521
        vtt6(43) = 35.1449

c - type 44 is under type 39


c -Hydrofluoroethers F-C-O-C vtorso=a0+a1*cosx+a2*cos(2*x)=a3*cos(3*x)
       vtt0(45) = 804.608 
       vtt1(45) = -6.3210
       vtt2(45) = 9.1809
       vtt3(45) = 785.878

c -Hydrofluoroethers H-C-O-C vtorso=a0+a1*cosx+a2*cos(2*x)+a3*cos(3*x)
       vtt0(46) = 327.282d0
       vtt1(46) = 5.29603d0
       vtt2(46) = 9.29972d0
       vtt3(46) = 324.084d0

c - Hydrofluoroethers F-C-C-O vtorso = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x) + a5*c
c os(5*x)
        vtt0(47) = 1738.42d0 
        vtt1(47) = -462.352d0
        vtt2(47) =  9.39616d0
        vtt3(47) =  1086.9d0
        vtt4(47) = 238.459d0
        vtt5(47) =  40.9771d0

c - Hydrofluorethers C-O-C-C vtorso = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x) + a4*cos(4*x) + a5
c *cos(5*x) + a6*cos(6*x) + a7*cos(7*x
        vtt0(48) = 1207.59d0
        vtt1(48) = 1146.14d0
        vtt2(48) = 90.5438d0
        vtt3(48) = 252.856d0
        vtt4(48) = 306.492d0
        vtt5(48) = 101.542d0
        vtt6(48) = 14.9379d0
        vtt7(48) = -104.586d0

c - Hydrofluorethers H-C-C-O vtorso = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x) + a4*cos(4*x) + a5
c *cos(5*x) + a6*cos(6*x) + a7*cos(7*x)
        vtt0(49) = 1434.36d0
        vtt1(49) = -56.6214d0
        vtt2(49) = 76.6241d0
        vtt3(49) = 767.995d0
        vtt4(49) = 360.543d0
        vtt5(49) = 307.889d0
        vtt6(49) = 555.559d0
        vtt7(49) = 354.076d0

c *** Ketones Jeff's Email

c --- Ch3-C(=O)-CH2-CH2
        vtt0(51) = -17.26d0
        vtt1(51) = 752.6d0
        vtt2(51) = 14.89d0
        vtt3(51) = 282.1d0
 
c --- O=CH2-CH2-CHx: Fit by me (Jeff) to HF/3-21G ab initio data

        vtt0(52) = 2035.5876d0
        vtt1(52) = -736.8992d0
        vtt2(52) = 57.8440d0
        vtt3(52) = -293.229d0


c -- From TraPPE-7 Added by Neeraj for amines and nitro
c *** 60 through 70 fit by the BEST TORSIONAL FITTING PROGRAM EVER!

c - torsional parameters for nitro H-C-C-N fit to OPLS
      vtt0(60) = 165.24d0
      vtt1(60) = -219.263d0
      vtt2(60) = 63.667d0
      vtt3(60) = 4.98368d0
      vtt4(60) = 8.18974d0
      vtt5(60) = -2.63063d0
      vtt6(60) = 0.78009d0

c - torsional parameters for nitro C-C-N-O fit to OPLS
      vtt0(61) = 69.1666d0
      vtt1(61) = -41.3563d0
      vtt2(61) = -14.5474d0
      vtt3(61) = -19.1091d0
      vtt4(61) = 8.02837d0
      vtt5(61) = -2.91134d0
      vtt6(61) = 0.954035d0

c - torsional parameters for nitro H-C-N-O fit to OPLS
      vtt0(62) = 75.4217d0
      vtt1(62) = -40.797d0
      vtt2(62) = 80.445d0
      vtt3(62) = -41.0586d0
      vtt4(62) = 16.0442d0
      vtt5(62) = -5.44169d0
      vtt6(62) = 1.67867d0

c - torsional parameters for ether H-C-O-H fit to OPLS
      vtt0(63) = 192.557d0
      vtt1(63) = -88.3325d0
      vtt2(63) = 10.2361d0
      vtt3(63) = -114.617d0
      vtt4(63) = 0.177971d0
      vtt5(63) = -0.0247285d0
      vtt6(63) = 0.00349949d0

c - torsional parameters for alkanol and ether H-C-C-O OPLS
      vtt0(64) = 215.758d0
      vtt1(64) = 94.6829d0
      vtt2(64) = 40.9651d0
      vtt3(64) = -144.295d0
      vtt4(64) = 10.7712d0
      vtt5(64) = -3.6513d0
      vtt6(64) = 1.11172d0

c - torsional parameters for alkanol H-C-O-C fit to OPLS
      vtt0(65) = 351.912d0
      vtt1(65) = -289.934d0
      vtt2(65) = 195.209d0
      vtt3(65) = -284.436d0
      vtt4(65) = 37.1459d0
      vtt5(65) = -13.173d0
      vtt6(65) = 4.28713d0

c - torisonal parameters for amine H-C-N-H
      vtt0(66) = 198.768d0
      vtt1(66) = -109.123d0
      vtt2(66) = 12.4603d0
      vtt3(66) = -102.29d0
      vtt4(66) = 0.210352d0
      vtt5(66) = -0.0287978d0
      vtt6(66) = 0.00401613

c - torisonal parameters for amine H-C-N-C
      vtt0(67) = 173.871d0
      vtt1(67) = -36.9908d0
      vtt2(67) = 4.69016d0
      vtt3(67) = -141.655d0
      vtt4(67) = 0.0976006d0
      vtt5(67) = -0.0148358d0
      vtt6(67) = 0.0022968d0

c - torisonal parameters for amine H-N-C-C
      vtt0(68) = 189.877d0
      vtt1(68) =47.8376d0
      vtt2(68) =104.991d0
      vtt3(68) =-105.243d0
      vtt4(68) = 0.0d0
      vtt5(68) = 0.0d0
      vtt6(68) = 0.0d0

c - torisonal parameters for amine C-N-C-C
      vtt0(69) = 1466.12d0
      vtt1(69) = -2188.07d0
      vtt2(69) = 1380.77d0
      vtt3(69) = -889.694d0
      vtt4(69) = 329.24d0
      vtt5(69) = -136.897d0
      vtt6(69) = 52.6532d0

c      vtt0(69) = 864.411d0
c      vtt1(69) = -1029.11d0
c      vtt2(69) = 718.434d0
c      vtt3(69) = -43.7331d0
c      vtt4(69) = 8.24626d0
c      vtt5(69) = -1.59901d0
c      vtt6(69) = 0.315767d0

c - torisonal parameters for amine H-C-C-N same as C-C-C-N
      vtt0(70) = 438.11d0
      vtt1(70) = 480.681d0
      vtt2(70) = 150.364d0
      vtt3(70) = -115.192d0
      vtt4(70) = -0.566972d0
      vtt5(70) = 0.0847927d0
      vtt6(70) = -0.0129149d0

c -- Starting methyl, dimethyl, diethyl acetamide torsions
c -- Mp2/6-311+g**//HF/6-311+g**

c --------------Beging All atom alkane potentials---------------

c -- All atom alkane. This torsion is fit to the C5H12.

c -- torsion H-C3-C2-C2 (Linear Alkanes) vtorso=a0 + a3*cos(3x)       
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

c -- Ethane H-C3-C3-H vtorso=a0+a1*(cosx)+a2*(cos2x)+a3*(cos3x)
c --- Using the same for H-C3-C2-H
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

c -- torsion H-C3-C0-C3 (NeoPentane) V = a0+a1*cosx+a2*cos(2x)+a3*cos(3*x)
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

c -- torsion for C2-C2-C2-C3 MP2/6-311+G** for C5H12
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

c Type 104   C*-C1-C1-C*


c Type 105   C*-C1-C0-C*


c Type 106   C*-C0-C0-C*



c --HF/6-311+g** C3-C0-C2-C1 fitted for 224trimethylhexane
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
c --HF/6-311+g** C3-C2-C1-C2 fitted for 224trimethylhexane
c y = a0 + a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)
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

c --MP2/6-311+G** C3-C1-C2-C2 2Methylhexane
c y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)
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

c--MP2/6-311+G** C2-C2-C1-C2 4Methyl Hexane
c y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)
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
 
c -- MP2/6-311+G** H-C2-C2-H
c -- y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)
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
 
c MP2/6-311+G** H-C2-C2-C2
c y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)
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

c Type 113   H-C2-C1-H
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


c Type 114   H-C3-C1-H
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

c Type 115   H-C1-C1-H 


c Type 116   H-C3-C1-C*
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
 

c Type 117   H-C2-C1-C*
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


c Type 118   H-C2-C0-C*
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


c Type 119   H-C1-C2-C*  (Hc1c2c3)
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

c Type 120   H-C1-C1-C*


c Type 121   H-C1-C0-C*







c ----------End All atom alkane potentials-----------------------------

c ---------Begin Amide potential------------------------------------

c -- Adding for Formamide
c -- MP2/6-311+G** H-C(=O)-N-H
c y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)
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
c MP2/6-311+G** O-C(=O)-N-H
c y = a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7*cos(7*x)+a8*cos(8*x)
Computed values:
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

c 1) H-C(=O)-N-C

c Fitting with formula: y =
c a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7*
c cos(7*x)+a8*cos(8*x)
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

c  H-C-N-C(=O)

c y=a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7*
c cos(7*x)+a8*cos(8*x)+a9*cos(9*x)
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

c H-C-N-C
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
 
c C-N-C(=O)-O

c y= a0+a1*cos(x)+a2*cos(2*x)+a3*cos(3*x)+a4*cos(4*x)+a5*cos(5*x)+a6*cos(6*x)+a7* cos(7*x)+a8*cos(8*x)+a9*cos(9*x)
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
c ---------------------------------------------------------------------
C Starting Alkanol
c -- Ethanol H-O-C-C vtorso=a0+a1*(1-cosx)+a2*(1+cos2x)+a3*(1-cos3x)+a4*(1+cos4x)
        vtt0(144) = 639.492d0
        vtt1(144) = -101.095d0
        vtt2(144) = 10.2389d0
        vtt3(144) = -321.075d0
        vtt4(144) = 89.8948d0

c -- Ethanol H-O-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x+a7*cos7x
        vtt0(145) = 262.743
        vtt1(145) = -72.2022
        vtt2(145) = 25.3956
        vtt3(145) = 261.653
        vtt4(145) = -38.3658
        vtt5(145) =  42.3685
        vtt6(145) = 7.93367
        vtt7(145) = 15.1805

c-- Ethanol O-C-C-H vtorso=a0+a1*cosx+a2*cos2x+a3*cos3x+a4*cos4x+a5*cos5x+a6*cos6x
        vtt0(146) = 853.463
        vtt1(146) = 11.4499
        vtt2(146) = -12.8932
        vtt3(146) = 887.455
        vtt4(146) = 12.9193
        vtt5(146) = -10.5521
        vtt6(146) = 35.1449

c Starting TATB part

c    H--N--C--C
       vtt0(200) = 0.5d0*17.0d0*503.25d0

c    O--N--C--C
       vtt0(201) = 0.5d0*5.6d0*503.25d0

c    C--C--C--C, C--C--C--N, N--C--C--N
       vtt0(202) = 0.5d0*25.0d0*503.25d0

      return
      end





