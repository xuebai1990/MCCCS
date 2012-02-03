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

c -- C-N in amines
      brvib(41) = 1.474d0
      brvibk(41) = 0.0d0

c -- N-H in amines
      brvib(42) = 1.041d0
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

c * unused 
c
c * methanol O-H from H2 parameter set from Monica's dissertation *
c      brvib() = 1.0285d0
c      brvibk() = 0.0d0
c
c * methanol C-O from H2 parameter set from Monica's dissertation *
c      brvib() = 1.4175d0
c      brvibk() = 0.0d0

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c - TraPPE-UA bond angle for alkane segment centered at methylene (CH2 sp3)-
      brben(1) = 114.0d0
c     write(6,*) '***** brben', brben(1) * raddeg
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
      brbenk(39) = 0.0d0
c      brbenk(39) = 18883.2d0

c - MMFF for c-c-h in alkanes
c      brben(39) = 110.549d0
c      brbenk(39) = 0.0d0
c     brbenk(39) = 23032.5d0

c - OPLS AA for h-c-h in alkanes
      brben(40) = 107.8d0
      brbenk(40) = 0.0d0
c - Methane h-c-h
c      brben(40) = 109.4712206344907d0
c      brbenk(40) = 16617.2d0
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

c -- C-N-H in amines
      brben(58) = 112.9d0
      brbenk(58) = 31250.0d0

c -- H-N-H in amines
      brben(59) = 105.9d0
      brbenk(59) = 31250.0d0
 
c -- C-N-C in tertiary amines
      brben(60) = 108.0d0
      brbenk(60) = 31250.0d0

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

      do i = 1, 69
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

      return
      end





