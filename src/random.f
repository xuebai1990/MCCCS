      FUNCTION RANDOM()

c random
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

      IMPLICIT NONE
      real(8)::RANDOM , RCARRY
      RANDOM = RCARRY()
      RETURN
C ----------------------------------------------------C
      END

      FUNCTION RANDx(ISEED)
C----------------------------------------------------------------------C
C  Random number generator, fast and rough, machine independent.
C  Returns an uniformly distributed deviate in the 0 to 1 interval.
C  This random number generator is portable, machine-independent and
C  reproducible, for any machine with at least 32 bits / real::number.
C  REF: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes (1986)
C----------------------------------------------------------------------C
      IMPLICIT NONE
      integer::IA , IC , ISEED , M1
      real(8)::RANDx , RM
      PARAMETER (M1=714025,IA=1366,IC=150889,RM=1.D+0/M1)
c
      ISEED = MOD(IA*ISEED+IC,M1)
      RANDx= ISEED*RM
      IF ( RANDx.LT.0.D+0 ) call cleanup('*** Random number is negative ***')
c
      RETURN
      END

      SUBROUTINE RANSET(ISEED)
      IMPLICIT NONE
      integer::ISEED

      CALL RSTART(ISEED)
      RETURN
      END

      SUBROUTINE RSTART(ISEEDA)
C----------------------------------------------------------------------C
C       Initialize Marsaglia list of 24 random numbers.
C----------------------------------------------------------------------C
      IMPLICIT NONE
      real(8)::Carry , ran , RANDx , Seed
      integer::i , I24 , Iseed , ISEEDA , J24
      COMMON /RANDM/ Seed(24) , Carry , I24 , J24 , Iseed

      I24 = 24
      J24 = 10
      Carry = 0.D+0
      Iseed = ISEEDA
c
c       get rid of initial correlations in rand by throwing
c       away the first 100 random numbers generated.
c
      DO 10 i = 1 , 100
        ran = RANDx(Iseed)
   10 CONTINUE
c
c       initialize the 24 elements of seed
c

      DO 20 i = 1 , 24
        Seed(i) = RANDx(Iseed)
   20 CONTINUE

      RETURN
      END


      FUNCTION RCARRY()
C----------------------------------------------------------------------C
C       Random number generator from Marsaglia.
C----------------------------------------------------------------------C
      IMPLICIT NONE
      real(8)::Carry , RCARRY , Seed , TWOM24 , TWOP24 , uni
      integer::I24 , Iseed , J24
      PARAMETER (TWOP24=16777216.D+0,TWOM24=1.D+0/TWOP24)
      COMMON /RANDM/ Seed(24) , Carry , I24 , J24 , Iseed
c
c       f.james Comp. Phys. Comm. 60, 329  (1990)
c       algorithm by G. Marsaglia and A. Zaman
c       base b = 2**24  lags r=24 and s=10
c
      uni = Seed(I24) - Seed(J24) - Carry
      IF ( uni.LT.0.D+0 ) THEN
        uni = uni + 1.D+0
        Carry = TWOM24
      ELSE
        Carry = 0.D+0
      ENDIF
      Seed(I24) = uni
      I24 = I24 - 1
      IF ( I24.EQ.0 ) I24 = 24
      J24 = J24 - 1
      IF ( J24.EQ.0 ) J24 = 24
      RCARRY = uni

      RETURN
      END
