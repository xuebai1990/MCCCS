      FUNCTION RANDOM()

! random
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

      IMPLICIT NONE
      real(8)::RANDOM , RCARRY
      RANDOM = RCARRY()
      RETURN
! ----------------------------------------------------C
      END

      FUNCTION RANDx(ISEED)
!----------------------------------------------------------------------C
!  Random number generator, fast and rough, machine independent.
!  Returns an uniformly distributed deviate in the 0 to 1 interval.
!  This random number generator is portable, machine-independent and
!  reproducible, for any machine with at least 32 bits / real::number.
!  REF: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes (1986)
!----------------------------------------------------------------------C
      IMPLICIT NONE
      integer::IA , IC , ISEED , M1
      real(8)::RANDx , RM
      PARAMETER (M1=714025,IA=1366,IC=150889,RM=1.D+0/M1)
!
      ISEED = MOD(IA*ISEED+IC,M1)
      RANDx= ISEED*RM
      IF ( RANDx.LT.0.D+0 ) call cleanup('*** Random number
     &                              is negative ***')
!
      RETURN
      END

      SUBROUTINE RANSET(ISEED)
      IMPLICIT NONE
      integer::ISEED

      CALL RSTART(ISEED)
      RETURN
      END

      SUBROUTINE RSTART(ISEEDA)
!----------------------------------------------------------------------C
!       Initialize Marsaglia list of 24 random numbers.
!----------------------------------------------------------------------C
      IMPLICIT NONE
      real(8)::Carry , ran , RANDx , Seed
      integer::i , I24 , Iseed , ISEEDA , J24
      COMMON /RANDM/ Seed(24) , Carry , I24 , J24 , Iseed

      I24 = 24
      J24 = 10
      Carry = 0.D+0
      Iseed = ISEEDA
!
!       get rid of initial correlations in rand by throwing
!       away the first 100 random numbers generated.
!
      DO 10 i = 1 , 100
        ran = RANDx(Iseed)
   10 CONTINUE
!
!       initialize the 24 elements of seed
!

      DO 20 i = 1 , 24
        Seed(i) = RANDx(Iseed)
   20 CONTINUE

      RETURN
      END


      FUNCTION RCARRY()
!----------------------------------------------------------------------C
!       Random number generator from Marsaglia.
!----------------------------------------------------------------------C
      IMPLICIT NONE
      real(8)::Carry , RCARRY , Seed , TWOM24 , TWOP24 , uni
      integer::I24 , Iseed , J24
      PARAMETER (TWOP24=16777216.D+0,TWOM24=1.D+0/TWOP24)
      COMMON /RANDM/ Seed(24) , Carry , I24 , J24 , Iseed
!
!       f.james Comp. Phys. Comm. 60, 329  (1990)
!       algorithm by G. Marsaglia and A. Zaman
!       base b = 2**24  lags r=24 and s=10
!
      uni = Seed(I24) - Seed(J24) - Carry
      IF ( uni.LT.0.D+0 ) THEN
        uni = uni + 1.D+0
        Carry = TWOM24
      ELSE
        Carry = 0.D+0
      end if
      Seed(I24) = uni
      I24 = I24 - 1
      IF ( I24.EQ.0 ) I24 = 24
      J24 = J24 - 1
      IF ( J24.EQ.0 ) J24 = 24
      RCARRY = uni

      RETURN
      END
