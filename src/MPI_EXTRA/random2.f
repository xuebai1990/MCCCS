      FUNCTION RANDOM2()

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
      real(8)::RANDOM2 , nRCARRY
      RANDOM2 = nRCARRY()
      RETURN
! ----------------------------------------------------C
      END

      FUNCTION nRANDx(nISEED)
!----------------------------------------------------------------------C
!  Random number generator, fast and rough, machine independent.
!  Returns an uniformly distributed deviate in the 0 to 1 interval.
!  This random number generator is portable, machine-independent and
!  reproducible, for any machine with at least 32 bits / real::number.
!  REF: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes (1986)
!----------------------------------------------------------------------C
      IMPLICIT NONE
! RP added for debugging
!      include 'mpif.h'
!      include 'mpi.inc'

      integer::nIA , nIC , nISEED , nM1
      real(8)::nRANDx , nRM
      PARAMETER (nM1=714025,nIA=1366,nIC=150889,nRM=1.D+0/nM1)
!
      nISEED = MOD(nIA*nISEED+nIC,nM1)
! ---RP added for debug
!      write(6,*)'in random2 nISEED=',nISEED,',myid=',myid

      nRANDx= nISEED*nRM
      IF ( nRANDx.LT.0.D+0 ) call cleanup('*** Random number is negative ***')
!
      RETURN
      END

      SUBROUTINE nRANSET(nISEED)
      IMPLICIT NONE
      integer::nISEED

      CALL nRSTART(nISEED)
      RETURN
      END

      SUBROUTINE nRSTART(nISEEDA)
!----------------------------------------------------------------------C
!       Initialize Marsaglia list of 24 random numbers.
!----------------------------------------------------------------------C

      IMPLICIT NONE
! ----RP added for debugging
!      include 'mpif.h'
!      include 'mpi.inc'

      real(8)::nCarry , nran , nRANDx , nSeed
      integer::ni , nI24 , nIseed , nISEEDA , nJ24
      COMMON /nRANDM/ nSeed(24) , nCarry , nI24 , nJ24 , nIseed

      nI24 = 24
      nJ24 = 10
      nCarry = 0.D+0
      nIseed = nISEEDA
!
!       get rid of initial correlations in rand by throwing
!       away the first 100 random numbers generated.
!
      DO 10 ni = 1 , 100
        nran = nRANDx(nIseed)
   10 CONTINUE
!
!       initialize the 24 elements of seed
!

      DO 20 ni = 1 , 24
        nSeed(ni) = nRANDx(nIseed)
   20 CONTINUE

      RETURN
      END


      FUNCTION nRCARRY()
!----------------------------------------------------------------------C
!       Random number generator from Marsaglia.
!----------------------------------------------------------------------C
      IMPLICIT NONE
      real(8)::nCarry , nRCARRY , nSeed , nTWOM24,nTWOP24,nuni
      integer::nI24 , nIseed , nJ24
      PARAMETER (nTWOP24=16777216.D+0,nTWOM24=1.D+0/nTWOP24)
      COMMON /nRANDM/ nSeed(24) , nCarry , nI24 , nJ24 , nIseed
!
!       f.james Comp. Phys. Comm. 60, 329  (1990)
!       algorithm by G. Marsaglia and A. Zaman
!       base b = 2**24  lags r=24 and s=10
!
      nuni = nSeed(nI24) - nSeed(nJ24) - nCarry
      IF ( nuni.LT.0.D+0 ) THEN
        nuni = nuni + 1.D+0
        nCarry = nTWOM24
      ELSE
        nCarry = 0.D+0
      end if
      nSeed(nI24) = nuni
      nI24 = nI24 - 1
      IF ( nI24.EQ.0 ) nI24 = 24
      nJ24 = nJ24 - 1
      IF ( nJ24.EQ.0 ) nJ24 = 24
      nRCARRY = nuni

      RETURN
      END
