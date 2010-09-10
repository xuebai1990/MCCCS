      FUNCTION RANDOM2()

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
      real(8)::RANDOM2 , nRCARRY
      RANDOM2 = nRCARRY()
      RETURN
C ----------------------------------------------------C
      END

      FUNCTION nRANDx(nISEED)
C----------------------------------------------------------------------C
C  Random number generator, fast and rough, machine independent.
C  Returns an uniformly distributed deviate in the 0 to 1 interval.
C  This random number generator is portable, machine-independent and
C  reproducible, for any machine with at least 32 bits / real number.
C  REF: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes (1986)
C----------------------------------------------------------------------C
      IMPLICIT NONE
c RP added for debugging
c      include 'mpif.h'
c      include 'mpi.inc'

      integer::nIA , nIC , nISEED , nM1
      real(8)::nRANDx , nRM
      PARAMETER (nM1=714025,nIA=1366,nIC=150889,nRM=1.D+0/nM1)
c
      nISEED = MOD(nIA*nISEED+nIC,nM1)
c ---RP added for debug
c      write(6,*)'in random2 nISEED=',nISEED,',myid=',myid

      nRANDx= nISEED*nRM
      IF ( nRANDx.LT.0.D+0 ) STOP '*** Random number is negative ***'
c
      RETURN
      END

      SUBROUTINE nRANSET(nISEED)
      IMPLICIT NONE
      integer::nISEED

      CALL nRSTART(nISEED)
      RETURN
      END

      SUBROUTINE nRSTART(nISEEDA)
C----------------------------------------------------------------------C
C       Initialize Marsaglia list of 24 random numbers.
C----------------------------------------------------------------------C

      IMPLICIT NONE
c ----RP added for debugging
c      include 'mpif.h'
c      include 'mpi.inc'

      real(8)::nCarry , nran , nRANDx , nSeed
      integer::ni , nI24 , nIseed , nISEEDA , nJ24
      COMMON /nRANDM/ nSeed(24) , nCarry , nI24 , nJ24 , nIseed

      nI24 = 24
      nJ24 = 10
      nCarry = 0.D+0
      nIseed = nISEEDA
c
c       get rid of initial correlations in rand by throwing
c       away the first 100 random numbers generated.
c
      DO 10 ni = 1 , 100
        nran = nRANDx(nIseed)
   10 CONTINUE
c
c       initialize the 24 elements of seed
c

      DO 20 ni = 1 , 24
        nSeed(ni) = nRANDx(nIseed)
   20 CONTINUE

      RETURN
      END


      FUNCTION nRCARRY()
C----------------------------------------------------------------------C
C       Random number generator from Marsaglia.
C----------------------------------------------------------------------C
      IMPLICIT NONE
      real(8)::nCarry , nRCARRY , nSeed , nTWOM24,nTWOP24,nuni
      integer::nI24 , nIseed , nJ24
      PARAMETER (nTWOP24=16777216.D+0,nTWOM24=1.D+0/nTWOP24)
      COMMON /nRANDM/ nSeed(24) , nCarry , nI24 , nJ24 , nIseed
c
c       f.james Comp. Phys. Comm. 60, 329  (1990)
c       algorithm by G. Marsaglia and A. Zaman
c       base b = 2**24  lags r=24 and s=10
c
      nuni = nSeed(nI24) - nSeed(nJ24) - nCarry
      IF ( nuni.LT.0.D+0 ) THEN
        nuni = nuni + 1.D+0
        nCarry = nTWOM24
      ELSE
        nCarry = 0.D+0
      ENDIF
      nSeed(nI24) = nuni
      nI24 = nI24 - 1
      IF ( nI24.EQ.0 ) nI24 = 24
      nJ24 = nJ24 - 1
      IF ( nJ24.EQ.0 ) nJ24 = 24
      nRCARRY = nuni

      RETURN
      END
