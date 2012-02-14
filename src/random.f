!> \brief (Pseudo)random number generator (RNG)
!>
!> -- Mersenne-Twister (MT19937) with initialization improved 2002-01-26
!> Originally coded in C by Takuji Nishimura and Makoto Matsumoto
!> Translated to Fortran 77 by Tsuyoshi TADA. (2005-12-19)
!> Modified for use in topmon by Peng Bai. (2012-02-11)

!-----------------------------------------------------------------------
!    initialize large number (over 32-bit constant number)
!-----------------------------------------------------------------------
      subroutine mt_initln
      implicit none
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
!C    TOPBIT_MASK = Z'80000000'
!C    ALLBIT_MASK = Z'ffffffff'
!C    UPPER_MASK  = Z'80000000'
!C    LOWER_MASK  = Z'7fffffff'
!C    MATRIX_A    = Z'9908b0df'
!C    T1_MASK     = Z'9d2c5680'
!C    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end

!-----------------------------------------------------------------------
!> \brief Initialize mt(0:N-1) with a seed
!>
!> \param s seed, must be non-negative
!>
!> This subroutine should be called once before using the RNG, otherwise
!> a default seed will be used.
!-----------------------------------------------------------------------
      subroutine RANSET(s)
      implicit none
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
!     
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
         mt(mti)=1812433253*ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
         mt(mti)=iand(mt(mti),ALLBIT_MASK)
 100  continue
      initialized=DONE
!     
      return
      end

!-----------------------------------------------------------------------
!> \brief Generates a random number on [0,1)-real-interval
!-----------------------------------------------------------------------
      double precision function random()
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
!     
      if(initialized.ne.DONE)then
         call RANSET(21641)
      endif
!     First generates a random number on [0,0xffffffff]-interval
      if(mti.ge.N)then
         do 100 kk=0,N-M-1
            y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 100     continue
         do 200 kk=N-M,N-1-1
            y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 200     continue
         y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
         mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
         mti=0
      endif
!     
      y=mt(mti)
      mti=mti+1
!     
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
!     
      if(y.lt.0.d0) then
         random=dble(y)+2.d0**32
      else
         random=dble(y)
      end if
! Divide by 2^32-1, for [0,1]-real-interval
!      random=random/4294967295.d0
! Divide by 2^32, for [0,1)-real-interval
      random=random/4294967296.d0
      return
      end
