      program topmon

! topmon
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

      integer(KIND=normal_int)::i,seed
      real(KIND=double_precision)::random,rtest
      dimension rtest(10)
      character::*50 fileout
! --- RP added for MPI for random2.f
      integer(KIND=normal_int)::nseed
      real(KIND=double_precision)::random2,nrtest
      dimension nrtest(10)
      character::*50 nfileout

! ----------------------------------------------------------------

      read(4,*)
      read(4,*) seed
      read(4,*)
      read(4,*) nseed
      fileout = 'Nrandomtest.dat'

      open(unit=71,FILE=fileout,status="unknown") 

! --- RP added for random2()
      nfileout = 'nNrandomtest.dat'
      open(unit=72,FILE=nfileout,status="unknown")
! -------------------
!      seed = 0

! --- initialize random number generator 
      call ranset(seed)

! ---RP added for new random2.f
      call nranset(nseed)
! ---------

! *** set up random number generator ***
!      call g05ccf
!      call g05cbf(54581)
!      idum = 5481
!      xini = ran1(idum)
 
! -------------------------------------------------------------------

! --- print 10 random numbers for control ---
      do i=1,10
         rtest(i) = random()
! ---RP added for random2()
         nrtest(i)= random2()
! -----------
      end do
      write(71,1000) (rtest(i),i=1,5)
      write(71,1000) (rtest(i),i=6,10)
 1000 format(2x,5f10.6)

! -- RP added for random2()
      write(72,2000) (nrtest(i),i=1,5)
      write(72,2000) (nrtest(i),i=6,10)
 2000 format(2x,5f10.6)
      close(72)
! ----------

      close(71)
! --- call main program
      call monola

! ----------------------------------------------------------------

      call cleanup('')
      end
