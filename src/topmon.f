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

      implicit none

      integer::i,seed
      real(8)::random,rtest
      dimension rtest(10)
      character*50 fileout


! ----------------------------------------------------------------

      read(4,*)
      read(4,*) seed

      fileout = 'Nrandomtest.dat'

      open(unit=71,FILE=fileout,status="unknown") 

!      seed = 0

! --- initialize random number generator 
      call ranset(seed)
! *** set up random number generator ***
!      call g05ccf
!      call g05cbf(54581)
!      idum = 5481
!      xini = ran1(idum)
 
! -------------------------------------------------------------------

! --- print 10 random numbers for control ---
      do i=1,10
         rtest(i) = random()
      end do
      write(71,1000) (rtest(i),i=1,5)
      write(71,1000) (rtest(i),i=6,10)
 1000 format(2x,5f10.6)

      close(71)
! --- call main program
      call monola

! ----------------------------------------------------------------

      end
