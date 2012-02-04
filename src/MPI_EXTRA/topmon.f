      program topmon

c topmon
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

      integer i,seed
      double precision random,rtest
      dimension rtest(10)
      character *50 fileout
c --- RP added for MPI for random2.f
      integer nseed
      double precision random2,nrtest
      dimension nrtest(10)
      character *50 nfileout

c ----------------------------------------------------------------

      read(4,*)
      read(4,*) seed
      read(4,*)
      read(4,*) nseed
      fileout = 'Nrandomtest.dat'

      open(unit=71,FILE=fileout,status="unknown") 

c --- RP added for random2()
      nfileout = 'nNrandomtest.dat'
      open(unit=72,FILE=nfileout,status="unknown")
c -------------------
c      seed = 0

c --- initialize random number generator 
      call ranset(seed)

c ---RP added for new random2.f
      call nranset(nseed)
c ---------

c *** set up random number generator ***
c      call g05ccf
c      call g05cbf(54581)
c      idum = 5481
c      xini = ran1(idum)
 
c -------------------------------------------------------------------

c --- print 10 random numbers for control ---
      do i=1,10
         rtest(i) = random()
c ---RP added for random2()
         nrtest(i)= random2()
c -----------
      enddo
      write(71,1000) (rtest(i),i=1,5)
      write(71,1000) (rtest(i),i=6,10)
 1000 format(2x,5f10.6)

c -- RP added for random2()
      write(72,2000) (nrtest(i),i=1,5)
      write(72,2000) (nrtest(i),i=6,10)
 2000 format(2x,5f10.6)
      close(72)
c ----------

      close(71)
c --- call main program
      call monola

c ----------------------------------------------------------------

      stop
      end
