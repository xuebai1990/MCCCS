      subroutine ztest

c ztest
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
      include 'grid.inc'
      include 'zeolite.inc'
      integer i,idi,tel
      double precision errm,errt,err,res,rest,erra
      integer samp
      double precision xi,yi,zi,random,exzeof,exzeo
     +     ,ebolt,bolt,boltf,eboltf
c     --- test accuracy
      idi=1
      errm=0
      errt=0
      erra=0
      samp=100 
      tel=0
c---  test accuracy
      write(16,*) ' Program uses a grid of the zeolite lattice '
      write(16,*) ' --- accuracy test '
      bolt=0
      ebolt=0
      boltf=0
      eboltf=0
      do i=1,samp
         if (mod(i,100).eq.0) print*,' test table ',i,
     +        ' out ',samp
         xi=(4-8*random())*zeorx
         yi=(4-8*random())*zeory
         zi=(4-8*random())*zeorz
         res=exzeo(xi,yi,zi,idi)
         ebolt=ebolt+res*dexp(-res/300.)
         bolt=bolt+dexp(-res/300.)
         rest=exzeof(xi,yi,zi,idi)
         boltf=boltf+dexp(-rest/300.)
         eboltf=eboltf+rest*dexp(-rest/300.)
         if (dexp(-res/300.).gt. 1.e-5) then
            tel=tel+1
            err=abs((res-rest)/res)
            if (err.gt.0.03) then
               print*,' warning err. interp.larger than 3% '
               print*,sngl(err),sngl(res),sngl(rest)
            endif
            if (erra.lt.abs(res-rest)) erra=abs(res-rest)
            errt=errt+err
            if (errm.lt.err) errm=err
         endif
      enddo
      write(6,*) ' test over : ',samp,tel, ' random positions '
      write(6,*) ' average error : ',sngl(errt/tel)
      write(6,*) ' maximum error : ',sngl(errm)
      write(6,*) ' maximum ab.er.: ',sngl(erra)
      write(6,*) ' Boltz avera energy (table): ',sngl(ebolt/bolt)
      write(6,*) ' Boltz avera energy (full ): ',sngl(eboltf/boltf)
      write(6,*) ' ---------------------------- '
      write(6,*) 

      return
      end
