      subroutine ztest

! ztest
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
      include 'grid.inc'
      include 'zeolite.inc'
      include 'control.inc'
      integer::i,idi,tel
      real(8)::errm,errt,err,res,rest,erra
      integer::samp
      real(8)::xi,yi,zi,random,exzeof,exzeo
     &     ,ebolt,bolt,boltf,eboltf
!     --- test accuracy
      idi=1
      errm=0
      errt=0
      erra=0
      samp=100 
      tel=0
!---  test accuracy
      write(16,*) ' Program uses a grid of the zeolite lattice '
      write(16,*) ' --- accuracy test '
      bolt=0
      ebolt=0
      boltf=0
      eboltf=0
      do i=1,samp
         if (mod(i,100).eq.0) print*,' test table ',i,
     &        ' out ',samp
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
               print*,err,res,rest
            end if
            if (erra.lt.abs(res-rest)) erra=abs(res-rest)
            errt=errt+err
            if (errm.lt.err) errm=err
         end if
      end do
      write(iou,*) ' test over : ',samp,tel, ' random positions '
      write(iou,*) ' average error : ',errt/tel
      write(iou,*) ' maximum error : ',errm
      write(iou,*) ' maximum ab.er.: ',erra
      write(iou,*) ' Boltz avera energy (table): ',ebolt/bolt
      write(iou,*) ' Boltz avera energy (full ): ',eboltf/boltf
      write(iou,*) ' ---------------------------- '
      write(iou,*) 

      return
      end
