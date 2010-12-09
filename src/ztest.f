      subroutine ztest(idi)

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

      use grid
      implicit none
      include 'zeolite.inc'
      include 'control.inc'
      include 'mpi.inc'
      integer::i,idi,tel,samp
      real(8)::errm,errt,err,res,rest,erra
     &     ,xi,yi,zi,random,exzeof,exzeo
     &     ,ebolt,bolt,boltf,eboltf
     &     ,halfpx,halfmx,halfpy,halfmy,halfpz,halfmz

!--- test accuracy
      if (myid.eq.0) write(16,*) ' Program uses a grid of the
     & zeolite lattice / --- accuracy test '
      samp=20 
      tel=0
      errm=0
      errt=0
      erra=0
      bolt=0
      ebolt=0
      boltf=0
      eboltf=0
      do i=1,samp
         if (myid.eq.0 .and. mod(i,5).eq.0) write(iou,*) 'test table '
     &        ,i,' out of ',samp
         xi=random()*zunitx
         yi=random()*zunity
         zi=random()*zunitz
         res=exzeo(xi,yi,zi,idi)
         rest=exzeof(xi,yi,zi,idi)
         halfpx=xi+zunitx/2
         if (halfpx>zunitx) halfpx=halfpx-zunitx
         halfmx=zunitx/2-xi
         if (halfmx<0) halfmx=halfmx+zunitx
         halfpy=yi+zunity/2
         if (halfpy>zunity) halfpy=halfpy-zunity
         halfmy=zunity/2-yi
         if (halfmy<0) halfmy=halfmy+zunity
         halfpz=zi+zunitz/2
         if (halfpz>zunitz) halfpz=halfpz-zunitz
         halfmz=zunitz/2-zi
         if (halfmz<0) halfmz=halfmz+zunitz
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',xi,yi,zi,')=',rest
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',zunitx-xi,halfpy,
     &        zunitz-zi,')=',exzeof(zunitx-xi,halfpy,zunitz-zi,idi)
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',halfpx,yi
     &        ,halfmz,')=',exzeof(halfpx,yi,halfmz,idi)
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',halfmx,halfpy,
     &        halfpz,')=',exzeof(halfmx,halfpy,halfpz,idi)
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',halfmx,zunity-yi,
     &        halfpz,')=',exzeof(halfmx,zunity-yi,halfpz,idi)
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',halfpx,halfmy,
     &        halfmz,')=',exzeof(halfpx,halfmy,halfmz,idi)
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',zunitx-xi,zunity-yi
     &        ,zunitz-zi,')=',exzeof(zunitx-xi,zunity-yi,zunitz-zi,idi)
         write(iou,'(A,3(F8.5,1X),A,G20.7)') 'U(',xi,halfmy,
     &        zi,')=',exzeof(xi,halfmy,zi,idi)
         write(iou,*) '------------------------'
         if (dexp(-res/300.).gt. 1.e-5) then
            tel=tel+1
            ebolt=ebolt+res*dexp(-res/300.)
            bolt=bolt+dexp(-res/300.)
            boltf=boltf+dexp(-rest/300.)
            eboltf=eboltf+rest*dexp(-rest/300.)
            err=abs((res-rest)/res)
            if (myid.eq.0 .and. err.gt.0.03) then
               write(iou,*) ' warning err. interp.larger than 3% '
               write(iou,*) err,res,rest
            end if
            errt=errt+err
            if (erra.lt.abs(res-rest)) erra=abs(res-rest)
            if (errm.lt.err) errm=err
         end if
      end do

      if (myid.eq.0) then
         write(iou,*) ' test over : ',samp,tel, ' random positions '
         write(iou,*) ' average error : ',errt/tel
         write(iou,*) ' maximum error : ',errm
         write(iou,*) ' maximum ab.er.: ',erra
         write(iou,*) ' Boltz avera energy (table): ',ebolt/bolt
         write(iou,*) ' Boltz avera energy (full ): ',eboltf/boltf
         write(iou,*) ' ---------------------------- '
         write(iou,*) 
      end if

      return
      end
