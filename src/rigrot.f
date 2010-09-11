      subroutine rigrot(lnew,lterm,iskip,imol,imolty,ibox,wadd )

! rigrot
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


!     *******************************************************************
!     **      performs a rotational configurational bias move          **
!     *******************************************************************


!     &*&*&*&*&*&*&*&*&*&*&*&* IMPORTANT *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&
!             All Rigid Molecules Should Come After Riutry
!     &*&*&*&*&*&*&*&*&*&*&*&* IMPORTANT *&*&*&*&*&*&*&*&*&*&*&*&*&*&*&*&

      implicit none

      include 'control.inc'
      include 'cbmc.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'rosen.inc'
      include 'system.inc'
      include 'conver.inc'
      include 'ipswpar.inc'

      logical::lnew,ovrlap,lovrr,lterm,ltors

      integer::i,ibox,igrow,iii,j,ip,iwalk,iunit,imolty,iu,iskip
      integer::ichoi,imol,istt,iend,stch,glist,ntogrow,count

      real(8)::rx,ry,rz,rxorig,ryorig,rzorig,rxnw,rynw,rznw
      real(8)::xdgamma,ydgamma,zdgamma,xcosdg,xsindg,ycosdg
     &     ,ysindg,zcosdg,zsindg,rbf,rxur,ryur,rzur,delo,length
      real(8)::delen,v,vtor,w,bsum,vdum,vintra,vinter
     &     ,velect,vewald,wadd,vext,random,maxlen,bs
  
      dimension lovrr(nchmax),delen(nchmax),rxur(numax)
      dimension rzur(numax),ryur(numax),glist(numax)
     

! ----------------------------------------------------------------------
      
!      write(iou,*) 'start RIGROT'
     
!     --- initialize conformation energies and weight

      wadd = 1.0d0
      w = 0.0d0
      ichoi = nchoir(imolty)
      ltors = .false.
      iunit = nunit(imolty)
      igrow = riutry(imolty,1) 
      istt = igrow + 1 
      iend = iunit

        
      if (lnew) then
         rxorig = rxnew(igrow)
         ryorig = rynew(igrow)
         rzorig = rznew(igrow)

         do j = igrow+1, iunit
            rxur(j) = rxnew(j)
            ryur(j) = rynew(j)
            rzur(j) = rznew(j)

         end do
      else
         rxorig = rxu(imol,igrow)
         ryorig = ryu(imol,igrow)
         rzorig = rzu(imol,igrow)
         
         do j = igrow+1, iunit
            rxur(j) = rxu(imol,j)
            ryur(j) = ryu(imol,j)
            rzur(j) = rzu(imol,j)
         end do
      end if

!     --- find maxlen
      maxlen = 0.0d0
     
      do j = igrow+1, iunit
         length = (rxur(j)-rxorig)**2+(ryur(j)-ryorig)**2
     &        +(rzur(j)-rzorig)**2

         if (length.gt.maxlen) maxlen = length
      end do
      
      maxlen = dsqrt(maxlen)


      do ip = 1, ichoi 
 
         if (lnew.or.ip.ne.1) then

!            xdgamma = 8.0d0*atan(1.0d0)*random()
!            ydgamma = 8.0d0*atan(1.0d0)*random()
!            zdgamma = 8.0d0*atan(1.0d0)*random()
            xdgamma = twopi*random()
            ydgamma = twopi*random()
            zdgamma = twopi*random()
!     --- set up rotation matrix
            xcosdg = dcos(xdgamma)
            xsindg = dsin(xdgamma)
            ycosdg = dcos(ydgamma)
            ysindg = dsin(ydgamma)
            zcosdg = dcos(zdgamma)
            zsindg = dsin(zdgamma)
!     --- set molecule to rotate around
             
!     --- rotate around all axis  
            count = 0            
            do j = igrow+1, iunit 
               count = count + 1
               ry = ryur(j) - ryorig
               rz = rzur(j) - rzorig                     
               rynw = xcosdg * ry + xsindg * rz
               rznw = xcosdg * rz - xsindg * ry
               ryp(count,ip) = ryorig + rynw
               rzp(count,ip) = rzorig + rznw
            
               rx = rxur(j) - rxorig
               rz = rzp(count,ip) - rzorig
               rxnw = ycosdg * rx - ysindg * rz
               rznw = ycosdg * rz + ysindg * rx
               rxp(count,ip) = rxorig + rxnw
               rzp(count,ip) = rzorig + rznw
            
               rx = rxp(count,ip) - rxorig
               ry = ryp(count,ip) - ryorig
               rxnw = zcosdg * rx + zsindg * ry
               rynw = zcosdg * ry - zsindg * rx
               rxp(count,ip) = rxorig + rxnw
               ryp(count,ip) = ryorig + rynw
            end do

         else
            count = 0
            do j = igrow+1, iunit
               count = count + 1
               rxp(count,ip) = rxu(imol,j)
               ryp(count,ip) = ryu(imol,j)
               rzp(count,ip) = rzu(imol,j)
            end do
         end if
      end do         
   
      ntogrow = iunit-igrow
      
      count = 0
      do j=1,igrow-1
          lexist(j) = .false.
      end do

!      write(iou,*) 'igrow',igrow
 
      lexist(igrow) = .true.
     
      do j = igrow+1, iunit
         count = count + 1
         glist(count) = j
         lexist(j) = .false.
      end do
      
      call boltz( lnew,.false.,ovrlap,iskip,imol,imolty,ibox
     &     ,ichoi,igrow,ntogrow,glist,maxlen )

      if (ovrlap) then
         lterm = .true.
         return
      end if
      
      bsum = 0.0d0
      do ip = 1, ichoi
         bsum = bsum + bfac(ip)
      end do

      wadd = wadd * bsum

      if (lnew) then
         if ( wadd .lt. softlog ) then
            lterm = .true.
            return
         end if

!     --- select one position at random ---
         rbf = bsum * random()
         bs = 0.0d0
         do ip = 1, ichoi
            if (.not. lovr(ip) ) then
               bs = bs + bfac(ip)
               if ( rbf .lt. bs ) then
!     --- select ip position
                  iwalk = ip
                  goto 15
               end if
            end if
         end do
         write(iou,*) 'screwup in rigrot'
         call cleanup('screwup in rigrot')
 15      continue
            
         
      else
         iwalk = 1
         if ( wadd .lt. softlog ) then
            write(iou,*) '###old rigrot weight too low'
         end if
      end if
      
      if ( lnew ) then
         vnewt = vnewt + vtry(iwalk)
         vnewext = vnewext + vtrext(iwalk)
         vnewinter = vnewinter + vtrinter(iwalk)
!         vnewintra = vnewintra + vtrintra(iwalk)
         vnewelect = vnewelect + vtrelect(iwalk)
         vnewewald = vnewewald + vtrewald(iwalk)
         vipswn = vipswn+vipswnt(iwalk)
         vwellipswn = vwellipswn+vwellipswnt(iwalk)
      else
         voldt = voldt + vtry(iwalk)
         voldext = voldext + vtrext(iwalk)
         voldinter = voldinter + vtrinter(iwalk)
!         voldintra = voldintra + vtrintra(iwalk)
         voldelect = voldelect + vtrelect(iwalk)
         voldewald = voldewald + vtrewald(iwalk)
         vipswo = vipswo+vipswot(iwalk)
         vwellipswo = vwellipswo+vwellipswot(iwalk)
      end if
     
!	write(iou,*) 'vtry', vtry(iwalk),iwalk
    
      do j = 1, ntogrow
         iu = glist(j)
         if (lnew) then
            rxnew(iu) = rxp(j,iwalk)
            rynew(iu) = ryp(j,iwalk)
            rznew(iu) = rzp(j,iwalk)
         end if
      end do


!      write(iou,*) 'end RIGROT'
 
      return
      end











