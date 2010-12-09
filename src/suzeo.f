      subroutine suzeo(rczeo,maxlayer)

! suzeo
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
      include 'control.inc'
      include 'coord.inc'
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'external.inc'
      include 'mpi.inc'
      integer::imol,iunit,igtype,idi,jerr,i,j,k,
     &     ngrxt,ngryt,ngrzt,maxlayer
      character(len=8)::filename
      real(8)::rczeo,zunitxt,zunityt,zunitzt,exzeof

! === load force field
!      call forcefield(rczeo)
      ztype(1)=177
      ztype(2)=178

! === load positions of zeolite atoms
      call zeocoord()
      zeorxi=1./zeorx
      zeoryi=1./zeory
      zeorzi=1./zeorz
!
! === tabulation of the zeolite potential
      if (lzgrid) then
! find all bead types and store them in an array
         gntype=0
         do imol=1,nmolty
            do iunit=1,nunit(imol)
               do igtype=1,gntype
                  if (gtable(igtype).eq.ntype(imol,iunit)) exit
               end do
               if (igtype.gt.gntype) then
                  gntype=igtype
                  gtable(gntype)=ntype(imol,iunit)
               end if
            end do
         end do

         zunitx = zeorx/nx
         zunity = zeory/ny
         zunitz = zeorz/nz
         zunitxi=1./zunitx
         zunityi=1./zunity 
         zunitzi=1./zunitz 
!     --- volume minimal box
!     ideal grid size: dgr
         dgr = 0.2
!     dgr=(dgr)**(1.d00/3.d00)
!     --- find closest values for x and y
         ngrx=int(zunitx/dgr)
         ngry=int(zunity/dgr)
         ngrz=int(zunitz/dgr)
!     if (ngrid.gt.maxtab) call cleanup('ngrid.gt.maxtab')
         dgrx=zunitx/ngrx
         dgry=zunity/ngry
         dgrz=zunitz/ngrz
         factx=ngrx/zunitx
         facty=ngry/zunity
         factz=ngrz/zunitz
         if (myid.eq.0) write(iou,1000) zunitx,zunity,zunitz,
     &        ngrx,dgrx,ngry,dgry,ngrz,dgrz


         allocate(egrid(0:ngrx-1,0:ngry-1,0:ngrz-1,gntype),
     &        xzz(-maxp:ngrx+maxp-1),yzz(-maxp:ngry+maxp-1),
     &        zzz(-maxp:ngrz+maxp-1),stat=jerr)
         if (jerr.ne.0) call cleanup('allocate failed')

! ---  set up hulp array: (allow for going beyond unit cell
!      for polynom fitting)
         do i=-maxp,ngrx+maxp-1
            xzz(i)=i*dgrx
         end do
         do i=-maxp,ngry+maxp-1
            yzz(i)=i*dgry
         end do
         do i=-maxp,ngrz+maxp-1
            zzz(i)=i*dgrz
         end do

         nlayermax=0

         do igtype=1,gntype
            idi=gtable(igtype)
            write(filename,'(I3.3,A)'),idi,'.ztb'
            open(91,file=filename,form='binary',iostat=jerr,
     &           status='old')
            if (jerr.eq.0) then
! ---    read zeolite table from disk
               if (myid.eq.0) write(iou,*) 'read in ',filename
               read(91) zunitxt,zunityt,zunitzt,ngrxt,ngryt,ngrzt
               if (abs(zunitxt-zunitx).gt.eps .or. abs(zunityt-zunity)
     &              .gt.eps .or. abs(zunitzt-zunitz).gt.eps .or.
     &              ngrxt.ne.ngrx .or. ngryt.ne.ngry .or. ngrzt.ne.ngrz)
     &              call cleanup('problem in zeolite potential table')
               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        read(91) egrid(i,j,k,igtype)
                     end do
                  end do
               end do
            else
! make a tabulated potential of the zeolite
               if (myid.eq.0) write(iou,*) 'make new ',filename
               if (myid.eq.0) open(91,file=filename,form='binary')
               if (myid.eq.0) write(91) zunitx,zunity,zunitz,ngrx,
     &              ngry,ngrz
               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        egrid(i,j,k,igtype)=exzeof(dgrx*i,dgry*j,
     &                       dgrz*k,idi)
                        if (myid.eq.0) write(91) egrid(i,j,k,igtype)
                     end do
                  end do
               end do
               if (myid.eq.0) write(iou,*) 'maxlayer = ',nlayermax
!               call ztest(idi)
            end if
            if (myid.eq.0) then
               close(91)
               close(16)
            end if
         end do        
      end if

      return

 1000 format(' Tabulated zeolite potential: ',/,
     &     ' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/
     &     '         x-dir           : ',i12,'  size: ',f7.4,/,
     &     '         y-dir           : ',i12,'  size: ',f7.4,/,
     &     '         z-dir           : ',i12,'  size: ',f7.4,/)

      end  
