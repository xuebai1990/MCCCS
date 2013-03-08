      subroutine suzeo(rczeo)

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

      implicit none 
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'external.inc'
      include 'grid.inc'
      include 'control.inc'
      integer::imol,iunit,igtype,idi,ierr,i,j,k,ngrxt,ngryt,ngrzt
      character(len=8)::filename
      real(8)::rczeo,zunitxt,zunityt,zunitzt

! === load force field
!      call forcefield(rczeo)

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
         write(iou,1000) zunitx,zunity,zunitz,
     &        ngrx,dgrx,ngry,dgry,ngrz,dgrz

         allocate(egrid(0:ngrx-1,0:ngry-1,0:ngrz-1,gntype),
     &        xzz(-maxp,ngrx+maxp-1),yzz(-maxp,ngry+maxp-1),
     &        zzz(-maxp,ngrz+maxp-1),stat=ierr)
         if (ierr.ne.0) call cleanup('allocate failed')

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

         do igtype=1,gntype
            idi=gtable(igtype)
            write(filename,'(I3.3,A)'),idi,'.ztb'
            print*,filename
            open(91,file=filename,iostat=ierr,status=old)
            if (ierr.eq.0) then
! ---    read zeolite table from disk
               read(91) zunitxt,zunityt,zunitzt,ngrxt,ngryt,ngrzt
               if (abs(zunitxt-zunitx).gt.eps .or. abs(zunityt-zunity)
     &              .gt.eps .or. abs(zunitzt-zunitz).gt.eps .or.
     &              ngrxt.ne.ngrx .or. ngryt.ne.ngry .or. ngrzt.ne.ngrz)
     &              call cleanup('problem in zeolite potential table')
               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        read(91) egrid(i,j,k,idi)
                     end do
                  end do
               end do
            else
! make a tabulated potential of the zeolite
               open(91,file=filename,status=unknown)
               write(91) zunitx,zunity,zunitz,ngrx,ngry,ngrz
               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        egrid(i,j,k,idi)=exzeof(dgrx*i,dgry*j,
     &                       dgrz*k,idi)
                        write(91) egrid(i,j,k,idi)
                     end do
                  end do
               end do
            end if
            close(91)
         end do
         call ztest
      end if

      return

 1000 format(' Tabulated zeolite potential: ',/,
     &     ' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/)
     &     '         x-dir           : ',i12,'  size: ',f7.4,/,
     &     '         y-dir           : ',i12,'  size: ',f7.4,/,
     &     '         z-dir           : ',i12,'  size: ',f7.4,/)

      end  
