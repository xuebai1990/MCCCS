      function exzeof(xi,yi,zi,idi)

! exzeof
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
      real(8)::exzeof,xi,yi,zi,r2,rcutsq,xr,yr,zr,r2i,r6,velect
      integer::j,idi,idj,ntij,layer,ii,jj,kk
      logical::precise
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'control.inc'
      include 'external.inc'
      include 'system.inc'
      include 'poten.inc'      

      exzeof=0.
!     calculate the Lennard-Jones part
      rcutsq = rcut(1)**2
      do j=1,nzeo
         idj=idzeo(j)
         ntij = (idi - 1) * nntype + idj
         xr=xi-zeox(j)
         xr=xr-zeorx*anint(xr*zeorxi)
         yr=yi-zeoy(j)
         yr=yr-zeory*anint(yr*zeoryi)
         zr=zi-zeoz(j)
         zr=zr-zeorz*anint(zr*zeorzi)
         r2=xr*xr+yr*yr+zr*zr
         if (r2.lt.rcutsq) then
            r2i=sig2ij(ntij)/r2
            r6=r2i*r2i*r2i
!            if (lshift) then     
!              exzeof=exzeof+4.*zeps(idi,idj)*(r6-1.0)*r6-zencut(idi,idj)
!            else
            exzeof=exzeof+4.*epsij(ntij)*(r6-1.0)*r6
!            end if
         end if
      end do

!     calculate the Coulombic interaction, include as many layers of
!     neighboring unit cells as needed for the specified precision
      layer=0
      precise=.false.
      do while (.not. precise)
         velect=0.
         do j=1,nzeo
            if (zeox(j).lt.zunitx .and. zeoy(j).lt.zunity .and.
     &           zeoz(j).lt.zunitz) then
               do ii=-layer,layer
                  do jj=-layer,layer
                     do kk=-layer,layer
                        if (abs(ii).eq.layer .or. abs(jj).eq.layer .or.
     &                       abs(kk).eq.layer) then
                           xr=zeox(j)+ii*zunitx-xi
                           yr=zeoy(j)+jj*zunity-yi
                           zr=zeoz(j)+kk*zunitz-zi
                           r2=xr*xr+yr*yr+zr*zr
                           velect=velect+qelect(idi)*qelect(idj)/
     &                          dsqrt(r2)
                        end if
                     end do
                  end do
               end do
            end if
         end do
         exzeof=exzeof+velect
         if (velect.le.eps) precise=.true.
         layer=layer+1
      end do
      return
      end
 
