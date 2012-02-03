      subroutine ztab

c ztab
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
      include 'control.inc'
      integer i,j,k,pgrid,idi,ngrid
      double precision xi,yi,zi,exzeof,dgr
c make a tabulated potential of the zeolite
c
c     ideal grid size: dgr

      zunitx = zeorx/nx
      zunity = zeory/ny
      zunitz = zeorz/nz
      zunitxi=1./zunitx
      zunityi=1./zunity 
      zunitzi=1./zunitz 
c --- volume minimal box
      dgr = zunitx * zunity * zunitz/maxtab
      dgr=(dgr)**(1.d00/3.d00)
c --- find closest values for x and y
      ngrx=int(zunitx/dgr)
      ngry=int(zunity/dgr)
      ngrz=int(zunitz/dgr)
      ngrid=ngrx*ngry*ngrz
      if (ngrid.gt.maxtab) stop  'ngrid.gt.maxtab'
      dgrx=zunitx/ngrx
      dgry=zunity/ngry
      dgrz=zunitz/ngrz
      factx=ngrx/zunitx
      facty=ngry/zunity
      factz=ngrz/zunitz
c ---  set up hulp array: (allow for going beyond unit cell
c      for polynom fitting)
      do i=-maxp,ngrx+maxp-1
         xzz(i)=i*dgrx
      enddo
      do i=-maxp,ngry+maxp-1
         yzz(i)=i*dgry
      enddo
      do i=-maxp,ngrz+maxp-1
         zzz(i)=i*dgrz
      enddo
      write(iou,1000) ngrid,ngrx,dgrx,ngry,dgry,ngrz,dgrz
c --- calculate interactions on grid
      idi=1
      xi=-dgrx
      do i=0,ngrx-1
         xi=xi+dgrx
	 yi=-dgry
         do j=0,ngry-1
            yi=yi+dgry
	    zi=-dgrz
            do k=0,ngrz-1
               zi=zi+dgrz
	       egrid(pgrid(i,j,k,ngrx,ngry))=sngl(exzeof(xi,yi,zi,idi))
	    enddo
	 enddo
         write(iou,*) 'Ztable: done ',i+1,' out of ',ngrx
      enddo
c--   write table to disk
      call rwztab(1) 

      return
     
 1000 format(' Making the zeolite grid: ',/,
     +       '   total number of cells : ',i12,/,
     +       '         x-dir           : ',i12,'  size: ',f7.3,/,
     +       '         y-dir           : ',i12,'  size: ',f7.3,/,
     +       '         z-dir           : ',i12,'  size: ',f7.3,/)
      end

