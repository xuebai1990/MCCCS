      function exzeof(xi,yi,zi,idi)

c exzeof
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
      real(8)::exzeof,xi,yi,zi,r2,rcutsq
     +     ,xr,yr,zr,r2i,r6
      integer::j,idi,idj,ntij
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'control.inc'
      
      rcutsq = rcut(1)**2
      exzeof=0.
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
         if (r2.lt.zrc2(idi,idj)) then
            r2i=zsig2(idi,idj)/r2
            r6=r2i*r2i*r2i
            if (lshift) then     
              exzeof=exzeof+4.*zeps(idi,idj)*(r6-1.0)*r6-zencut(idi,idj)
            else
              exzeof=exzeof+4.*zeps(idi,idj)*(r6-1.0)*r6
            endif
         endif
      enddo
      return
      end
 
