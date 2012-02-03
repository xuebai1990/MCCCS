      function exzeot(xi,yi,zi,idi)

c exzeot
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
      double precision exzeot,xi,yi,zi,xr,yr,zr,dy
      integer m,idi,j0,j,jp,k,k0,kp,l,l0,lp,mt,mp,
     +        pgrid
      parameter (m=2,mt=2*m+1)
      double precision yjtmp(mt),yktmp(mt),yltmp(mt)
      double precision xt(mt),yt(mt),zt(mt)
      include 'grid.inc'
      include 'zeolite.inc' 
c
c --- determine cell parameters
c
      xr=xi-zunitx*anint(xi*zunitxi)
      if (xr.lt.0) xr=xr+zunitx
      j=int(xr*factx)
      yr=yi-zunity*anint(yi*zunityi)
      if (yr.lt.0) yr=yr+zunity
      k=int(yr*facty)
      zr=zi-zunitz*anint(zi*zunitzi)
      if (zr.lt.0) zr=zr+zunitz
      l=int(zr*factz)

c -- block m*m*m centered around: j,k,l
      mp=m+1
c --- store x,y,z values around xi,yi,zi in arrays
      do j0=-m,m
	 xt(j0+mp)=xzz(j+j0)
         yt(j0+mp)=yzz(k+j0)
         zt(j0+mp)=zzz(l+j0)
      enddo
      do j0=-m,m
         jp=j+j0
	 if (jp.lt.0)    jp=jp+ngrx
	 if (jp.ge.ngrx) jp=jp-ngrx
         do k0=-m,m
	    kp=k+k0
	    if (kp.lt.0)    kp=kp+ngry
	    if (kp.ge.ngry) kp=kp-ngry
	    do l0=-m,m
	       lp=l+l0
	       if (lp.lt.0)    lp=lp+ngrz
	       if (lp.ge.ngrz) lp=lp-ngrz
               yltmp(l0+mp)=egrid(pgrid(jp,kp,lp,ngrx,ngry))
           enddo
           call polint(zt,yltmp,mt,zr,yktmp(k0+mp),dy)
	 enddo
         call polint(yt,yktmp,mt,yr,yjtmp(j0+mp),dy)
      enddo
      call polint(xt,yjtmp,mt,xr,exzeot,dy)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +3Y.

      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      integer n,nmax
      double precision dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      double precision den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do  i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo
      return
      END



