      function exzeot(xi,yi,zi,idi)

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'grid.inc'
!$$$      include 'zeolite.inc' 
      real(KIND=double_precision)::exzeot,xi,yi,zi,xr,yr,zr,dy
      integer(KIND=int)::m,idi,j0,j,jp,k,k0,kp,l,l0,lp,mt,mp,
     &        pgrid
      parameter (m=2,mt=2*m+1)
      real(KIND=double_precision)::yjtmp(mt),yktmp(mt),yltmp(mt)
      real(KIND=double_precision)::xt(mt),yt(mt),zt(mt)
!
! --- determine cell parameters
!
      xr=xi-zunitx*anint(xi*zunitxi)
      if (xr.lt.0) xr=xr+zunitx
      j=int(xr*factx)
      yr=yi-zunity*anint(yi*zunityi)
      if (yr.lt.0) yr=yr+zunity
      k=int(yr*facty)
      zr=zi-zunitz*anint(zi*zunitzi)
      if (zr.lt.0) zr=zr+zunitz
      l=int(zr*factz)

! -- block m*m*m centered around: j,k,l
      mp=m+1
! --- store x,y,z values around xi,yi,zi in arrays
      do j0=-m,m
	 xt(j0+mp)=xzz(j+j0)
         yt(j0+mp)=yzz(k+j0)
         zt(j0+mp)=zzz(l+j0)
      end do
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
           end do
           call polint(zt,yltmp,mt,zr,yktmp(k0+mp),dy)
	 end do
         call polint(yt,yktmp,mt,yr,yjtmp(j0+mp),dy)
      end do
      call polint(xt,yjtmp,mt,xr,exzeot,dy)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software +3Y.

      subroutine polint(xa,ya,n,x,y,dy)
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      integer(KIND=int)::n,nmax
      real(KIND=double_precision)::dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer(KIND=int)::i,m,ns
      real(KIND=double_precision)::den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do  i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        end if
        c(i)=ya(i)
        d(i)=ya(i)
      end do
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
        end do
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        end if
        y=y+dy
      end do
      return
      END



