      function exzeo(xi,yi,zi,idi)

      use grid
      implicit none 
      real(8)::exzeo,xi,yi,zi,r2,rcutsq
     &     ,xr,yr,zr,r2i,r6
      integer::m,mt,mp,idi,idj,ntij,igtype
      integer::j,j0,jp,k,k0,kp,l,l0,lp
      parameter (m=2,mt=2*m+1)
      real(8)::yjtmp(mt),yktmp(mt),yltmp(mt)
      real(8)::xt(mt),yt(mt),zt(mt),dy
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'control.inc'
      include 'external.inc'
      include 'system.inc'
      include 'poten.inc'

      rcutsq = rcut(1)**2
      if (.not.lzgrid) then
         exzeo=0.
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
!            if (r2.lt.zrc2(idi,idj)) then
            if (r2 .lt. rcutsq) then
              r2i=sig2ij(ntij)/r2
              r6=r2i*r2i*r2i
!              if (lshift) then     
!                 exzeo=exzeo+4.*epsij(ntij)*(r6-1.0)*r6-
!     &                 zencut(idi,idj)
!              else
                 exzeo=exzeo+4.*epsij(ntij)*(r6-1.0)*r6+
     &             qelect(idi)*qelect(idj)/dsqrt(r2)
!              end if
           end if
         end do
      else
!     calculation using a grid
!         write(iou,*) 'entering exzeo. xi yi zi idi',xi,yi,zi,idi

         do igtype=1,gntype
            if (gtable(igtype).eq.idi) exit
         end do
         if (igtype.gt.gntype) then
            call cleanup('exzeo: no such bead type')
         end if
         
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
! ---    test if in the reosanble regime
         exzeo=1.0d+6
         if ( egrid(j,k,l,igtype).gt.exzeo) return
! --     block m*m*m centered around: j,k,l
         mp=m+1
! ---    store x,y,z values around xi,yi,zi in arrays
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
                  yltmp(l0+mp)=egrid(jp,kp,lp,igtype)
               end do
               call polint(zt,yltmp,mt,zr,yktmp(k0+mp),dy)
	    end do
            call polint(yt,yktmp,mt,yr,yjtmp(j0+mp),dy)
         end do
         call polint(xt,yjtmp,mt,xr,exzeo,dy)
      end if
      return
      end
 
      subroutine polint(xa,ya,n,x,y,dy)
!  (C) Copr. 1986-92 Numerical Recipes Software +3Y.

      implicit none
      integer::n,nmax
      real(8)::dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer::i,m,ns
      real(8)::den,dif,dift,ho,hp,w,c(nmax),d(nmax)
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

