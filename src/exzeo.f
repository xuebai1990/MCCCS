      function exzeo(xi,yi,zi,idi)

      use grid
      implicit none 
      real(8)::exzeo,xi,yi,zi,r2,xr,yr,zr,r2i,r6
      integer::m,mt,mp,idi,idj,ntij,igtype,ibox=1
      integer::j,j0,jp,k,k0,kp,l,l0,lp
      parameter (m=2,mt=2*m+1,mp=m+1)
      real(8)::yjtmp(mt),yktmp(mt),yltmp(mt)
      real(8)::xt(mt),yt(mt),zt(mt),sx,sy,sz
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'control.inc'
      include 'external.inc'
      include 'system.inc'
      include 'cell.inc'
      include 'poten.inc'

      if (.not.lzgrid) then
! to be implemented
      else
!     calculation using a grid
!         write(iou,*) 'entering exzeo. xi yi zi idi',xi,yi,zi,idi

         do igtype=1,gntype
            if (gtable(igtype).eq.idi) exit
         end do
         if (igtype.gt.gntype) then
            call cleanup('exzeo: no such bead type')
         end if
         
! --- determine cell parameters
         sx=(xi*hmati(ibox,1)+yi*hmati(ibox,4)+zi*hmati(ibox,7))*nx
         sy=(xi*hmati(ibox,2)+yi*hmati(ibox,5)+zi*hmati(ibox,8))*ny
         sz=(xi*hmati(ibox,3)+yi*hmati(ibox,6)+zi*hmati(ibox,9))*nz
         sx=sx-floor(sx)
         sy=sy-floor(sy)
         sz=sz-floor(sz)
         xr=sx*hmat(ibox,1)/nx+sy*hmat(ibox,4)/ny+sz *hmat(ibox,7)/nz
         yr=sy*hmat(ibox,5)/ny+sz*hmat(ibox,8)/nz
         zr=sz*hmat(ibox,9)/nz
         j = sx*ngrx
         k = sy*ngry
         l = sz*ngrz

! ---    test if in the reasonable regime
         exzeo=1.0d+6
         if ( egrid(j,k,l,igtype).gt.exzeo) return
! --     block m*m*m centered around: j,k,l
! ---  set up hulp array: (allow for going beyond unit cell
!      for polynom fitting)
         do l0=-m,m
            lp=l+l0
            sz=dble(lp)/ngrz/nz
! ---    store x,y,z values around xi,yi,zi in arrays
            zt(l0+mp)=sz*hmat(ibox,9)
            if (lp.lt.0)    lp=lp+ngrz
            if (lp.ge.ngrz) lp=lp-ngrz
            do k0=-m,m
	       kp=k+k0
               sy=dble(kp)/ngry/ny
               yt(k0+mp)=sy*hmat(ibox,5)+sz*hmat(ibox,8)
	       if (kp.lt.0)    kp=kp+ngry
	       if (kp.ge.ngry) kp=kp-ngry
               do j0=-m,m
                  jp=j+j0
                  sx=dble(jp)/ngrx/nx
                  xt(j0+mp)=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox
     &                 ,7)
                  if (jp.lt.0)    jp=jp+ngrx
                  if (jp.ge.ngrx) jp=jp-ngrx
                  yjtmp(j0+mp)=egrid(jp,kp,lp,igtype)
               end do
               call polint(xt,yjtmp,mt,xr,yktmp(k0+mp))
	    end do
            call polint(yt,yktmp,mt,yr,yltmp(l0+mp))
         end do
         call polint(zt,yltmp,mt,zr,exzeo)
      end if
      return
      end

      subroutine polint(xa,ya,n,x,y)
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

