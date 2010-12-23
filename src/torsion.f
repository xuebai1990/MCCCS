      subroutine lininter(theta,spltor,ttyp)

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
!$$$      include 'torsion.inc'
!$$$      include 'conver.inc'
!$$$      include 'control.inc'

      integer(KIND=normal_int)::ttyp,klo,khi,k,xa,bin,addl
      real(KIND=double_precision)::theta,spltor,thetarem,tordiff,torstep
     & ,left

! Routine to calculate torsion potential using linear interpolation between 2 points
! Requires a file (fort.40) running from -180 to 180 in 1/4 degree intervals

      theta=raddeg*theta

      klo=1
      khi=splpnts(ttyp)

      xa=idint(theta)
      left=(theta-xa)
!      write(iou,*) 'left',left

!     selecting bin of correct degree
      addl = idint(left*4)
      bin = (xa-deg(1,ttyp))*4 + 1 + addl
      if (theta .lt. 0.0d0.and.theta.ge.-180.0d0) then
!         write(iou,*) 'negative bin',bin
         khi = bin
         klo = khi-1
      elseif (theta .ge. 0.0d0.and.theta.le.180.0d0) then
!         write(iou,*) 'positive bin',bin
         klo = bin
         khi = klo+1
      else
         write(iou,*) 'Error in lininter.f - theta',theta
      end if

!      write(iou,*) 'klo,khi',klo,deg(klo,ttyp),khi,deg(khi,ttyp)
! check
      if(deg(klo,ttyp).gt.theta.or.deg(khi,ttyp).lt.theta) then
          write(iou,*) 'problem below'
          write(iou,*) 'theta',theta,' ttyp',ttyp
          write(iou,*) 'klo',klo,deg(klo,ttyp),'khi',khi,deg(khi,ttyp)
          write(iou,*)
       end if

      thetarem=theta-deg(klo,ttyp)
      thetarem=thetarem*4.0d0

!      tordiff=tabtorso(khi,ttyp)-tabtorso(klo,ttyp)
      spltor=thetarem*tordif(klo,ttyp)
      spltor=spltor+tabtorso(klo,ttyp)

      return
      end



      subroutine spline(yp1,ypn,tortyp)

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'torsion.inc'
!$$$      include 'control.inc'
      
!      integer(KIND=normal_int)::n,nmax
!      real(KIND=double_precision)::yp1,ypn,x(n),y(n),y2(n)
!      parameter (nmax=500)

!     Numerical recipes, 2nd ed, 1992. - cubic spline
!     Given arrays x(1:n) and y(1:n) containing a tabulated function, ie 
!     y(i) = f(x(i)) with x ascending in order, and given values yp1 and ypn
!     for the first derivative of the interpolating function at the point 1
!     and n, respectively, this routine returns an array y2(1:n) of length n
!     which contains the second derivatives of the interpolating function at 
!     the tabulated points x(i). If yp1 and/or ypn are equal to 1E30 or more,
!     the routine is signaled to set the corresponding boundary condition for 
!     a natural spline, with zero second derivative on that boundary
!     nmax = largest value of n
!     x=deg y=tabtorso y2=torderiv2 n=points

      integer(KIND=normal_int)::i,k,points,tortyp
      real(KIND=double_precision)::p,qn,sig,un,u(500,10),yp1,ypn
      
! Routine sets up derivatives for use with spline interpolations
! Requires file (fort.40) running from -195 to 195 in degree steps
! (Extra 15 degrees on each side required so that second derivatives are reasonable 
!     by the time the degrees of interest are reached.)
      points=splpnts(tortyp)

      write(iou,*) 'beginning of spline',splpnts(tortyp),yp1,ypn

      if (yp1.gt.0.99d30) then
         torderiv2(1,tortyp) = 0.0d0
         u(1,tortyp) = 0.0d0

      else
         torderiv2(1,tortyp) = -0.5d0
         u(1,tortyp) = (3.0d0/(deg(2,tortyp)-deg(1,tortyp)))*
     &        ((tabtorso(2,tortyp)-tabtorso(1,tortyp))/
     &        (deg(2,tortyp)-deg(1,tortyp))-yp1)

      end if

      do i = 2,points-1
         sig = (deg(i,tortyp)-deg(i-1,tortyp))/
     &        (deg(i+1,tortyp)-deg(i-1,tortyp))
         p = sig*torderiv2(i-1,tortyp)+2.0d0
         torderiv2(i,tortyp) = (sig-1.0d0)/p
         u(i,tortyp) = (6.0d0*((tabtorso(i+1,tortyp)-tabtorso(i,tortyp))
     &        /(deg(i+1,tortyp)-deg(i,tortyp))-(tabtorso(i,tortyp)-
     &        tabtorso(i-1,tortyp))/(deg(i,tortyp)-deg(i-1,tortyp)))/
     &        (deg(i+1,tortyp)-deg(i-1,tortyp))-sig*u(i-1,tortyp))/p
      end do


      if (ypn .gt. 0.99d30) then
         qn = 0.0d0
         un = 0.0d0

      else
         qn = 0.5d0
         un = (3.0d0/(deg(points,tortyp)-deg(points-1,tortyp)))*
     &        (ypn-(tabtorso(points,tortyp)-tabtorso(points-1,tortyp))
     &        /(deg(points,tortyp)-deg(points-1,tortyp)))

      end if

      torderiv2(points,tortyp) = (un-qn*u(points-1,tortyp))/
     &     (qn*torderiv2(points-1,tortyp)+1.0d0)
      do k = points-1, 1, -1
         torderiv2(k,tortyp) = torderiv2(k,tortyp)*torderiv2(k+1,tortyp)
     &        +u(k,tortyp)
      end do

      do i=1,points
         write(55,*) deg(i,tortyp),tabtorso(i,tortyp),
     &        torderiv2(i,tortyp)
      end do

      return
      end


      subroutine splint(x,y,tortyp)
!     From Numerical recipes, 2nd ed, 1992. - spline interpolation
!     Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
!     function and given the array y2a(1:n), which is from the ouput of 
!     spline, and given a value of x, this routine returns a cubic-spline
!     interpolated value of y.
!     xa = deg, ya = tabtorso, y2a = torderiv2, n = points
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'torsion.inc'
!$$$      include 'conver.inc'
!$$$      include 'control.inc'

!      integer(KIND=normal_int)::n
!      real(KIND=double_precision)::x,y,xa(n),y2a(n),ya(n)
      integer(KIND=normal_int)::k, khi,klo,points,jttor,tortyp,xa
      real(KIND=double_precision)::a,bb,h,x,y,vtorso,vtorsoa

      points = splpnts(tortyp)

      klo = 1
      khi = points

      x=x*raddeg

! Below is correct for tabulated data from -195 to 195 degrees
      xa=int(x)+197

      if (x.lt.0.0d0) then
         klo=xa-1
         khi=xa
      elseif (xa.ge.0.0d0) then
         klo=xa
         khi=xa+1
      end if

!      write(iou,*) 'klo,khi',deg(klo,tortyp),deg(khi,tortyp)

      h = deg(khi,tortyp)-deg(klo,tortyp)
     
      if (dabs(h).lt.1.0d-8) write(iou,*) 'bad deg input in splint',h,
     &     khi,deg(khi,tortyp),klo,deg(klo,tortyp)
      a = (deg(khi,tortyp)-x)/h
      bb = (x-deg(klo,tortyp))/h
      y = a*tabtorso(klo,tortyp)+bb*tabtorso(khi,tortyp)+
     &     ((a**3-a)*torderiv2(klo,tortyp)+(bb**3-bb)*
     &     torderiv2(khi,tortyp))*(h**2)/6.0d0

!      write(iou,*) x,y,tortyp

      return
      end
