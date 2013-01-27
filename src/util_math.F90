module util_math
  use const_math,only:onepi
  implicit none
  private
  public::erfunc,mbessel,polint,spline,splint
contains
!> \brief complementary error function
  pure function erfunc(x)
    real::erfunc
    real,intent(in)::x

    real,parameter::p=0.3275911d0,a1=0.254829592d0,a2=-0.284496736d0,a3=1.421413741d0,a4=-1.453152027d0,a5=1.061405429d0
    real::tt,eee

    eee = dexp(-x*x)
    tt = 1.0d0/(1.0d0 + p*x)
    erfunc = ((((a5*tt+a4)*tt+a3)*tt+a2)*tt+a1)*tt*eee
    return
  end function erfunc

  pure function mbessel(z,nu)
    real::mbessel
    real,intent(in)::z,nu

! -- simple form
    mbessel = sqrt(onepi/(2.0d0*z))*exp(-z)
!     &         *(1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     &         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))
  end function mbessel

!  (C) Copr. 1986-92 Numerical Recipes Software +3Y.
  pure subroutine polint(xa,ya,n,x,y)
    real,intent(in)::xa(n),ya(n),x
    integer,intent(in)::n
    real,intent(out)::y

    integer,parameter::nmax=10
    integer::i,m,ns
    real::dy,den,dif,dift,ho,hp,w,c(nmax),d(nmax)
    ns=1
    dif=abs(x-xa(1))
    do i=1,n
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
  end subroutine polint

!> \brief Set up cubic spline derivative array
!>
!>   Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.,
!> y(i) = f(x(i)) with x ascending in order, and given values yp1 and ypn
!> for the first derivative of the interpolating function at the point 1
!> and n, respectively, this routine returns an array y2(1:n) of length n
!> which contains the second derivatives of the interpolating function at
!> the tabulated points x(i). If yp1 and/or ypn are equal to 1E30 or larger,
!> the routine is signaled to set the corresponding boundary condition for
!> a natural spline, with zero second derivative on that boundary
!> Parameter: NMAX is the largest anticipated value of n
!>
!>   Numerical recipes, 2nd ed, 1992.
  subroutine spline(x,y,n,yp1,ypn,y2)
    integer::n,NMAX
    real::yp1,ypn,x(n),y(n),y2(n)
    parameter (NMAX=500)

    integer::i,k
    real::p,qn,sig,un,u(NMAX)

    if (yp1.gt.0.99d30) then
       y2(1) = 0.0d0
       u(1) = 0.0d0
    else
       y2(1) = -0.5d0
       u(1) = (3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    end if
    do i = 2,n-1
       sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p = sig*y2(i-1)+2.0d0
       y2(i) = (sig-1.0d0)/p
       u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    if (ypn .gt. 0.99d30) then
       qn = 0.0d0
       un = 0.0d0
    else
       qn = 0.5d0
       un = (3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
    do k = n-1, 1, -1
       y2(k) = y2(k)*y2(k+1)+u(k)
    end do
    return
  end subroutine spline

!> \brief Spline interpolation
!>
!>   Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
!> function (with the xa(i) in order), and given the array y2a(1:n), which
!> is from the ouput from spline, and given a value of x, this routine
!> returns a cubic-spline interpolated value y.
  subroutine splint(xa,ya,y2a,n,x,y)
    use util_runtime,only:err_exit
    use util_search,only:locate
    integer::n
    real::x,y,xa(n),y2a(n),ya(n)
    integer::khi,klo
    real::a,b,h

    klo = locate(xa,n,x,2)
    khi = klo+1
    h = xa(khi)-xa(klo)
    if (h.eq.0.) call err_exit('bad xa input in splint')
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
    return
  end subroutine splint

!      subroutine coordinate_transform(x,y,z,invh,sx,sy,sz)
!      real,intent(in)::x,y,z,invhmat
!      real,intent(out)::sx,sy,sz

!      end subroutine coordinate_transform
end module util_math
