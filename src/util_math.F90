module util_math
  use const_math,only:onepi
  implicit none
  private
  public::erfunc,mbessel,polint
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

!      subroutine coordinate_transform(x,y,z,invh,sx,sy,sz)
!      real,intent(in)::x,y,z,invhmat
!      real,intent(out)::sx,sy,sz

!      end subroutine coordinate_transform
end module util_math
