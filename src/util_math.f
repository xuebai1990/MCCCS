module util_math
  use var_type,only:normal_int,double_precision
  use const_math,only:onepi
  implicit none
  private
  public::erfunc,mbessel,polint
contains
! complementary error function
  function erfunc(x)
    real(KIND=double_precision)::erfunc
    real(KIND=double_precision),intent(in)::x

    real(KIND=double_precision),parameter::p=0.3275911d0,a1=0.254829592d0,a2=-0.284496736d0,a3=1.421413741d0,a4=-1.453152027d0,a5=1.061405429d0
    real(KIND=double_precision)::tt,eee

    eee = dexp(-x*x)
    tt = 1.0d0/(1.0d0 + p*x)
    erfunc = ((((a5*tt+a4)*tt+a3)*tt+a2)*tt+a1)*tt*eee
    return
  end function erfunc

  function mbessel(z,nu)
    real(KIND=double_precision)::mbessel
    real(KIND=double_precision),intent(in)::z,nu
        
!       mbessel = sqrt(onepi/(2.0d0*z))*exp(-z)*
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))

! -- simple form     
    mbessel = sqrt(onepi/(2.0d0*z))*exp(-z)
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))
  end function mbessel

!  (C) Copr. 1986-92 Numerical Recipes Software +3Y.
  subroutine polint(xa,ya,n,x,y)
    real(KIND=double_precision),intent(in)::xa(n),ya(n),x
    integer(KIND=normal_int),intent(in)::n
    real(KIND=double_precision),intent(out)::y

    integer(KIND=normal_int),parameter::nmax=10
    integer(KIND=normal_int)::i,m,ns
    real(KIND=double_precision)::dy,den,dif,dift,ho,hp,w,c(nmax),d(nmax)
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
    end subroutine polint


!      subroutine coordinate_transform(x,y,z,invh,sx,sy,sz)
!      real(KIND=double_precision),intent(in)::x,y,z,invhmat
!      real(KIND=double_precision),intent(out)::sx,sy,sz
    
!      end subroutine coordinate_transform
end module util_math
