      function mbessel(z,nu)

        real(KIND=double_precision)::z,nu
        real(KIND=double_precision)::mbessel
        real(KIND=double_precision)::pi
        parameter (pi = 3.14159265359d0)
        
!       mbessel = sqrt(pi/(2.0d0*z))*exp(-z)*
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))

! -- simple form     
        mbessel = sqrt(pi/(2.0d0*z))*exp(-z)
!     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
!     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))
        end function mbessel

