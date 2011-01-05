      module const_math
      use var_type, only: double_precision
      implicit none
      save
      real(KIND=double_precision),parameter::onepi =3.14159265358979323846264338_double_precision ,twopi=6.28318530717958647692528677_double_precision

      real(KIND=double_precision),parameter::raddeg =180.0_double_precision/onepi,degrad=onepi/180.0_double_precision

      real(KIND=double_precision),parameter::eps=1.0_double_precision
      end module const_math
