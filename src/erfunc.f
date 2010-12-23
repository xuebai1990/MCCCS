      function erfunc(x)
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none

! complementary error function

      real(KIND=double_precision)::p,a1,a2,a3,a4,a5,x,erfunc,tt,eee
      parameter (p=0.3275911d0,a1=0.254829592d0,
     &     a2=-0.284496736d0,a3=1.421413741d0,
     &     a4=-1.453152027d0,a5=1.061405429d0)

      eee = dexp(-x*x)
      tt = 1.0d0/(1.0d0 + p*x)
      erfunc = ((((a5*tt+a4)*tt+a3)*tt+a2)*tt+a1)*tt*eee
      return
      end
