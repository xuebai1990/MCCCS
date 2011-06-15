module const_math
  use var_type, only: double_precision,single_precision
  implicit none
  private
  real(KIND=double_precision),parameter,public::onepi=3.1415926535897932384626433832795028841971693993751058209749_double_precision,twopi=6.2831853071795864769252867665590057683943387987502116419498_double_precision,fourpi=12.566370614359172953850573533118011536788677597500423283899_double_precision

  real(KIND=double_precision),parameter,public::raddeg=180.0_double_precision/onepi,degrad=onepi/180.0_double_precision
  
  real(KIND=double_precision),parameter,public::eps=epsilon(1.0_single_precision)
end module const_math
