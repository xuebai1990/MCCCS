module const_phys
  use var_type, only: double_precision
  implicit none

  real(KIND=double_precision),parameter::eXV_to_K=11600.0_double_precision ! e times V to Kelvin

  ! conversion factor for coulomb interactions (incl e-->C, m-->A, J-->K, 4pi epsi naught)
  real(KIND=double_precision),parameter::qqfact=1.67125E+5_double_precision
  
  ! qqfact for bohr/hartree - from pot_KAng.f code
  !         qqfact = 0.99865377d0

end module const_phys
