module const_phys
  use var_type, only: double_precision
  implicit none

  real,parameter::cal2joule=4.184_double_precision,joule2cal=1.0_double_precision/cal2joule

  real,parameter::eXV_to_K=11600.0_double_precision ! e times V to Kelvin

  real,parameter::MPa2SimUnits=7.2429e-02_double_precision ! Conversion factor for Mpa to simulation unit

  ! conversion factor for coulomb interactions (incl e-->C, m-->A, J-->K, 4pi epsi naught)
  real,parameter::qqfact=1.67125E+5_double_precision

  ! qqfact for bohr/hartree - from pot_KAng.f code
  !         qqfact = 0.99865377d0

end module const_phys
