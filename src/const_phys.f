      module const_phys
      use var_type, only: double_precision
      implicit none
      save
      real(KIND=double_precision),parameter::eXV_to_K
     & =11600.0_double_precision
      real(KIND=double_precision),parameter::qqfact=1.67E
     & +5_double_precision

      end module const_phys
