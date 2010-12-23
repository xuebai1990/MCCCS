      module var_type
! Define variable types(kinds) based on required precision, rather than
! the default values that rely on the computer architecture and may
! cause inconsistent behavior across platforms

      implicit none
      save

      integer,parameter::single_precision=selected_real_kind(6,30)
     & ,single_precision_size=4
      integer,parameter::double_precision=selected_real_kind(14,200)
     & ,double_precision_size=8

      integer,parameter::normal_int=selected_int_kind(5),int_size=4
      integer,parameter::long_int=selected_int_kind(10),long_int_size=8

      integer,parameter::default_string_length=128
      integer,parameter::default_path_length=256
      
      end module var_type
