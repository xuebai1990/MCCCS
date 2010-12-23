      function v_elect_field(i, j, rzfield,E)
! **********************************************************************
! **  calculates interaction of molecule i with an external field E  ***
! **  added 06/24/07 by KM                                           ***
! **********************************************************************
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'external.inc'
!$$$      include 'coord.inc'

      real(KIND=double_precision)::v_elect_field, convert, rzfield, E
      integer(KIND=normal_int)::i, j,ibox


! ********************************************
! **  units
! **  E in V/A, q in e, rz in A
! **  E*q*rz = V*e
! **  1 V*e = 11600 K
! ********************************************

      v_elect_field = -E*rzfield*qqu(i,j)  
      
!      write(6,*) 'E ', E, ' exfield ', exfield

      return
      end
