      module util_files
      implicit none

! *****************************************************************************
!> \brief returns the first fortran unit that is not preconnected
!> \note
!>       -1 if no free unit exists
!> taken from CP2K
! *****************************************************************************
  FUNCTION get_unit_number() RESULT(unit_number)
    INTEGER                                  :: unit_number

    INTEGER                                  :: istat
    LOGICAL                                  :: exists, opened

    DO unit_number=1,max_unit_number
       IF (ANY(unit_number == reserved_unit_numbers)) CYCLE
       INQUIRE (UNIT=unit_number,EXIST=exists,OPENED=opened,IOSTAT=istat)
       IF (exists.AND.(.NOT.opened).AND.(istat == 0)) RETURN
    END DO

    unit_number = -1

  END FUNCTION get_unit_number
      end module util_files
