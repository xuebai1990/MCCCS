module util_files
  implicit none
  private
  public::get_iounit

  INTEGER,PARAMETER::max_unit_number=999,reserved_unit_numbers(2)=(/5,6/)

CONTAINS

! *****************************************************************************
!> \brief returns the first fortran unit that is not preconnected
!> \note
!>       -1 if no free unit exists
!> taken from CP2K
! *****************************************************************************
  INTEGER FUNCTION get_iounit()
    INTEGER                                  :: istat
    LOGICAL                                  :: lexist, lopen

    DO get_iounit=1,max_unit_number
       IF (ANY(get_iounit == reserved_unit_numbers)) CYCLE
       INQUIRE(UNIT=get_iounit,EXIST=lexist,OPENED=lopen,IOSTAT=istat)
       IF (lexist.AND.(.NOT.lopen).AND.(istat == 0)) RETURN
    END DO

    get_iounit = -1

  END FUNCTION get_iounit
end module util_files
