module util_files
  use util_string, only: str_trim
  implicit none
  private
  public::get_iounit

  INTEGER,PARAMETER::max_unit_number=999,reserved_unit_numbers(2)=(/5,6/)
  CHARACTER(LEN=1),PARAMETER::commentChar(3)=(/"#","!","%"/)
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

! *****************************************************************************
!> \brief returns .true. if the line is a comment line or an empty line
!> \par History
!>      03.2009 [tlaino] - Teodoro Laino
! *****************************************************************************
  LOGICAL FUNCTION is_comment_line(line)
    CHARACTER(LEN=*), INTENT(IN)             :: line
    INTEGER                                  :: ia,ib

    is_comment_line=.false.
    call str_trim(line,ia,ib)
    if (any(commentChar==line(ia:ia))) is_comment_line=.true.

  END FUNCTION is_comment_line

! read the following line, skipping comment or blank lines
  subroutine readLine(iounit,iostat,line)
    integer,intent(in)::iounit
    character(len=*),intent(out)::line
    integer,intent(out)::iostat

    do
       read(unit=iounit,FMT='(A)',iostat=iostat) line
       if  (iostat.eq.0) then
          if (.not.is_comment_line(line)) exit
       else
          exit
       end if
    end do

  end subroutine readLine

! read the n-th line (i.e., skipping n-1 lines and read the following line), skipping
! comment or blank lines
  subroutine readNthLine(iounit,n,iostat,line)
    integer,intent(in)::iounit,n
    character(len=*),intent(out)::line
    integer,intent(out)::iostat
    integer::i

    do i=1,n
       call readLine(iounit,iostat,line)
       if (iostat.ne.0) exit
    end do

  end subroutine readNthLine
end module util_files
