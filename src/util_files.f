module util_files
  use util_string, only: is_blank_line
  implicit none
  private
  public::get_iounit,readLine,readNthLine

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

! read the following line, skipping comment or blank lines
  subroutine readLine(iounit,line,skipComment,iostat)
    integer,intent(in)::iounit
    logical,intent(in)::skipComment
    character(len=*),intent(out)::line
    integer,intent(out)::iostat

    do
       line=''
       read(unit=iounit,FMT='(A)',iostat=iostat) line
       if  (iostat.eq.0) then
          if (.not.is_blank_line(line,skipComment)) exit
       else
          exit
       end if
    end do

  end subroutine readLine

! read the n-th line (i.e., skipping n-1 lines and read the following line), skipping
! comment or blank lines
  subroutine readNthLine(iounit,n,line,skipComment,iostat)
    integer,intent(in)::iounit,n
    logical,intent(in)::skipComment
    character(len=*),intent(out)::line
    integer,intent(out)::iostat
    integer::i

    do i=1,n
       call readLine(iounit,line,skipComment,iostat)
       if (iostat.ne.0) exit
    end do

  end subroutine readNthLine
end module util_files
