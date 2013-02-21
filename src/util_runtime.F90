module util_runtime
  implicit none
  include 'common.inc'
  private
  public::err_exit
contains
  subroutine err_exit(msg)
    character(LEN=*),intent(in)::msg
    integer::ierr

    write(*,*) msg

    call MPI_ABORT(ierr)

    stop

  end subroutine err_exit
end module util_runtime
