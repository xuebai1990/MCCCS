subroutine cleanup(msg)
  use global_data,only:iou,ierr
  implicit none
  include 'common.inc'
  !$$$      include 'control.inc'
  !$$$      include 'mpi.inc'
  character(LEN=*),intent(in)::msg

  call MPI_FINALIZE(ierr)

  write(iou,*) msg

  stop

end subroutine cleanup
