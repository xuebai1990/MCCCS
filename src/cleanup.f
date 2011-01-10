      subroutine cleanup(msg)
      use global_data,only:iou,ierr
      use util_timings
      implicit none
      include 'common.inc'
!$$$      include 'control.inc'
!$$$      include 'mpi.inc'
      character(LEN=*)::msg

      call MPI_FINALIZE(ierr)

      write(iou,*) msg

      stop

      end
