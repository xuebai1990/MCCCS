      subroutine cleanup(msg)
      implicit none
      include 'control.inc'
      include 'mpi.inc'
      character(len=*)::msg

      call MPI_FINALIZE(ierr)

      write(iou,*) msg

      stop

      end
