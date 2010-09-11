      subroutine cleanup(msg)
      implicit none
      include 'mpi.inc'
      character(len=*)::msg

      call MPI_FINALIZE(ierr)

      print *,msg

      stop -1

      end
