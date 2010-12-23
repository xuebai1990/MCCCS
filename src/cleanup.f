      subroutine cleanup(msg)
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
!$$$      include 'mpi.inc'
      character(LEN=*)::msg

      call MPI_FINALIZE(ierr)

      write(iou,*) msg

      stop

      end
