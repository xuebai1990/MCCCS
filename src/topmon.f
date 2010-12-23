      program topmon

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'mpi.inc'
!$$$      include 'mpif.h'
! ----------------------------------------------------------------
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

! --- call main program
      call monola

      call MPI_FINALIZE(ierr)
! ----------------------------------------------------------------
      end program topmon
