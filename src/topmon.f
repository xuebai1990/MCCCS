      program topmon

      implicit none
      include 'common.inc'
! ----------------------------------------------------------------
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

! --- call main program
      call monola

      call MPI_FINALIZE(ierr)
! ----------------------------------------------------------------
      end program topmon
