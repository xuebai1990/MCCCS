    INTEGER, INTENT(IN) :: mydata
    INTEGER, INTENT(OUT) :: alldata(:)
    INTEGER, OPTIONAL, INTENT(IN) :: gid
#ifndef ALLGATHER
    INTEGER, INTENT(IN) :: root
#endif
    INTEGER :: group, ierr

#ifdef __MPI__
    group = MPI_COMM_WORLD
    IF( PRESENT( gid ) ) group = gid

#ifdef ALLGATHER
    CALL MPI_ALLGATHER(mydata, 1, MPI_INTEGER, alldata, 1, MPI_INTEGER, group, IERR)
#else
    CALL MPI_GATHER(mydata, 1, MPI_INTEGER, alldata, 1, MPI_INTEGER, root, group, IERR)
#endif

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    alldata(1) = mydata
#endif
    RETURN
