    INTEGER, INTENT(IN) :: mydata
    INTEGER, INTENT(OUT) :: alldata(:)
    INTEGER, INTENT(IN) :: comm
#ifndef ALLGATHER
    INTEGER, INTENT(IN) :: root
#endif
    INTEGER :: ierr

#ifdef __MPI__

#ifdef ALLGATHER
    CALL MPI_ALLGATHER(mydata, 1, MPI_INTEGER, alldata, 1, MPI_INTEGER, comm, IERR)
#else
    CALL MPI_GATHER(mydata, 1, MPI_INTEGER, alldata, 1, MPI_INTEGER, root, comm, IERR)
#endif

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    alldata(1) = mydata
#endif
    RETURN
