    DATA_TYPE, INTENT(IN) :: mydata(:)
    DATA_TYPE, INTENT(OUT) :: alldata(:)
    INTEGER, OPTIONAL, INTENT(IN) :: gid
#ifndef ALLGATHER
    INTEGER, INTENT(IN)::root
#endif
    INTEGER :: group
    INTEGER :: msglen, ierr

    msglen = SIZE(mydata)
    !IF( msglen .NE. SIZE(alldata, 1) ) CALL mp_stop(__LINE__)

#ifdef __MPI__
    group = MPI_COMM_WORLD
    IF( PRESENT( gid ) ) group = gid

#ifdef ALLGATHER
    CALL MPI_ALLGATHER(mydata, msglen, MP_TYPE, alldata, msglen, MP_TYPE, group, IERR)
#else
    CALL MPI_GATHER(mydata, msglen, MP_TYPE, alldata, msglen, MP_TYPE, root, group, IERR)
#endif

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    alldata(1:msglen) = mydata
#endif
    RETURN
