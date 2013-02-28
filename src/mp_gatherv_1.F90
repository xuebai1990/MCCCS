    DATA_TYPE:: mydata(:)
    DATA_TYPE:: alldata(:)
    INTEGER, INTENT(IN) :: mycount, recvcount(:), displs(:)
    INTEGER, OPTIONAL, INTENT(IN) :: gid
#ifndef ALLGATHER
    INTEGER, INTENT(IN) :: root
#endif
    INTEGER :: group
    INTEGER :: ierr, npe, myid

#ifdef __MPI__
    group = MPI_COMM_WORLD
    IF( PRESENT( gid ) ) group = gid

    CALL MPI_comm_size( group, npe, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)

    CALL MPI_comm_rank( group, myid, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)

    IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop(__LINE__)

#ifndef ALLGATHER
    IF ( myid == root ) THEN
#endif

       IF ( SIZE( alldata ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop(__LINE__)

#ifndef ALLGATHER
    END IF
#endif

    IF ( SIZE( mydata ) < mycount ) CALL mp_stop(__LINE__)

#ifdef ALLGATHER
    CALL MPI_ALLGATHERV( mydata, mycount, MP_TYPE, &
     alldata, recvcount, displs, MP_TYPE, group, ierr )
#else
    CALL MPI_GATHERV( mydata, mycount, MP_TYPE, &
     alldata, recvcount, displs, MP_TYPE, root, group, ierr )
#endif

    IF (ierr/=0) CALL mp_stop(__LINE__)

#else
    IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop(__LINE__)
    IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop(__LINE__)
    alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
#endif
    RETURN
