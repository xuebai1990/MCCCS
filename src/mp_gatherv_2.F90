    DATA_TYPE :: mydata(:,:) !< \warning first dimension is supposed constant!
    DATA_TYPE :: alldata(:,:)
    INTEGER, INTENT(IN) :: recvcount(:), displs(:)
    INTEGER, INTENT(IN) :: comm
#ifndef ALLGATHER
    INTEGER, INTENT(IN) :: root
#endif
    INTEGER :: ierr, npe, myid
    INTEGER, ALLOCATABLE :: nrecv(:), ndisp(:)

#ifdef __MPI__

    CALL MPI_comm_size( comm, npe, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)

    CALL MPI_comm_rank( comm, myid, ierr )
    IF (ierr/=0) CALL mp_stop(__LINE__)

    IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop(__LINE__)

#ifndef ALLGATHER
    IF ( myid == root ) THEN
#endif

       IF ( SIZE( alldata, 2 ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop(__LINE__)
       IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop(__LINE__)

#ifndef ALLGATHER
    END IF
#endif

    IF ( SIZE( mydata, 2 ) < recvcount( myid + 1 ) ) CALL mp_stop(__LINE__)

    ALLOCATE( nrecv( npe ), ndisp( npe ) )

    nrecv( 1:npe ) = recvcount( 1:npe ) * SIZE( mydata, 1 )
    ndisp( 1:npe ) = displs( 1:npe ) * SIZE( mydata, 1 )

#ifdef ALLGATHER
    CALL MPI_ALLGATHERV( mydata, nrecv( myid + 1 ), MP_TYPE, &
     alldata, nrecv, ndisp, MP_TYPE, comm, ierr )
#else
    CALL MPI_GATHERV( mydata, nrecv( myid + 1 ), MP_TYPE, &
     alldata, nrecv, ndisp, MP_TYPE, root, comm, ierr )
#endif

    IF (ierr/=0) CALL mp_stop(__LINE__)

    DEALLOCATE( nrecv, ndisp )

#else
    IF ( SIZE( alldata, 1 ) /= SIZE( mydata, 1 ) ) CALL mp_stop(__LINE__)
    IF ( SIZE( alldata, 2 ) < recvcount( 1 ) ) CALL mp_stop(__LINE__)
    IF ( SIZE( mydata, 2  ) < recvcount( 1 ) ) CALL mp_stop(__LINE__)
    alldata( :, 1:recvcount( 1 ) ) = mydata( :, 1:recvcount( 1 ) )
#endif
    RETURN
