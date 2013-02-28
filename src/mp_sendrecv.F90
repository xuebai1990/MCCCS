    INTEGER, INTENT(IN) :: dest, src, tag
    INTEGER, INTENT(INOUT) :: msglen
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: myid
#ifdef __MPI__
    INTEGER :: group
    INTEGER :: istatus(MPI_STATUS_SIZE)
    INTEGER :: ierr, nrcv

    group = MPI_COMM_WORLD
    IF( PRESENT( gid ) ) group = gid
    CALL MPI_comm_rank( group, myid, ierr )
    IF( ierr /= 0 ) CALL mp_stop(__LINE__)
#else
    myid = 0
#endif

    IF (src .NE. dest) THEN
#ifdef __MPI__
       IF(myid .EQ. src) THEN
          CALL MPI_SEND( msg_src, msglen, MP_TYPE, dest, tag, group, ierr)
          IF (ierr/=0) CALL mp_stop(__LINE__)
       ELSE IF(myid .EQ. dest) THEN
          CALL MPI_RECV( msg_dest, msglen, MP_TYPE, src, tag, group, istatus, ierr )
          IF (ierr/=0) CALL mp_stop(__LINE__)
          CALL MPI_GET_COUNT(istatus, MP_TYPE, nrcv, ierr)
          IF (ierr/=0) CALL mp_stop(__LINE__)
          msglen=nrcv
       ELSE
          ! processors not taking part in the communication have 0 length message
          msglen = 0
       END IF
#endif
    ELSE IF (myid .EQ. src)THEN
       msg_dest = msg_src
    END IF

#ifdef __MPI__
    CALL MPI_BARRIER(group, IERR)
    IF (ierr/=0) CALL mp_stop(__LINE__)
#endif

    RETURN
