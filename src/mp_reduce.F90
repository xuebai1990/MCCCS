    INTEGER, INTENT(IN) :: msglen
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
#ifdef __MPI__
    group = MPI_COMM_WORLD
    IF( PRESENT( gid ) ) group = gid
    CALL MP_AUX_FUNC( msglen, msg, group, -1 )
#endif
