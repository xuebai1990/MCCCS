    INTEGER, INTENT(IN) :: source, msglen
    INTEGER, OPTIONAL, INTENT(IN) :: gid
    INTEGER :: group
#ifdef __MPI__
    group = MPI_COMM_WORLD
    IF( PRESENT( gid ) ) group = gid
    CALL MP_AUX_FUNC( msg, msglen, source, group )
#endif
