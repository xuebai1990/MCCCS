! NSIX.INC

! ***************************************************
! * stores interaction parameters for 9-6 potential *
! ***************************************************

      integer::nxatom
      parameter (nxatom=5)
      real(8)::rzero(nxatom*nxatom),epsnx(nxatom*nxatom)
      real(8)::shiftnsix(nxatom*nxatom)

      common /formic/ rzero,epsnx,shiftnsix
