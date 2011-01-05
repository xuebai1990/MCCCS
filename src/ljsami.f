       function ljsami ( rijsq, ntij )
 
!    *********************************************************
!    ** calculates SAMI's LJ and 12+3 energy for a bead.    **
!    *********************************************************
 
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'external.inc'
!$$$      include 'ljsamipara.inc'

      real(KIND=double_precision)::ljsami, rijsq, rij, sr
      integer(KIND=normal_int)::ntij

! --------------------------------------------------------------------

      rij = dsqrt( rijsq )
      sr = sij(ntij) / rij

      if ( ntij .eq. 1 ) then
! *** head-head interaction ( repulsive 12+3 interaction )
         ljsami = ( eij(1) * sr**3 * ( 1.0d0 + sr**9 ) ) - vsh(1) + ( rij * vsha(1) ) 
      else
! *** head-tail or tail-tail interaction ( LJ 12-6 interaction )
         ljsami = ( 4.0d0 * eij(ntij) * sr**6 * ( sr**6 - 1.0d0 ) ) - vsh(ntij) + ( rij * vsha(ntij) ) 
      end if

      return

! ----------------------------------------------------------------------------

      end
