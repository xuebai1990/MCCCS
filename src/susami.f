      subroutine susami
 
!    *******************************
!    ** Set-Up SAMI's potentials. **
!    *******************************
 
      use sim_system
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

      real(KIND=double_precision)::hsig,heps,tsig,teps

      parameter ( hsig=4.220d0, heps=110.68816d0,  tsig=3.527d0, teps=79.982210d0 )

      real(KIND=double_precision)::rcsami
      integer(KIND=normal_int)::ij

! --------------------------------------------------------------------

      rcsami = 2.5d0*tsig

      sij(1)=hsig
      sij(2)=0.5d0*(hsig+tsig)
      sij(3)=sij(2)
      sij(4)=sij(2)
      sij(5)=tsig
      sij(6)=tsig
      sij(7)=sij(2)
      sij(8)=tsig
      sij(9)=tsig

      eij(1)=heps
      eij(2)=dsqrt(heps*teps)
      eij(3)=eij(2)
      eij(4)=eij(2)
      eij(5)=teps
      eij(6)=teps
      eij(7)=eij(2)
      eij(8)=teps
      eij(9)=teps

      vsh(1)  = eij(1) * ( ( 13.0d0 * (sij(1)/rcsami)**12 ) + (  4.0d0 * (sij(1)/rcsami)**3  ) )
      vsha(1) = eij(1) * ( ( 12.0d0 * sij(1)**12 / rcsami**13 ) + (  3.0d0 * sij(1)**3  / rcsami**4  ) )

      do ij = 2, 9
         vsh(ij)  = 4.0d0 * eij(ij) *  ( ( 13.0d0 * (sij(ij)/rcsami)**12 ) - (  7.0d0 * (sij(ij)/rcsami)**6  ) )
         vsha(ij) = 4.0d0 * eij(ij) * ( ( 12.0d0 * sij(ij)**12 / rcsami**13 ) - (  6.0d0 * sij(ij)**6  / rcsami**7  ) )
      end do

!      do ij = 1,9
!         write(iou,*) 'ij',ij,'vsh',(vsh(ij)/80.0d0),
!     +                      'vsha',(vsha(ij)/80.0d0),
!     +                      'eij',(eij(ij)/80.0d0)
!      end do

      return

! ----------------------------------------------------------------------------

      end
