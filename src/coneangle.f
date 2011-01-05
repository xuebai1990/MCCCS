      subroutine coneangle( thetaone, phione, thetatwo, phitwo, angle )

!     ******************************************************************
!     ** takes two unit vectors in spherical coordinates and computes **
!     ** the angle between them.                                      **
!     ** thetaone is the angle between vector one and the z-axis      **
!     ** phione is the angle between vector one and the x-axis        **
!     ** thetatwo is the angle between vector two and the z-axis      **
!     ** phitwo is the angle between vector two and the x-axis        **
!     ** angle is the angle between the two vectors                   **
!     ** x = r sin (theta) cos (phi)                                  **
!     ** y = r sin (theta) sin (phi)                                  **
!     ** z = r cos (theta)                                            **
!     ** M.G. Martin 2-4-98                                           **
!     ******************************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'control.inc'
      real(KIND=double_precision)::thetaone, thetatwo, phione, phitwo, angle,sintheone,costheone,sinthetwo,costhetwo,sinphione,cosphione ,sinphitwo,cosphitwo,cosangle

      sintheone = dsin(thetaone)
      costheone = dcos(thetaone)
      sinthetwo = dsin(thetatwo)
      costhetwo = dcos(thetatwo)
      sinphione = dsin(phione)
      cosphione = dcos(phione)
      sinphitwo = dsin(phitwo)
      cosphitwo = dcos(phitwo)
      

      cosangle = sintheone*cosphione*sinthetwo*cosphitwo + sintheone*sinphione*sinthetwo*sinphitwo + costheone*costhetwo

      angle = dacos(cosangle)

      return
      end



