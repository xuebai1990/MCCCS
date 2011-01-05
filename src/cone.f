      subroutine cone(iinit,x,y,z,alpha,gamma,ux,uy,uz)
!     ******************************************************************
!     * if iinit = 1 then it sets up the rotation matrix for the cone  *
!     * using x,y,z as a unit vector pointing in the +z direction      *
!     ******************************************************************
!     * if iinit = 2 then it creates a unit vector that has an angle   *
!     * of alpha from the +z direction (previous vector) and an angle  *
!     * of gamma (0,2Pi) around the cone circle and returns this as    *
!     * ux,uy,uz                                                       *
!     ******************************************************************
!     * if iinit = 3 then this computes gamma (0,2Pi) for a unit       *
!     * vector ux,uy,uz with angle alpha from the -z direction         *
!     ******************************************************************
!     * note that the rotation matrix is saved after each call so it   *
!     * needs to be reset when you wish to use another cone            *
!     *                                                                *
!     * originally written prior to 1995                               *
!     * last modified 02-12-2001 by M.G. Martin                        *
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
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'conver.inc'

!     --- variables passed to/from the subroutine
      integer(KIND=normal_int)::iinit
      real(KIND=double_precision)::x,y,z,alpha,gamma,ux,uy,uz

!     --- local variables
      real(KIND=double_precision)::a11,a12,a13,a21,a22,a31,a32,a33 ,sinthe,costhe,sinpsi,cospsi,singamma,cosgamma ,uxtemp,uytemp,uztemp,sinalph,cosalph
      real(KIND=double_precision)::determ,arccos
      save a11,a12,a13,a21,a22,a31,a32,a33,cosalph,sinalph

!      write(iou,*) 'start CONE'

      if ( iinit .eq. 1 ) then      
!       --- setup the unit cone
!       --- assume x,y,z is a rotation of vector (0,0,1)
!       --- conversions from xyz to theta,psi
!           x = -dsin(theta)*dcos(psi)
!           y = dsin(theta)*dsin(psi)
!           z = dcos(theta)
!       --- rotation matrix (note that sin(phi) = 0, cos(phi) = 1)
!       a11 = cos(psi) cos(theta) cos(phi) - sin(psi) sin(phi)
!           = cos(psi) cos(theta)
!       a12 = -sin(psi) cos(theta) cos(phi) - cos(psi) sin(phi)
!           = -sin(psi) cos(theta)
!       a13 = sin(theta) cos(phi) 
!           = sin(theta)
!       a21 = cos(psi) cos(theta) sin(phi) + sin (psi) cos(phi)
!           = sin(psi)
!       a22 = -sin(psi) cos(theta) sin(phi) + cos (psi) cos(phi)
!           = cos(psi)
!       a23 = sin(theta) sin(phi) 
!           = 0.0
!       a31 = -cos(psi) sin(theta)
!           = x
!       a32 = sin(psi)  sin (theta)
!           = y
!       a33 = cos(theta)
!           = z

        costhe = z
        sinthe = dsqrt(1.0d0 - costhe*costhe)
        if (sinthe .ne. 0.0d0) then
           sinpsi=(y)/sinthe
           cospsi=(-x)/sinthe
        else
           sinpsi = 0.0d0
           cospsi = 1.0d0
        end if

! rotation matrix
        a11 = cospsi*costhe
        a12 = -(sinpsi*costhe)
        a13 = sinthe
        a21 = sinpsi
        a22 = cospsi
        a31 = x
        a32 = y
        a33 = z

      elseif ( iinit .eq. 2 ) then
!        --- find the vector on a cone:
!        --- alpha is the angle with the cone axis (theta)
         sinalph = dsin(alpha)
         cosalph = dcos(alpha)
!        --- gamma is now passed to cone - (phi)
         uztemp = -cosalph
         uytemp = sinalph*dcos(gamma)
         uxtemp = sinalph*dsin(gamma)
         ux = a11*uxtemp + a21*uytemp + a31*uztemp
         uy = a12*uxtemp + a22*uytemp + a32*uztemp
         uz = a13*uxtemp +              a33*uztemp

      elseif ( iinit .eq. 3 ) then
!        --- find gamma for a given unit vector on the cone
!        --- use cramer's rule where a23 has already been replaced with 0
!        --- we don't need the uztemp so it is not computed

         determ =  ( a11*a22*a33 + a21*(a32*a13 - a12*a33) - a31*a22*a13 )

         uxtemp =   ( ux*a22*a33 + a21*(a32*uz - uy*a33) - a31*a22*uz ) / determ

         uytemp =  ( a11*(uy*a33 - a32*uz) + ux*(a32*a13 - a12*a33) + a31*(a12*uz - uy*a13) ) / determ

         sinalph = dsin(alpha)
         singamma = uxtemp/sinalph
         cosgamma = uytemp/sinalph

!        --- now need to find the gamma on [-Pi,Pi] that satisfies cos and sin
         gamma = arccos(cosgamma)
         if ( singamma .lt. 0.0d0 ) gamma = -gamma

      else
         write(iou,*) 'iinit ',iinit
         call cleanup('non valid iinit in cone.f')
      end if
      
!      write(iou,*) 'finish CONE'

      return
      end

      function arccos( value )
!     ******************************************************************
!     * computes the arc cosine of a value, with safety checks to make *
!     * sure that value is between -1.0 and 1.0.  If value is outside  *
!     * of that range, it is set to the nearest extreme of the range   *
!     *                                                                *
!     * originally written 02-12-2001 by M.G. Martin                   *
!     * last modified 02-12-2001 by M.G. Martin                        *
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
      include 'common.inc'

!$$$      include 'conver.inc'

!#include "constant.inc"      

      real(KIND=double_precision)::arccos,value

!      real(KIND=double_precision)::arccos,value,onepi
!      onepi = 3.141592653d0

      if ( value .gt. 1.0d0 ) then
         arccos = 0.0d0
      elseif ( value .lt. -1.0d0 ) then
         arccos = onepi
      else
         arccos = dacos(value)
      end if

      return
      end
