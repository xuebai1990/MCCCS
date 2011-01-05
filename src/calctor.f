      subroutine calctor(iu1,iu2,iu3,iu4,jttor,vtor)

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
!$$$      include 'fix.inc'

      integer(KIND=normal_int)::iu1,iu2,iu3,iu4,jttor
      
      real(KIND=double_precision)::thetac,xaa1,yaa1,zaa1,xa1a2,ya1a2 ,za1a2,daa1,da1a2,dot,vtor,vtorso,tcc,xcc,ycc,zcc,theta,spltor
      
!     --- calculate cross products d_a x d_a-1 
      xaa1 = yvec(iu2,iu1) * zvec(iu3,iu2) + zvec(iu2,iu1) * yvec(iu2,iu3)
      yaa1 = zvec(iu2,iu1) * xvec(iu3,iu2)  + xvec(iu2,iu1) * zvec(iu2,iu3)
      zaa1 = xvec(iu2,iu1) * yvec(iu3,iu2)  + yvec(iu2,iu1) * xvec(iu2,iu3)
      
!     --- calculate cross products d_a-1 x d_a-2
      xa1a2 = yvec(iu2,iu3) * zvec(iu3,iu4) - zvec(iu2,iu3) * yvec(iu3,iu4)
      ya1a2 = zvec(iu2,iu3) * xvec(iu3,iu4) - xvec(iu2,iu3) * zvec(iu3,iu4)
      za1a2 = xvec(iu2,iu3) * yvec(iu3,iu4) - yvec(iu2,iu3) * xvec(iu3,iu4)
   
!     --- calculate lengths of cross products ***
      daa1 = dsqrt ( xaa1**2 + yaa1**2 + zaa1**2 )
      da1a2 = dsqrt ( xa1a2**2 + ya1a2**2  + za1a2**2 )
      
! ----Addition for table look up for Torsion potential
!     --- calculate dot product of cross products ***
      dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
      thetac = - (dot / ( daa1 * da1a2 ))

      if (thetac.gt.1.0d0) thetac=1.0d0
      if (thetac.lt.-1.0d0) thetac=-1.0d0
!     KEA -- added for extending range to +/- 180 and additional defns of torsions
!     if torsion type is greater than 50, call spline program to use table of torsion
!     potentials and fit from these. Especially useful for asymmetric potentials

      if (jttor .ge. 50) then
!     *** calculate cross product of cross products ***
         xcc = yaa1*za1a2 - zaa1*ya1a2
         ycc = zaa1*xa1a2 - xaa1*za1a2
         zcc = xaa1*ya1a2 - yaa1*xa1a2
!     *** calculate scalar triple product ***
         tcc = xcc*xvec(iu2,iu3) + ycc*yvec(iu2,iu3) + zcc*zvec(iu2,iu3)
!     determine angle between -180 and 180, not 0 to 180
         theta = dacos(thetac)
         if (tcc .lt. 0.0d0) theta = -theta

         if (jttor.lt.100) then
            call splint(theta,spltor,jttor)
         else
            call lininter(theta,spltor,jttor)
         end if
         vtor = spltor

      else
         vtor = vtorso(thetac,jttor)
      end if
 
      return
      end
