      subroutine calctor(iu1,iu2,iu3,iu4,jttor,vtor)

c calctor
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
c John Stubbs, and Collin Wick and Ilja Siepmann  
c                     
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to 
c
c Free Software Foundation, Inc. 
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'control.inc'
      include 'fix.inc'

      integer iu1,iu2,iu3,iu4,jttor
      
      double precision thetac,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2
     &     ,daa1,da1a2,dot,vtor,vtorso,tcc,xcc,ycc,zcc,theta,spltor
      
c     --- calculate cross products d_a x d_a-1 
      xaa1 = yvec(iu2,iu1) * zvec(iu3,iu2)
     &     + zvec(iu2,iu1) * yvec(iu2,iu3)
      yaa1 = zvec(iu2,iu1) * xvec(iu3,iu2) 
     &     + xvec(iu2,iu1) * zvec(iu2,iu3)
      zaa1 = xvec(iu2,iu1) * yvec(iu3,iu2) 
     &     + yvec(iu2,iu1) * xvec(iu2,iu3)
      
c     --- calculate cross products d_a-1 x d_a-2
      xa1a2 = yvec(iu2,iu3) * zvec(iu3,iu4) -
     +     zvec(iu2,iu3) * yvec(iu3,iu4)
      ya1a2 = zvec(iu2,iu3) * xvec(iu3,iu4) -
     +     xvec(iu2,iu3) * zvec(iu3,iu4)
      za1a2 = xvec(iu2,iu3) * yvec(iu3,iu4) -
     +     yvec(iu2,iu3) * xvec(iu3,iu4)
   
c     --- calculate lengths of cross products ***
      daa1 = dsqrt ( xaa1**2 + yaa1**2 + zaa1**2 )
      da1a2 = dsqrt ( xa1a2**2 + ya1a2**2 
     &     + za1a2**2 )
      
c ----Addition for table look up for Torsion potential
c     --- calculate dot product of cross products ***
      dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
      thetac = - (dot / ( daa1 * da1a2 ))

      if (thetac.gt.1.0d0) thetac=1.0d0
      if (thetac.lt.-1.0d0) thetac=-1.0d0
c     KEA -- added for extending range to +/- 180 and additional defns of torsions
c     if torsion type is greater than 50, call spline program to use table of torsion
c     potentials and fit from these. Especially useful for asymmetric potentials

      if (jttor .ge. 50) then
c     *** calculate cross product of cross products ***
         xcc = yaa1*za1a2 - zaa1*ya1a2
         ycc = zaa1*xa1a2 - xaa1*za1a2
         zcc = xaa1*ya1a2 - yaa1*xa1a2
c     *** calculate scalar triple product ***
         tcc = xcc*xvec(iu2,iu3) + ycc*yvec(iu2,iu3)
     &        + zcc*zvec(iu2,iu3)
c     determine angle between -180 and 180, not 0 to 180
         theta = dacos(thetac)
         if (tcc .lt. 0.0d0) theta = -theta

         if (jttor.lt.100) then
            call splint(theta,spltor,jttor)
         else
            call lininter(theta,spltor,jttor)
         endif
         vtor = spltor

      else
         vtor = vtorso(thetac,jttor)
      endif
 
      return
      end
