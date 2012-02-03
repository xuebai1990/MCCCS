      subroutine coneangle( thetaone, phione, thetatwo, phitwo, angle )

c coneangle
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

c     ******************************************************************
c     ** takes two unit vectors in spherical coordinates and computes **
c     ** the angle between them.                                      **
c     ** thetaone is the angle between vector one and the z-axis      **
c     ** phione is the angle between vector one and the x-axis        **
c     ** thetatwo is the angle between vector two and the z-axis      **
c     ** phitwo is the angle between vector two and the x-axis        **
c     ** angle is the angle between the two vectors                   **
c     ** x = r sin (theta) cos (phi)                                  **
c     ** y = r sin (theta) sin (phi)                                  **
c     ** z = r cos (theta)                                            **
c     ** M.G. Martin 2-4-98                                           **
c     ******************************************************************

      implicit none
      double precision thetaone, thetatwo, phione, phitwo, angle
     &     ,sintheone,costheone,sinthetwo,costhetwo
     &     ,sinphione,cosphione,sinphitwo,cosphitwo,cosangle

      sintheone = dsin(thetaone)
      costheone = dcos(thetaone)
      sinthetwo = dsin(thetatwo)
      costhetwo = dcos(thetatwo)
      sinphione = dsin(phione)
      cosphione = dcos(phione)
      sinphitwo = dsin(phitwo)
      cosphitwo = dcos(phitwo)
      

      cosangle = sintheone*cosphione*sinthetwo*cosphitwo
     &     + sintheone*sinphione*sinthetwo*sinphitwo
     &     + costheone*costhetwo

      angle = dacos(cosangle)

      return
      end



