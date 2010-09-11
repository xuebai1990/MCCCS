      subroutine coneangle( thetaone, phione, thetatwo, phitwo, angle )

! coneangle
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
! John Stubbs, and Collin Wick and Ilja Siepmann  
!                     
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to 
!
! Free Software Foundation, Inc. 
! 59 Temple Place - Suite 330
! Boston, MA  02111-1307, USA.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

      implicit none
      include 'control.inc'
      real(8)::thetaone, thetatwo, phione, phitwo, angle
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



