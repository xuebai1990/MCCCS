      subroutine close(iinit,rx,ry,rz,bondl,angle,lterm)

! close
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

!     **************************************************************
!     **    Takes three or four points and determines a point     ** 
!     **       that is a certain length from all of them          **
!     **            - ALL LENGTHS MUST BE THE SAME -              **
!     **************************************************************
!     **   This wonderful subroutine can also find a unit vector  **
!     **  connected to two other unit vectors with all the angles **
!     **                    between them given.                   **
!     **************************************************************
!     **    This mess was unfortunately created by Collin Wick    **
!     **            on December 1999, BUT IT DOES WORK            **
!     **************************************************************

      implicit none

      include 'control.inc'

      logical::lterm

      integer::n,iinit

      real(8)::x,y,z,rx,ry,rz,length,lengtha,lengthb
     &     ,xa,ya,za,theta,thetac,ux,uy,uz,bondl,avar,bvar,cvar
     &     ,rxa,rya,rza,lengthc,angle,rxf,ryf,rzf,var,dvar,a,b,c

      dimension rx(6),ry(6),rz(6),x(4),y(4),z(4),ux(3),uy(3),uz(3)
      dimension angle(3)

! ---------------------------------------------------------------------
!     Determines a point from three others with equal bond lengths


!     *******************************************************
!     if iinit=1 it finds two possibilities to close 3 beads
!     if iinit=2 it finds one possibility to close 4 beads
!     if iinit=3 it finds a vector connected two others given 
!     *******************************************************


!      write(iou,*) 'START CLOSE iinit=',iinit


      if (iinit.ne.3) then
    
         if (iinit.eq.2) then
            rxf = rx(4)
            ryf = ry(4)
            rzf = rz(4)
         end if


!     --- determine lengths from all coordinates
         lengtha = dsqrt( (rx(2)-rx(1))**2 + (ry(2)-ry(1))**2
     &        + (rz(2)-rz(1))**2)
         lengthb = dsqrt( (rx(2)-rx(3))**2 + (ry(2)-ry(3))**2
     &        + (rz(2)-rz(3))**2)
         lengthc = dsqrt( (rx(3)-rx(1))**2 + (ry(3)-ry(1))**2
     &        + (rz(3)-rz(1))**2)
      
!     --- detemine the angle at 1
         thetac = - 0.5d0 * ( lengthb**2 - lengtha**2 - lengthc**2 )
     &        / ( lengtha * lengthc )
         theta = dacos(thetac)

!     --- now lets determine the in-plane coordinates
         x(1) = 0
         y(1) = 0
         x(2) = lengtha
         y(2) = 0
         x(3) = lengthc * thetac
         y(3) = lengthc * dsin( theta )
         
!     --- determine the in-plane coordinates with the same length
         ya = 0.5d0 * ( ( x(2)**2 - x(3)**2 - y(3)**2 )
     &        / ( - y(3) ) + ( x(2)**2 * ( x(3) - x(2) ) ) 
     &        / ( ( - y(3) ) * ( x(2) ) ) )
         
         xa = 0.5d0 * x(2)
      
!     --- determine length to new position
         length = dsqrt( (x(1)-xa)**2 + (y(1)-ya)**2 )
         
         if (length.gt.bondl) then
!     --- the distances are too far to be able to close
            lterm = .true. 
            return
         end if
         
!     --- determine perpendicular distance from here to final point
         lengthc = dsqrt( bondl**2 - length**2 )
         
!     --- find vectors from 1 to 2, 2 to 3, and 3 to 1 in real::space
         x(1) = (rx(2) - rx(1)) 
         y(1) = (ry(2) - ry(1))
         z(1) = (rz(2) - rz(1))
         
         x(2) = (rx(3) - rx(2))
         y(2) = (ry(3) - ry(2))
         z(2) = (rz(3) - rz(2))
         
         x(3) = (rx(1) - rx(3))
         y(3) = (ry(1) - ry(3)) 
         z(3) = (rz(1) - rz(3)) 
         
!     --- find middle point from 1 to 2
         
         rx(4) = 0.5d0 * x(1) + rx(1)
         ry(4) = 0.5d0 * y(1) + ry(1)
         rz(4) = 0.5d0 * z(1) + rz(1)
         
!     --- find distance from this point to a
         
         length = dsqrt( length**2 - (0.5d0*lengtha)**2 )
      
!     --- cross 1 with 2 to find perpendicular vector
         
         ux(1) = y(1)*z(2) - z(1)*y(2)
         uy(1) = z(1)*x(2) - x(1)*z(2)
         uz(1) = x(1)*y(2) - y(1)*x(2)
         
!     --- cross this vector with 1 to find vector to final spot
         
         ux(2) = y(1)*uz(1) - z(1)*uy(1)
         uy(2) = z(1)*ux(1) - x(1)*uz(1)
         uz(2) = x(1)*uy(1) - y(1)*ux(1)
         
!     --- normalize this vector and give it the appropriate distance
         
         lengtha = dsqrt( ux(2)**2 + uy(2)**2 + uz(2)**2)
         
         ux(2) = length * ux(2) / lengtha
         uy(2) = length * uy(2) / lengtha
         uz(2) = length * uz(2) / lengtha
         
!     --- now determine the two possible points here
         
         rx(5) = rx(4) + ux(2)
         ry(5) = ry(4) + uy(2)
         rz(5) = rz(4) + uz(2)
         
         rx(6) = rx(4) - ux(2)
         ry(6) = ry(4) - uy(2)
         rz(6) = rz(4) - uz(2)
         
!     --- with the two possibilities, determine which is closer to 3
         
         lengtha = dsqrt((rx(3)-rx(5))**2 + (ry(3)-ry(5))**2
     &        + (rz(3) - rz(5))**2)
         lengthb = dsqrt((rx(3)-rx(6))**2 + (ry(3)-ry(6))**2
     &        + (rz(3) - rz(6))**2)
         
         if (lengtha.lt.lengthb) then
!     --- number 5 is right
            rxa = rx(5)
            rya = ry(5)
            rza = rz(5)
         else
!     --- number 6 is right
            rxa = rx(6)
            rya = ry(6)
            rza = rz(6)
         end if
         
!     --- find vector from 1 to a
      
         xa = rxa - rx(1)
         ya = rya - ry(1)
         za = rza - rz(1)
         
!     --- cross vector a with 1
         
         ux(1) = ya*z(1) - za*y(1)
         uy(1) = za*x(1) - xa*z(1) 
         uz(1) = xa*y(1) - ya*x(1)
         
!     --- normalize this and give it the length to our final point
         
         length = dsqrt(ux(1)**2 + uy(1)**2 + uz(1)**2)
         
         ux(1) = ux(1) * lengthc / length
         uy(1) = uy(1) * lengthc / length
         uz(1) = uz(1) * lengthc / length
         
!     --- finally lets find our final points and send them back

         rx(1) = rxa + ux(1)
         ry(1) = rya + uy(1)
         rz(1) = rza + uz(1)
         
         rx(2) = rxa - ux(1)
         ry(2) = rya - uy(1)
         rz(2) = rza - uz(1)

         if (iinit.eq.2) then
            
!     --- for four bead closes, only one of these will work
            
            x(1) = rxf - rx(1)
            y(1) = ryf - ry(1)
            z(1) = rzf - rz(1)
            lengtha = dsqrt(x(1)**2 + y(1)**2 + z(1)**2)
            
            x(2) = rxf - rx(2)
            y(2) = ryf - ry(2)
            z(2) = rzf - rz(2)
            lengthb = dsqrt(x(2)**2 + y(2)**2 + z(2)**2)

            write(iou,*) lengtha,lengthb,bondl

            
            if ((lengtha-bondl).lt.(lengthb-bondl)) then
               rx(1) = rx(1)
               ry(1) = ry(1)
               rz(1) = rz(1)
               if (abs(lengtha-bondl).gt.0.0001) then
                  lterm = .true.
                  return
               end if
            else
               if (abs(lengthb-bondl).gt.0.0001) then
                  lterm = .true.
                  return
               end if

               rx(1) = rx(2)
               ry(1) = ry(2)
               rz(1) = rz(2)
            end if
         end if

      else
!     --- lets find our open vector with angles from prev and closed
  
!     angle 1 is from vector 1 to the new
!     angle 2 is from vector 2 to the new
!     angle 3 is between vectors 1 and 2

         avar = dcos(angle(2)) - rx(2)*dcos(angle(1))/rx(1) 
         bvar = rx(2) * rz(1) / rx(1) - rz(2)
         cvar = ry(2) - rx(2) * ry(1) / rx(1)

         avar = avar / cvar
         bvar = bvar / cvar

         cvar = (ry(1) * bvar + rz(1)) / rx(1)
         dvar = (dcos(angle(1)) - ry(1) * avar) / rx(1)

         a = cvar**2 + bvar**2 + 1.0d0
         b = 2.0d0 * (avar*bvar - cvar*dvar)
         c = dvar**2 + avar**2 - 1.0d0

         var = b**2 - 4.0d0 * a * c

         if (var.lt.0) then
            lterm = .true.
            return
         end if

         rz(3) = (-b + dsqrt(var)) / (2.0d0 * a )
         rz(4) = (-b - dsqrt(var)) / (2.0d0 * a )

         rx(3) = dvar - rz(3) * cvar
         rx(4) = dvar - rz(4) * cvar
         
         ry(3) = rz(3) * bvar + avar
         ry(4) = rz(4) * bvar + avar

         return


         avar = dcos(angle(1))*(rx(2)*(rx(2)*rz(1) - rx(1)*rz(2))
     &        + ry(2)*(ry(2)*rz(1) - ry(1)*rz(2)))
     &        + dcos(angle(2))*(rx(1)*(rx(1)*rz(2) - rx(2)*rz(1))
     &        + ry(1)*(ry(1)*rz(2) - ry(2)*rz(1)))
      
         var = rx(2)**2*(ry(1)**2 + rz(1)**2)
     &        + ry(2)**2*(rx(1)**2 + rz(1)**2)
     &        + rz(2)**2*(rx(1)**2 + ry(1)**2)
     &        - 2.0d0*(rx(1)*rx(2)*ry(1)*ry(2)
     &        + rx(1)*rx(2)*rz(1)*rz(2)
     &        + ry(1)*ry(2)*rz(1)*rz(2))
     &        - (rx(2)**2 + ry(2)**2 + rz(2)**2)*dcos(angle(1))**2
     &        - (rx(1)**2 + ry(1)**2 + rz(1)**2)*dcos(angle(2))**2
     &        + 2.0d0*(rx(1)*rx(2) + ry(1)*ry(2) + rz(1)*rz(2))
     &        *dcos(angle(1))*dcos(angle(2))

         if (var.lt.0) then
            var = abs(var)
            lterm = .true.

!            return
         end if

         bvar = (rx(1)*ry(2) - rx(2)*ry(1))*dsqrt(var)

         cvar = rx(2)**2*(ry(1)**2 + rz(1)**2)
     &        + (ry(2)*rz(1) - ry(1)*rz(2))**2
     &        - 2.0d0*rx(1)*rx(2)*(ry(1)*ry(2) + rz(1)*rz(2))
     &        + rx(1)**2*(ry(2)**2 + rz(2)**2)

         rz(3) = (avar + bvar) / cvar
         rz(4) = (avar - bvar) / cvar

         avar = rx(2)*rz(1) - rx(1)*rz(2)
         bvar = rx(1)*dcos(angle(2)) - rx(2)*dcos(angle(1))
         cvar = rx(1)*ry(2) - rx(2)*ry(1)

         ry(3) = (rz(3)*avar + bvar) / cvar
         ry(4) = (rz(4)*avar + bvar) / cvar

         rx(3) = (dcos(angle(1)) - ry(1)*ry(3) - rz(1)*rz(3))
     &        / rx(1)
         rx(4) = (dcos(angle(1)) - ry(1)*ry(4) - rz(1)*rz(4))
     &        / rx(1)
                

      end if

!      write(iou,*) 'END CLOSE iinit=',iinit
      
!     ---------------------------------------------------------------


      return
      end
