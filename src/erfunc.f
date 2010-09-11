      function erfunc(x)
      implicit none

! erfunc
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

      real(8)::p,a1,a2,a3,a4,a5,x,erfunc,tt,eee
      parameter (p=0.3275911d0,a1=0.254829592d0,
     &     a2=-0.284496736d0,a3=1.421413741d0,
     &     a4=-1.453152027d0,a5=1.061405429d0)

      eee = dexp(-x*x)
      tt = 1.0d0/(1.0d0 + p*x)
      erfunc = ((((a5*tt+a4)*tt+a3)*tt+a2)*tt+a1)*tt*eee
      return
      end
