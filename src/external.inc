! EXTERNAL.INC
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

      real(8)::alpha1,alpha2,beta1,beta2
      real(8)::a1,delta,rsol
      integer::ntsubst
      integer::tau1,tau2
      logical::ljoe,lsami,lexzeo,lpsurf,lzgrid,lmuir,lgraphite
      logical::lslit
      logical::lcorreg,lelect_field
!  external electric field
      parameter (lelect_field=.false.)
! surfaces:
      parameter (ljoe=.false.,lsami=.false.,lmuir=.false.)
      parameter (lslit = .false.)
      parameter (lpsurf=.false.)
      parameter (lgraphite=.false.)
      parameter (lcorreg = .false.)
!      parameter (lcorreg = .false.)
! Graphite surface
! -- a1, delta are in Angstroms
! -- rsol [=] 1/A^3
      parameter (a1 = 2.460d0)
      parameter (delta = 3.40d0)
      parameter (ntsubst = 190)
      parameter (rsol = 0.114d0)
! zeolites: 
      parameter (lexzeo=.true.,lzgrid=.false.)
! Joe Hautman's parameters are read from standard input:
!               to be used with ljoe = .true. 
! Sami's parameters: to be used with lsami = .true. AND lmuir = .true.
      parameter ( alpha1=21.162d0, alpha2=-21.162d0,
     &              beta1=661.6d0, beta2=6616.0d0,
     &              tau1=-32, tau2=-16 )
! parameters for zeolites are read in subroutine SUZEO
! AT PRESENT no parameters for polymeric surfactants (see LJPSUR)
