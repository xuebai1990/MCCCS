C EXTERNAL.INC
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

      double precision alpha1,alpha2,beta1,beta2
      double precision a1,delta,rsol
      integer ntsubst
      integer tau1,tau2
      logical ljoe,lsami,lexzeo,lpsurf,lzgrid,lmuir,lgraphite
      logical lslit
      logical lcorreg,lelect_field
c  external electric field
      parameter (lelect_field=.false.)
c surfaces:
      parameter (ljoe=.false.,lsami=.false.,lmuir=.false.)
      parameter (lslit = .false.)
      parameter (lpsurf=.false.)
      parameter (lgraphite=.false.)
      parameter (lcorreg = .false.)
c      parameter (lcorreg = .false.)
c Graphite surface
c -- a1, delta are in Angstroms
c -- rsol [=] 1/A^3
      parameter (a1 = 2.460d0)
      parameter (delta = 3.40d0)
      parameter (ntsubst = 190)
      parameter (rsol = 0.114d0)
c zeolites: 
      parameter (lexzeo=.false.,lzgrid=.false.)
c Joe Hautman's parameters are read from standard input:
c               to be used with ljoe = .true. 
c Sami's parameters: to be used with lsami = .true. AND lmuir = .true.
      parameter ( alpha1=21.162d0, alpha2=-21.162d0,
     &              beta1=661.6d0, beta2=6616.0d0,
     &              tau1=-32, tau2=-16 )
c parameters for zeolites are read in subroutine SUZEO
c AT PRESENT no parameters for polymeric surfactants (see LJPSUR)
