! EWALDSUM.INC
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

      logical::L_Ewald_Auto
      integer::kmax,numvect
! - kalp, a parameter to control the real::space sum
! - kmax, a parameter to control the total number of reciprocal vectors
! - numvect, the total number of reciprocal vectors
      dimension numvect(nbxmax)
      real(8)::kx,ky,kz,prefact,ssumr,ssumi,ssumrn,ssumin,
     &                 ssumro,ssumio
      real(8)::sself,correct,kalp,calp,dipolex,dipoley,dipolez,
     &                 dipolexo,dipoleyo,dipolezo
! - calp = kalp / boxlen
      dimension kx(vectormax,nbxmax),ky(vectormax,nbxmax),
     &          kz(vectormax,nbxmax),prefact(vectormax,nbxmax),
     &          ssumr(vectormax,nbxmax),ssumi(vectormax,nbxmax),
     &          ssumrn(vectormax,nbxmax),ssumin(vectormax,nbxmax),
     &          calp(nbxmax),kmax(nbxmax),kalp(nbxmax),
     &          ssumro(vectormax,nbxmax),ssumio(vectormax,nbxmax)
      dimension dipolex(nbxmax),dipoley(nbxmax),dipolez(nbxmax)

      common /ewaldl/ L_Ewald_Auto 
      common /numvect/ numvect,kmax 
      common /vector/ kx,ky,kz
      common /ssum/ calp,kalp,ssumr,ssumi,prefact,ssumrn,
     &                  ssumin,ssumro,ssumio
      common /selfterm/ sself,correct
      common /dipoterm/ dipolex,dipoley,dipolez,dipolexo,dipoleyo,
     &                  dipolezo
