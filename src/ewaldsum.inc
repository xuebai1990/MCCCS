C EWALDSUM.INC
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

      logical::L_Ewald_Auto
      integer::kmax,numvect
c - kalp, a parameter to control the real::space sum
c - kmax, a parameter to control the total number of reciprocal vectors
c - numvect, the total number of reciprocal vectors
      dimension numvect(nbxmax)
      real(8)::kx,ky,kz,prefact,ssumr,ssumi,ssumrn,ssumin,
     &                 ssumro,ssumio
      real(8)::sself,correct,kalp,calp,dipolex,dipoley,dipolez,
     &                 dipolexo,dipoleyo,dipolezo
c - calp = kalp / boxlen
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
