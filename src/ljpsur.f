       function ljpsur ( rijsq, ntij )
 
c ljpsur
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

c    *********************************************************
c    ** calculates energy for a polymeric surfactant bead.  **
c    *********************************************************
 
      implicit none

c *** common blocks ***
      include 'external.inc'

      real(8)::ljpsur, rijsq, sr, sr6
      integer::ntij

C --------------------------------------------------------------------
c AT PRESENT: all sigma = 1.0
c             epsilon   = 1.0
C --------------------------------------------------------------------

      sr = 1.0d0 / rijsq
      sr6 = sr**3

      if ( ntij .eq. 1 ) then
c *** nonpolar-nonpolar interaction ( LJ interaction )
         ljpsur = sr6*sr6 - sr6
      else
c *** polar-polar or polar-nonpolar interaction ( repulsive LJ interaction )
         if ( rijsq .le. 1.259921d0 ) then
            ljpsur = sr6*sr6 - sr6 + 0.25d0
         else
            ljpsur = 0.0d0
         endif
      endif

      return

C ----------------------------------------------------------------------------

      end
