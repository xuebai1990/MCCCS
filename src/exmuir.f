       function exmuir ( z, ntj )

c exmuir
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
c    ** calculates the lmuir external energy for a bead.    **
c    *********************************************************
 
      implicit none

c *** common blocks ***
      include 'external.inc'
      include 'externalmuir.inc'

      real(8)::exmuir, z
      integer::ntj

C --------------------------------------------------------------------

      if ( ntj .eq. 1 ) then
c --- HEADgroup potential ---
         if ( z .le. alpha2 ) then
            exmuir = 0.0d0
         else
            exmuir = beta2 / ( 1.0d0 + ( (z/alpha2) - 1.0d0 )**tau2 )
         endif
      else
c --- TAILgroup potential ---
         if ( z .ge. alpha1 ) then
            exmuir = 0.0d0
         else
            if ( z .lt. zprmin ) then
               if ( ntj .eq. 2 ) then
                  exmuir = betac2 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) +
     +                     v2prmin
               else
                  exmuir = betac3 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) +
     +                     v3prmin
               endif
            else
               if ( ntj .eq. 2 ) then
                  exmuir = betac2 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) +
     +                     c9ch2 / z**9 -  c3ch2 / z**3
               else
                  exmuir = betac3 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) +
     +            c9ch3 / z**9 -  c3ch3 / z**3
               endif
            endif
         endif
      endif
      return

C ----------------------------------------------------------------------------

      end
