       function ljsami ( rijsq, ntij )
 
c ljsami
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
c    ** calculates SAMI's LJ and 12+3 energy for a bead.    **
c    *********************************************************
 
      implicit none

c *** common blocks ***
      include 'external.inc'
      include 'ljsamipara.inc'

      real(8)::ljsami, rijsq, rij, sr
      integer::ntij

C --------------------------------------------------------------------

      rij = dsqrt( rijsq )
      sr = sij(ntij) / rij

      if ( ntij .eq. 1 ) then
c *** head-head interaction ( repulsive 12+3 interaction )
         ljsami = ( eij(1) * sr**3 * ( 1.0d0 + sr**9 ) )
     &        - vsh(1) + ( rij * vsha(1) ) 
      else
c *** head-tail or tail-tail interaction ( LJ 12-6 interaction )
         ljsami = ( 4.0d0 * eij(ntij) * sr**6 * ( sr**6 - 1.0d0 ) )
     &        - vsh(ntij) + ( rij * vsha(ntij) ) 
      endif

      return

C ----------------------------------------------------------------------------

      end
