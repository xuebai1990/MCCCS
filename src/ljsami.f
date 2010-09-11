       function ljsami ( rijsq, ntij )
 
! ljsami
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
 
!    *********************************************************
!    ** calculates SAMI's LJ and 12+3 energy for a bead.    **
!    *********************************************************
 
      implicit none

! *** common blocks ***
      include 'external.inc'
      include 'ljsamipara.inc'

      real(8)::ljsami, rijsq, rij, sr
      integer::ntij

! --------------------------------------------------------------------

      rij = dsqrt( rijsq )
      sr = sij(ntij) / rij

      if ( ntij .eq. 1 ) then
! *** head-head interaction ( repulsive 12+3 interaction )
         ljsami = ( eij(1) * sr**3 * ( 1.0d0 + sr**9 ) )
     &        - vsh(1) + ( rij * vsha(1) ) 
      else
! *** head-tail or tail-tail interaction ( LJ 12-6 interaction )
         ljsami = ( 4.0d0 * eij(ntij) * sr**6 * ( sr**6 - 1.0d0 ) )
     &        - vsh(ntij) + ( rij * vsha(ntij) ) 
      end if

      return

! ----------------------------------------------------------------------------

      end
