       function ljmuir ( rijsq, ntij )

! ljmuir
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
 
!    *************************************************
!    ** calculates SAMI's 12+3 energy for headgroup **
!    **            and normal LJ energy for tail    **
!    *************************************************
 
      implicit none

! *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'poten.inc'

      real(8)::ljmuir, rijsq, sr, sr2, sr6, epshead, sighead
      integer::ntij

! --- attention: eps_hh / 4 used, since later multiplied by 4 --- 
!      parameter (epshead=27.67204d0,sighead=4.22d0)
      parameter (epshead=27.7204d0,sighead=6.5d0)

! --------------------------------------------------------------------

!       write(iou,*) 'sig2ij',sig2ij
!       write(iou,*) 'epsij',epsij

      if ( ntij .eq. 1 ) then
         sr = sighead / dsqrt( rijsq )
!       write(iou,*) 'sr',sr,'v',4.0d0*epshead*sr**3*(sr**9+1.0d0)
         ljmuir = epshead * sr**3 * ( sr**9 + 1.0d0 )
      else
         sr2 = sig2ij(ntij) / rijsq
         sr6 = sr2 * sr2 * sr2
         ljmuir = epsij(ntij) * sr6 * ( sr6 - 1.0d0)
!         if (ljmuir .gt. 100.0d0)
!     &       write(18,*) sig2ij(ntij),rijsq,'sr',dsqrt(sr2),'v',
!     &            4.0d0 * epsij(ntij) * sr6 * ( sr6 - 1.0d0)
      end if

      return

! ----------------------------------------------------------------------------

      end
