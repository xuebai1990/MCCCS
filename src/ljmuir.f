       function ljmuir ( rijsq, ntij )

c ljmuir
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
 
c    *************************************************
c    ** calculates SAMI's 12+3 energy for headgroup **
c    **            and normal LJ energy for tail    **
c    *************************************************
 
      implicit none

c *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'poten.inc'

      double precision ljmuir, rijsq, sr, sr2, sr6, epshead, sighead
      integer ntij

c --- attention: eps_hh / 4 used, since later multiplied by 4 --- 
c      parameter (epshead=27.67204d0,sighead=4.22d0)
      parameter (epshead=27.7204d0,sighead=6.5d0)

C --------------------------------------------------------------------

c       write(iou,*) 'sig2ij',sig2ij
c       write(iou,*) 'epsij',epsij

      if ( ntij .eq. 1 ) then
         sr = sighead / dsqrt( rijsq )
c       write(iou,*) 'sr',sr,'v',4.0d0*epshead*sr**3*(sr**9+1.0d0)
         ljmuir = epshead * sr**3 * ( sr**9 + 1.0d0 )
      else
         sr2 = sig2ij(ntij) / rijsq
         sr6 = sr2 * sr2 * sr2
         ljmuir = epsij(ntij) * sr6 * ( sr6 - 1.0d0)
c         if (ljmuir .gt. 100.0d0)
c     &       write(18,*) sig2ij(ntij),rijsq,'sr',dsqrt(sr2),'v',
c     &            4.0d0 * epsij(ntij) * sr6 * ( sr6 - 1.0d0)
      endif

      return

C ----------------------------------------------------------------------------

      end
