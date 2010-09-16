! SWTCMOVE.INC
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

      integer::nswaty,nswatb,gswatc,nswach
      integer::nsampos,splist,nswtcb,box3,box4,ncut

      real(8)::bnswat,bniswat,pmswtcb
      real(8)::bnswat_empty	

      dimension nswach(npamax),nsampos(npamax)
      dimension nswatb(npamax,npabmax)
      dimension bnswat(npamax,npabmax),ncut(npamax,2)
      dimension gswatc(npamax,2,2*npamax),splist(npamax,numax,2)
      dimension nswtcb(npamax)
      dimension pmswtcb(npamax,npabmax)
      dimension box3(npamax,npabmax)
      dimension box4(npamax,npabmax)
      dimension bnswat_empty(npamax,npabmax)

      common /swtcma/ nswaty,nswtcb
      common /swtcmb/ nswach,nsampos
      common /swtcmc/ box3,box4,bnswat,bnswat_empty,nswatb
      common /swtcmd/ gswatc,splist,ncut
      common /swtcme/ pmswtcb
