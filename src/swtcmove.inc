c SWTCMOVE.INC
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
      common /swtcmc/ nswatb,bnswat,box3,box4,bnswat_empty
      common /swtcmd/ gswatc,splist,ncut
      common /swtcme/ pmswtcb
