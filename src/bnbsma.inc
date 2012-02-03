C BNBSMA.INC
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

      double precision bsvol,bnvol,bntrax,bntray,bntraz,bshmat
     &  ,bstrax,bstray,bstraz,bnrotx,bnroty,bnrotz,bnhmat
     &  ,bsrotx,bsroty,bsrotz,bnflcq,bsflcq,bnexpc,bsexpc
     &  ,bnflcq2,bsflcq2

      dimension bsvol(nbxmax), bnvol(nbxmax),bshmat(nbxmax,9),
     &	        bnhmat(nbxmax,9)
      dimension bntrax(ntmax,nbxmax),bntray(ntmax,nbxmax),
     &          bntraz(ntmax,nbxmax)
      dimension bstrax(ntmax,nbxmax),bstray(ntmax,nbxmax),
     &          bstraz(ntmax,nbxmax)
      dimension bnrotx(ntmax,nbxmax),bnroty(ntmax,nbxmax),
     &          bnrotz(ntmax,nbxmax)
      dimension bsrotx(ntmax,nbxmax),bsroty(ntmax,nbxmax),
     &          bsrotz(ntmax,nbxmax)
      dimension bnexpc(ntmax,nbxmax),bsexpc(ntmax,nbxmax)
      dimension bnflcq(ntmax,nbxmax),bsflcq(ntmax,nbxmax),
     &          bnflcq2(ntmax,nbxmax),bsflcq2(ntmax,nbxmax)

      common /bvolum/ bnvol,bsvol,bshmat,bnhmat
      common /btrans/ bntrax,bntray,bntraz,bstrax,bstray,bstraz
      common /brotat/ bnrotx,bnroty,bnrotz,bsrotx,bsroty,bsrotz
      common /bflucq/ bnflcq,bsflcq,bnexpc,bsexpc,bnflcq2,bsflcq2
