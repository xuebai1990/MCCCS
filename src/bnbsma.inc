! BNBSMA.INC
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

      real(8)::bsvol,bnvol,bntrax,bntray,bntraz,bshmat
     &  ,bstrax,bstray,bstraz,bnrotx,bnroty,bnrotz,bnhmat
     &  ,bsrotx,bsroty,bsrotz,bnflcq,bsflcq,bnexpc,bsexpc
     &  ,bnflcq2,bsflcq2

       real(8)::Abntrax,Abntray,Abntraz,   
     &                  Abstrax,Abstray,Abstraz

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
      common /atomtra/ Abntrax,Abntray,Abntraz,Abstrax,
     &          Abstray,Abstraz  
