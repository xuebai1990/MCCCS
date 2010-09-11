! MERCK.INC
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

      integer::natomtyp
      parameter (natomtyp = 3)
      real(8)::epsimmff,sigimmff,smmff,sigisq,ammff,nmmff
      real(8)::gmmff,alphammff,coru_cons,corp_cons 
      dimension epsimmff(natomtyp),sigimmff(natomtyp),smmff(natomtyp)
      dimension ammff(natomtyp),nmmff(natomtyp),gmmff(natomtyp)
      dimension sigisq(natomtyp),alphammff(natomtyp)
      dimension coru_cons(natomtyp),corp_cons(natomtyp) 
      common /merk/ epsimmff,sigimmff,smmff,sigisq,ammff,nmmff
     &   ,alphammff,gmmff,coru_cons,corp_cons 

 
