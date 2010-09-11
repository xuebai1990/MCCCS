      module grid
! GRID.INC
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
      save
      integer::ngrx,ngry,ngrz,maxl,maxp
      real(8),allocatable::egrid(:,:,:,:),xzz(:),yzz(:),zzz(:)
      real(8)::factx,facty,factz,dgrx,dgry,dgrz,eps=1d-5,dgr
! - parameters currently set low since zeolites not used
      parameter (maxl=20,maxp=10)
!      common/gridd/factx,facty,factz,dgrx,dgry,dgrz,
!     &             xzz,yzz,zzz,egrid
!     &            ,ngrx,ngry,ngrz,eps
      end module grid
