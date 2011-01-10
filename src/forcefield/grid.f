      module grid
      implicit none
      save
      integer::ngrx,ngry,ngrz,maxp,nlayermax
      real(8),allocatable::egrid(:,:,:,:),xzz(:),yzz(:),zzz(:)
      real(8)::factx,facty,factz,dgrx,dgry,dgrz,eps=1d-5,dgr
      parameter (maxp=10)
!      common/gridd/factx,facty,factz,dgrx,dgry,dgrz,
!     &             xzz,yzz,zzz,egrid
!     &            ,ngrx,ngry,ngrz,eps
      end module grid
