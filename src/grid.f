      module grid

      implicit none
      save
      integer::ngrx,ngry,ngrz,maxp,nlayermax
      real(8),allocatable::egrid(:,:,:,:)
      real(8)::dgrx,dgry,dgrz,eps=1d-2,dgr
      parameter (maxp=10,dgr=0.2)
!      common/gridd/dgrx,dgry,dgrz,egrid,ngrx,ngry,ngrz,eps
      end module grid
