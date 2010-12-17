! ZEOLITE.INC

      integer::idzeo,nzeomax,nzeo,nx,ny,nz
      real(8)::zeox,zeoy,zeoz,zeorx,zeory,zeorz,zunitx,zunity,zunitz
      logical::lunitcell

      parameter       (nzeomax=57600)

      common/zeolit/zeox(nzeomax),zeoy(nzeomax),zeoz(nzeomax),
     &              idzeo(nzeomax),lunitcell(nzeomax),nzeo
      common/box  /zeorx,zeory,zeorz,zunitx,zunity,zunitz,nx,ny,nz

