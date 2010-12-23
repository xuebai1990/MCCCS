! NRTAB.INC

      integer::nrtcon,nrtmax,nrtbin
! - parameters currently set low to minimize storage
      parameter(nrtcon=50,nrtmax=10,nrtbin=10)

      real(8)::dmrtab(nrtmax),dbrtab(nrtmax),
     &                 wnrtab(nrtmax,nrtbin)

      common /nrtabc/ dmrtab,dbrtab,wnrtab
