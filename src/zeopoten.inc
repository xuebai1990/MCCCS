! ZEOPOTEN.INC

!      real(8)::zeps,zsig2,zrc2,zencut
      integer,parameter::zntype=2,gntypemax=32
      integer::gntype,gtable,ztype(zntype),znum(zntype)
      common/forces/gtable(gntypemax),znum,ztype,gntype
!     &     zeps(zntype,zntype),zsig2(zntype,zntype)
!     +             ,zrc2(zntype,zntype),zencut(zntype,zntype)
