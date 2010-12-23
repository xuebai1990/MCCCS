! GAROFALINI.INC

      integer::pair_max
      parameter(pair_max=50000)

      real(8)::galpha(6),grho(6),gbeta(6),ga(6,3),
     &	 gb(6,3),gc(6,3),glambda(4),grij(4,2),ggamma(4,2),
     &   gtheta(4),grijsq(4,2)

      real(8)::dxij(pair_max),dyij(pair_max),
     &	 dzij(pair_max),dij(pair_max),dik(pair_max),dxik(pair_max)
     &   ,dyik(pair_max),dzik(pair_max)

      integer::itr1(pair_max),itr2(pair_max),itr3(pair_max),
     &   ntr,tagged

      real(8)::v3garo,garofalini

      common /garo/ ga,gb,gc,ggamma,grij,galpha,grho,gbeta,
     &              glambda,gtheta,grijsq,tagged

      common /Ddist/ dij,dxij,dyij,dzij,dik,dxik,dyik,dzik
      common /triads/ itr1,itr2,itr3,ntr

      common /garoenrgy/ v3garo
