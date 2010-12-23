! BLKAVG.INC

      integer::nener,nprop,blockm
      real(8)::naccu,nccold,accum,bccold,aver,baver
     &     ,acsvol,acnvol,acshmat,acnhmat
! ** Neeraj adding for solubility parameter and heat of vaporization

      real(8)::naccu1,nccold1,accum1,bccold1,aver1,baver1
      integer::nprop1
      parameter (nener=18,nprop=4+nener+(4*ntmax)+3,blockm=100)
      parameter (nprop1=11)
      dimension naccu(nprop,nbxmax),nccold(nprop,nbxmax)
      dimension accum(nprop,nbxmax),bccold(nprop,nbxmax)
      dimension aver(nprop,nbxmax), baver(nprop,nbxmax,blockm)
      dimension acsvol(nbxmax),acnvol(nbxmax)
      dimension acshmat(nbxmax,9),acnhmat(nbxmax,9)
      
      dimension naccu1(nprop1,nbxmax,nbxmax), 
     &  nccold1(nprop1,nbxmax,nbxmax), aver1(nprop1,nbxmax,nbxmax)
      dimension accum1(nprop1,nbxmax,nbxmax),
     &  bccold1(nprop1,nbxmax,nbxmax), 
     &  baver1(nprop1,nbxmax,nbxmax,blockm)

      double precision, dimension (ntmax,nbxmax) :: acntray,acntrax,
     &     acntraz,acnrotx,acnroty,acnrotz,acstrax,acstray,acstraz,
     &     acsrotx,acsroty,acsrotz

      common /bokavg/ baver,naccu,nccold,accum,bccold,aver,acsvol
     &                ,acnvol,acshmat,acnhmat
      common /blkavg/ baver1,naccu1,nccold1,accum1,bccold1,aver1

      common /trarotavg/ acntrax,acntray,acntraz,acnrotx,acnroty,
     &     acnrotz,acstrax,acstray,acstraz,acsrotx,acsroty,
     &     acsrotz
