      module coord
      implicit none
      save

      logical,allocatable::lplace(:,:)
      logical::lelect,lflucq,lqtrans,lexpand,lpresim,lring,lrigid
     &  ,lrig,lrigi,lchiral,licell
      integer::nbox,nchain,nmolty,moltyp,nunit,nugrow,ntype,nchbox,ncmt
     &  ,ncmt2,iring,nrig,irig,frig,counthist,counttot,boxlink
      integer::nboxi, eetype, parbox,iurot,iupdatefix,nrigmin,nrigmax
      integer::parall,temtyp,riutry,rindex,prior,maxgrow,nrotbd,irotbd
      real(8)::rxu, ryu, rzu,masst, beta,xcm,ycm,zcm
     &  ,exp_c,eta,eta2,sxcm,sycm,szcm
      real(8)::pmrotbd
      real(8),allocatable::mass(:)
      real(8)::qqu,rcmu,rintramax,hist,probf,Elect_field
      character(len=18),allocatable::chname(:)
      character(len=4),allocatable::chemid(:)
      dimension lelect(ntmax),lflucq(ntmax),lqtrans(ntmax)
     &  ,lexpand(ntmax)
      dimension maxgrow(ntmax),probf(30,30,maxbin)
     &  ,hist(30,30,maxbin),iring(ntmax),nrig(ntmax),irig(ntmax,6)
     &  ,prior(ntmax,numax),frig(ntmax,6)
      dimension moltyp(nmax),nrigmin(ntmax),nrigmax(ntmax)
      dimension nunit(ntmax),nugrow(ntmax)
      dimension ntype(ntmax,numax)
      dimension nboxi(nmax),eetype(ntmax),rcmu(nmax)
      dimension iurot(ntmax)
      dimension nrotbd(ntmax),irotbd(numax,ntmax),pmrotbd(numax,ntmax)
      dimension rxu(nmax,numax), ryu(nmax,numax), rzu(nmax,numax)
      dimension qqu(nmax,numax), Elect_field(nbxmax)
      dimension nchbox(nbxmax),ncmt(nbxmax,ntmax)
     &  ,ncmt2(nbxmax,ntmax,20),eta(nbxmax,ntmax,20)
      dimension parall(ntmax,nmax),temtyp(ntmax)
      dimension masst(ntmax),parbox(nmax,nbxmax,ntmax)
      dimension xcm(nmax),ycm(nmax),zcm(nmax),exp_c(nmax)
     &  ,eta2(nbxmax,ntmax),sxcm(nmax),sycm(nmax),szcm(nmax)
      dimension riutry(ntmax,3),rindex(ntmax),lrigid(ntmax),lrig(ntmax)
     &  ,lrigi(ntmax,numax),lchiral(ntmax,numax),lring(ntmax)

!$$$      common /ncnunt/ moltyp,nunit,nugrow,ntype,parbox,boxlink
!$$$      common /ccord1/ probf, hist,rxu, ryu, rzu, qqu, rcmu, Elect_field
!$$$      common /ccord2/ nboxi,eetype,iurot,nchbox,ncmt,ncmt2
!$$$     &  ,parall,temtyp,licell,rintramax
!$$$      common /ccord3/ nbox,nchain,nmolty,nrigmin,nrigmax
!$$$      common /electa/ lelect,lflucq,lqtrans,lexpand,lring,lpresim
!$$$     &  ,lplace,lrigid,lrig,lrigi,lchiral
!$$$      common /masses/ mass, masst, xcm, ycm, zcm, exp_c, beta
!$$$     &  ,eta,eta2,sxcm,sycm,szcm
!$$$      common /histfix/ counthist,counttot,iupdatefix,maxgrow,iring
!$$$      common /rigidstuff/ riutry,rindex,nrig,irig,frig
!$$$      common /verbose/ chname, chemid
!$$$      common /ccord4/ pmrotbd
!$$$      common /ccord5/ irotbd,nrotbd

      end module coord