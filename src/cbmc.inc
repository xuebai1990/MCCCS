! CBMC.INC

      logical::lovr,lexist,lexshed,llplace,lpnow,llrig,lsave,lrplc

! *** adding this for iswatch move ***
      logical::liswinc, liswatch
      dimension liswinc(numax,ntmax),llplace(ntmax)
      integer::other,iparty,iplace, pfrom, pnum, pprev,nplace
! *** end add ***

      real(8)::bncb,bscb,rxp,ryp,rzp,rxnew,rynew,rznew
     &   ,bnregr,bsregr,bfac,bsswat,bsiswat,fbncb,fbscb
     &   ,vtry,vtrext,vtrintra,vtrinter,vtrelect,vtrewald
     &   ,vtrorient,vtrelect_intra,vtrelect_inter
      integer::nchoi1,nchoi,icbdir,icbsta,nmaxcbmc,nchoih,nchtor
     &  ,nchbna,nchbnb,nrigi,nchoir,rlist,rfrom,rprev,rnum
      real(8)::acchem, bnchem
      dimension lexist(numax),lexshed(numax)
      common /pbias/  bncb(ntmax,numax), bscb(ntmax,2,numax)
     &   ,bnregr(ntmax),bsregr(ntmax,2)
     &   ,bsswat(npamax,npabmax),bsiswat(npamax,npabmax)
     &   ,fbncb(ntmax,numax), fbscb(ntmax,2,numax)
      common /cbpara/ nchoi1(ntmax),nchoi(ntmax),icbdir(ntmax)
     &   ,icbsta(ntmax),nmaxcbmc(ntmax),nchoih(ntmax)
     &   ,nchtor(ntmax),nchbna(ntmax),nchbnb(ntmax)
     &   ,nchoir(ntmax),lrplc(ntmax)
      common /bbias/  bfac(nchmax), vtry(nchmax),vtrinter(nchmax)
     &   ,vtrext(nchmax),vtrintra(nchmax)
     &   ,vtrelect(nchmax),vtrewald(nchmax),lovr(nchmax)
     &   ,vtrorient(nchmax),vtrelect_intra(nchmax),
     &    vtrelect_inter(nchmax)
      common /xyzp/   lexist,lexshed,rxp(numax,nchmax)
     &   ,ryp(numax,nchmax), rzp(numax,nchmax)
     &   ,liswinc,liswatch,iparty,other,llplace
     &   ,llrig,lpnow(numax),lsave(numax)
      common /xyznew/ rxnew(numax), rynew(numax), rznew(numax)
      common /achepot/ acchem(nbxmax,ntmax),bnchem(nbxmax,ntmax)
      common /arig/ rlist(numax,numax),rfrom(numax),rprev(numax)
     &   ,rnum(numax),nrigi
      common /aplace/ nplace,iplace(numax,numax),pfrom(numax)
     &   ,pnum(numax),pprev(numax)

