C CBMC.INC
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff,
c John Stubbs, and Collin Wick and Ilja Siepmann
c
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to
c
c Free Software Foundation, Inc.
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical lovr,lexist,lexshed,llplace,lpnow,llrig,lsave,lrplc

c *** adding this for iswatch move ***
      logical liswinc, liswatch
      dimension liswinc(numax,ntmax),llplace(ntmax)
      integer other,iparty,iplace, pfrom, pnum, pprev,nplace
c *** end add ***

      double precision bncb,bscb,rxp,ryp,rzp,rxnew,rynew,rznew
     +   ,bnregr,bsregr,bfac,bsswat,bsiswat,fbncb,fbscb
     +   ,vtry,vtrext,vtrintra,vtrinter,vtrelect,vtrewald
     +   ,vtrorient,vtrelect_intra,vtrelect_inter
      integer nchoi1,nchoi,icbdir,icbsta,nmaxcbmc,nchoih,nchtor
     &  ,nchbna,nchbnb,nrigi,nchoir,rlist,rfrom,rprev,rnum
      double precision acchem, bnchem
      dimension lexist(numax),lexshed(numax)
      common /pbias/  bncb(ntmax,numax), bscb(ntmax,2,numax)
     +   ,bnregr(ntmax),bsregr(ntmax,2)
     +   ,bsswat(npamax,npabmax),bsiswat(npamax,npabmax)
     +   ,fbncb(ntmax,numax), fbscb(ntmax,2,numax)
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

