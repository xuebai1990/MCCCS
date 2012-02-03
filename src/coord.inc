C COORD.INC
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

      logical lelect,lflucq,lqtrans,lexpand,lpresim,lring,lrigid,lplace
     &  ,lrig,lrigi,lchiral,licell
      integer nbox,nchain,nmolty,moltyp,nunit,nugrow,ntype,nchbox,ncmt
     &  ,ncmt2,iring,nrig,irig,frig,counthist,counttot,boxlink
      integer nboxi, eetype, parbox,iurot,iupdatefix,nrigmin,nrigmax
      integer parall,temtyp,riutry,rindex,prior,maxgrow,nrotbd,irotbd
      double precision rxu, ryu, rzu, mass, masst, beta,xcm,ycm,zcm
     &  ,exp_c,eta,eta2,sxcm,sycm,szcm
      double precision pmrotbd
      double precision qqu,rcmu,rintramax,hist,probf,Elect_field
      character*18 chname
      character*4 chemid
      dimension lelect(ntmax),lflucq(ntmax),lqtrans(ntmax)
     &  ,lexpand(ntmax),lplace(ntmax,nntype)
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
      dimension mass(nntype),masst(ntmax),parbox(nmax,nbxmax,ntmax)
      dimension xcm(nmax),ycm(nmax),zcm(nmax),exp_c(nmax)
     &  ,eta2(nbxmax,ntmax),sxcm(nmax),sycm(nmax),szcm(nmax)
      dimension riutry(ntmax,3),rindex(ntmax),lrigid(ntmax),lrig(ntmax)
     &  ,lrigi(ntmax,numax),lchiral(ntmax,numax),lring(ntmax)
      dimension chname(nntype), chemid(nntype)

      common /ncnunt/ moltyp,nunit,nugrow,ntype,parbox,boxlink
      common /ccord1/ probf, hist,rxu, ryu, rzu, qqu, rcmu, Elect_field
      common /ccord2/ nboxi,eetype,iurot,nchbox,ncmt,ncmt2
     &  ,parall,temtyp,licell,rintramax
      common /ccord3/ nbox,nchain,nmolty,nrigmin,nrigmax
      common /electa/ lelect,lflucq,lqtrans,lexpand,lring,lpresim
     &  ,lplace,lrigid,lrig,lrigi,lchiral

      common /masses/ mass, masst, xcm, ycm, zcm, exp_c, beta
     &  ,eta,eta2,sxcm,sycm,szcm

      common /histfix/ counthist,counttot,iupdatefix,maxgrow,iring
      common /rigidstuff/ riutry,rindex,nrig,irig,frig

      common /verbose/ chname, chemid

      common /ccord4/ pmrotbd
      common /ccord5/ irotbd,nrotbd

