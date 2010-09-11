! INPUTDATA.INC
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff,
! John Stubbs, and Collin Wick and Ilja Siepmann
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to
!
! Free Software Foundation, Inc.
! 59 Temple Place - Suite 330
! Boston, MA  02111-1307, USA.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical::lstop, ldie
      logical::lrdf,lintra,lstretch,lgvst,lbend,lete, lrhoz
      logical::L_movie_xyz
      character*1 suffix
      integer::run_num
      logical::L_add, L_sub
      integer::N_add,N_box2add,N_moltyp2add
      integer::N_sub,N_box2sub,N_moltyp2sub
      integer::ghost_particles
      integer::ianalyze,nbin,tmcc
      integer::nequil,ninstf,ninsth,ndumph
      integer::nstep, iprint, imv, iratio, iblock, idiele, isolute
     &       ,iratv,iratp,irsave,nnstep,nchoiq,numcoeff,rmexpc
     &       ,nswapb,box1,box2,nvolb,box5,box6,ntemp, iheatcapacity

      real(8)::Armtrax,Armtray,Armtraz
      real(8)::Num_cell_a, Num_cell_b, Num_cell_c

      real(8)::pmcb, pmtra, temp, express, ewald_precision
     &                ,rmtrax, rmtray, rmtraz, rmrotx, rmroty, rmrotz
     &                ,rmflcq,taflcq,fqbeta,tatra, tarot,pmvolx,pmvoly
     &                ,pmvol, pmswap,rmvol, rmhmat, tavol,pmfix
     &                ,upnn, upnnsq, upnndg, B, pmswmt,pmiswat,pmisatc
     &                ,pmswat,pmsatc,pmcbmt,pmtrmt,pmromt
     &                ,pmvlmt,pmall,pmflcq, pmfqmt,pmexpc,pmeemt
     &                ,pmswapb,pmvolb,virtemp,bin_width,pm_atom_tra
     &                ,pmexpc1
    
      dimension pmswmt(ntmax),pmcbmt(ntmax)
     & ,pmtrmt(ntmax),pmromt(ntmax)
     & ,pmall(ntmax),pmfqmt(ntmax),pmeemt(ntmax)
     & ,rmexpc(ntmax),pmswapb(ntmax,npabmax)
     & ,box1(ntmax,npabmax),box2(ntmax,npabmax)
     & ,box5(nbxmax),box6(nbxmax),pmvolb(nbxmax)
      dimension nchoiq(nbxmax),nswapb(ntmax)
      dimension pmsatc(npamax),pmisatc(npamax)
      dimension pmvlmt(nbxmax),virtemp(maxntemp)
      dimension numcoeff(ntmax),pmfix(ntmax)
      dimension B(ntmax),isolute(ntmax),ewald_precision(nbxmax)
      dimension express(nbxmax)
      dimension ghost_particles(nbxmax)   
 
      dimension rmvol(nbxmax),rmhmat(nbxmax,9)
      dimension rmtrax(ntmax,nbxmax),rmtray(ntmax,nbxmax)
     &   ,rmtraz(ntmax,nbxmax),rmrotx(ntmax,nbxmax)
     &   ,rmroty(ntmax,nbxmax),rmrotz(ntmax,nbxmax)
      dimension rmflcq(ntmax,nbxmax)
      common /indata/ nstep,iprint,imv,iratio,iblock, idiele, isolute
     &   ,iratv,iratp,irsave,nnstep,lstop,nchoiq,numcoeff,ldie
     &   ,rmexpc,nswapb,box1,box2,nvolb,box5,box6,ntemp,tmcc
      common /indatc/ pmswmt,pmcbmt,pmtrmt,pmromt
     &   ,pmfqmt,pmexpc,pmeemt,pmfix,pmexpc1
     &   ,pmall,pmsatc,pmisatc,pmswat,pmiswat,pmcb,pmflcq
     &   ,pmvlmt,pmtra,pmvol,pmswap,pmswapb,pmvolb,virtemp
     &   ,ewald_precision

      common /indatb/ rmvol,rmtrax,rmtray,rmtraz,rmrotx,rmroty,rmrotz
     &         ,rmflcq,temp,express,tatra, tarot,tavol,taflcq,fqbeta
     &         ,upnn, upnnsq, upnndg, B,rmhmat,pmvolx,pmvoly,bin_width,
     &          pm_atom_tra
      common /indatd/ ianalyze,lrdf,lintra,lstretch,lgvst,lbend
     &         ,nbin,lete,lrhoz,L_movie_xyz
     
      common /indatd2/ nequil,ninstf,ninsth,ndumph,run_num,suffix      
      common /addsol/  N_add,N_box2add,N_moltyp2add,L_add
      common /subsol/  N_sub,N_box2sub,N_moltyp2sub,L_sub
      common /atomdisp/ Armtrax,Armtray,Armtraz,Num_cell_a, 
     &                  Num_cell_b, Num_cell_c
      common /ghost/ ghost_particles
      common /heatcap/ iheatcapacity
