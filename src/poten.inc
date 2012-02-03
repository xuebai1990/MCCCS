c POTEN.INC
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

      logical lspecial,lqchg,lij,lpl,lq14scale
      logical L_Coul_CBMC
      logical L_tor_table, L_spline, L_linear
      integer a15type	
      double precision sig2ij,epsij,epsilon,sigma,epsi,sigi,ecut
     +   ,extc12, extc3, extz0,aspecial,bspecial,q1
      double precision qqfact,xiq,jayq,jayself,fqegp,a15
     &   ,qscale
      double precision ljscale,qscale2
      
      dimension fqegp(ntmax),a15(2),a15type(ntmax,numax,numax)
      dimension extc12(nntype), extc3(nntype), extz0(nntype)
     &   ,xiq(nntype),jayself(nntype),lqchg(nntype)
     &   ,lij(nntype)

      dimension sig2ij(nntype*nntype), epsij(nntype*nntype), 
     +   ecut(nntype*nntype), jayq(nntype*nntype),
     +   epsilon(2,numax),sigma(2,numax),epsi(nntype),
     +   sigi(nntype),q1(nntype)
      dimension aspecial(nntype*nntype),bspecial(nntype*nntype)
     +   ,lspecial(nntype*nntype),lpl(nntype)
      
      dimension ljscale(ntmax,numax,numax),qscale2(ntmax,numax,numax)
      dimension lq14scale(ntmax),qscale(ntmax)
      double precision n0,n1

c ***********************************************************
c *** parameters for Generalized Lennard Jones Potential  ***
c repulsive part
      parameter(n0=12.0d0)
c   attractive part
      parameter(n1=6.0d0)  
c ****  Ref:  J. Chem. Phys. 120, 4994 (2004)         ****
c ***********************************************************
 
c *** common blocks *** 
      common /ljpot/ sig2ij,epsij,ecut,epsilon,sigma,epsi,sigi
     &   ,q1,qscale,a15,a15type
      common /ljepot/ ljscale,qscale2,lq14scale
      common /extpot/ extc12, extc3, extz0
      common /specail/ aspecial,bspecial,lspecial,lpl
      common /electr/ qqfact,xiq,jayq,jayself,fqegp,lqchg,lij
      common /etects/ L_Coul_CBMC,L_tor_table, L_spline, L_linear
