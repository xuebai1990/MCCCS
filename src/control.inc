c CONTROL.INC
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

      logical lgibbs,lgrand,lnpt,lanes,llj
      logical lpbc,lpbcx,lpbcy,lpbcz,lfold,lijall,lewald
      logical lchgall,ldielect
      logical lexpsix,lmmff 
      logical lninesix
      logical ltailc,lshift
      logical lneigh,lvirial
      logical lcutcm,ldual
      integer nmax,numax,ntmax,nntype,nvmax,nbxmax,npamax,npabmax
     & ,maxbin,cmax,cmaxa,tor_bin_max,nbinmax_ete
      integer vectormax,maxvir,maxntemp,nchmax,nchtor_max,nchbn_max
      integer ntdifmx,nbinmx,angle_max,ang_bin_max,tor_max	
c --------------------------------------------------------------
c *******************************
c *** PARAMETERS FOR ENSEMBLE ***
c *******************************
c if LNPT=.TRUE.. 
c    then a NPT volume move is used to equilibrate with a pressure bath
c    (implies cubic simulation boxes)
c    else an NVT simulation is performed 
      parameter (lnpt = .false.)
c if LGIBBS=.TRUE.. 
c    then a Gibbs-ensemble simulation is performed
c    (implies cubic simulation boxes)
      parameter (lgibbs = .true.)
c if LGRAND=.TRUE. 
c    then simulation is performed in the grand-canonical ensemble
      parameter (lgrand = .false.)
c if LANES=.TRUE. 
c    then simulation is performed in the adiabatic nuclear and electronic
c    sampling technique for polarizable force fields
      parameter (lanes = .false.)
c if LVIRIAL=.TRUE.
c    then one chain will be simulated in each box independently and
c    the second viral coefficient will be calculated for their 
c    interactions at a series of distances along the x-axis
      parameter (lvirial = .false.)
c ******************************************
c *** PARAMETERS FOR BOUNDARY CONDITIONS ***
c ******************************************
c if LPBC = .TRUE.
c    then periodic boundaries are used
      parameter (lpbc = .true.)
c if LPBCX = .TRUE.
c    then periodic boundary in x-directions is used
      parameter (lpbcx = .true.)
c if LPBCY = .TRUE.
c    then periodic boundary in y-directions is used
      parameter (lpbcy = .true.)
c if LPBCZ = .TRUE.
c    then periodic boundary in z-directions is used
      parameter (lpbcz = .true.)
c if LFOLD = .TRUE.
c    then coordinates are always folded into central box
      parameter (lfold = .true.)
c ********************************************
c *** PARAMETERS FOR INTERACTION CUT-OFFS  ***
c ********************************************
c if LIJALL = .TRUE.
c    then all i-j-interactions are considered (no potential cut-off)
c    must set lallcm and ltailc to false for lijall = true
c    * check top of sumup.f if lijall = true and lchgall is false
      parameter (lijall = .false.)
c if LCHGALL= .TRUE.
c    then all the electrostatic interaction are considered 
      parameter (lchgall = .true.)
C if LCUTCM=.TRUE. then a cutoff of the centers of mass will be used
c with a value of rcmu as calculated in ctrmas
      parameter (lcutcm = .true.)
c if LEWALD=.TRUE. then ewald-sum will be used to calculate the
c electrostatic interactions.
      parameter (lewald = .true.)
c if LDIELECT=.TRUE. then dielectric constant will be calculated and
c LEWALD must be .TRUE.
      parameter (ldielect = .false.)
c ********************************************
c *** PARAMETERS FOR CBMC-PSEUDO-POTENTIAL ***
c ********************************************
c if LDUAL=.TRUE. then the external potential during a CBMC growth will
c only go out to a radius of rcutin and then will be corrected to the 
c full rcut at the end.  This is Dual Cutoff Configurational-bias Monte 
c Carlo (DC-CBMC)
      parameter (ldual = .true.) 

c ********************************************
c ***   PARAMETER FOR LENNARD-JONES 12-6   ***
c*********************************************
      parameter (llj = .true.)
c ********************************************
c *** PARAMETERS FOR EXP-6 POTENTIAL       ***
c ********************************************
c if LEXPSIX=.TRUE. exp-6 potential will be used instead of leonard jones
      parameter (lexpsix = .false.)
c ***************************************************
c *** PARAMETERS FOR BUFFERED 14-7 POTENTIAL      ***
c ***************************************************
c if LMMFF=.TRUE. buffered 14-7 potential will be used.
      parameter (lmmff = .false.)
c ************************************
c *** PARAMETERS FOR 9-6 POTENTIAL ***
c ************************************
c if LNINESIX=.TRUE. then the 9-6 potential will be used.
      parameter (lninesix = .false.)
c ***************************************
c *** PARAMETERS FOR TAIL CORRECTIONS ***
c ***************************************
c if LTAILC=.TRUE. tail corrections are added
c    (WARNING:  .lsami. in external.inc switches an intrinsic
c               tailcorrection on)
      parameter (ltailc = .true.)
c truncated and shifted potentials
      parameter (lshift = .false.)
c ***************************************
c *** PARAMETERS OF NEIGHBOR LIST     ***
C ***************************************
c if LNIEGH=.TRUE. the nearest neighbor list will be used with a 
c value of rcutnn specified by fort.4
      parameter (lneigh = .false.)
c ***************************************
c ** DIMENSIONS FOR ARRAYS             **
c ***************************************
c - nmax = maximum number of chains +2
      parameter (nmax = 700)
c - numax = maximum number of units
      parameter (numax = 30)
c - ntmax = maximum number of types of chains
      parameter (ntmax = 3)
c - nbxmax = maximum number of boxes
      parameter (nbxmax = 2) 
c - nntype = number of types of beads
      parameter (nntype = 200)
c - nvmax = maximum number of bond vibration and bending constants
      parameter (nvmax = 100)
c - npamax = maximum number of pairs to switch or swatch
      parameter (npamax = 20)
c - npabmax = maximum number of box pairs (for swatch and swap) 
      parameter (npabmax = 3)
c - maxvir = maximum number of bins for the 2nd virial coefficient
      parameter (maxvir = 1,maxntemp=3)
c - nchmax = maximum number of choices for trial sites in CBMC growth
      parameter (nchmax=100)
c - nchbn_max = maxium number of choices for bending CBMC growth
      parameter (nchbn_max=2000)
c - nchtor_max = maximum number of choices for torsion CBMC growth
      parameter (nchtor_max = 200)
c - vectormax = the maximum number of reciprocal vectors for Ewald sum
      parameter (vectormax=10000) 

c ****************************************
c *** SAFECBMC max number of bins        ***
c ****************************************
      parameter (maxbin = 201)

c *************************************************************
c *** cmax is the max number of cells for the linkcell list ***
c *** cmaxa is the max number of molecules per cell         ***
c *************************************************************
      parameter (cmax = 1000)
      parameter (cmaxa = 100)

c*************************************************************
C***
C***      FOR ANALYSIS PART
C***
C*************************************************************

C NTDIFMX = max number of diff beads in simulation

	parameter(ntdifmx=10)	

C NBINMX = max number of bins for histograms

	parameter(nbinmx=200)


C THESE ARE FOR ANGLE AND TORSIONAL DISTRIBUTION

        parameter(angle_max=50)
	parameter(ang_bin_max=360)
	parameter(tor_bin_max=360)
        parameter(tor_max=4)

C  For end to end vector calculation 

       parameter(nbinmax_ete=2000)
