! CONTROL.INC

      logical::lgibbs,lgrand,lnpt,lanes,llj,lmipsw,lexpee
      logical::lpbc,lpbcx,lpbcy,lpbcz,lfold,lijall,lewald
      logical::lchgall,ldielect
      logical::lexpsix,lmmff 
      logical::lninesix,lgaro,lionic
      logical::lgenlj
      logical::ltailc,lshift
      logical::lneigh,lvirial
      logical::lcutcm,ldual
      integer::nmax,numax,ntmax,nntype,nvmax,nbxmax,npamax,npabmax
     & ,maxbin,cmax,cmaxa,tor_bin_max,nbinmax_ete,smax
      integer::vectormax,maxvir,maxntemp,nchmax,nchtor_max,nchbn_max
      integer::ntdifmx,nbinmx,angle_max,ang_bin_max,tor_max,iou

      common /iounit/ iou	
! --------------------------------------------------------------
! *******************************
! *** PARAMETERS FOR ENSEMBLE ***
! *******************************
! if LNPT=.TRUE.. 
!    then a NPT volume move is used to equilibrate with a pressure bath
!    (implies cubic simulation boxes)
!    else an NVT simulation is performed 
      parameter (lnpt = .true.)
! if LGIBBS=.TRUE.. 
!    then a Gibbs-ensemble simulation is performed
!    (implies cubic simulation boxes)
      parameter (lgibbs = .true.)
! if LGRAND=.TRUE. 
!    then simulation is performed in the grand-canonical ensemble
      parameter (lgrand = .false.)
! if LANES=.TRUE. 
!    then simulation is performed in the adiabatic nuclear and electronic
!    sampling technique for polarizable force fields
      parameter (lanes = .false.)
! if LVIRIAL=.TRUE.
!    then one chain will be simulated in each box independently and
!    the second viral coefficient will be calculated for their 
!    interactions at a series of distances along the x-axis
      parameter (lvirial = .false.)
! if lmipsw is true, then thermodynamic integration is performed
! for the phases, fort.35 must be supplied if so
      parameter (lmipsw = .false.)
! if lexpee is true, then expanded esnemble is performed
! for the phases, fort.44 must be supplied if so
      parameter (lexpee = .false.)
! ******************************************
! *** PARAMETERS FOR BOUNDARY CONDITIONS ***
! ******************************************
! if LPBC = .TRUE.
!    then periodic boundaries are used
      parameter (lpbc = .true.)
! if LPBCX = .TRUE.
!    then periodic boundary in x-directions is used
      parameter (lpbcx = .true.)
! if LPBCY = .TRUE.
!    then periodic boundary in y-directions is used
      parameter (lpbcy = .true.)
! if LPBCZ = .TRUE.
!    then periodic boundary in z-directions is used
      parameter (lpbcz = .true.)
! if LFOLD = .TRUE.
!    then coordinates are always folded into central box
      parameter (lfold = .true.)
! ********************************************
! *** PARAMETERS FOR INTERACTION CUT-OFFS  ***
! ********************************************
! if LIJALL = .TRUE.
!    then all i-j-interactions are considered (no potential cut-off)
!    must set lcutcm and ltailc to false for lijall = true
!    * check top of sumup.f if lijall = true and lchgall is false
      parameter (lijall = .false.)
! if LCHGALL= .TRUE.
!    then all the electrostatic interaction are considered 
      parameter (lchgall = .false.)
! if LCUTCM=.TRUE. then a cutoff of the centers of mass will be used
! with a value of rcmu as calculated in ctrmas
      parameter (lcutcm = .true.)
! if LEWALD=.TRUE. then ewald-sum will be used to calculate the
! electrostatic interactions.
      parameter (lewald = .true.)
! if LDIELECT=.TRUE. then dielectric constant will be calculated and
! LEWALD must be .TRUE.
! Correct only in NVT ensemble
      parameter (ldielect = .false.)
! ********************************************
! *** PARAMETERS FOR CBMC-PSEUDO-POTENTIAL ***
! ********************************************
! if LDUAL=.TRUE. then the external potential during a CBMC growth will
! only go out to a radius of rcutin and then will be corrected to the 
! full rcut at the end.  This is Dual Cutoff Configurational-bias Monte 
! Carlo (DC-CBMC)
      parameter (ldual = .true.) 

! ********************************************
! ***   PARAMETER FOR LENNARD-JONES 12-6   ***
!*********************************************
! This is only used if you are using regular lj potential ie without shifting or modifying 
! it in any other way. 
      parameter (llj = .false.)
! ********************************************
! *** PARAMETERS FOR EXP-6 POTENTIAL       ***
! ********************************************
! if LEXPSIX=.TRUE. exp-6 potential will be used instead of lennard jones
      parameter (lexpsix = .false.)
! ***************************************************
! *** PARAMETERS FOR BUFFERED 14-7 POTENTIAL      ***
! ***************************************************
! if LMMFF=.TRUE. buffered 14-7 potential will be used.
      parameter (lmmff = .false.)
! ************************************
! *** PARAMETERS FOR 9-6 POTENTIAL ***
! ************************************
! if LNINESIX=.TRUE. then the 9-6 potential will be used.
      parameter (lninesix = .false.)
! *****************************************************
! *** PARAMETERS FOR GENERALIZED LJ POTENTIAL      ***
! *****************************************************
! if LGENLJ=.TRUE. then the the Generalized lennard jones
! potential will be used.
      parameter (lgenlj=.false.)
! ********************************************
! *** PARAMETERS FOR GAROFALINI POTENTIAL  ***
! ********************************************
! if LGARO=.TRUE. garofalini potential will be used instead of lennard jones
      parameter (lgaro = .false.)
! ********************************************
! *** PARAMETER FOR IONIC SYSTEMS  ***
! ********************************************
! if LIONIC=.TRUE. System contains charged species, so system may not neutral
      parameter (lionic = .false.)
! ***************************************
! *** PARAMETERS FOR TAIL CORRECTIONS ***
! ***************************************
! if LTAILC=.TRUE. tail corrections are added
!    (WARNING:  .lsami. in external.inc switches an intrinsic
!               tailcorrection on)
      parameter (ltailc = .false.)
! truncated and shifted potentials
      parameter (lshift = .false.)
! ***************************************
! *** PARAMETERS OF NEIGHBOR LIST     ***
! ***************************************
! if LNIEGH=.TRUE. the nearest neighbor list will be used with a 
! value of rcutnn specified by fort.4
      parameter (lneigh = .false.)
! ***************************************
! ** DIMENSIONS FOR ARRAYS             **
! ***************************************
! - nmax = maximum number of chains +2
      parameter (nmax = 5000)
! - numax = maximum number of units
      parameter (numax = 12)
! - ntmax = maximum number of types of chains
      parameter (ntmax = 10)
! - smax = max no. of mstates
      parameter (smax = 35)
! - nbxmax = maximum number of boxes
      parameter (nbxmax = 3) 
! - nntype = number of types of beads
      parameter (nntype = 450)
! - nvmax = maximum number of bond vibration and bending constants
      parameter (nvmax = 300)
! - npamax = maximum number of pairs to switch or swatch
      parameter (npamax = 10)
! - npabmax = maximum number of box pairs (for swatch and swap) 
      parameter (npabmax = 3)
! - maxvir = maximum number of bins for the 2nd virial coefficient
      parameter (maxvir = 1,maxntemp=3)
! - nchmax = maximum number of choices for trial sites in CBMC growth
      parameter (nchmax=100)
! - nchbn_max = maxium number of choices for bending CBMC growth
      parameter (nchbn_max=2000)
! - nchtor_max = maximum number of choices for torsion CBMC growth
      parameter (nchtor_max = 200)
! - vectormax = the maximum number of reciprocal vectors for Ewald sum
      parameter (vectormax=100000) 

! ****************************************
! *** SAFECBMC max number of bins        ***
! ****************************************
      parameter (maxbin = 201)

! *************************************************************
! *** cmax is the max number of cells for the linkcell list ***
! *** cmaxa is the max number of molecules per cell         ***
! *************************************************************
      parameter (cmax = 1000)
      parameter (cmaxa = 100)

!*************************************************************
!***
!***      FOR ANALYSIS PART
!***
!*************************************************************

! NTDIFMX = max number of diff beads in simulation

	parameter(ntdifmx=2)	

! NBINMX = max number of bins for histograms

	parameter(nbinmx=2)

! THESE ARE FOR ANGLE AND TORSIONAL DISTRIBUTION

        parameter(angle_max=2)
	parameter(ang_bin_max=2)
	parameter(tor_bin_max=2)
        parameter(tor_max=2)

!  For end to end vector calculation 

       parameter(nbinmax_ete=2)
