module sim_system
  use var_type,only:dp,default_path_length
  use util_runtime,only:err_exit
  use util_string,only:integer_to_string
  use util_files,only:get_iounit
  use util_search,only:LookupTable
  implicit none
  save

! CONTROL.INC
! *******************************
! *** PARAMETERS FOR ENSEMBLE ***
! *******************************
  logical::lnpt& !< if LNPT=.TRUE. then a NPT volume move is used to equilibrate with a pressure bath (implies cubic simulation boxes) else an NVT simulation is performed
   ,lgibbs& !< if LGIBBS=.TRUE. then a Gibbs-ensemble simulation is performed (implies cubic simulation boxes)
   ,lgrand& !< if LGRAND=.TRUE. then simulation is performed in the grand-canonical ensemble
   ,lanes& !< if LANES=.TRUE. then simulation is performed in the adiabatic nuclear and electronic sampling technique for polarizable force fields
   ,lvirial& !<! if LVIRIAL=.TRUE. then one chain will be simulated in each box independently and the second viral coefficient will be calculated for their interactions at a series of distances along the x-axis
   ,lmipsw& !< if lmipsw is true, then thermodynamic integration is performed for the phases, fort.35 must be supplied if so
   ,lexpee !< if lexpee is true, then expanded esnemble is performed for the phases, fort.44 must be supplied if so
! ******************************************
! *** PARAMETERS FOR BOUNDARY CONDITIONS ***
! ******************************************
  logical,parameter::lpbc=.true.& !< if LPBC = .TRUE. then periodic boundaries are used
   ,lpbcx=.true.& !< if LPBCX = .TRUE. then periodic boundary in x-directions is used
   ,lpbcy=.true.& !< if LPBCY = .TRUE. then periodic boundary in y-directions is used
   ,lpbcz=.true.& !< if LPBCZ = .TRUE. then periodic boundary in z-directions is used
   ,lfold=.true. !< if LFOLD = .TRUE. then coordinates are always folded into central box
! ********************************************
! *** PARAMETERS FOR INTERACTION CUT-OFFS  ***
! ********************************************
  logical::lijall& !< if LIJALL = .TRUE. then all i-j-interactions are considered (no potential cut-off). Must set lcutcm and ltailc to false for lijall = true *** check top of sumup.f if lijall = true and lchgall is false
   ,lchgall& !< if LCHGALL= .TRUE. then all the electrostatic interaction are considered
   ,lcutcm& !< if LCUTCM=.TRUE. then a cutoff of the centers of mass will be used with a value of rcmu as calculated in ctrmas
   ,lewald& !< if LEWALD=.TRUE. then ewald-sum will be used to calculate the electrostatic interactions.
   ,ldielect !< if LDIELECT=.TRUE. then dielectric constant will be calculated and LEWALD must be .TRUE. Correct only in NVT ensemble
! ********************************************
! *** PARAMETERS FOR CBMC-PSEUDO-POTENTIAL ***
! ********************************************
  logical::ldual !< if LDUAL=.TRUE. then the external potential during a CBMC growth will only go out to a radius of rcutin and then will be corrected to the full rcut at the end. This is Dual Cutoff Configurational-bias Monte Carlo (DC-CBMC)
! ***************************************
! *** PARAMETERS FOR TAIL CORRECTIONS ***
! ***************************************
  logical::ltailc& !< if LTAILC=.TRUE. tail corrections are added (WARNING: .lsami. in external.inc switches an intrinsic tail correction on)
   ,lshift& !< truncated and shifted potentials
   ,ltailcZeo
! ***************************************
! *** PARAMETERS OF NEIGHBOR LIST     ***
! ***************************************
  logical::lneigh !< if LNIEGH=.TRUE. the nearest neighbor list will be used with a value of rcutnn specified by fort.4
! ********************************************
! *** PARAMETER FOR IONIC SYSTEMS  ***
! ********************************************
  logical::lionic !< if LIONIC=.TRUE. System contains charged species, so system may not neutral
! ********************************************
! ***   PARAMETER FOR LENNARD-JONES 12-6   ***
!*********************************************
  logical::llj !< This is only used if you are using regular lj potential, i.e. without shifting or modifying it in any other way.
! ********************************************
! *** PARAMETERS FOR EXP-6 POTENTIAL       ***
! ********************************************
  logical::lexpsix !< if LEXPSIX=.TRUE. exp-6 potential will be used instead of lennard jones
! ***************************************************
! *** PARAMETERS FOR BUFFERED 14-7 POTENTIAL      ***
! ***************************************************
  logical::lmmff !< if LMMFF=.TRUE. buffered 14-7 potential will be used.
! ************************************
! *** PARAMETERS FOR 9-6 POTENTIAL ***
! ************************************
  logical::lninesix !< if LNINESIX=.TRUE. then the 9-6 potential will be used.
! *****************************************************
! *** PARAMETERS FOR GENERALIZED LJ POTENTIAL      ***
! *****************************************************
  logical::lgenlj !< if LGENLJ=.TRUE. then the the Generalized lennard jones potential will be used.
! ********************************************
! *** PARAMETERS FOR GAROFALINI POTENTIAL  ***
! ********************************************
  logical::lgaro !< if LGARO=.TRUE. garofalini potential will be used instead of lennard jones
  logical::lexzeo& !< implicit rigid framework for zeolites and metal-organic frameworks. Parameters for zeolites are read in subroutine SUZEO
   ,lzgrid& !< whether to use tabulated potential for the interactions with the rigid framework
   ,lelect_field !< external electric field
! surfaces:
  logical::lsami& !<
   ,lmuir& !<
   ,lslit& !<
   ,lgraphite& !< Graphite surface
   ,ljoe& !<
   ,lpsurf& !<
   ,lcorreg !<
  logical,parameter::lfepsi=.false. !< If lfepsi is true, the fluctuation of epsilon is used instead of the fluctuation of the sigma.
! ***************************************
! ** DIMENSIONS FOR ARRAYS             **
! ***************************************
  integer::io_output=6,checkpoint_interval=1800
  character(LEN=default_path_length)::file_input='fort.4',file_restart='fort.77',file_struct='input_struc.xyz',file_run='run1a.dat',file_movie='movie1a.dat',file_dipole='dipole1a.dat'
  namelist /io/ file_input,file_restart,file_struct,file_run,file_movie,file_dipole,io_output,checkpoint_interval
  integer::nntype& !< number of types of beads
   ,nbox=1,nbxmax& !< maximum number of boxes
   ,npabmax& !< maximum number of box pairs (for swatch and swap)
   ,nchain=4,nmax& !< maximum number of chains + 2
   ,nmolty=1,ntmax& !< maximum number of types of chains
   ,npamax& !< maximum number of pairs to switch or swatch
   ,numax& !< maximum number of units
   ,nchmax& !< maximum number of choices for trial sites in CBMC growth
   ,nchtor_max& !< maximum number of choices for torsion CBMC growth
   ,nchbn_max& !< maxium number of choices for bending CBMC growth
   ,nprop

  integer,parameter::smax=35& !< max no. of mstates
   ,maxvir=1,maxntemp=3& !< maximum number of bins for the 2nd virial coefficient
   ,maxbin=201& !< SAFECBMC max number of bins
   ,maxneigh=20,nener=18,nprop1=11,blockm=100

! OPENMP.inc
  integer::thread_id,thread_num,thread_num_max,thread_num_proc

! MPI.INC
  integer::myid,numprocs,groupid
  integer,parameter::rootid=0

  real,allocatable,target::boxlx(:),boxly(:),boxlz(:),rcut(:),rcutnn(:),kalp(:)
  logical,allocatable,target::lsolid(:),lrect(:)

  real,allocatable::sigi(:),epsi(:),qelect(:),mass(:),sig2ij(:),epsij(:),q1(:),xiq(:),jayself(:),rminee(:),ecut(:),jayq(:)
  logical,allocatable::lij(:),lqchg(:),lpl(:)
  character(len=4),allocatable::chemid(:)
  real,allocatable::brvib(:),brvibk(:),brben(:),brbenk(:)

  real(kind=dp),allocatable::vbox(:),vinterb(:),vtailb(:),vintrab(:),vvibb(:),vbendb(:),vtgb(:),vextb(:),velectb(:),vflucqb(:),v3garob(:),vipswb(:),vwellipswb(:)
  real,allocatable::pmrotbd(:,:),vtry(:),vtrext(:),vtrintra(:),vtrinter(:),vtrelect(:),vtrewald(:),vtrorient(:),vtrelect_intra(:),vtrelect_inter(:),bfac(:),rxp(:,:),ryp(:,:),rzp(:,:),vwellipswot(:),vwellipswnt(:),vipswnt(:),vipswot(:),epsilon_f(:,:),sigma_f(:,:),ljscale(:,:,:),qscale2(:,:,:),ee_qqu(:,:),rxnew(:),rynew(:),rznew(:),rxu(:,:),ryu(:,:),rzu(:,:),xcm(:),ycm(:),zcm(:),qqu(:,:),rxuion(:,:),ryuion(:,:),rzuion(:,:),qquion(:,:),xvec(:,:),yvec(:,:),zvec(:,:),awell(:,:,:),pmsatc(:),pmswtcb(:,:),pmisatc(:),B(:),fqegp(:),eta2(:,:),qscale(:),pmbias(:),pmbsmt(:),pmbias2(:),rmtrax(:,:),rmtray(:,:),rmtraz(:,:),rmrotx(:,:),rmroty(:,:),rmrotz(:,:),rmflcq(:,:),pmswmt(:),pmswapb(:,:),pmcbmt(:),pmall(:),pmfix(:),pmfqmt(:),pmeemt(:),pmtrmt(:),pmromt(:),masst(:),avsolinter(:,:),avsolintra(:,:),avsolbend(:,:),avsoltor(:,:),avsolelc(:,:),rxwell(:,:),rywell(:,:),rzwell(:,:),sxwell(:,:),sywell(:,:),szwell(:,:),naccu(:,:),nccold(:,:),accum(:,:),bccold(:,:),aver(:,:),baver(:,:,:),naccu1(:,:,:),nccold1(:,:,:),rcmu(:),exp_c(:),sxcm(:),sycm(:),szcm(:),ndij(:,:),nxij(:,:),nyij(:,:),nxijo(:,:),nzij(:,:),ndijo(:,:),nyijo(:,:),nzijo(:,:),favor(:),favor2(:),express(:),Elect_field(:),zshift(:),dshift(:),rmvol(:),pmvlmt(:),pmvolb(:),rmhmat(:,:),dipolex(:),dipoley(:),dipolez(:),wbox(:)
  real,allocatable::bnflcq(:,:),bsflcq(:,:),bnflcq2(:,:),bsflcq2(:,:) ! *** temporary accumulators for max. displacement updates ***
  real,allocatable::aver1(:,:,:),accum1(:,:,:),bccold1(:,:,:),baver1(:,:,:,:) ! ** Neeraj adding for solubility parameter and heat of vaporization
  integer,allocatable::ucheck(:),ntype(:,:),leaderq(:,:),invib(:,:),itvib(:,:,:),ijvib(:,:,:),inben(:,:),itben(:,:,:),ijben2(:,:,:),ijben3(:,:,:),intor(:,:),ittor(:,:,:),ijtor2(:,:,:),ijtor3(:,:,:),ijtor4(:,:,:),nrotbd(:),irotbd(:,:),splist(:,:,:),rlist(:,:),rfrom(:),rprev(:),rnum(:),iplace(:,:),pfrom(:),pnum(:),pprev(:),a15type(:,:,:),prior(:,:),wsched(:),growfrom(:),growprev(:),grownum(:),growlist(:,:),nswatb(:,:),nsampos(:),ncut(:,:),gswatc(:,:,:),nswtcb(:),box3(:,:),box4(:,:),temtyp(:),nunit(:),nugrow(:),nmaxcbmc(:),iurot(:),maxgrow(:),isolute(:),iring(:),nrig(:),irig(:,:),frig(:,:),nrigmin(:),nrigmax(:),rindex(:),riutry(:,:),ininch(:,:),nswapb(:),box1(:,:),box2(:,:),nchoi1(:),nchoi(:),nchoir(:),nchoih(:),nchtor(:),nchbna(:),nchbnb(:),icbdir(:),icbsta(:),rmexpc(:),eetype(:),ncmt(:,:),ncmt2(:,:,:),parall(:,:),parbox(:,:,:),solcount(:,:),nwell(:),moltyp(:),nboxi(:),neighbor(:,:),neigh_cnt(:),neighboro(:,:),neigh_o(:),ghost_particles(:),numberDimensionIsIsotropic(:),inix(:),iniy(:),iniz(:),inirot(:),inimix(:),nchoiq(:),box5(:),box6(:),nchbox(:)
  logical,allocatable::lovr(:),lexist(:),lexclu(:,:,:,:),lrigi(:,:),liswinc(:,:),lchiral(:,:),linclu(:,:,:),lqinclu(:,:,:),lainclu(:,:,:),wschvib(:,:),wschben(:,:),wschtor(:,:),lelect(:),lflucq(:),lqtrans(:),lexpand(:),lavbmc1(:),lavbmc2(:),lavbmc3(:),lbias(:),lring(:),lrigid(:),lrig(:),lq14scale(:),lbranch(:),lrplc(:),lplace(:,:),lwell(:),lideal(:),ltwice(:)

  logical::ldebug=.true.

! INPUTDATA.INC
  integer::run_num
  character(LEN=1)::suffix
  logical::L_add,L_sub
  integer::N_add,N_box2add,N_moltyp2add
  integer::N_sub,N_box2sub,N_moltyp2sub
  logical::L_Coul_CBMC
  real::Num_cell_a,Num_cell_b,Num_cell_c
  integer::nstep,iupdatefix
  logical::lstop,lpresim
  logical::L_tor_table,L_spline,L_linear,L_vib_table,L_bend_table,L_vdW_table,L_elect_table
  real::temp
  integer::ianalyze,nbin
  logical::lrdf,lintra,lstretch,lgvst,lbend,lete,lrhoz
  real::bin_width
  integer::iprint,imv,iratio,iblock,idiele,iheatcapacity
  logical::L_movie_xyz
  integer::nequil,ninstf,ninsth,ndumph
  real::rmin,softcut,softlog,rcutin,rbsmax,rbsmin,vol_eff
  logical::licell
  real::rintramax
  integer::boxlink
  real::Armtrax,Armtray,Armtraz,tatra,tarot
  real::tavol,taflcq
  integer::iratv,iratp
  real::pmvol,pmvolx,pmvoly
  integer::nvolb
  real::pmswat
  integer::nswaty
  real::pmswap,pmcb
  real::pmflcq
  real::pmexpc,pmexpc1
  real::pm_atom_tra
  real::pmtra
  integer::counttot
  real::probf(30,30,maxbin)
  integer::ntemp
  real::virtemp(maxntemp)

! global
  real::beta,fqbeta,hist(30,30,maxbin)
  integer::tmcc
  logical::lneighbor
  integer::counthist,nrigi
  integer::nnstep
  real::pmiswat

! *** adding this for iswatch move ***
  logical::liswatch
  integer::other,iparty,nplace
! *** end add ***

! POTEN.INC
  real::dipolexo,dipoleyo,dipolezo
  real(kind=dp)::v3garo

! EEPAR.INC
  logical::leemove,lmstate,leeacc
  integer::fmstate,sstate1,sstate2,nstate,box_state(smax),eepointp,eeirem,boxrem1,boxins1,ee_prob(smax)
  real::wee_ratio,psi(smax),um_markov(smax,smax),eeratio
  integer::mstate,ee_moltyp(smax),nmolty1

! COORD2.INC
  integer::moltion(2)
  real::exp_cion(2)

! NEIGH.INC
  integer::neigh_icnt,neighi(maxneigh)
  real::ndiji(maxneigh),nxiji(maxneigh),nyiji(maxneigh),nziji(maxneigh)

! BNBSMA.INC
  real::Abntrax,Abntray,Abntraz,Abstrax,Abstray,Abstraz

! ROSEN.INC
  real::weight,weiold,voldt,voldbb,voldtg,voldext,voldintra,voldinter,voldbvib,voldelect,voldewald,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib,vnewinter,vnewelect,vnewewald,vneworient,vnewinterr,vnewextr,vnewelectr,voldinterr ,voldextr,voldelectr,voldorient

! IPSWPAR.INC
  integer::nw,iratipsw
  parameter (nw = 4000)
  real::dvdl,acdvdl,acipsw,vipsw,pipsw ,vwellipsw,pwellipsw,etais,lambdais,bwell,vipswo,vipswn,vwellipswo,vwellipswn,lena,lenc,pwellips(3,3),pips(3,3),dhmat(3,3)
  logical::lstagea,lstageb,lstagec

! BOLTZMANN.INC
  real::rxi1,ryi1,rzi1,vi1,wi1,vext1,velect1

! FEPSI.INC
  real::aslope,bslope,ashift,bshift,a0,b0

! *** slope=0.3 a=3.05
!	parameter (a0 = 0.2003E0_dp)
!	parameter (b0 = 1.3946E0_dp)
!	parameter (aslope = 8.85E5_dp)
!	parameter (bslope = 158.25E0_dp)
!	parameter (ashift = 7.227E5_dp)
!	parameter (bshift = 505.97E0_dp)
! *** slope=0.3 a=2.90
!	parameter (a0 = 0.15561E0_dp)
!	parameter (b0 = 1.2960E0_dp)
!	parameter (aslope = 5.5475E5_dp)
!	parameter (bslope = 131.23E0_dp)
!	parameter (ashift = 4.121E5_dp)
!	parameter (bshift = 381.96E0_dp)
! *** slope=0.3 a=2.95
!	parameter (a0 = 0.16833E0_dp)
!	parameter (b0 = 1.3299E0_dp)
!	parameter (aslope = 6.5125E5_dp)
!	parameter (bslope = 139.75E0_dp)
!	parameter (ashift = 4.9975E5_dp)
!	parameter (bshift = 419.81E0_dp)
! *** slope=0.3 a=2.82
!	parameter (a0 = 0.13150E0_dp)
!	parameter (b0 = 1.2574E0_dp)
!	parameter (aslope = 4.2863E5_dp)
!	parameter (bslope = 117.5E0_dp)
!	parameter (ashift = 3.0208E5_dp)
!	parameter (bshift = 323.6E0_dp)
! *** slope=0.3 a=2.75
!	parameter (a0 = 0.11037E0_dp)
!	parameter (b0 = 1.1942E0_dp)
!	parameter (aslope = 3.4013E5_dp)
!	parameter (bslope = 108E0_dp)
!	parameter (ashift = 2.2868E5_dp)
!	parameter (bshift = 285.04E0_dp)
! *** slope=0.3 a=2.85
!	parameter (a0 = 0.14035E0_dp)
!	parameter (b0 = 1.2596E0_dp)
!	parameter (aslope = 4.7263E5_dp)
!	parameter (bslope = 123.25E0_dp)
!	parameter (ashift = 3.3978E5_dp)
!	parameter (bshift = 347.63E0_dp)
! *** slope=0.3 a=2.88
!	parameter (a0 = 0.14820E0_dp)
!	parameter (b0 = 1.2689E0_dp)
!	parameter (aslope = 5.2125E5_dp)
!	parameter (bslope = 128.75E0_dp)
!	parameter (ashift = 3.8220E5_dp)
!	parameter (bshift = 371.29E0_dp)
! *** slope=0.5 a=2.71
  parameter (a0 = -0.22818E0_dp)
  parameter (b0 = 0.44662E0_dp)
  parameter (aslope = 14.0738E5_dp)
  parameter (bslope = 351.25E0_dp)
  parameter (ashift = 3.5852E5_dp)
  parameter (bshift = 358.94E0_dp)
! *** slope=0.5 a=2.68
!	parameter (a0 = -0.23379E0_dp)
!	parameter (b0 = 0.43630E0_dp)
!	parameter (aslope = 12.78125E5_dp)
!	parameter (bslope = 337.5E0_dp)
!	parameter (ashift = 3.1904E5_dp)
!	parameter (bshift = 337.86E0_dp)
! *** slope=0.5 a=2.655
!	parameter (a0 = -0.23808E0_dp)
!	parameter (b0 = 0.425507E0_dp)
!	parameter (aslope = 11.7775E5_dp)
!	parameter (bslope = 322.04E0_dp)
!	parameter (ashift = 2.8902E5_dp)
!	parameter (bshift = 326.875E0_dp)
! *** slope=0.5 a=2.665
!	parameter (a0 = -0.23641E0_dp)
!	parameter (b0 = 0.43409E0_dp)
!	parameter (aslope = 12.171E5_dp)
!	parameter (bslope = 330E0_dp)
!	parameter (ashift = 3.0075E5_dp)
!	parameter (bshift = 326.52E0_dp)
! *** slope=0.5 a=2.65
!	parameter (a0 = -0.23903E0_dp)
!	parameter (b0 = 0.42231E0_dp)
!	parameter (aslope = 11.5875E5_dp)
!	parameter (bslope = 325E0_dp)
!	parameter (ashift = 2.8339E5_dp)
!	parameter (bshift = 319.44E0_dp)
! *** slope=0.5 a=2.60
!	parameter (a0 = -0.24796E0_dp)
!	parameter (b0 = 0.40345E0_dp)
!	parameter (aslope = 9.8238E5_dp)
!	parameter (bslope = 304E0_dp)
!	parameter (ashift = 2.3206E5_dp)
!	parameter (bshift = 288.72E0_dp)
! *** slope=0.5 a=2.55
!	parameter (a0 = -0.25703E0_dp)
!	parameter (b0 = 0.38299E0_dp)
!	parameter (aslope = 8.3075E5_dp)
!	parameter (bslope = 284.38E0_dp)
!	parameter (ashift = 1.8945E5_dp)
!	parameter (bshift = 261.09E0_dp)
! *** slope = 0.7 a=2.52
!	parameter (a0 = -0.39106E0_dp)
!	parameter (b0 = 0.08431E0_dp)
!	parameter (aslope = 25.22875E5_dp)
!	parameter (bslope = 662.75E0_dp)
!	parameter (ashift = 3.0689E5_dp)
!	parameter (bshift = 335.43E0_dp)
! *** slope = 0.7 a=2.51
!	parameter (a0 = -0.39233E0_dp)
!	parameter (b0 = 0.08203E0_dp)
!	parameter (aslope = 24.425E5_dp)
!	parameter (bslope = 653.75E0_dp)
!	parameter (ashift = 2.9499E5_dp)
!	parameter (bshift = 328.60E0_dp)
! *** slope = 0.7 a=2.50
!	parameter (a0 = -0.39357E0_dp)
!	parameter (b0 = 0.07946E0_dp)
!	parameter (aslope = 23.641E5_dp)
!	parameter (bslope = 645E0_dp)
!	parameter (ashift = 2.835E5_dp)
!	parameter (bshift = 322.13E0_dp)
! *** slope = 0.7 a=2.49
!	parameter (a0 = -0.394906E0_dp)
!	parameter (b0 = 0.075678E0_dp)
!	parameter (aslope = 22.88625E5_dp)
!	parameter (bslope = 637.735E0_dp)
!	parameter (ashift = 2.72469E5_dp)
!	parameter (bshift = 316.22E0_dp)
! *** slope = 0.7 a=2.48
!	parameter (a0 = -0.396175E0_dp)
!	parameter (b0 = 0.072973E0_dp)
!	parameter (aslope = 22.149125E5_dp)
!	parameter (bslope = 629E0_dp)
!	parameter (ashift = 2.61798E5_dp)
!	parameter (bshift = 309.95E0_dp)
! *** slope = 0.7 a=2.47
!	parameter (a0 = -0.39745E0_dp)
!	parameter (b0 = 0.070157E0_dp)
!	parameter (aslope = 21.4335E5_dp)
!	parameter (bslope = 620.75E0_dp)
!	parameter (ashift = 2.5151E5_dp)
!	parameter (bshift = 306.89E0_dp)
! *** slope = 0.7 a=2.45
!	parameter (a0 = -0.399983E0_dp)
!	parameter (b0 = 0.064159E0_dp)
!	parameter (aslope = 20.06425E5_dp)
!	parameter (bslope = 604.75E0_dp)
!	parameter (ashift = 2.32037E5_dp)
!	parameter (bshift = 292.09E0_dp)
! *** slope = 0.7 a=2.40
!	parameter (a0 = -0.40629E0_dp)
!	parameter (b0 = 0.050088E0_dp)
!	parameter (aslope = 16.9775E5_dp)
!	parameter (bslope = 565.50E0_dp)
!	parameter (ashift = 1.89213E5_dp)
!	parameter (bshift = 263.93E0_dp)
! *** slope = 0.9 a=2.30
!	parameter (a0 = -0.48332E0_dp)
!	parameter (b0 = -0.12334E0_dp)
!	parameter (aslope = 34.64125E5_dp)
!	parameter (bslope = 1014.25E0_dp)
!	parameter (ashift = 2.2815E5_dp)
!	parameter (bshift = 294.26E0_dp)
! *** slope = 0.9 a=2.32
!	parameter (a0 = -0.48137E0_dp)
!	parameter (b0 = -0.11887E0_dp)
!	parameter (aslope = 36.99125E5_dp)
!	parameter (bslope = 1041.25E0_dp)
!	parameter (ashift = 2.4744E5_dp)
!	parameter (bshift = 306.23E0_dp)
! *** slope = 0.3 a=2.85
!	parameter (a0 = 0.0E0_dp)
!	parameter (b0 = 0.0E0_dp)
!	parameter (aslope = 3.0E5_dp)
!	parameter (bslope = 0.0E0_dp)
!	parameter (ashift = 8.0E5_dp)
!	parameter (bshift = 1200.0E0_dp)

  namelist /system/nchain,nmolty,nbox

  type(LookupTable)::atoms

CONTAINS
  subroutine read_system(io_input)
    integer,intent(in)::io_input
    integer::jerr

    read(UNIT=io_input,NML=io,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) then
       call err_exit(__FILE__,__LINE__,'reading namelist: io',jerr)
    end if

    ! rewind(io_input)
    ! read(UNIT=io_input,NML=system,iostat=jerr)
    ! if (jerr.ne.0.and.jerr.ne.-1) then
    !    call err_exit(__FILE__,__LINE__,'reading namelist: system',jerr)
    ! end if
    ! nchain=nchain+2

    if (io_output.eq.5) then
       call err_exit(__FILE__,__LINE__,'unit 5 is for standard input',myid+1)
    else if(io_output.ne.6.and.io_output.ne.0.and.myid.eq.0) then
       io_output=get_iounit()
       open(unit=io_output,access='sequential',action='write',file=file_run,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) then
          call err_exit(__FILE__,__LINE__,'cannot open output file '//file_run,myid+1)
       end if
    end if
  end subroutine read_system

  subroutine allocate_cell()
    integer::jerr
    allocate(boxlx(nbxmax),boxly(nbxmax),boxlz(nbxmax),rcut(nbxmax),rcutnn(nbxmax),kalp(nbxmax),lsolid(nbxmax),lrect(nbxmax),express(nbxmax),Elect_field(nbxmax),ghost_particles(nbxmax),numberDimensionIsIsotropic(nbxmax),inix(nbxmax),iniy(nbxmax),iniz(nbxmax),inirot(nbxmax),inimix(nbxmax),nchoiq(nbxmax),box5(nbxmax),box6(nbxmax),zshift(nbxmax),dshift(nbxmax),rmvol(nbxmax),pmvlmt(nbxmax),pmvolb(nbxmax),lideal(nbxmax),ltwice(nbxmax),rmhmat(nbxmax,9),dipolex(nbxmax),dipoley(nbxmax),dipolez(nbxmax),nchbox(nbxmax),vbox(nbxmax),wbox(nbxmax),vinterb(nbxmax),vtailb(nbxmax),vintrab(nbxmax),vvibb(nbxmax),vbendb(nbxmax),vtgb(nbxmax),vextb(nbxmax),velectb(nbxmax),vflucqb(nbxmax),v3garob(nbxmax),vipswb(nbxmax),vwellipswb(nbxmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_system: allocation failed',jerr)
    end if
  end subroutine allocate_cell

  subroutine allocate_system()
    integer,parameter::initial_size=15
    integer::jerr
    allocate(ucheck(ntmax),nrotbd(ntmax),xcm(nmax),ycm(nmax),zcm(nmax),pmsatc(npamax),pmswtcb(npamax,npabmax),nswatb(npamax,2),nsampos(npamax),ncut(npamax,2),gswatc(npamax,2,2*npamax),nswtcb(npamax),box3(npamax,npabmax),box4(npamax,npabmax),pmisatc(npamax),temtyp(ntmax),B(ntmax),nunit(ntmax),nugrow(ntmax),nmaxcbmc(ntmax),iurot(ntmax),maxgrow(ntmax),isolute(ntmax),iring(ntmax),nrig(ntmax),irig(ntmax,6),frig(ntmax,6),nrigmin(ntmax),nrigmax(ntmax),rindex(ntmax),riutry(ntmax,3),lelect(ntmax),lflucq(ntmax),lqtrans(ntmax),lexpand(ntmax),lavbmc1(ntmax),lavbmc2(ntmax),lavbmc3(ntmax),lbias(ntmax),lring(ntmax),lrigid(ntmax),lrig(ntmax),lq14scale(ntmax),fqegp(ntmax),eta2(nbxmax,ntmax),qscale(ntmax),pmbias(ntmax),pmbsmt(ntmax),pmbias2(ntmax),rmtrax(ntmax,nbxmax),rmtray(ntmax,nbxmax),rmtraz(ntmax,nbxmax),rmrotx(ntmax,nbxmax),rmroty(ntmax,nbxmax),rmrotz(ntmax,nbxmax),lbranch(ntmax),ininch(ntmax,nbxmax),rmflcq(ntmax,nbxmax),pmswmt(ntmax),pmswapb(ntmax,npabmax),pmcbmt(ntmax),pmall(ntmax),pmfix(ntmax),pmfqmt(ntmax),pmeemt(ntmax),pmtrmt(ntmax),pmromt(ntmax),nswapb(ntmax),box1(ntmax,npabmax),box2(ntmax,npabmax),nchoi1(ntmax),nchoi(ntmax),nchoir(ntmax),nchoih(ntmax),nchtor(ntmax),nchbna(ntmax),nchbnb(ntmax),icbdir(ntmax),icbsta(ntmax),lrplc(ntmax),masst(ntmax),rmexpc(ntmax),eetype(ntmax),ncmt(nbxmax,ntmax),ncmt2(nbxmax,ntmax,20),parall(ntmax,nmax),parbox(nmax,nbxmax,ntmax),solcount(nbxmax,ntmax),avsolinter(nbxmax,ntmax),avsolintra(nbxmax,ntmax),avsolbend(nbxmax,ntmax),avsoltor(nbxmax,ntmax),avsolelc(nbxmax,ntmax),bnflcq(ntmax,nbxmax),bsflcq(ntmax,nbxmax),bnflcq2(ntmax,nbxmax),bsflcq2(ntmax,nbxmax),rxwell(nw,ntmax),rywell(nw,ntmax),rzwell(nw,ntmax),sxwell(nw,ntmax),sywell(nw,ntmax),szwell(nw,ntmax),nwell(ntmax),lwell(ntmax),naccu(nprop,nbxmax),nccold(nprop,nbxmax),accum(nprop,nbxmax),bccold(nprop,nbxmax),aver(nprop,nbxmax),baver(nprop,nbxmax,blockm),naccu1(nprop1,nbxmax,nbxmax),nccold1(nprop1,nbxmax,nbxmax),aver1(nprop1,nbxmax,nbxmax),accum1(nprop1,nbxmax,nbxmax),bccold1(nprop1,nbxmax,nbxmax),baver1(nprop1,nbxmax,nbxmax,blockm),moltyp(nmax),rcmu(nmax),exp_c(nmax),sxcm(nmax),sycm(nmax),szcm(nmax),nboxi(nmax),neighbor(maxneigh,nmax),neigh_cnt(nmax),neighboro(maxneigh,nmax),neigh_o(nmax),ndij(maxneigh,nmax),nxij(maxneigh,nmax),nyij(maxneigh,nmax),nxijo(maxneigh,nmax),nzij(maxneigh,nmax),ndijo(maxneigh,nmax),nyijo(maxneigh,nmax),nzijo(maxneigh,nmax),favor(nmax),favor2(nmax),ntype(ntmax,initial_size),leaderq(ntmax,initial_size),lplace(ntmax,initial_size),lrigi(ntmax,initial_size),invib(ntmax,initial_size),itvib(ntmax,initial_size,6),ijvib(ntmax,initial_size,6),inben(ntmax,initial_size),itben(ntmax,initial_size,12),ijben2(ntmax,initial_size,12),ijben3(ntmax,initial_size,12),intor(ntmax,initial_size),ittor(ntmax,initial_size,12),ijtor2(ntmax,initial_size,12),ijtor3(ntmax,initial_size,12),ijtor4(ntmax,initial_size,12),irotbd(initial_size,ntmax),pmrotbd(initial_size,ntmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_system: allocation failed',jerr)
    end if
  end subroutine allocate_system

  subroutine allocate_molecule()
    integer::jerr
    allocate(splist(npamax,numax,2),lexist(numax),lexclu(ntmax,numax,ntmax,numax),rlist(numax,numax),rfrom(numax),rprev(numax),rnum(numax),liswinc(numax,ntmax),iplace(numax,numax),pfrom(numax),pnum(numax),pprev(numax),a15type(ntmax,numax,numax),epsilon_f(2,numax),sigma_f(2,numax),ljscale(ntmax,numax,numax),qscale2(ntmax,numax,numax),ee_qqu(numax,smax),rxnew(numax),rynew(numax),rznew(numax),prior(ntmax,numax),rxu(nmax,numax),ryu(nmax,numax),rzu(nmax,numax),qqu(nmax,numax),lchiral(ntmax,numax),rxuion(numax,2),ryuion(numax,2),rzuion(numax,2),qquion(numax,2),linclu(ntmax,numax,numax),lqinclu(ntmax,numax,numax),lainclu(ntmax,numax,numax),wschvib(numax,6),wschben(numax,12),wschtor(numax,12),wsched(numax),xvec(numax,numax),yvec(numax,numax),zvec(numax,numax),growfrom(numax),growprev(numax),grownum(numax),growlist(numax,numax),awell(numax,numax,ntmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_molecule: allocation failed',jerr)
    end if
  end subroutine allocate_molecule

  subroutine checkAtom()
    if (.not.allocated(atoms%list)) call err_exit(__FILE__,__LINE__,": ATOMS section has not been defined!",myid+1)
  end subroutine checkAtom
end module sim_system
