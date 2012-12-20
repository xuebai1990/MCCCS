module sim_system
  use var_type
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
  integer::io_output=6
  character(LEN=default_path_length)::file_input='fort.4',file_restart='fort.77',file_struct='input_struc.xyz',file_run='run1a.dat',file_movie='movie1a.dat',file_dipole='dipole1a.dat'
  namelist /io/ file_input,file_restart,file_struct,file_run,file_movie,file_dipole,io_output
  integer,parameter::nmax = 7000& !< maximum number of chains +2
   ,numax = 35& !< - numax = maximum number of units
   ,ntmax = 35& !< ntmax = maximum number of types of chains
   ,smax = 35& !< smax = max no. of mstates
   ,nbxmax = 3& !< nbxmax = maximum number of boxes
   ,nntype = 450& !< nntype = number of types of beads
   ,natom=3& !< for exsix
   ,nxatom=5& !< for ninesix
   ,nvmax = 300& !< nvmax = maximum number of bond vibration and bending constants
   ,npamax = 40& !< npamax = maximum number of pairs to switch or swatch
   ,npabmax = 3& !< npabmax = maximum number of box pairs (for swatch and swap)
   ,maxvir = 1,maxntemp=3& !< maxvir = maximum number of bins for the 2nd virial coefficient
   ,nchmax=100& !< nchmax = maximum number of choices for trial sites in CBMC growth
   ,nchbn_max=2000& !< nchbn_max = maxium number of choices for bending CBMC growth
   ,nchtor_max = 200& !< nchtor_max = maximum number of choices for torsion CBMC growth
   ,maxbin = 201& !< SAFECBMC max number of bins
   ,cmax = 1000& !< cmax is the max number of cells for the linkcell list
   ,cmaxa = 100& !< cmaxa is the max number of molecules per cell
   ,maxneigh = 20
!*************************************************************
!***
!***      FOR ANALYSIS PART
!***
!*************************************************************
   integer,parameter::ntdifmx=2& !< max number of diff beads in simulation
    ,nbinmx=2& !< NBINMX = max number of bins for histograms
    ,angle_max=2,ang_bin_max=2,tor_bin_max=2,tor_max=2& !< THESE ARE FOR ANGLE AND TORSIONAL DISTRIBUTION
    ,nbinmax_ete=2 !< For end to end vector calculation

! OPENMP.inc
  integer::thread_id,thread_num,thread_num_max,thread_num_proc

! MPI.INC
  integer::myid,numprocs,ierr
  integer,parameter::numprocmax=32

! INPUTDATA.INC
  logical::ldie
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
  integer::nbox=1
  real::express(nbxmax)
  integer::ghost_particles(nbxmax)
  real::temp,Elect_field(nbxmax)
  integer::ianalyze,nbin
  logical::lrdf,lintra,lstretch,lgvst,lbend,lete,lrhoz
  real::bin_width
  integer::iprint,imv,iratio,iblock,idiele,iheatcapacity
  logical::L_movie_xyz
  integer::nequil,ninstf,ninsth,ndumph
  real,target::boxlx(nbxmax),boxly(nbxmax),boxlz(nbxmax),rcut(nbxmax),rcutnn(nbxmax),kalp(nbxmax)
  logical,target::lsolid(nbxmax),lrect(nbxmax)
  integer::numberDimensionIsIsotropic(nbxmax)
  integer::nchain=4,nmolty=1,temtyp(ntmax)
  real::B(ntmax)
  real::rmin,softcut,softlog,rcutin,rbsmax,rbsmin,vol_eff
  integer::nunit(ntmax),nugrow(ntmax),nmaxcbmc(ntmax),iurot(ntmax),maxgrow(ntmax),isolute(ntmax),iring(ntmax),nrig(ntmax),irig(ntmax,6),frig(ntmax,6),nrigmin(ntmax),nrigmax(ntmax),rindex(ntmax),riutry(ntmax,3),ntype(ntmax,numax),leaderq(ntmax,numax),invib(ntmax,numax),itvib(ntmax,numax,6),ijvib(ntmax,numax,6),inben(ntmax,numax),itben(ntmax,numax,12),ijben2(ntmax,numax,12),ijben3(ntmax,numax,12),intor(ntmax,numax) ,ittor(ntmax,numax,12),ijtor2(ntmax,numax,12),ijtor3(ntmax,numax ,12),ijtor4(ntmax,numax,12),nrotbd(ntmax),irotbd(numax,ntmax)
  logical::lelect(ntmax),lflucq(ntmax),lqtrans(ntmax),lexpand(ntmax),lavbmc1(ntmax),lavbmc2(ntmax),lavbmc3(ntmax),lbias(ntmax),lring(ntmax),lrigid(ntmax),lrig(ntmax),lq14scale(ntmax)
  real::fqegp(ntmax),eta2(nbxmax,ntmax),qscale(ntmax),pmbias(ntmax),pmbsmt(ntmax),pmbias2(ntmax),pmrotbd(numax,ntmax)
  logical::licell
  real::rintramax
  integer::boxlink
  real::Armtrax,Armtray,Armtraz,rmtrax(ntmax,nbxmax),rmtray(ntmax,nbxmax),rmtraz(ntmax,nbxmax),rmrotx(ntmax,nbxmax),rmroty(ntmax,nbxmax),rmrotz(ntmax,nbxmax),tatra, tarot
  logical::lbranch(ntmax)
  integer::ininch(ntmax,nbxmax),inix(nbxmax),iniy(nbxmax),iniz(nbxmax),inirot(nbxmax),inimix(nbxmax),nchoiq(nbxmax)
  real::zshift(nbxmax),dshift(nbxmax)
  real::rmvol(nbxmax),tavol,rmflcq(ntmax,nbxmax),taflcq
  integer::iratv,iratp
  real::pmvol,pmvlmt(nbxmax),pmvolb(nbxmax),pmvolx,pmvoly
  integer::nvolb,box5(nbxmax),box6(nbxmax)
  real::pmswat,pmsatc(npamax),pmswtcb(npamax,npabmax)
  integer::nswaty,nswatb(npamax,npabmax),nsampos(npamax),ncut(npamax,2),splist(npamax,numax,2),gswatc(npamax,2,2*npamax),nswtcb(npamax),box3(npamax,npabmax),box4(npamax,npabmax)
  real::pmswap,pmswmt(ntmax),pmswapb(ntmax,npabmax)
  integer::nswapb(ntmax),box1(ntmax,npabmax),box2(ntmax,npabmax)
  real::pmcb,pmcbmt(ntmax),pmall(ntmax),pmfix(ntmax)
  real::pmflcq,pmfqmt(ntmax)
  real::pmexpc,pmeemt(ntmax),pmexpc1
  real::pm_atom_tra
  real::pmtra,pmtrmt(ntmax)
  real::pmromt(ntmax)
  integer::nchoi1(ntmax),nchoi(ntmax),nchoir(ntmax),nchoih(ntmax),nchtor(ntmax)
  integer::counttot
  real::probf(30,30,maxbin)
  integer::nchbna(ntmax),nchbnb(ntmax),icbdir(ntmax),icbsta(ntmax)
  logical::lexclu(ntmax,numax,ntmax,numax)
  integer::ntemp
  real::virtemp(maxntemp)
  logical::lideal(nbxmax),ltwice(nbxmax),lrplc(ntmax)

! global
  real::beta,fqbeta,mass(0:nntype),masst(ntmax),hist(30,30,maxbin)
  integer::moltyp(nmax),tmcc
  logical::lrigi(ntmax,numax),lpl(nntype),lneighbor

  integer::nchbox,ncmt,ncmt2,counthist,nrigi,rlist(numax ,numax),rfrom(numax),rprev(numax),rnum(numax)
  integer::irsave,nnstep,rmexpc
  real::rmhmat,upnn,upnnsq,upnndg,pmiswat,pmisatc
  dimension rmexpc(ntmax)
  dimension pmisatc(npamax)
  dimension rmhmat(nbxmax,9)

! *** adding this for iswatch move ***
  logical::liswinc(numax,ntmax),liswatch
  integer::other,iparty,iplace(numax,numax),pfrom(numax),pnum(numax),pprev(numax),nplace
! *** end add ***

! POTEN.INC
  logical::lspecial,lqchg,lij,lovr(nchmax),lexist(numax)
  real::dipolex(nbxmax),dipoley(nbxmax),dipolez(nbxmax),dipolexo,dipoleyo,dipolezo
  integer::a15type(ntmax,numax,numax)
  real::sig2ij,epsij,epsilon,sigma,epsi,sigi,ecut,extc12, extc3, extz0,aspecial,bspecial,q1
  real::xiq,jayq,jayself,a15(2),qelect,vtry(nchmax),vtrext(nchmax),vtrintra(nchmax),vtrinter(nchmax),vtrelect(nchmax),vtrewald(nchmax),vtrorient(nchmax),vtrelect_intra(nchmax),vtrelect_inter(nchmax),bfac(nchmax)
  real::ljscale,qscale2
  real::v3garo,rminee(nntype*nntype),rxp(numax,nchmax),ryp(numax,nchmax),rzp(numax,nchmax)
  dimension extc12(nntype),extc3(nntype),extz0(nntype),xiq(0:nntype),jayself(nntype),lqchg(0:nntype),lij(0:nntype)
  dimension sig2ij(nntype*nntype),epsij(nntype*nntype),ecut(nntype*nntype),jayq(nntype*nntype),epsilon(2,numax),sigma(2,numax),epsi(0:nntype),sigi(0:nntype),q1(nntype),qelect(0:nntype)
  dimension aspecial(nntype*nntype),bspecial(nntype*nntype),lspecial(nntype*nntype)
  dimension ljscale(ntmax,numax,numax),qscale2(ntmax,numax,numax)

! EEPAR.INC
  logical::leemove,lmstate,leeacc
  integer::fmstate,sstate1,sstate2,nstate,box_state,eepointp,eeirem,boxrem1,boxins1,ee_prob
  real::wee_ratio,psi,ee_qqu,um_markov,eeratio
  dimension psi(smax),ee_qqu(numax,smax),box_state(smax),um_markov(smax,smax),ee_prob(smax)
  integer::mstate,ee_moltyp(smax),nmolty1

! COORD.INC
  logical::lplace,lchiral
  integer::nboxi,eetype,parbox
  integer::parall,prior
  real::rxu, ryu, rzu, xcm ,ycm,zcm,exp_c,sxcm,sycm,szcm
  real::rxnew(numax),rynew(numax),rznew(numax)
  real::qqu,rcmu
  character(LEN=18)::chname
  character(LEN=4)::chemid
  dimension lplace(ntmax,nntype)
  dimension prior(ntmax,numax)
  dimension nboxi(nmax),eetype(ntmax),rcmu(nmax)
  dimension rxu(nmax,numax),ryu(nmax,numax),rzu(nmax,numax)
  dimension qqu(nmax,numax)
  dimension nchbox(nbxmax),ncmt(nbxmax,ntmax),ncmt2(nbxmax,ntmax,20)
  dimension parall(ntmax,nmax)
  dimension parbox(nmax,nbxmax,ntmax)
  dimension xcm(nmax),ycm(nmax),zcm(nmax),exp_c(nmax),sxcm(nmax),sycm(nmax),szcm(nmax)
  dimension lchiral(ntmax,numax)
  dimension chname(0:nntype), chemid(nntype)

! COORD2.INC
  integer::moltion(2)
  real::rxuion(numax,2),ryuion(numax,2),rzuion(numax,2),qquion(numax,2),exp_cion(2)

! CONNECT.INC
  logical::linclu(ntmax,numax,numax),lqinclu(ntmax,numax,numax),lainclu(ntmax,numax,numax),wschvib(numax,6),wschben(numax,12),wschtor(numax,12)
  integer::wsched(numax)
  real::brvib(nvmax),brvibk(nvmax),brvibmin(60),brvibmax(60),brben(nvmax),brbenk(nvmax)

! SYSTEM.INC
  integer::icell,nicell,iucell,solcount
  real::avsolinter,avsolintra,avsolbend,avsoltor,avsolelc
  dimension icell(nmax),nicell(cmax),iucell(cmax,cmaxa) ,solcount(nbxmax,ntmax)
  dimension avsolinter(nbxmax,ntmax),avsolintra(nbxmax,ntmax) ,avsolbend(nbxmax,ntmax),avsoltor(nbxmax,ntmax) ,avsolelc(nbxmax,ntmax)

! NEIGH.INC
! lnn(1,1):replace the 1s with nmax to use neighbor list
  logical::lnn(1,1)
  integer::neighbor(maxneigh,nmax),neigh_cnt(nmax) ,neigh_icnt,neighi(maxneigh),neighboro(maxneigh,nmax) ,neigh_o(nmax)
  real::ndij(maxneigh,nmax),nxij(maxneigh,nmax),nyij(maxneigh,nmax),nxijo(maxneigh,nmax),nzij(maxneigh ,nmax),ndiji(maxneigh),nxiji(maxneigh),nyiji(maxneigh),nziji(maxneigh),ndijo(maxneigh,nmax),nyijo(maxneigh,nmax),nzijo(maxneigh,nmax)

! NEIGH2.INC
! disvec(2,1,3):replace the 1 with nmax to use neighbor lists
  real::disvec(2,1,3)

! ENSEMBLE.INC
  real::vbox(nbxmax),wbox(nbxmax) ,vinterb(nbxmax),vtailb(nbxmax),vintrab(nbxmax),vvibb(nbxmax) ,vbendb(nbxmax),vtgb(nbxmax),vextb(nbxmax), velectb(nbxmax) ,vflucqb(nbxmax),v3garob(nbxmax)
!      real::velectb_intra(nbxmax), velectb_inter(nbxmax)

! FIX.INC
  real::xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)

! BNBSMA.INC
! *** temporary accumulators for max. displacement updates ***
  real::bnflcq,bsflcq,bnflcq2,bsflcq2
  real::Abntrax,Abntray,Abntraz,Abstrax,Abstray,Abstraz
  dimension bnflcq(ntmax,nbxmax),bsflcq(ntmax,nbxmax),bnflcq2(ntmax,nbxmax),bsflcq2(ntmax,nbxmax)

! ROSEN.INC
  integer::growfrom,growprev,growlist,grownum
  dimension growfrom(numax),growprev(numax),grownum(numax)
  dimension growlist(numax,numax)
  real::weight,weiold,voldt,voldbb,voldtg,voldext,voldintra,voldinter,voldbvib,voldelect,voldewald,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib,vnewinter,vnewelect,vnewewald,vneworient,vnewinterr,vnewextr,vnewelectr,voldinterr ,voldextr,voldelectr,voldorient

! IPSWPAR.INC
  integer::nw,nwell,iratipsw
  parameter (nw = 4000)
  real::dvdl,acdvdl,acipsw,vipsw,pipsw ,vwellipsw,pwellipsw,vwellipswb,vipswb,etais,lambdais,rxwell ,rywell,rzwell,awell,bwell,vipswo,vipswn,vwellipswo,vwellipswn ,lena,lenc,pwellips,pips,dhmat,sxwell,sywell,szwell,vwellipswot ,vwellipswnt,vipswnt,vipswot
  logical::lwell,lstagea,lstageb,lstagec
  dimension vipswb(nbxmax),vwellipswb(nbxmax),rxwell(nw,ntmax) ,rywell(nw,ntmax) ,rzwell(nw,ntmax),pwellips(3,3),pips(3,3),nwell(ntmax) ,lwell(ntmax),awell(numax,numax,ntmax),dhmat(3,3),sxwell(nw,ntmax),sywell(nw,ntmax),szwell(nw,ntmax) ,vwellipswot(nchmax),vwellipswnt(nchmax) ,vipswot(nchmax),vipswnt(nchmax)

! BLKAVG.INC
  integer::nener,nprop,blockm
  real::naccu,nccold,accum,bccold,aver,baver
! ** Neeraj adding for solubility parameter and heat of vaporization
  real::naccu1,nccold1,accum1,bccold1,aver1,baver1
  integer::nprop1
  parameter (nener=18,nprop=4+nener+(4*ntmax)+3,blockm=100)
  parameter (nprop1=11)
  dimension naccu(nprop,nbxmax),nccold(nprop,nbxmax)
  dimension accum(nprop,nbxmax),bccold(nprop,nbxmax)
  dimension aver(nprop,nbxmax), baver(nprop,nbxmax,blockm)
  dimension naccu1(nprop1,nbxmax,nbxmax), nccold1(nprop1,nbxmax,nbxmax), aver1(nprop1,nbxmax,nbxmax)
  dimension accum1(nprop1,nbxmax,nbxmax), bccold1(nprop1,nbxmax,nbxmax), baver1(nprop1,nbxmax,nbxmax,blockm)

! BOLTZMANN.INC
  real::rxi1,ryi1,rzi1,vi1,wi1,vext1,velect1

! FEPSI.INC
  real::aslope,bslope,ashift,bshift,a0,b0
  real::favor(nmax),favor2(nmax)

! *** slope=0.3 a=3.05
!	parameter (a0 = 0.2003d0)
!	parameter (b0 = 1.3946d0)
!	parameter (aslope = 8.85d5)
!	parameter (bslope = 158.25d0)
!	parameter (ashift = 7.227d5)
!	parameter (bshift = 505.97d0)
! *** slope=0.3 a=2.90
!	parameter (a0 = 0.15561d0)
!	parameter (b0 = 1.2960d0)
!	parameter (aslope = 5.5475d5)
!	parameter (bslope = 131.23d0)
!	parameter (ashift = 4.121d5)
!	parameter (bshift = 381.96d0)
! *** slope=0.3 a=2.95
!	parameter (a0 = 0.16833d0)
!	parameter (b0 = 1.3299d0)
!	parameter (aslope = 6.5125d5)
!	parameter (bslope = 139.75d0)
!	parameter (ashift = 4.9975d5)
!	parameter (bshift = 419.81d0)
! *** slope=0.3 a=2.82
!	parameter (a0 = 0.13150d0)
!	parameter (b0 = 1.2574d0)
!	parameter (aslope = 4.2863d5)
!	parameter (bslope = 117.5d0)
!	parameter (ashift = 3.0208d5)
!	parameter (bshift = 323.6d0)
! *** slope=0.3 a=2.75
!	parameter (a0 = 0.11037d0)
!	parameter (b0 = 1.1942d0)
!	parameter (aslope = 3.4013d5)
!	parameter (bslope = 108d0)
!	parameter (ashift = 2.2868d5)
!	parameter (bshift = 285.04d0)
! *** slope=0.3 a=2.85
!	parameter (a0 = 0.14035d0)
!	parameter (b0 = 1.2596d0)
!	parameter (aslope = 4.7263d5)
!	parameter (bslope = 123.25d0)
!	parameter (ashift = 3.3978d5)
!	parameter (bshift = 347.63d0)
! *** slope=0.3 a=2.88
!	parameter (a0 = 0.14820d0)
!	parameter (b0 = 1.2689d0)
!	parameter (aslope = 5.2125d5)
!	parameter (bslope = 128.75d0)
!	parameter (ashift = 3.8220d5)
!	parameter (bshift = 371.29d0)
! *** slope=0.5 a=2.71
  parameter (a0 = -0.22818d0)
  parameter (b0 = 0.44662d0)
  parameter (aslope = 14.0738d5)
  parameter (bslope = 351.25d0)
  parameter (ashift = 3.5852d5)
  parameter (bshift = 358.94d0)
! *** slope=0.5 a=2.68
!	parameter (a0 = -0.23379d0)
!	parameter (b0 = 0.43630d0)
!	parameter (aslope = 12.78125d5)
!	parameter (bslope = 337.5d0)
!	parameter (ashift = 3.1904d5)
!	parameter (bshift = 337.86d0)
! *** slope=0.5 a=2.655
!	parameter (a0 = -0.23808d0)
!	parameter (b0 = 0.425507d0)
!	parameter (aslope = 11.7775d5)
!	parameter (bslope = 322.04d0)
!	parameter (ashift = 2.8902d5)
!	parameter (bshift = 326.875d0)
! *** slope=0.5 a=2.665
!	parameter (a0 = -0.23641d0)
!	parameter (b0 = 0.43409d0)
!	parameter (aslope = 12.171d5)
!	parameter (bslope = 330d0)
!	parameter (ashift = 3.0075d5)
!	parameter (bshift = 326.52d0)
! *** slope=0.5 a=2.65
!	parameter (a0 = -0.23903d0)
!	parameter (b0 = 0.42231d0)
!	parameter (aslope = 11.5875d5)
!	parameter (bslope = 325d0)
!	parameter (ashift = 2.8339d5)
!	parameter (bshift = 319.44d0)
! *** slope=0.5 a=2.60
!	parameter (a0 = -0.24796d0)
!	parameter (b0 = 0.40345d0)
!	parameter (aslope = 9.8238d5)
!	parameter (bslope = 304d0)
!	parameter (ashift = 2.3206d5)
!	parameter (bshift = 288.72d0)
! *** slope=0.5 a=2.55
!	parameter (a0 = -0.25703d0)
!	parameter (b0 = 0.38299d0)
!	parameter (aslope = 8.3075d5)
!	parameter (bslope = 284.38d0)
!	parameter (ashift = 1.8945d5)
!	parameter (bshift = 261.09d0)
! *** slope = 0.7 a=2.52
!	parameter (a0 = -0.39106d0)
!	parameter (b0 = 0.08431d0)
!	parameter (aslope = 25.22875d5)
!	parameter (bslope = 662.75d0)
!	parameter (ashift = 3.0689d5)
!	parameter (bshift = 335.43d0)
! *** slope = 0.7 a=2.51
!	parameter (a0 = -0.39233d0)
!	parameter (b0 = 0.08203d0)
!	parameter (aslope = 24.425d5)
!	parameter (bslope = 653.75d0)
!	parameter (ashift = 2.9499d5)
!	parameter (bshift = 328.60d0)
! *** slope = 0.7 a=2.50
!	parameter (a0 = -0.39357d0)
!	parameter (b0 = 0.07946d0)
!	parameter (aslope = 23.641d5)
!	parameter (bslope = 645d0)
!	parameter (ashift = 2.835d5)
!	parameter (bshift = 322.13d0)
! *** slope = 0.7 a=2.49
!	parameter (a0 = -0.394906d0)
!	parameter (b0 = 0.075678d0)
!	parameter (aslope = 22.88625d5)
!	parameter (bslope = 637.735d0)
!	parameter (ashift = 2.72469d5)
!	parameter (bshift = 316.22d0)
! *** slope = 0.7 a=2.48
!	parameter (a0 = -0.396175d0)
!	parameter (b0 = 0.072973d0)
!	parameter (aslope = 22.149125d5)
!	parameter (bslope = 629d0)
!	parameter (ashift = 2.61798d5)
!	parameter (bshift = 309.95d0)
! *** slope = 0.7 a=2.47
!	parameter (a0 = -0.39745d0)
!	parameter (b0 = 0.070157d0)
!	parameter (aslope = 21.4335d5)
!	parameter (bslope = 620.75d0)
!	parameter (ashift = 2.5151d5)
!	parameter (bshift = 306.89d0)
! *** slope = 0.7 a=2.45
!	parameter (a0 = -0.399983d0)
!	parameter (b0 = 0.064159d0)
!	parameter (aslope = 20.06425d5)
!	parameter (bslope = 604.75d0)
!	parameter (ashift = 2.32037d5)
!	parameter (bshift = 292.09d0)
! *** slope = 0.7 a=2.40
!	parameter (a0 = -0.40629d0)
!	parameter (b0 = 0.050088d0)
!	parameter (aslope = 16.9775d5)
!	parameter (bslope = 565.50d0)
!	parameter (ashift = 1.89213d5)
!	parameter (bshift = 263.93d0)
! *** slope = 0.9 a=2.30
!	parameter (a0 = -0.48332d0)
!	parameter (b0 = -0.12334d0)
!	parameter (aslope = 34.64125d5)
!	parameter (bslope = 1014.25d0)
!	parameter (ashift = 2.2815d5)
!	parameter (bshift = 294.26d0)
! *** slope = 0.9 a=2.32
!	parameter (a0 = -0.48137d0)
!	parameter (b0 = -0.11887d0)
!	parameter (aslope = 36.99125d5)
!	parameter (bslope = 1041.25d0)
!	parameter (ashift = 2.4744d5)
!	parameter (bshift = 306.23d0)
! *** slope = 0.3 a=2.85
!	parameter (a0 = 0.0d0)
!	parameter (b0 = 0.0d0)
!	parameter (aslope = 3.0d5)
!	parameter (bslope = 0.0d0)
!	parameter (ashift = 8.0d5)
!	parameter (bshift = 1200.0d0)

  namelist /system/nchain,nmolty,nbox

  type(LookupTable)::atoms !,bonds,angles,dihedrals

CONTAINS
  subroutine read_system(io_input)
    integer,intent(in)::io_input
    integer::jerr

    read(UNIT=io_input,NML=io,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) then
       write(io_output,*) 'ERROR ',jerr,' in ',TRIM(__FILE__),':',__LINE__
       call err_exit('reading namelist: io')
    end if

    rewind(io_input)
    read(UNIT=io_input,NML=system,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) then
       write(io_output,*) 'ERROR ',jerr,' in ',TRIM(__FILE__),':',__LINE__
       call err_exit('reading namelist: system')
    end if
    nchain=nchain+2

    if (io_output.eq.5) then
       call err_exit('unit 5 is for standard input')
    else if(io_output.ne.6.and.io_output.ne.0.and.myid.eq.0) then
       io_output=get_iounit()
       open(unit=io_output,access='sequential',action='write',file=file_run,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) then
          call err_exit('cannot open output file '//file_run)
       end if
    end if

  end subroutine read_system

  subroutine checkAtom()
    if (.not.allocated(atoms%list)) call err_exit(TRIM(__FILE__)//integer_to_string(__LINE__)//": ATOMS section has not been defined!")
  end subroutine checkAtom
end module sim_system
