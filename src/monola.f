! -----------------------------------------------------------------
! subroutine monola
! reads the control-data from unit 4
! reads a starting configuration from unit 7
! calculates interaction table
! starts and controls the simulation
! -----------------------------------------------------------------
subroutine monola(file_in)
  use var_type,only:dp,default_string_length
  use const_math,only:onepi,twopi
  use const_phys,only:debroglie_factor
  use util_random,only:random
  use util_string,only:format_n
  use util_runtime,only:err_exit
  use util_timings,only:time_date_str,time_now
  use util_mp,only:mp_barrier
  use sim_particle,only:init_neighbor_list
  use sim_system
  use sim_cell
  use energy_pairwise,only:sumup
  use moves_simple,only:traxyz,rotxyz,Atom_traxyz,output_translation_rotation_stats
  use moves_volume,only:volume,prvolume,output_volume_stats
  use moves_cbmc,only:config,schedule,output_cbmc_stats
  use moves_ee,only:eesetup,eemove,ee_index_swap,expand,numcoeff,output_ee_stats
  use transfer_shared,only:opt_bias,lopt_bias,freq_opt_bias
  use transfer_swap,only:swap,cnt,output_swap_stats,acchem,bnchem
  use transfer_swatch,only:swatch,output_swatch_stats
  use prop_pressure,only:pressure
  implicit none

  character(LEN=*),intent(in)::file_in

  real,allocatable::pres(:),acvkjmol(:,:),acsolpar(:,:,:),acEnthalpy(:),acEnthalpy1(:),stdev1(:,:,:),sterr1(:,:,:)&
   ,errme1(:,:,:),vstart(:),vend(:),avv(:,:),acv(:,:),acvsq(:,:),aflv(:),acboxl(:,:),acboxa(:,:),acpres(:),acsurf(:)&
   ,acvolume(:),acnbox(:,:),acnbox2(:,:,:),asetel(:,:),acdens(:,:),molfra(:,:),dsq(:,:),stdev(:,:),dsq1(:,:,:),sterr(:,:)&
   ,errme(:,:),molvol(:),speden(:),flucmom(:),flucmom2(:),acvol(:),acvolsq(:)
  integer,allocatable::nminp(:),nmaxp(:),ncmt_list(:,:),ndist(:,:),mnbox(:,:)
  logical,allocatable::lratfix(:),lsolute(:)
  character(LEN=default_path_length),allocatable::file_ndis(:)

  character(LEN=default_path_length)::file_flt,file_hist
  character(LEN=default_path_length)::file_config,file_cell
  character(LEN=default_path_length)::fileout

  ! dimension statements for block averages ---
  character(LEN=default_string_length)::vname(nener)
  character(LEN=default_string_length)::enth,enth1

  ! variables added for GCMC histogram reweighting
  integer,parameter::fmax=1E6_dp

  integer::nummol,jerr
  integer::ntii
  integer::imax,itmax
  integer::n,nconfig,nentry
  real::time_prev,time_cur,vhist,eng_list(fmax)

  real::Temp_Energy,Temp_Mol_Vol

  integer::point_of_start,point_to_end

  integer::i,j,nblock,ibox,jbox,nnn,ii,itemp,itype,itype2,intg,imolty,ilunit,nbl,itel,ig,il,k,histtot,Temp_nmol
  integer::nvirial,zzz,steps,igrow
  real::starvir,stepvir,starviro
  real::acnp,acmove,v(nEnergy),press1,surf
  real::rm,temvol,setx,sety,setz,setel,temacd,temspd,dblock,dbl1
  real::ostwald,stdost,dummy,debroglie

  real:: enthalpy,enthalpy2,sigma2Hsimulation
  real:: inst_enth, inst_enth2, tmp,sigma2H ,Cp

  real:: ennergy,ennergy2,sigma2Esimulation
  real:: inst_energy, inst_energy2, sigma2E ,Cv

  real::binvir(maxvir,maxntemp),binvir2(maxvir,maxntemp),inside,bvirial
  real::HSP_T, HSP_LJ, HSP_COUL
  real::CED_T, CED_LJ, CED_COUL
  real::Heat_vapor_T,Heat_vapor_LJ ,Heat_vapor_COUL
  real::pdV

  real::molfrac,gconst,dielect

  logical::ovrlap,lratio,lratv,lprint,lmv,lrsave,lblock,lucall,ltfix,ltsolute

  integer::iprop,jblock

  integer::bin
  real::profile(1000)

!****************************************************************
! initializes some variables that otherwise cause errors      ***
! Matt McGrath, October 1, 2009                               ***
!****************************************************************
  nmolty1=0
  leemove=.false.
  iratipsw=0
  acipsw=0.0E0_dp
  do bin = 1,1000
     profile(bin) = 0.0E0_dp
  end do
  ! binstep = 0.05E0_dp
  ! lvirial2 = .false.

  ! KM for MPI
  ! only one processor at a time reads and writes data from files
  do i=1,numprocs
     if (myid.eq.i-1) then
        call readdat(file_in,lucall,nvirial,starvir,stepvir)
     end if
     call mp_barrier(groupid)
  end do

  allocate(nminp(ntmax),nmaxp(ntmax),ncmt_list(fmax,ntmax),ndist(0:nmax,ntmax),pres(nbxmax),acvkjmol(nener,nbxmax)&
   ,acsolpar(nprop1,nbxmax,nbxmax),acEnthalpy(nbxmax),acEnthalpy1(nbxmax),stdev1(nprop1,nbxmax,nbxmax)&
   ,sterr1(nprop1,nbxmax,nbxmax),errme1(nprop1,nbxmax,nbxmax),vstart(nbxmax),vend(nbxmax),avv(nener,nbxmax)&
   ,acv(nener,nbxmax),acvsq(nener,nbxmax),aflv(nbxmax),acboxl(nbxmax,3),acboxa(nbxmax,3),acpres(nbxmax),acsurf(nbxmax)&
   ,acvolume(nbxmax),acnbox(nbxmax,ntmax),acnbox2(nbxmax,ntmax,20),mnbox(nbxmax,ntmax),asetel(nbxmax,ntmax)&
   ,acdens(nbxmax,ntmax),molfra(nbxmax,ntmax),dsq(nprop,nbxmax),stdev(nprop,nbxmax),dsq1(nprop1,nbxmax,nbxmax)&
   ,sterr(nprop,nbxmax),errme(nprop,nbxmax),lratfix(ntmax),lsolute(ntmax),molvol(nbxmax),speden(nbxmax),flucmom(nbxmax)&
   ,flucmom2(nbxmax),acvol(nbxmax),acvolsq(nbxmax),file_ndis(ntmax),stat=jerr)
  if (jerr.ne.0) then
     call err_exit(__FILE__,__LINE__,'monola: allocation failed',jerr)
  end if

  DO iprop=1,nprop
     DO ibox=1,nbxmax
        DO jblock=1,blockm
           baver(iprop,ibox,jblock)=0.0E0_dp
        end do
     end do
  end do

  ! KM for MPI
  ! check that if lneigh or lgaro numprocs .eq. 1
  if (lneigh.and.numprocs.ne.1) then
     call err_exit(__FILE__,__LINE__,'Cannot run on more than 1 processor with neighbor list!!',myid+1)
  end if
  if (lgaro.and.numprocs.ne.1) then
     call err_exit(__FILE__,__LINE__,'Cannot run on more than 1 processor with lgaro = .true.!!',myid+1)
  end if

  ! kea don't stop for lgaro
  ! if (lchgall .and. (.not. lewald).and.(.not.lgaro)) then
  ! call err_exit(__FILE__,__LINE__,'lchgall is true and lewald is false. Not checked for accuracy!',myid+1)
  ! end if

  vname(1:10)=(/' Total energy',' Inter LJ    ',' Bond bending',' Torsion     ',' Intra LJ    ',' External pot',' Stretch     ',' Coulomb     ',' Tail  LJ    ',' Fluc Q      '/)
  enth      = ' Enthalpy Inst. Press '
  enth1     = ' Enthalpy Ext.  Press '

  ! SETTING UP ARRAYS FOR ANALYSYS PURPOSE

  ! JLR 11-11-09
  ! do not call analysis if you set ianalyze to be greater than number of cycles
  ! KM 01/10 remove analysis
  ! if (ianalyze.le.nstep) then
  ! call analysis(0)
  ! end if
  ! END JLR 11-11-09

  ! set up initial linkcell
  if (licell) then
     call build_linked_cell()
  end if

  if ( lneigh ) then
     ! call for initial set-up of the near-neighbour bitmap ***
     call init_neighbor_list()
  end if

  ! set up thermodynamic integration stuff
  if (lmipsw) call ipswsetup

  if (.not.lmipsw) then
     lstagea = .false.
     lstageb = .false.
     lstagec = .false.
  end if

  ! set up expanded ensemble stuff
  if (lexpee) then
     call eesetup
     if (lmipsw) call err_exit(__FILE__,__LINE__,'not for BOTH lexpee AND lmipsw',myid+1)
  end if

  ! KM for MPI
  ! only processor 0 opens files
  if (myid.eq.0) then
     write(file_cell,'("cell_param",I1.1,A,".dat")') run_num,suffix
     open(13,file=file_cell, status='unknown')
     close(13)
  end if

  ! setup files for histogram reweighting
  ! KM fom MPI
  ! will need to check this file I/O if want to run grand canonical in parallel
  if(lgrand) then
     if (myid.eq.0) then
        write(file_flt,'("nfl",I1.1,A,".dat")') run_num,suffix
        write(file_hist,'("his",I1.1,A,".dat")') run_num,suffix

        do i=1,nmolty
           write(file_ndis(i),'("n",I2.2,"dis",I1.1,A,".dat")') i,run_num,suffix
        end do

        open(50, file = file_flt, status='unknown')
        close(50)

        open(51,file=file_hist,status='unknown')
        write(51,'(f8.4,2x,i5,2x,g15.5,3f12.3)')  temp, nmolty, (temp*log(B(i)),i=1,nmolty), boxlx(1), boxly(1), boxlz(1)
        close(51)
     end if

     ! extra zero accumulator for grand-canonical ensemble
     nentry = 0
     nconfig = 0
     do itmax = 1,ntmax
        nminp(itmax) = 1E6_dp
        nmaxp(itmax) = -1E6_dp
        do imax=0,nmax
           ndist(imax,itmax) = 0
        end do
     end do
  end if

  ! zero accumulators ***
  do i = 1,11
     do ibox = 1,nbox-1
        do jbox = ibox+1, nbox
           acsolpar(i,ibox,jbox)=0.0E0_dp
        end do
     end do
  end do

  do i=1,nbox
     do j=1,nener
        acv(j,i) = 0.0E0_dp
        acvsq(j,i) = 0.0E0_dp
        acvkjmol(j,i) = 0.0E0_dp
     end do
     aflv(i) = 0.0E0_dp
     do j = 1, nmolty
        solcount(i,j) = 0
        avsolinter(i,j) = 0.0E0_dp
        avsolintra(i,j) = 0.0E0_dp
        avsolbend(i,j) = 0.0E0_dp
        avsoltor(i,j) = 0.0E0_dp
        avsolelc(i,j) = 0.0E0_dp
     end do
     acpres(i) = 0.0E0_dp
     acsurf(i) = 0.0E0_dp
     acEnthalpy(i) = 0.0E0_dp
     acEnthalpy1(i) = 0.0E0_dp
  end do

  inst_enth = 0.0E0_dp
  inst_enth2 = 0.0E0_dp
  inst_energy = 0.0E0_dp
  inst_energy2 = 0.0E0_dp

  acnp = 0.0E0_dp

  do ibox = 1, nbox
     acvol(ibox) = 0.0E0_dp
     acvolsq(ibox) = 0.0E0_dp
  end do

  acmove = 0.0E0_dp

  ! 2nd viral coefficient
  if (lvirial) then
     do i=1,maxvir
        do j = 1, ntemp
           binvir(i,j) = 0.0E0_dp
           binvir2(i,j) = 0.0E0_dp
        end do
     end do
     ! if ( lvirial2 ) then
     ! call virial2(binvir,binvir2,nvirial,starvir,stepvir)
     ! goto 2000
     ! end if
  end if
  ! permanent accumulators for box properties ***
  do ibox = 1, nbox
     do i = 1,3
        acboxl(ibox,i) = 0.0E0_dp
        acboxa(ibox,i) = 0.0E0_dp
     end do
  end do

  do i = 1, nmolty
     do ibox = 1, nbox
        asetel(ibox,i) = 0.0E0_dp
        mnbox(ibox,i)  = 0
        acdens(ibox,i) = 0.0E0_dp
        molfra(ibox,i) = 0.0E0_dp
        acnbox(ibox,i) = 0.0E0_dp
     end do
     if ( lexpand(i) ) then
        do itype = 1, numcoeff(i)
           do ibox = 1,2
              acnbox2(1,i,itype) = 0.0E0_dp
              acnbox2(2,i,itype) = 0.0E0_dp
           end do
        end do
     end if
  end do

  ! For the atom displacements
  Abstrax = 0.0E0_dp
  Abstray = 0.0E0_dp
  Abstraz = 0.0E0_dp
  Abntrax = 0.0E0_dp
  Abntray = 0.0E0_dp
  Abntraz = 0.0E0_dp

  acvolume = 0.0E0_dp

  ! accumulators for fluctuating charge performance
  do i = 1, nmolty
     do j = 1,nbox
        bnflcq(i,j) = 0.0E0_dp
        bnflcq2(i,j) = 0.0E0_dp
        bsflcq(i,j) = 0.0E0_dp
        bsflcq2(i,j) = 0.0E0_dp
     end do
  end do
  ! accumulators for block averages ---
  do i = 1, nprop
     do j = 1, nbox
        naccu(i,j) = 0.0E0_dp
        accum(i,j) = 0.0E0_dp
        nccold(i,j) = 0.0E0_dp
        bccold(i,j) = 0.0E0_dp
        dsq(i,j) = 0.0E0_dp
     end do
  end do
  do i = 1, nprop1
     do ibox = 1,nbox-1
        do jbox = ibox+1,nbox
           naccu1(i,ibox,jbox) = 0.0E0_dp
           accum1(i,ibox,jbox) = 0.0E0_dp
           nccold1(i,ibox,jbox) = 0.0E0_dp
           bccold1(i,ibox,jbox) = 0.0E0_dp
           dsq1(i,ibox,jbox) = 0.0E0_dp
        end do
     end do
  end do

  nblock = 0

! -----------------------------------------------------------------

  if (lneighbor) then
     do i = 1,maxneigh
        do ii = 1,nmax
           neighbor(i,ii) = 0
        end do
     end do
  end if

  ! calculate initial energy and check for overlaps ***
  do ibox=1,nbox
     call sumup(ovrlap,v,ibox,.false.)

     vbox(:,ibox) = v
     vbox(2,ibox)  = v(2)
     vbox(3,ibox)   = v(3)
     vbox(4,ibox)  = v(4)
     vbox(5,ibox)    = v(5)
     vbox(6,ibox)   = v(6)
     vbox(7,ibox)     = v(7)
     vbox(9,ibox)    = v(9)
     vbox(8,ibox)  = v(8)
     vbox(11,ibox)  = v(11)
     vbox(12,ibox) = vipsw
     vbox(13,ibox) = vwellipsw
     ! kea
     vbox(10,ibox) = v3garo

     if( ovrlap ) then
        call err_exit(__FILE__,__LINE__,'overlap in initial configuration',myid+1)
     end if
     vstart(ibox) = vbox(1,ibox)
     if (myid.eq.0) then
        write(io_output,*)
        write(io_output,*) 'box  ',ibox,' initial v   = ', vbox(1,ibox)
     end if
     ! calculate initial pressure ***
     call pressure( press1, surf, ibox )
     if (myid.eq.0) then
        write(io_output,"(' surf. tension :   box',i2,' =',f14.5)") ibox, surf
        write(io_output,"(' pressure check:   box',i2,' =',f14.2)") ibox, press1
     end if
  end do

  if (myid.eq.0) then
     write(io_output,*)
     write(io_output,*) '+++++ start of markov chain +++++'
     write(io_output,*)
     write(io_output,*)  'Cycle   Total   Energy    Boxlength   Pressure  Molecules'
     ! set up info at beginning of fort.12 for analysis
     write(12,*) nstep,nmolty,(masst(i),i=1,nmolty)
  end if
!************************************************************
! loops over all cycles and all molecules                  **
!************************************************************

  if (lneighbor) then
     write(21,*) 'ii:',nnstep,(neigh_cnt(i),i=1,nchain)
  end if

  time_prev = time_now()

  do nnn = 1, nstep
     tmcc = nnstep + nnn
     do ii = 1, nchain
        ! write(io_output,*) 'nstep',(nnn-1)*nchain+ii
        ! select a move-type at random ***
        rm = random(-1)
        ! write(io_output,*) 'tmcc, random: ',tmcc,rm
        ! special ensemble dependent moves ###
        if  (rm .le. pmvol) then
           ! volume move ---
           if ( lnpt ) then
              call prvolume
           else
              call volume
           end if
        else if (rm .le. pmswat) then
           ! CBMC switch move ---
           call swatch
        else if ( rm .le. pmswap ) then
           ! swap move for linear and branched molecules ---
           call swap
        else if ( rm .le. pmcb ) then
           ! configurational bias move ---
           call config
        else if ( rm .le. pmflcq ) then
           ! displacement of fluctuating charges ---
           call flucq(2,0)
        else if (rm .le. pmexpc ) then
           ! expanded-ensemble move ---
           call expand
        else if (rm .le. pmexpc1 ) then
           ! new expanded-ensemble move ---
           ! call expand
           if (random(-1).le.eeratio) then
              call ee_index_swap
           else
              call eemove
           end if
        else if ( rm .le. pm_atom_tra) then
           rm = 3.0E0_dp * random(-1)
           if ( rm .le. 1.0E0_dp ) then
              call Atom_traxyz (.true.,.false.,.false.)
           else if ( rm .le. 2.0E0_dp ) then
              call Atom_traxyz (.false.,.true.,.false.)
           else
              call Atom_traxyz (.false.,.false.,.true.)
           end if
        else if ( rm .le. pmtra ) then
           ! translational move in x,y, or z direction ---
           rm = 3.0E0_dp * random(-1)
           if ( rm .le. 1.0E0_dp ) then
              call traxyz(.true.,.false.,.false.)
           else if ( rm .le. 2.0E0_dp ) then
              call traxyz(.false.,.true.,.false.)
           else
              call traxyz(.false.,.false.,.true.)
           end if
        else
           ! rotation around x,y, or z axis move --
           rm = 3.0E0_dp * random(-1)
           if ( rm .le. 1.0E0_dp ) then
              call rotxyz(.true.,.false.,.false.)
           else if ( rm .le. 2.0E0_dp ) then
              call rotxyz(.false.,.true.,.false.)
           else
              call rotxyz(.false.,.false.,.true.)
           end if
        end if

        acmove = acmove + 1.0E0_dp
        ! accumulate probability of being in an expanded ensemble
        ! state
        if (lexpee) then
           ee_prob(mstate) = ee_prob(mstate)+1
        end if

        ! calculate instantaneous values ***
        ! accumulate averages ***
        do ibox=1,nbox
           do itype = 1, nmolty
              acnbox(ibox,itype) = acnbox(ibox,itype) +  real(ncmt(ibox,itype),dp)
              if ( lexpand(itype) ) then
                 do itype2 = 1, numcoeff(itype)
                    acnbox2(ibox,itype,itype2) =  acnbox2(ibox,itype,itype2) + real(ncmt2(ibox,itype,itype2),dp)
                    ! write(io_output,*) '1:',acnbox2(ibox,itype,itype2)
                 end do
                 ! write(io_output,*) '2:', acnbox(ibox,itype)
              end if
           end do

           acv(1,ibox)    = acv(1,ibox)   + vbox(1,ibox)
           acvsq(1,ibox)  = acvsq(1,ibox) + vbox(1,ibox)**2
           acv(2,ibox)    = acv(2,ibox)   + vbox(2,ibox)
           acvsq(2,ibox)  = acvsq(2,ibox) + vbox(2,ibox)**2
           acv(3,ibox)    = acv(3,ibox)   + vbox(6,ibox)
           acvsq(3,ibox)  = acvsq(3,ibox) + vbox(6,ibox)**2
           acv(4,ibox)    = acv(4,ibox)   + vbox(7,ibox)
           acvsq(4,ibox)  = acvsq(4,ibox) + vbox(7,ibox)**2
           acv(5,ibox)    = acv(5,ibox)   + vbox(4,ibox)
           acvsq(5,ibox)  = acvsq(5,ibox) + vbox(4,ibox)**2
           acv(6,ibox)    = acv(6,ibox)   + vbox(9,ibox)
           acvsq(6,ibox)  = acvsq(6,ibox) + vbox(9,ibox)**2
           acv(7,ibox)    = acv(7,ibox)   + vbox(5,ibox)
           acvsq(7,ibox)  = acvsq(7,ibox) + vbox(5,ibox)**2
           acv(8,ibox)    = acv(8,ibox)   + vbox(8,ibox)
           acvsq(8,ibox)  = acvsq(8,ibox) + vbox(8,ibox)**2
           acv(9,ibox)    = acv(9,ibox)   + vbox(3,ibox)
           acvsq(9,ibox)  = acvsq(9,ibox) + vbox(3,ibox)**2
           acv(10,ibox)    = acv(10,ibox)   + vbox(11,ibox)
           acvsq(10,ibox)  = acvsq(10,ibox) + vbox(11,ibox)**2
           ! KEA added 17 for v3garo
           acv(17,ibox)    = acv(17,ibox) + vbox(10,ibox)
           acvsq(17,ibox)  = acvsq(17,ibox) + vbox(10,ibox)**2

           ! leftover from Bin, not currently used
           if ( ldielect ) then
              acv(11,ibox) = acv(11,ibox)+dipolex(ibox)
              acvsq(11,ibox) = acvsq(11,ibox)+dipolex(ibox)**2
              acv(12,ibox) = acv(12,ibox)+dipoley(ibox)
              acvsq(12,ibox) = acvsq(12,ibox)+dipoley(ibox)**2
              acv(13,ibox) = acv(13,ibox)+dipolez(ibox)
              acvsq(13,ibox) = acvsq(13,ibox)+dipolez(ibox)**2
              acvsq(14,ibox) = acvsq(11,ibox) + acvsq(12,ibox) + acvsq(13,ibox)
              acv(15,ibox) = acv(15,ibox)+sqrt(dipolex(ibox)* dipolex(ibox)+dipoley(ibox)*dipoley(ibox)+ dipolez(ibox)*dipolez(ibox))
           end if

           if (lsolid(ibox) .and. .not. lrect(ibox)) then
              temvol = cell_vol(ibox)
           else
              if ( lpbcz ) then
                 temvol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
              else
                 temvol = boxlx(ibox)*boxly(ibox)
              end if
           end if

           ! KMB/KEA Energy in kJ/mol
           Temp_nmol = 0
           do itype=1,nmolty
              Temp_nmol =   Temp_nmol + ncmt(ibox,itype)
           end do
           Temp_Energy  = (vbox(1,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(1,ibox) = acvkjmol(1,ibox) + Temp_Energy
           Temp_Energy  = (vbox(2,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(2,ibox) = acvkjmol(2,ibox) + Temp_Energy
           Temp_Energy  = (vbox(6,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(3,ibox) = acvkjmol(3,ibox) + Temp_Energy
           Temp_Energy  = (vbox(7,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(4,ibox) = acvkjmol(4,ibox) + Temp_Energy
           Temp_Energy  = (vbox(4,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(5,ibox) = acvkjmol(5,ibox) + Temp_Energy
           Temp_Energy  = (vbox(9,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(6,ibox) = acvkjmol(6,ibox) + Temp_Energy
           Temp_Energy  = (vbox(5,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(7,ibox) = acvkjmol(7,ibox) + Temp_Energy
           Temp_Energy  = (vbox(8,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(8,ibox) = acvkjmol(8,ibox) + Temp_Energy
           Temp_Energy  = (vbox(3,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(9,ibox) = acvkjmol(9,ibox) + Temp_Energy
           Temp_Energy  = (vbox(11,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(10,ibox) = acvkjmol(10,ibox) + Temp_Energy
           Temp_Energy  = (vbox(10,ibox)/Temp_nmol)*0.00831451E0_dp
           acvkjmol(17,ibox) = acvkjmol(17,ibox) + Temp_Energy

           if ( lnpt ) then
              ! acv(16,ibox) = acv(16,ibox) + vbox(1,ibox)*boxlx(ibox)*boxly(ibox)*boxlz(ibox)
              acv(16,ibox) = acv(16,ibox) + vbox(1,ibox)*temvol
           end if

           if (lsolid(ibox) .and. .not. lrect(ibox)) then
              boxlx(ibox) = cell_length(ibox,1)
              boxly(ibox) = cell_length(ibox,2)
              boxlz(ibox) = cell_length(ibox,3)

              acboxa(ibox,1) = acboxa(ibox,1) + cell_ang(ibox,1)
              acboxa(ibox,2) = acboxa(ibox,2) + cell_ang(ibox,2)
              acboxa(ibox,3) = acboxa(ibox,3) + cell_ang(ibox,3)
           end if

           acboxl(ibox,1) = acboxl(ibox,1) + boxlx(ibox)
           acboxl(ibox,2) = acboxl(ibox,2) + boxly(ibox)
           acboxl(ibox,3) = acboxl(ibox,3) + boxlz(ibox)

           acvol(ibox) = acvol(ibox) + temvol
           acvolsq(ibox) = acvolsq(ibox) + temvol*temvol
           acvolume(ibox) = acvolume(ibox) + temvol

           do itype = 1, nmolty
              acdens(ibox,itype) = acdens(ibox,itype) + ncmt(ibox,itype) / temvol
              if ( nchbox(ibox) .gt. 0 ) then
                 molfra(ibox,itype) = molfra(ibox,itype) + real(ncmt(ibox,itype),dp) / real(nchbox(ibox),dp)
              end if
           end do
	    end do


        if (lstop) then
	       if ( int(acmove) .ge. nstep ) goto 101
	    end if

        If  ( lnpt.and..not.lgibbs ) then
           ibox = 1
           tmp= vbox(1,ibox) + express(ibox) * ( boxlx(ibox) *boxly(ibox)*boxlz(ibox))
           ! write(49,*) nnn, tmp
           inst_enth=inst_enth+ tmp
           inst_enth2=inst_enth2+(tmp*tmp)
        end if

        If  (.not. lnpt .and..not.lgibbs ) then
           ibox = 1
           tmp= vbox(1,ibox)
           ! write(49,*) nnn, tmp
           inst_energy=inst_energy + tmp
           inst_energy2=inst_energy2+ (tmp*tmp)
        end if

        ! collect histogram data (added 8/30/99)
        if(lgrand) then
           ibox = 1
           vhist = vbox(2,ibox) + vbox(8,ibox) + vbox(11,ibox)
           nconfig = nconfig + 1
           ! KM for MPI
           ! check this if want to run grand canonical
           if (mod(nconfig,ninstf).eq.0.and.myid.eq.0) then
              open(50, file = file_flt, status='old', position='append')
              !	write(50, '(i10,5x,i7,5x,g15.6)') nconfig,(ncmt(ibox,i), i=1,nmolty),vhist

              ! Formatting removed until I can figure out how to get <n>i7 to work
              write(50, * ) nconfig,(ncmt(ibox,i), i=1,nmolty),vhist
              close(50)
           end if

           if(mod(nconfig,ninsth).eq.0.and.nconfig.gt.nequil) then
              do imolty = 1, nmolty
                 nminp(imolty) = min(nminp(imolty),ncmt(ibox,imolty))
                 nmaxp(imolty) = max(nmaxp(imolty),ncmt(ibox,imolty))
              end do
              nentry = nentry + 1

              do imolty=1,nmolty
                 ncmt_list(nentry,imolty) = ncmt(ibox,imolty)
                 ndist(ncmt(ibox,imolty),imolty) = ndist(ncmt(ibox,imolty),imolty) + 1
              end do

              eng_list(nentry) = vhist
           end if

           if (mod(nconfig,ndumph).eq.0.and.myid.eq.0) then
              open(51,file = file_hist,status='old', position='append')
              do i=1,nentry
                 write(51, * ) (ncmt_list(i,imolty), imolty=1,nmolty), eng_list(i)
              end do
              close(51)
              nentry = 0

              do imolty=1,nmolty
                 open(52, file=file_ndis(imolty), status='unknown')
                 do n=nminp(imolty),nmaxp(imolty)
                    write(52,*) n,ndist(n,imolty)
                 end do
                 close(52)
              end do
           end if
        end if
     end do
!*********************************************
! ends loop over chains                     **
!*********************************************

     ! perform periodic operations  ***
     time_cur = time_now()
     if ( time_cur - time_prev .gt. checkpoint_interval ) then
        lrsave = .true.
        time_prev = time_cur
     else
        lrsave = .false.
     end if

     if (any(lopt_bias).and.mod(nnn,freq_opt_bias).eq.0) then
        call opt_bias
     end if

     if ( lgibbs .or. lnpt .and. (.not. lvirial) ) then
        do ibox = 1,nbox
           if ( lpbcz ) then
              if (lsolid(ibox) .and. .not. lrect(ibox).and. myid.eq.0) then
                 write(12,FMT='(7E13.5,'//format_n(nmolty,"i5")//')') hmat(ibox,1) ,hmat(ibox,4),hmat(ibox,5) ,hmat(ibox,7),hmat(ibox,8) ,hmat(ibox,9),vbox(1,ibox), (ncmt(ibox,itype),itype=1,nmolty)

                 open(13,file = file_cell,status='old', position='append')
                 write(13,'(i8,6f12.4)') nnn+nnstep, cell_length(ibox,1)/Num_cell_a, cell_length(ibox,2)/Num_cell_b, cell_length(ibox,3)/Num_cell_c, cell_ang(ibox,1)*180.0E0_dp/onepi, cell_ang(ibox,2)*180.0E0_dp/onepi, cell_ang(ibox,3)*180.0E0_dp/onepi
                 close(13)

                 ! write(13,'(i8,3f12.4)') nnn,cell_ang(ibox,1),cell_ang(ibox,2),cell_ang(ibox,3)
              else
                 ! do ibox = 1, nbox
                 if (myid.eq.0) then
                    write(12,'(4E13.5,'//format_n(nmolty,"i5")//')')boxlx(ibox),boxly(ibox) ,boxlz(ibox),vbox(1,ibox), (ncmt(ibox,itype),itype=1,nmolty)
                 end if
                 ! end do
              end if
           else
              ! do ibox = 1, nbox
              if (myid.eq.0) then
                 write(12,'(2E12.5,'//format_n(nmolty,"i4")//')') boxlx(ibox)*boxly(ibox) ,vbox(1,ibox),(ncmt(ibox,itype),itype=1,nmolty)
              end if
              ! end do
           end if
        end do
     end if

     if (lucall) then
        call err_exit(__FILE__,__LINE__,'not recently checked for accuracy',myid+1)
        ! do j = 1,nmolty
        !    if ( ucheck(j) .gt. 0 ) then
        !       call chempt(bsswap,j,ucheck(j))
        !    end if
        ! end do
     end if

     if ( mod(nnn,iratp) .eq. 0 ) then
        ! calculate pressure ***
        acnp = acnp + 1.0E0_dp
        do ibox = 1, nbox
           call pressure ( press1, surf, ibox )
           ! write(io_output,*) 'control pressure', press1
           pres(ibox) = press1
           acpres(ibox) = acpres(ibox) + press1
           acsurf(ibox) = acsurf(ibox) + surf
        end do

        ! Enthalpy calculation
        do ibox = 1,nbox
           Temp_nmol = 0
           do itype=1,nmolty
              Temp_nmol =   Temp_nmol + ncmt(ibox,itype)
           end do
           ! Volume in m3/mol, energies in kJ/mol
           Temp_Mol_Vol = temvol/Temp_nmol*0.6022E-06_dp
           Temp_Energy  = (vbox(1,ibox)/Temp_nmol)*0.00831451E0_dp
           acEnthalpy(ibox) = acEnthalpy(ibox) + Temp_Energy + pres(ibox)*Temp_Mol_Vol
           acEnthalpy1(ibox) = acEnthalpy1(ibox) + Temp_Energy + (express(ibox)/7.2429E-5_dp)*Temp_Mol_Vol
        end do

        ! cannot calculate a heat of vaporization for only one box,
        ! and some compilers choke because Heat_vapor_T will not be
        ! defined if nbox == 1
        if(lgibbs) then
           do ibox = 1,nbox-1
              do jbox = ibox+1,nbox
                 ! WRITE(io_output,*) 'ieouwfe ',ibox,jbox
                 call calcsolpar(pres,Heat_vapor_T,Heat_vapor_LJ, Heat_vapor_COUL,pdV, CED_T,CED_LJ,CED_COUL, HSP_T,HSP_LJ, HSP_COUL,ibox,jbox)

                 ! Heat of vaporization
                 acsolpar(1,ibox,jbox)= acsolpar(1,ibox,jbox)+Heat_vapor_T
                 acsolpar(2,ibox,jbox)= acsolpar(2,ibox,jbox)+Heat_vapor_LJ
                 acsolpar(3,ibox,jbox)= acsolpar(3,ibox,jbox)+Heat_vapor_COUL
                 acsolpar(4,ibox,jbox)= acsolpar(4,ibox,jbox)+CED_T
                 acsolpar(5,ibox,jbox)= acsolpar(5,ibox,jbox)+CED_LJ
                 acsolpar(6,ibox,jbox)= acsolpar(6,ibox,jbox)+CED_COUL
                 acsolpar(7,ibox,jbox)= acsolpar(7,ibox,jbox)+HSP_T
                 acsolpar(8,ibox,jbox)= acsolpar(8,ibox,jbox)+HSP_LJ
                 acsolpar(9,ibox,jbox)= acsolpar(9,ibox,jbox)+HSP_COUL
                 ! acsolpar(10,ibox,jbox) = acsolpar(10,ibox,jbox)+DeltaU_Ext
                 acsolpar(11,ibox,jbox) = acsolpar(11,ibox,jbox)+pdV
              end do
           end do
        end if
     end if

     ! calculate the integrand of thermosynamic integration
     if (lmipsw.and.(mod(nnn,iratipsw).eq.0)) then
        acipsw = acipsw+1.0E0_dp
        call deriv(1)
        acdvdl = acdvdl+dvdl
     end if

     ! Add a call for subroutine to compute the
     lratio = .false.
     lratv = .false.
     lprint = .false.
     lmv = .false.
     lblock = .false.

     if ( mod(nnn,iratio) .eq. 0 ) then
        lratio = .true.
     end if

     if ( lgibbs .or. lnpt ) then
        if ( mod(nnn,iratv) .eq. 0 .and. pmvol .gt. 0.0E0_dp ) then
           lratv = .true.
        end if
     end if

     if ( mod(nnn,iprint) .eq. 0 ) then
        lprint = .true.
     end if

     if ( mod(nnn,imv) .eq. 0 ) then
        if ( lvirial ) then
           call virial(binvir,binvir2,nvirial,starvir,stepvir)
        else
           lmv = .true.
        end if
     end if

     ! JLR 11-11-09
     ! do not call analysis if ianalyze is greater than number of cycles
     ! KM 01/10 remove analysis
     !	 if(mod(nnn,ianalyze).eq.0) then
     !	    call analysis(1)
     ! end if
     ! END JLR 11-11-09

     do intg = 1, nchain
        ibox = nboxi(intg)
        imolty = moltyp(intg)
        ! accumulate m-n-box and m-s-e-t-e-l ***
        ! only count the main chain - not the hydrogens
        ilunit = nugrow(imolty)
        setx = rxu(intg,1) - rxu(intg,ilunit)
        sety = ryu(intg,1) - ryu(intg,ilunit)
        setz = rzu(intg,1) - rzu(intg,ilunit)
        setel = setx*setx + sety*sety + setz*setz
        ! if ( imolty .eq. 2 ) then
        ! write(??,*) imolty,setel
        ! end if
        mnbox( ibox, imolty ) = mnbox( ibox, imolty ) + 1
        asetel( ibox, imolty ) = asetel( ibox, imolty ) + setel
     end do

     if ( mod(nnn,iblock) .eq. 0 ) then
        lblock = .true.
     end if

     ltsolute = .false.
     ltfix    = .false.

     do imolty = 1, nmolty
        if (pmfix(imolty).gt.0.0001E0_dp .and. mod(nnn,iupdatefix).eq.0) then
           ltfix = .true.
           lratfix(imolty) = .true.
        else
           lratfix(imolty) = .false.
        end if

        if (mod(nnn,isolute(imolty)).eq.0) then
           ltsolute = .true.
           lsolute(imolty) = .true.
        else
           lsolute(imolty) = .false.
        end if
     end do

     if (lratio .or. lratv .or. lprint .or. lmv .or. lrsave .or. lblock .or. ltfix .or. ltsolute) then
        call monper(acv,acpres,acsurf,acvolume,molfra,mnbox,asetel,acdens,acmove,acnp,pres,nbox,nnn,nblock,lratio,lratv,lprint,lmv,lrsave,lblock,lratfix,lsolute,acsolpar,acEnthalpy,acEnthalpy1)
     end if
     ! not currently used
     if (ldielect.and.(mod(nnn,idiele).eq.0).and.myid.eq.0) then
        dielect = acvsq(14,ibox)/acmove

        ! If you really want this quantity comment should be taken out**
        ! write(14,*) nnn,6.9994685465110493E5_dp*dielect*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))

        dielect = acvsq(14,ibox)/acmove   -(acv(11,ibox)/acmove)**2 -(acv(12,ibox)/acmove)**2 - (acv(13,ibox)/acmove)**2
        ! use fort.27 to calculate dielectric constant
        ! write(15,*) nnn,6.9994685465110493E5_dp*dielect*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
        ! write(16,*) nnn,acv(11,ibox)/acmove, acv(12,ibox)/acmove,acv(13,ibox)/acmove
     end if
     if ( lnpt .and. nmolty .eq. 1 ) then
        ! output the fluctuation information
        if (mod(nnn,idiele) .eq. 0.and.myid.eq.0) then
           write(14,*) nnn,acvol(ibox)/acmove
           write(15,*) nnn,acvolsq(ibox)/acmove
           write(16,*) nnn,acv(1,ibox)/acmove
           write(17,*) nnn,acvsq(1,ibox)/acmove
           write(18,*) nnn,acvol(ibox)*acv(1,ibox)/(acmove*acmove)
           write(19,*) nnn,acv(16,ibox)/acmove - acvol(ibox) *acv(1,ibox)/(acmove*acmove)
        end if
     end if

     ! set idiele = 1 to print every cycle
     if ( ldielect .and. mod(nnn,idiele).eq. 0.and.myid.eq.0) then
        do ibox = 1,nbox
           write(27,*) dipolex(ibox),dipoley(ibox),dipolez(ibox)
        end do
     end if

     ! if ( mod(nnn,idiele) .eq. 0 ) write(25,*) nnn+nnstep,vbox(1,1)

     ! ibox = 1
     ! imolty = 1
     ! do i = 1,nchain
     !    if ( nboxi(i) .eq. ibox ) then
     !       if ( moltyp(i) .eq. imolty ) then
     !          bin = aint(zcm(i)/binstep) + 1
     !          temvol = boxlx(ibox)*boxly(ibox)*binstep
     !          profile(bin) = profile(bin)+1.0E0_dp/temvol
     !       end if
     !    end if
     ! end do

     ! Residual Heat capacity ---
     if(mod(nnn,iheatcapacity) .eq. 0) then

        if(lnpt.and..not.lgibbs) then
           enthalpy= inst_enth/real(nchain*nnn,dp)
           enthalpy2= inst_enth2/real(nchain*nnn,dp)
           sigma2Hsimulation=(enthalpy2)-(enthalpy*enthalpy)
           sigma2H=sigma2Hsimulation*(6.022E23_dp)*((1.38066E-23_dp)**2) / real(nchain,dp) !(J2/mol)
           Cp=sigma2H/((1.38066E-23_dp)*(temp**2))
           if (myid.eq.0) then
              write(56,'(I12,F18.6,F18.2,F18.6)')nnn,Cp,enthalpy2, enthalpy
           end if

        else if( .not.lnpt .and..not.lgibbs) then
           ennergy= inst_energy/real(nchain*nnn,dp)
           ennergy2= inst_energy2/real(nchain*nnn,dp)
           sigma2Esimulation=(ennergy2)-(ennergy*ennergy)
           sigma2E=sigma2Esimulation*(6.022E23_dp)*((1.38066E-23_dp)**2) / real(nchain,dp) !(J2/mol)
           Cv=sigma2E/((1.38066E-23_dp)*(temp**2))
           if (myid.eq.0) then
              write(55,'(I12,F18.6,F18.2,F18.6)')nnn,Cv,ennergy2, ennergy
           end if
        end if
     end if
  end do
!********************************************
! ends the loop over cycles                **
!********************************************
  call cnt

  ! do bin = 1,1000
  ! write(26,*) binstep*(real(bin,dp)-0.5E0_dp),profile(bin)/nstep
  ! end do

101 continue

  if (lneighbor) then
     write(21,*) 'ii:',tmcc,(neigh_cnt(i),i=1,nchain)
  end if
  if (myid.eq.0) then
     write(io_output,*)
     write(io_output,*) '+++++ end of markov chain +++++'

     call output_translation_rotation_stats(io_output)

     if ( pmcb .gt. 0.0E0_dp ) then
        call output_cbmc_stats(io_output)
     end if
     ! write some information about volume performance ***
     if ( lgibbs .or. lnpt) then
        call output_volume_stats(io_output)
     end if

     call output_swap_stats(io_output)
     call output_swatch_stats(io_output)

     write(io_output,*)
     write(io_output,*)    '### Charge Fluctuation  ###'
     write(io_output,*)

     do i = 1, nmolty
        do j = 1,nbox
           bnflcq2(i,j) = bnflcq2(i,j) + bnflcq(i,j)
           bsflcq2(i,j) = bsflcq2(i,j) + bsflcq(i,j)
           if (bnflcq2(i,j) .gt. 0.5E0_dp) then
              write(io_output,*) 'molecule typ =',i,'  box =',j
              bsflcq2(i,j) = bsflcq2(i,j)/bnflcq2(i,j)
              write(io_output,"(' attempts =',f8.1,'   ratio =',f6.3, '   max.displ. =',e11.4)") bnflcq2(i,j),bsflcq2(i,j),rmflcq(i,j)
           end if
        end do
     end do

     call output_ee_stats(io_output)

     write(io_output,"(/,'New Biasing Potential')")
     do i=1,nmolty
        write(io_output,"(/,'molecule ',I2,': ',$)") i
        do j=1,nbox
           write(io_output,"(G16.9,1X,$)") eta2(j,i)
        end do
     end do
     write(io_output,*)
  end if

  do ibox=1,nbox
     if ( ldielect ) then
        ! store old dipole moment
        call dipole(ibox,2)
     end if

     ! checks final value of the potential energy is consistent ***
     call sumup(ovrlap,v,ibox,.false.)
     vend(ibox) = v(1)

     ! need to check
     if (myid.eq.0) then
        if ( abs(v(1) - vbox(1,ibox)) .gt. 0.0001) then
           write(io_output,*) '### problem with energy ###  box ',ibox
           write(io_output,*) ' Total energy: ',v(1),vbox(1,ibox),v(1)-vbox(1,ibox)
        end if
        if ( abs(v(2) - vbox(2,ibox)) .gt. 0.000001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Inter mol.en.: ',v(2),vbox(2,ibox)
           if (lsolid(ibox) .and. .not. lrect(ibox)) then
              write(io_output,*)'You might check the cutoff wrt box widths'
              write(io_output,*) 'Normal PBC might be failing'
           end if
        end if
        if ( abs(v(3) - vbox(3,ibox)) .gt. 0.000001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Tail corr.en.: ',v(3),vbox(3,ibox)
        end if
        if ( abs(v(4) - vbox(4,ibox)) .gt. 0.000001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Intra mol.en.: ',v(4),vbox(4,ibox)
        end if
        if ( abs(v(5) - vbox(5,ibox)) .gt. 0.001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' bond vib. en.: ',v(5),vbox(5,ibox)
        end if
        if ( abs(v(6) - vbox(6,ibox)) .gt. 0.001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Bond ben.en.: ',v(6),vbox(6,ibox)
        end if
        if ( abs(v(7) - vbox(7,ibox)) .gt. 0.001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Torsion.en.: ',v(7),vbox(7,ibox)
        end if
        if ( abs(v(9) - vbox(9,ibox)) .gt. 0.0001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Externa.en.: ',v(9),vbox(9,ibox)
        end if
        if ( abs(v(8) - vbox(8,ibox)) .gt. 0.000001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Coulomb.en.: ',v(8),vbox(8,ibox)
        end if
        if ( abs(v(11) - vbox(11,ibox)) .gt. 0.0001) then
           write(io_output,*) '### problem  ###'
           write(io_output,*) ' Fluc Q en.: ',v(11),vbox(11,ibox)
        end if
        if ( abs(v3garo - vbox(10,ibox) ) .gt.0.001) then
           write(io_output,*) '### problem ###'
           write(io_output,*) ' 3-body en.: ',v3garo,vbox(10,ibox)
        end if
        if ( ldielect ) then
           if ( abs(dipolexo - dipolex(ibox)) .gt. 0.0001) then
              write(io_output,*) '### problem  ###'
              write(io_output,*) ' Dipole X: ',dipolexo,dipolex(ibox)
           end if
        end if
        if (lmipsw) then
           if (abs(vwellipsw-vbox(13,ibox)).gt.0.001) then
              write(io_output,*) '### problem  ###'
              write(io_output,*) ' well en.: ',vwellipsw,vbox(13,ibox)
           end if
        end if
     end if
  end do

  ! KM for MPI
  ! only processor 0 needs to calculate and write out averages, final config, etc
  if (myid.eq.0) then
     write(io_output,*)
     write(io_output,"(' vstart       =',3f24.10)") (vstart(i) ,i=1,nbox)
     write(io_output,"(' vend         =',3f24.10)") (vend(i)   ,i=1,nbox)
     write(io_output,"(' vbox         =',3f24.10)") (vbox(1,i)   ,i=1,nbox)
     write(io_output,*)

     ! normalize and write out presim results in fort.22 **
     if (lpresim) then
        if (counttot.eq.0) then
           write(21,*) counthist
        else
           write(21,*) counttot
        end if

        do j = 1, iring(1)
           do k = 1, iring(1)
              if (j.eq.k) goto 150
              histtot = 0
              do bin = 1, maxbin
                 hist(j,k,bin) = hist(j,k,bin) + 1.0E0_dp
                 histtot = histtot + hist(j,k,bin)
              end do

              do bin = 1, maxbin
                 hist(j,k,bin) = hist(j,k,bin) / histtot
                 write(21,*) bin,hist(j,k,bin)
              end do
150           continue
           end do
        end do
     end if

     ! put new distribution back into a file
     do imolty = 1, nmolty
        if (pmfix(imolty).gt.0) then
           if (counttot.eq.0) then
              write(21,*) counthist
           else
              write(21,*) counttot
           end if
           do j = 1, iring(1)
              do k = 1, iring(1)
                 if (j.eq.k) goto 160
                 do bin = 1, maxbin
                    write(21,*) bin,probf(j,k,bin)
                 end do
160              continue
              end do
           end do
        end if
     end do

     ! write out the final configuration for each box, Added by Neeraj 06/26/2006 3M ***
     do ibox = 1,nbox
        write(fileout,'("box",I1.1,"config",I1.1,A,".xyz")') ibox ,run_num,suffix
        open (200+ibox,FILE=fileout,status="unknown")

        nummol = 0
        do i = 1,nchain
           if (nboxi(i).eq.ibox) then
              nummol = nummol + nunit(moltyp(i))
           end if
        end do
        write(200+ibox,*) nummol
        write(200+ibox,*)
        do i = 1,nchain
           if(nboxi(i).eq.ibox) then
              imolty = moltyp(i)
              do ii = 1,nunit(imolty)
                 ntii = ntype(imolty,ii)
                 write(200+ibox,'(a4,5x,3f15.4)') chemid(ntii), rxu(i,ii), ryu(i,ii), rzu(i,ii)
              end do
           end if
        end do
        close(200+ibox)
     end do

     if (L_add) then
        do i = nchain+1, nchain+N_add
           moltyp(i)=N_moltyp2add
           nboxi(i) = N_box2add
           rxu(i,1) = random(-1)*boxlx(N_box2add)
           ryu(i,1) = random(-1)*boxly(N_box2add)
           rzu(i,1) = random(-1)*boxlz(N_box2add)
           qqu(i,1) = 0.0
        end do
        nchain = nchain+N_add
     else if (L_sub) then
        point_of_start = 0
        do i=1,N_moltyp2sub
           point_of_start=point_of_start+temtyp(i)
        end do
        point_of_start = point_of_start-N_sub+1

        point_to_end = nchain-N_sub

        do i = point_of_start,point_to_end
           moltyp(i)= moltyp(i+N_sub)
           nboxi(i) = nboxi(i+N_sub)
           imolty = moltyp(i)
           do  j = 1, nunit(imolty)
              rxu(i,j) = rxu(i+N_sub,j)
              ryu(i,j) = ryu(i+N_sub,j)
              rzu(i,j) = rzu(i+N_sub,j)
              qqu(i,j) = qqu(i+N_sub,j)
           end do
        end do
        nchain = nchain - N_sub
     end if

     ! write out the final configuration from the run
     write(file_config,'("config",I1.1,A,".dat")') run_num,suffix
     call dump(file_config)

     ! calculate and write out running averages ***
     do ibox=1,nbox
        ! energies
        do j=1,nener
           avv(j,ibox)   = acv(j,ibox) / acmove
           acvsq(j,ibox) = (acvsq(j,ibox)/acmove) - avv(j,ibox) ** 2
           acvkjmol(j,ibox) = acvkjmol(j,ibox)/acmove
        end do
        if ( ldielect ) then
           flucmom(ibox) = acvsq(14,ibox)-avv(15,ibox)*avv(15,ibox)
           ! momconst = 6.9994685465110493E5_dp
           flucmom(ibox) = 6.9994685465110493E5_dp*flucmom(ibox)*beta/ (boxlx(ibox)*boxly(ibox)*boxlz(ibox))
           flucmom2(ibox) = acvsq(14,ibox)-avv(11,ibox)*avv(11,ibox) -avv(12,ibox)*avv(12,ibox)  - avv(13,ibox)*avv(13,ibox)
           ! momconst = 6.9994685465110493E5_dp
           flucmom2(ibox) = 6.9994685465110493E5_dp*flucmom2(ibox)*beta /(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
        end if
        ! boxlength
        acboxl(ibox,1) = acboxl(ibox,1) / acmove
        acboxl(ibox,2) = acboxl(ibox,2) / acmove
        acboxl(ibox,3) = acboxl(ibox,3) / acmove

        if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
           acboxa(ibox,1) = acboxa(ibox,1) / acmove
           acboxa(ibox,2) = acboxa(ibox,2) / acmove
           acboxa(ibox,3) = acboxa(ibox,3) / acmove
        end if

        acvol(ibox) = acvol(ibox) / acmove
        acvolsq(ibox) = acvolsq(ibox) / acmove

        do itype = 1, nmolty
           ! number of molecules
           acnbox(ibox,itype) = acnbox(ibox,itype) / acmove
           ! molfraction
           molfra(ibox,itype) = molfra(ibox,itype) / acmove
           ! square end-to-end length
           if ( mnbox(ibox,itype) .gt. 0 ) then
              asetel(ibox,itype) =  asetel(ibox,itype) / real(mnbox(ibox,itype),dp)
           end if
        end do

        if ( lpbcz ) then
           do itype = 1, nmolty
              ! number density
              acdens(ibox,itype)=1000.0E0_dp*acdens(ibox,itype)/acmove
           end do
           ! sum over all types of molecules
           temacd = 0.0E0_dp
           do itype = 1, nmolty
              temacd = temacd + acdens(ibox,itype)
           end do

           ! molar volume
           molvol(ibox) = 602.2045E0_dp / temacd
           temspd = 0.0E0_dp
           do itype = 1, nmolty
              temspd = temspd +  ( acdens(ibox,itype) * masst(itype)/602.2045E0_dp)
           end do
           ! specific density
           speden(ibox) = temspd
        else
           do itype = 1, nmolty
              ! number density
              acdens(ibox,itype)=100.0E0_dp*acdens(ibox,itype)/acmove
           end do
           temacd = 0.0E0_dp
           do itype = 1, nmolty
              temacd = temacd + acdens(ibox,itype)
           end do
           ! molar volume
           molvol(ibox) = 100.0E0_dp / temacd
        end if

        ! system volume- convert to average
        acvolume(ibox) = acvolume(ibox) / acmove

        ! pressure and surface tension
        if ( acnp .gt. 0.5E0_dp ) then
           acpres(ibox) = acpres(ibox) / acnp
           acsurf(ibox) = acsurf(ibox) / acnp
        end if

        ! thermodynamic integration stuff
        if (acipsw.gt.0.5E0_dp) acdvdl = acdvdl/acipsw

        ! chemical potential
        ! This expression is suitable only for NPT ensemble
        do itype = 1, nmolty
           if( bnchem(ibox,itype) .gt. 0 ) then
              debroglie = debroglie_factor*sqrt(beta/masst(itype))
              ! determine how many steps it takes to grow molecule
              ! not counting the first inserted bead
              igrow = nugrow(itype)
              if (.not. lrigid(itype)) then
                 call schedule(igrow,itype,steps,1,0,2)
                 acchem(ibox,itype) = (-1.0E0_dp/beta) * log((acchem(ibox,itype)/bnchem(ibox,itype)) / ( real(nchoi1(itype)*(nchoi(itype)**steps)*nchoih(itype),dp) * (debroglie**3) ) )
              else
                 call schedule(igrow,itype,steps,1,0,4)
                 acchem(ibox,itype) = (-1.0E0_dp/beta) * log((acchem(ibox,itype)/bnchem(ibox,itype)) / ( real(nchoi1(itype)*nchoir(itype)*(nchoi(itype)**steps)*nchoih(itype),dp) * (debroglie**3) ) )
              end if
           end if
        end do
        if (acvsq(1,ibox).gt.0.0E0_dp) aflv(ibox)=sqrt(acvsq(1,ibox))
     end do

     write(io_output,"(' Averages and fluctuations',21x,4(a11,i1))") ('       Box ',i,i=1,nbox)
     write(io_output,*)
     write(io_output,"(' pressure                               [kPa] =',3f12.2)") (acpres(i) ,i=1,nbox)
     write(io_output,"(' pressure                  [simulation units] =',3f12.6)") ((acpres(i)*7.2429E-5_dp),i=1,nbox)
     write(io_output,"(' surface tension                       [mN/m] =',3f12.4)") (acsurf(i) ,i=1,nbox)
     do itype = 1, nmolty
        write(io_output,"(' chem. potential of type    ',i4,'          [K] =',3f12.3)") itype, (acchem(i,itype) ,i=1,nbox)
     end do
     write(io_output,*)

     do i = 1,3
        write(io_output,"(' boxlength                                [A] =',3f12.3)") (acboxl(ibox,i) ,ibox=1,nbox)
     end do

     do ibox = 1, nbox
        if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
           do i = 1,3
              write(io_output,"(' box angle                                deg =',3f12.3)") acboxa(ibox,i)*180.0E0_dp/onepi
           end do
        end if
     end do

     do itype = 1, nmolty
        write(io_output,"(' no. of chains of type      ',i4,'              =',3f12.3)") itype, (acnbox(i,itype) ,i=1,nbox)
     end do
     if ( lpbcz ) then
        write(io_output,"(' molar volume                      [cm^3/mol] =',3f12.3)") (molvol(i) ,i=1,nbox)
        write(io_output,"(' specific density                    [g/cm^3] =',3f12.6)") (speden(i) ,i=1,nbox)
        do itype = 1, nmolty
           write(io_output,"(' number density of type     ',i4,' [chain/nm^3] =',3f12.5)") itype, (acdens(i,itype) ,i=1,nbox)
           if ( lexpand(itype) ) then
              do itype2 = 1, numcoeff(itype)
                 write(io_output,"(' number density of type     ',i4,' eetype ',i4,'  =', 2f12.5)") itype, itype2,acdens(itype,itype)* acnbox2(itype,itype,itype2)/ (acnbox(itype,itype)*acmove), acnbox2(itype,itype,itype2)/ (acnbox(itype,itype)*acmove)
              end do
           end if
        end do
     else
        write(io_output,"(' area per chain                     [A^2/chain] =',3f12.4)") (molvol(i), i=1,nbox)
        do itype = 1, nmolty
           write(io_output,"(' number density of type     ',i4,' [chain/nm^2] =',3f12.6)") itype, (acdens(i,itype), i=1,nbox)
        end do
     end if
     do itype = 1, nmolty
        write(io_output,"(' molfraction of type        ',i4,'              =',3f12.7)") itype, (molfra(i,itype), i=1,nbox)
     end do
     do itype = 1, nmolty
        write(io_output,"(' mean sete length of type   ',i4,'        [A^2] =',3f12.3)") itype, (asetel(i,itype) ,i=1,nbox)
     end do
     write(io_output,*)
     do j=1,10
        ! only 1 to 10 is the energy information
        write(io_output,"(a15,'[K per system and kJ/mol per chain] =',3(f14.2,f12.2))") vname(j),avv(j,1:nbox),acvkjmol(j,1:nbox)
     end do

     write(io_output,*)
     write(io_output,"(' fluctuation in <vtot>                          =',3f12.2)") (aflv(i) ,i=1,nbox)
     write(io_output,*)

     ! Output 2nd virial coefficient data
2000 if (lvirial) then
        starviro = starvir
        dummy = real(nstep/imv,dp)
        do itemp = 1,ntemp
           starvir = starviro
           binvir(1,itemp) = binvir(1,itemp)/dummy
           ! write(45,*) starvir,binvir(1,itemp)
           inside = starvir*starvir*binvir(1,itemp)
           ! write(46,*) starvir,inside
           bvirial = 0.5E0_dp*inside
           ! write(47,*) starvir,bvirial
           starvir = starvir + stepvir

           do i = 2,nvirial-1
              binvir(i,itemp) = binvir(i,itemp)/dummy
              ! write(45,*) starvir,binvir(i,itemp)
              inside = starvir*starvir*binvir(i,itemp)
              ! write(46,*) starvir,inside
              bvirial = bvirial + inside
              ! write(47,*) starvir,bvirial
              starvir = starvir + stepvir
           end do

           binvir(nvirial,itemp) = binvir(nvirial,itemp)/dummy
           ! write(45,*) starvir,binvir(nvirial,itemp)
           inside = starvir*starvir*binvir(nvirial,itemp)
           ! write(46,*) starvir,inside
           ! write(47,*) starvir,bvirial
           starvir = starvir + stepvir
           bvirial = bvirial + 0.5E0_dp*inside

           write(io_output,*) 'At temperature of',virtemp(itemp)
           write(io_output,*) 'bvirial ', -(twopi*stepvir*bvirial),' [A^3 / molecule]'
           write(io_output,*) 'bvirial ',-0.602E0_dp*twopi* stepvir*bvirial,' [cm^3 / mole]'

           ! if ( lvirial2 ) then
           starvir = starviro + 0.5E0_dp*stepvir
           do i = 2, nvirial
              binvir2(i,itemp) =  binvir2(i,itemp)/dummy
              inside = starvir*starvir*binvir2(i,itemp)
              bvirial = bvirial + inside
              starvir = starvir + stepvir
           end do
           bvirial = -(twopi*stepvir*bvirial)
           write(io_output,*) 'With quantum correction:'
           write(io_output,*) 'bvirial ',bvirial,' [A^3 / molecule]'
           write(io_output,*) 'bvirial ',0.602E0_dp*bvirial,' [cm^3 / mole]'
        end do
     end if

     ! solute values
     write(io_output,*) 'type  box     vinter      vintra      vtor', '        vbend       vtail'

     do itype = 1, nmolty
        do ibox = 1, nbox
           if (solcount(ibox,itype).gt.0) then
              write(io_output,"(i5,i5,3f12.5,3f12.5,3f12.5,3f12.5,3f12.5)") itype,ibox,avsolinter(ibox,itype) /solcount(ibox,itype)&
               ,avsolintra(ibox,itype) /solcount(ibox,itype),avsoltor(ibox,itype) /solcount(ibox,itype)&
               ,avsolbend(ibox,itype) /solcount(ibox,itype),avsolelc(ibox,itype) /solcount(ibox,itype)
           else
              write(io_output,"(i5,i5,3f12.5,3f12.5,3f12.5,3f12.5,3f12.5)") itype,ibox,0.0,0.0,0.0,0.0,0.0
           end if
        end do
     end do

     ! calculate statistical errors ---
     if ( nblock .ge. 2 ) then
        dblock = real(nblock,dp)
        dbl1 = dblock - 1.0E0_dp
        ! global averages -
        do i = 1,nprop
           do j = 1,nbox
              if ( naccu(i,j) .lt. 0.5E-5_dp ) then
                 aver(i,j) = 0.0E0_dp
              else
                 aver(i,j) = accum(i,j) / naccu(i,j)
              end if
           end do
        end do
        do i = 1,nprop
           do j = 1,nbox
              do nbl = 1, nblock
                 dsq(i,j) = dsq(i,j) +  ( baver(i,j,nbl) - aver(i,j) )**2
              end do
              stdev(i,j) = sqrt( dsq(i,j) / dblock )
              sterr(i,j) = sqrt( dsq(i,j) / dbl1 )
              errme(i,j) = sterr(i,j) / sqrt(dblock)
           end do
        end do

        do i = 1,nprop1
           do ibox = 1,nbox-1
              do jbox = ibox+1,nbox
                 if ( naccu1(i,ibox,jbox) .lt. 0.5E-5_dp ) then
                    aver1(i,ibox,jbox) = 0.0E0_dp
                 else
                    aver1(i,ibox,jbox) = accum1(i,ibox,jbox) /  naccu1(i,ibox,jbox)
                 end if
              end do
           end do
        end do

        do i = 1,nprop1
           do ibox = 1,nbox-1
              do jbox = ibox+1,nbox
                 do nbl = 1, nblock
                    dsq1(i,ibox,jbox) = dsq1(i,ibox,jbox) + ( baver1(i,ibox,jbox,nbl) - aver1(i,ibox,jbox) )**2
                 end do
                 stdev1(i,ibox,jbox) = sqrt( dsq1(i,ibox,jbox) / dblock )
                 sterr1(i,ibox,jbox) = sqrt( dsq1(i,ibox,jbox) /  dbl1 )
                 errme1(i,ibox,jbox) = sterr1(i,ibox,jbox) /  sqrt(dblock)
              end do
           end do
        end do

        ! write out the heat of vaporization and solubility parameters
        do ibox = 1,nbox-1
           do jbox = ibox+1,nbox
              write(io_output,"(' H_vap      [kJ/mol] btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(1,ibox,jbox), stdev1(1,ibox,jbox),errme1(1,ibox,jbox)
              write(io_output,"(' H_vap LJ  [kJ/mol] btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(2,ibox,jbox), stdev1(2,ibox,jbox),errme1(2,ibox,jbox)
              write(io_output,"(' H_vap Coul [kJ/mol] btwn box  ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(3,ibox,jbox), stdev1(3,ibox,jbox),errme1(3,ibox,jbox)
              ! write(io_output,"(' DeltaU Ext [kJ/mol] btwn box   ',i4,' and',i4, ' =',3f15.4)") ibox,jbox,aver1(10,ibox,jbox),stdev1(10,ibox,jbox), errme1(10,ibox,jbox)
              write(io_output,"(' pdV        [kJ/mol] btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(11,ibox,jbox), stdev1(11,ibox,jbox), errme1(11,ibox,jbox)

              write(io_output,"(' CED [cal/cc]   btwn box        ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(4,ibox,jbox), stdev1(4,ibox,jbox),errme1(4,ibox,jbox)
              write(io_output,"(' CED_LJ[cal/cc] btwn box        ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(5,ibox,jbox), stdev1(5,ibox,jbox),errme1(5,ibox,jbox)
              write(io_output,"(' CED_Coul[cal/cc] btwn box      ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(6,ibox,jbox), stdev1(6,ibox,jbox),errme1(6,ibox,jbox)
              write(io_output,"(' HSP [(cal/cc)^1/2]  btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(7,ibox,jbox), stdev1(7,ibox,jbox),errme1(7,ibox,jbox)
              write(io_output,"(' HSP_LJ[(cal/cc)^1/2] btwn box  ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(8,ibox,jbox), stdev1(8,ibox,jbox),errme1(8,ibox,jbox)
              write(io_output,"(' HSP_Cou[(cal/cc)^1/2] btwn box ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,aver1(9,ibox,jbox), stdev1(9,ibox,jbox),errme1(9,ibox,jbox)
           end do
        end do

        ! specific density
        do ibox = 1, nbox
           write(io_output,"(' specific density    box ',i3, ' = ',3E12.5)") ibox,aver(1,ibox),stdev(1,ibox), errme(1,ibox)
        end do

        ! system volume
        itel = 4 + nener + 4*nmolty
        do ibox = 1, nbox
           write(io_output,"(' system volume       box ',i3, ' = ',3E12.5)") ibox, aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
        end do

        ! pressure
        do ibox = 1, nbox
           write(io_output,"(' pressure            box ',i3, ' = ',3g12.5)") ibox,aver(2,ibox),stdev(2,ibox), errme(2,ibox)
        end do

        ! surface tension
        itel = 2+nener+ 4*nmolty+1
        do ibox = 1, nbox
           write(io_output,"(' surface tension     box ',i3, ' = ',3f12.5)") ibox, aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
        end do

        write(io_output,*)
        ! energies
        ! write(io_output,*) 'average value', 'STD', 'SEM'
        do ibox = 1, nbox
           do j=3,2+10
              ! only 1 to 10 is the energy information
              write(io_output,"(a15 ,' box ',i3, ' = ',3E14.5)") vname(j-2),ibox,aver(j,ibox), stdev(j,ibox),errme(j,ibox)
           end do
        end do

        write(io_output,*)

        ! Enthalpy
        do ibox = 1,nbox
           j = 4+nener +4*nmolty + 1
           write(io_output, "(a15,'[kJ/mol] for box',i3,' =',3(f12.4))") enth, ibox, aver(j,ibox), stdev(j,ibox),errme(j,ibox)
           j = 4+nener + 4*nmolty + 2
           write(io_output,"(a15,'[kJ/mol] for box',i3,' =',3(f12.4))") enth1,ibox,aver(j,ibox), stdev(j,ibox),errme(j,ibox)
        end do
        write(io_output,*)

        ! Residual Heat capacity ---
        if(lnpt.and..not.lgibbs) then
           inst_enth= inst_enth/real(nchain*nstep,dp)
           inst_enth2= inst_enth2/real(nchain*nstep,dp)
           sigma2Hsimulation=((inst_enth2)- (inst_enth*inst_enth))
           sigma2H=sigma2Hsimulation*(6.022E23_dp)*((1.38066E-23_dp)**2) / real(nchain,dp) !(J2/mol)
           Cp=sigma2H/((1.38066E-23_dp)*(temp**2))
           write(io_output,*) 'Cp residual(J/Kmol) =', Cp
           write(io_output,*) ' H2=', inst_enth2
           write(io_output,*)  ' H=', inst_enth
        end if

        if( .not. lnpt .and..not.lgibbs) then
           inst_energy= inst_energy/real(nchain*nstep,dp)
           inst_energy2= inst_energy2/real(nchain*nstep,dp)
           sigma2Esimulation=((inst_energy2)-  (inst_energy*inst_energy))
           sigma2E=sigma2Esimulation*(6.022E23_dp)*((1.38066E-23_dp)**2) / real(nchain,dp) !(J2/mol)
           Cv=sigma2E/((1.38066E-23_dp)*(temp**2))
           write(io_output,*) 'Cv residual(J/Kmol) =', Cv
           write(io_output,*) ' E2=', inst_energy2
           write(io_output,*) ' E=', inst_energy
        end if

        ! chemical potential
        do itype = 1, nmolty
           itel = (2+nener) + itype
           do ibox = 1, nbox
              if ( aver(itel,ibox) .ne. 0.0E0_dp ) then
                 write(io_output,"(' chemical potential  itype ',i3,' box ',i3, ' = ',3f12.3)") itype,ibox, ((-1.0E0_dp)/beta)*log(aver(itel,ibox)), (1.0E0_dp/beta)*stdev(itel,ibox)/aver(itel,ibox), (1.0E0_dp/beta)*errme(itel,ibox)/aver(itel,ibox)
              else
                 write(io_output,"(' chemical potential  itype ',i3,' box ',i3, ' = ',3f12.3)") itype,ibox,0.0E0_dp,0.0E0_dp,0.0E0_dp
              end if
           end do
        end do

        ! square end-to-end length
        do itype = 1, nmolty
           itel = (2+nener) + nmolty + itype
           do ibox = 1, nbox
              write(io_output,"(' mean sete length    itype ',i3,' box ',i3, ' = ',3f12.3)") itype,ibox,aver(itel,ibox) ,stdev(itel,ibox),errme(itel,ibox)
           end do
        end do

        ! number density
        do itype = 1, nmolty
           itel = (2+nener) + 2 * nmolty + itype
           do ibox = 1, nbox
              if ( lpbcz ) then
                 write(io_output,"(' number density      itype ',i3,' box ',i3, ' = ',3E12.5)") itype,ibox,1.0E3_dp*aver(itel,ibox) ,1.0E3_dp*stdev(itel,ibox),1.0E3_dp*errme(itel,ibox)
              else
                 write(io_output,"(' number density      itype ',i3,' box ',i3, ' = ',3E12.5)") itype,ibox,1.0E2_dp*aver(itel,ibox) ,1.0E2_dp*stdev(itel,ibox),1.0E2_dp*errme(itel,ibox)
              end if
              if ( lexpand(itype) .and.  acnbox(ibox,itype) .gt. 0.5) then
                 do itype2 = 1, numcoeff(itype)
                    molfrac = acnbox2(ibox,itype,itype2) /(acmove*acnbox(ibox,itype))
                    write(io_output,"(' number density      itype ',i3,' typ ',i3, ' = ',2E12.5)") itype,itype2,1.0E3_dp* aver(itel,ibox)*molfrac,molfrac
                 end do
              end if
           end do
        end do

        ! molfraction
        do itype = 1, nmolty
           itel = (2+nener) + 3 * nmolty + itype
           do ibox = 1, nbox
              write(io_output,"(' mole fraction       itype ',i3,' box ',i3, ' = ',3f12.7)") itype,ibox,aver(itel,ibox), stdev(itel,ibox),errme(itel,ibox)
           end do
        end do

        if (lgibbs) then
           ! write density results in fitting format ---
           do ibox = 1, nbox-1
              do jbox = ibox+1,nbox
                 if (speden(ibox).lt.speden(jbox)) then
                    ig = ibox
                    il = jbox
                 else
                    ig = jbox
                    il = ibox
                 end if
                 ! write(41,"(2x,f6.1,4(1x,e12.5))") temp,aver(1,ig),stdev(1,ig),aver(1,il),stdev(1,il)
                 ! write ostwald values for each moltyp
                 gconst = 8.314/(1000*beta)
                 do itype = 1,nmolty
                    itel = (2+nener) + 2 * nmolty + itype
                    ostwald = aver(itel,il)/aver(itel,ig)
                    stdost  = ostwald * sqrt( (stdev(itel,il)/aver(itel,il))**2 + (stdev(itel,ig)/aver(itel,ig))**2 )
                    ! write(42,*) nunit(itype),ostwald,stdost
                    ! write(43,*) nunit(itype),-(gconst*log(ostwald)) + (eta2(ig,itype) - eta2(il,itype)) / 120.27167,gconst*stdost/ostwald
                    write(io_output,"('Ostwald Coefficient  itype ',i3,' between box ',i2, ' and ',i2,f18.6,f18.6)") itype,ig,il,ostwald,stdost
                    write(io_output,"('Free Enrgy of Transf itype ',i3,' between box ',i2, ' and ',i2,f18.6,f18.6,' kJ/mol')") itype,ig,il, -(gconst*log(ostwald)) + (eta2(ig,itype)  - eta2(il,itype)) / 120.27167 ,gconst*stdost/ostwald
                 end do
              end do
           end do
        end if

        write(io_output,*)

        ! write block averages  ---
        write(io_output,*)
        write(io_output,*) '-----block averages ------'
        do ibox=1,nbox
           write(io_output,"('  ------------ box: ' ,i4,/, ' block   energy    density   pressure   surf ten  mol fracs')") ibox
           do nbl = 1, nblock
              ! changed so output the same for all ensembles
              ! 06/08/09 KM
              write(io_output,"(2x,i2,15(2x,e10.3))") nbl,baver(3,ibox,nbl), baver(1,ibox ,nbl),baver(2,ibox,nbl), baver(3+nener+4*nmolty,ibox ,nbl), (baver(2+nener+3*nmolty+zzz,ibox,nbl),zzz=1 ,nmolty)
           end do
           if (lmipsw) then
              write(io_output,*) 'lambdais', lambdais
              write(io_output,*) 'maginn interphase switch integrand'
              do nbl = 1, nblock
                 write(io_output,*) nbl,baver(nprop,ibox,nbl)
              end do
           end if
        end do
     end if

     ! KM 01/10 remove analysis
     ! if (ianalyze.le.nstep) then
     ! call analysis(2)
     ! end if

     ! ee prob
     IF(lexpee) then
        write(io_output,*)
        write(io_output,*) 'probability of each mstate in ee'
        do nnn = 1, fmstate
           write(io_output,"(i5,1x,i10)") nnn,ee_prob(nnn)
        end do
     end if
     if (L_movie_xyz) then
        do ibox=1,nbox
           close (210+ibox)
        end do
     end if
     write(io_output,*) 'Program ended at ',time_date_str()
  end if

  deallocate(nminp,nmaxp,ncmt_list,ndist,pres,acvkjmol,acsolpar,acEnthalpy,acEnthalpy1,stdev1,sterr1,errme1,vstart,vend,avv,acv&
   ,acvsq,aflv,acboxl,acboxa,acpres,acsurf,acvolume,acnbox,acnbox2,mnbox,asetel,acdens,molfra,dsq,stdev,dsq1,sterr,errme,lratfix&
   ,lsolute,molvol,speden,flucmom,flucmom2,acvol,acvolsq,file_ndis,stat=jerr)
  if (jerr.ne.0) then
     call err_exit(__FILE__,__LINE__,'monola: deallocation failed',jerr)
  end if

  return
end subroutine monola
