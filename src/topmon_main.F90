MODULE topmon_main
  use sim_system
  use sim_cell
  implicit none
  private
  save
  public::monola

  ! variables added for GCMC histogram reweighting
  integer,parameter::fmax=1E6,nprop1=11
  logical::use_checkpoint
  integer::blockm,checkpoint_interval=1800,checkpoint_copies=1,nstep,nnstep,nnn,acmove,acnp,acipsw,nblock&
   ,io_movie,io_solute,io_cell,io_traj
  real::enthalpy,enthalpy2& !< enthalpy (NpT) or internal energy (NVT)
     ,acdvdl,binvir(maxvir,maxntemp),binvir2(maxvir,maxntemp)
  real,allocatable::acdens(:,:)& !< (ibox,itype): accumulators of box density
   ,molfra(:,:)& !< (ibox,itype): accumulators of mole fraction
   ,acnbox(:,:)& !< accumulators of ncmt
   ,acnbox2(:,:,:)& !< accumulators of ncmt2
   ,acv(:,:),acvsq(:,:),acvkjmol(:,:)& !< (j,ibox): accumulators of vbox
   ,acdipole(:,:),acdipolesq(:,:)& !< (j,ibox): accumulators of dipole values; j = 1: dipolex; 2: dipoley; 3: dipolez; 4: norm of dipole
   ,acboxa(:,:)& !< (ibox,j) accumulators of cell angles for ibox
   ,acboxl(:,:)& !< (ibox,j) accumulators of cell lengths for ibox
   ,acvol(:),acvolsq(:),acvolume(:)& !< accumulators of box volume
   ,acpres(:),acsurf(:)& !< accumulators of pressure and surface tension
   ,acEnthalpy(:),acEnthalpy1(:)& !< accumulators of enthalpies using calculated pressure (acEnthalpy) and specified enthalpy (acEnthalpy1)
   ,acsolpar(:,:,:),avsolinter(:,:),avsolintra(:,:),avsolbend(:,:),avsoltor(:,:),avsolelc(:,:)&
   ,asetel(:,:)& !< (ibox,itype): accumulators of square end-to-end distance; counter is in mnbox
   ,vstart(:),vend(:),pres(:),molvol(:),speden(:)&
   ,aver1(:,:,:),stdev1(:,:,:),errme1(:,:,:),bccold1(:,:,:),baver1(:,:,:,:)& !< accumulators for solubility parameter and heat of vaporization (Neeraj)
   ,aver(:,:),stdev(:,:),errme(:,:),bccold(:,:),baver(:,:,:) !< (j,ibox): properties of ibox; j =
  !< \verbatim
  !< ---------------------------------------------------------
  !< 1                                              = specific density
  !< 2                                              = pressure
  !< 3                      to (2+nEnergy)          = energies
  !< 1+(2+nEnergy)          to (2+nEnergy)+  nmolty = chemical potential
  !< 1+(2+nEnergy)+  nmolty to (2+nEnergy)+2*nmolty = square-end-to-end-length
  !< 1+(2+nEnergy)+2*nmolty to (2+nEnergy)+3*nmolty = number density
  !< 1+(2+nEnergy)+3*nmolty to (2+nEnergy)+4*nmolty = mole fraction
  !< 1+(2+nEnergy)+4*nmolty                         = surface tension
  !< 2+(2+nEnergy)+4*nmolty                         = system volume
  !< 3+(2+nEnergy)+4*nmolty                         = enthalpy inst
  !< 4+(2+nEnergy)+4*nmolty                         = enthalpy ext
  !< ---------------------------------------------------------
  !< \endverbatim
  integer,allocatable::io_box_movie(:),nminp(:),nmaxp(:),ncmt_list(:,:),ndist(:,:)& !< GCMC reweighting histograms
   ,mnbox(:,:),solcount(:,:),nccold1(:,:,:),nccold(:,:)

contains
!> \brief Main control logic of topmon
!>
!> reads the control-data from unit 4
!> starts and controls the simulation
  subroutine monola(file_in)
    use var_type,only:dp,default_string_length
    use const_math,only:onepi,twopi,raddeg
    use const_phys,only:debroglie_factor,N_Avogadro,R_gas,k_B,MPa2SimUnits
    use util_math,only:update_average,calculate_statistics
    use util_random,only:random
    use util_string,only:format_n,integer_to_string
    use util_runtime,only:err_exit
    use util_timings,only:time_date_str,time_now
    use util_files,only:get_iounit
    use util_mp,only:mp_barrier
    use sim_particle,only:init_neighbor_list
    use energy_pairwise,only:sumup
    use moves_simple,only:translation,rotation,Atom_translation,output_translation_rotation_stats
    use moves_volume,only:volume_1box,volume_2box,output_volume_stats
    use moves_cbmc,only:config,schedule,output_safecbmc,output_cbmc_stats
    use moves_ee,only:eesetup,eemove,ee_index_swap,expand,numcoeff,output_ee_stats
    use transfer_swap,only:swap,cnt,output_swap_stats,acchem,bnchem
    use transfer_swatch,only:swatch,output_swatch_stats
    use prop_pressure,only:pressure

    character(LEN=*),intent(in)::file_in

    character(LEN=default_path_length)::file_flt,file_hist,file_ndis,file_cnt,file_config

    ! descriptions of different kinds of energies
    character(LEN=default_string_length)::vname(nEnergy)=(/' Total energy',' Inter LJ',' Tail  LJ',' Intra LJ',' Stretch',' Bond bending',' Torsion',' Coulomb',' External pot',' 3-body Garo',' Fluc Q','','',''/)

    integer::io_flt,io_hist,io_cnt,io_ndis,io_config,i,jerr,ibox,itype,itype2,Temp_nmol,nentry,j,nummol,imolty,ii,ntii,point_of_start,point_to_end,igrow,steps,itemp,jbox,itel,ig,il,nbl,n,zzz,ichkpt,nnn_1st
    real::v(nEnergy),press1,surf,time_prev,time_cur,rm,temvol,tmp,vhist,eng_list(fmax),temacd,temspd,debroglie,starviro,dummy,inside,bvirial,gconst,ostwald,stdost,molfrac
    logical::ovrlap

    ! KM for MPI
    ! only one processor at a time reads and writes data from files
    do i=1,numprocs
       if (myid.eq.i-1) then
          call readdat(file_in)
       end if
       call mp_barrier(groupid)
    end do

    allocate(nminp(ntmax),nmaxp(ntmax),ncmt_list(fmax,ntmax),ndist(0:nmax,ntmax),vstart(nbxmax),vend(nbxmax),pres(nbxmax),molvol(nbxmax),speden(nbxmax),acdens(nbxmax,ntmax),molfra(nbxmax,ntmax),acnbox(nbxmax,ntmax),acnbox2(nbxmax,ntmax,20),acv(nEnergy,nbxmax),acvsq(nEnergy,nbxmax),acvkjmol(nEnergy,nbxmax),acdipole(4,nbxmax),acdipolesq(4,nbxmax),acboxa(nbxmax,3),acboxl(nbxmax,3),acvol(nbxmax),acvolsq(nbxmax),acvolume(nbxmax),acpres(nbxmax),acsurf(nbxmax),acEnthalpy(nbxmax),acEnthalpy1(nbxmax),solcount(nbxmax,ntmax),acsolpar(nprop1,nbxmax,nbxmax),avsolinter(nbxmax,ntmax),avsolintra(nbxmax,ntmax),avsolbend(nbxmax,ntmax),avsoltor(nbxmax,ntmax),avsolelc(nbxmax,ntmax),mnbox(nbxmax,ntmax),asetel(nbxmax,ntmax),nccold1(nprop1,nbxmax,nbxmax),bccold1(nprop1,nbxmax,nbxmax),baver1(nprop1,nbxmax,nbxmax,blockm),nccold(nprop,nbxmax),bccold(nprop,nbxmax),baver(nprop,nbxmax,blockm),aver1(nprop1,nbxmax,nbxmax),stdev1(nprop1,nbxmax,nbxmax),errme1(nprop1,nbxmax,nbxmax),aver(nprop,nbxmax),stdev(nprop,nbxmax),errme(nprop,nbxmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'monola: allocation failed',jerr)
    end if

    ! SETTING UP ARRAYS FOR ANALYSYS PURPOSE
    ! JLR 11-11-09
    ! do not call analysis if you set ianalyze to be greater than number of cycles
    ! KM 01/10 remove analysis
    ! if (ianalyze.le.nstep) then
    !    call analysis(0)
    ! end if
    ! END JLR 11-11-09

    ! set up initial linkcell
    if (licell) then
       call build_linked_cell()
    end if

    if (lneigh) then
       ! call for initial set-up of the near-neighbour bitmap ***
       call init_neighbor_list()
    end if

    ! set up thermodynamic integration stuff
    if (lmipsw) then
       call ipswsetup()
    else
       lstagea = .false.
       lstageb = .false.
       lstagec = .false.
    end if

    ! set up expanded ensemble stuff
    if (lexpee) then
       call eesetup
       if (lmipsw) call err_exit(__FILE__,__LINE__,'not for BOTH lexpee AND lmipsw',myid+1)
    else
       leemove=.false.
    end if

    if (lneighbor) then
       neighbor = 0
       if (myid.eq.rootid) then
          file_cnt='fort.21'
          io_cnt=get_iounit()
          open(unit=io_cnt,access='sequential',action='write',file=file_cnt,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open cnt file '//trim(file_cnt),jerr)
          write(io_cnt,*) 'ii:',nnstep,(neigh_cnt(i),i=1,nchain)
       end if
    end if

    ! setup files for histogram reweighting
    ! KM fom MPI
    ! will need to check this file I/O if want to run grand canonical in parallel
    if(lgrand) then
       if (myid.eq.rootid) then
          write(file_flt,'("nfl",I1.1,A,".dat")') run_num,suffix
          io_flt=get_iounit()
          open(unit=io_flt,access='sequential',action='write',file=file_flt,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open flt file '//trim(file_flt),jerr)

          write(file_hist,'("his",I1.1,A,".dat")') run_num,suffix
          io_hist=get_iounit()
          open(unit=io_hist,access='sequential',action='write',file=file_hist,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open hist file '//trim(file_hist),jerr)
          write(io_hist,'(f8.4,2x,i5,2x,g15.5,3f12.3)')  temp,nmolty,(temp*log(B(i)),i=1,nmolty),boxlx(1),boxly(1),boxlz(1)
       end if

       ! extra zero accumulator for grand-canonical ensemble
       nentry = 0
       nminp = 1E6_dp
       nmaxp = -1E6_dp
       ndist = 0
    end if

    ! 2nd viral coefficient
    if (lvirial) then
       binvir = 0.0E0_dp
       binvir2 = 0.0E0_dp
       ! if ( lvirial2 ) then
       !    call virial2(binvir,binvir2,nvirial,starvir,stepvir)
       !    goto 2000
       ! end if
       ! profile=0.0E0_dp
       ! binstep = 0.05E0_dp
       ! lvirial2 = .false.
    end if

    ! initialize variables
    acmove = 0
    acdens = 0.0E0_dp
    molfra = 0.0E0_dp
    acnbox = 0.0E0_dp
    acnbox2 = 0.0E0_dp
    acv = 0.0E0_dp
    acvsq = 0.0E0_dp
    acvkjmol = 0.0E0_dp
    acdipole = 0.0_dp
    acdipolesq = 0.0_dp
    acboxa = 0.0E0_dp
    acboxl = 0.0E0_dp
    acvol = 0.0E0_dp
    acvolsq = 0.0E0_dp
    acvolume = 0.0E0_dp
    enthalpy = 0.0E0_dp
    enthalpy2 = 0.0E0_dp
    acnp = 0
    acpres = 0.0E0_dp
    acsurf = 0.0E0_dp
    acEnthalpy = 0.0E0_dp
    acEnthalpy1 = 0.0E0_dp
    solcount = 0
    acsolpar=0.0E0_dp
    avsolinter = 0.0E0_dp
    avsolintra = 0.0E0_dp
    avsolbend = 0.0E0_dp
    avsoltor = 0.0E0_dp
    avsolelc = 0.0E0_dp
    mnbox  = 0
    asetel = 0.0E0_dp
    acipsw=0
    acdvdl=0.0E0_dp
    ! accumulators for block averages ---
    nblock = 0
    nccold1 = 0
    bccold1 = 0.0E0_dp
    baver1=0.0E0_dp
    nccold = 0
    bccold = 0.0E0_dp
    baver=0.0E0_dp
    ! accumulators for fluctuating charge performance
    bnflcq = 0.0E0_dp
    bnflcq2 = 0.0E0_dp
    bsflcq = 0.0E0_dp
    bsflcq2 = 0.0E0_dp
    ichkpt = 0

! -----------------------------------------------------------------
    if (use_checkpoint) then
       call read_checkpoint_main('save-stats')
       nnn_1st=nnn+1
    else
       nnn_1st=1
       ! calculate initial energy and check for overlaps ***
       do ibox=1,nbox
          call sumup(ovrlap,v,ibox,lvol=.false.)
          vbox(:,ibox) = v
          vbox(10,ibox) = v3garo
          vbox(12,ibox) = vipsw
          vbox(13,ibox) = vwellipsw
          if (ovrlap) then
             call err_exit(__FILE__,__LINE__,'overlap in initial configuration',myid+1)
          end if

          vstart(ibox) = vbox(1,ibox)
          if (myid.eq.rootid) then
             write(io_output,*)
             write(io_output,*) 'box  ',ibox,' initial v   = ', vbox(1,ibox)
          end if

          ! calculate initial pressure ***
          call pressure( press1, surf, ibox )
          if (myid.eq.rootid) then
             write(io_output,"(' surf. tension :   box',i2,' =',f14.5)") ibox, surf
             write(io_output,"(' pressure check:   box',i2,' =',f14.2)") ibox, press1
          end if
       end do

       if (myid.eq.rootid) then
          write(io_output,*)
          write(io_output,*) '+++++ start of markov chain +++++'
          write(io_output,*)
          write(io_output,*)  'Cycle   Total   Energy    Boxlength   Pressure  Molecules'
       end if
    end if

!************************************************************
! loops over all cycles and all molecules                  **
!************************************************************
    time_prev = time_now()
    MC_cycle: do nnn = nnn_1st, nstep
       tmcc = nnstep + nnn
       do ii = 1, nchain
          acmove = acmove + 1.0E0_dp
          ! select a move-type at random ***
          rm = random(-1)

          ! special ensemble dependent moves ###
          if  (rm .le. pmvol) then
             ! volume move ---
             if ( lnpt ) then
                call volume_1box()
             else
                call volume_2box()
             end if
          else if (rm .le. pmswat) then
             ! CBMC switch move ---
             call swatch()
          else if ( rm .le. pmswap ) then
             ! swap move for linear and branched molecules ---
             call swap()
          else if ( rm .le. pmcb ) then
             ! configurational bias move ---
             call config()
          else if ( rm .le. pmflcq ) then
             ! displacement of fluctuating charges ---
             call flucq(2,0)
          else if (rm .le. pmexpc ) then
             ! expanded-ensemble move ---
             call expand()
          else if (rm .le. pmexpc1 ) then
             ! new expanded-ensemble move ---
             ! call expand
             if (random(-1).le.eeratio) then
                call ee_index_swap()
             else
                call eemove()
             end if
          else if ( rm .le. pm_atom_tra) then
             call Atom_translation()
          else if ( rm .le. pmtra ) then
             ! translational move ---
             call translation()
          else
             ! rotation move --
             call rotation()
          end if

          ! accumulate probability of being in an expanded ensemble state
          if (lexpee) then
             ee_prob(mstate) = ee_prob(mstate)+1
          end if

          ! calculate instantaneous values ***
          ! accumulate averages ***
          do ibox=1,nbox
             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                temvol = cell_vol(ibox)
             else
                if ( lpbcz ) then
                   temvol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
                else
                   temvol = boxlx(ibox)*boxly(ibox)
                end if
             end if

             do itype = 1, nmolty
                call update_average(acdens(ibox,itype),ncmt(ibox,itype)/temvol,acmove)
                if (nchbox(ibox).gt.0) then
                   call update_average(molfra(ibox,itype),real(ncmt(ibox,itype),dp)/real(nchbox(ibox),dp),acmove)
                end if

                call update_average(acnbox(ibox,itype),real(ncmt(ibox,itype),dp),acmove)
                if (lexpand(itype)) then
                   do itype2 = 1, numcoeff(itype)
                      call update_average(acnbox2(ibox,itype,itype2),real(ncmt2(ibox,itype,itype2),dp),acmove)
                   end do
                end if
             end do

             call update_average(acv(1:11,ibox),vbox(1:11,ibox),acmove)
             call update_average(acvsq(1:11,ibox),vbox(1:11,ibox)**2,acmove)

             if (lnpt) then
                call update_average(acv(14,ibox),vbox(1,ibox)*temvol,acmove)
             end if

             ! KMB/KEA Energy in kJ/mol
             Temp_nmol = sum(ncmt(ibox,1:nmolty))
             call update_average(acvkjmol(1:11,ibox),vbox(1:11,ibox)/Temp_nmol,acmove)

             ! leftover from Bin, not currently used
             if ( ldielect ) then
                call update_average(acdipole(1,ibox),dipolex(ibox),acmove)
                call update_average(acdipolesq(1,ibox),dipolex(ibox)**2,acmove)
                call update_average(acdipole(2,ibox),dipoley(ibox),acmove)
                call update_average(acdipolesq(2,ibox),dipoley(ibox)**2,acmove)
                call update_average(acdipole(3,ibox),dipolez(ibox),acmove)
                call update_average(acdipolesq(3,ibox),dipolez(ibox)**2,acmove)
                call update_average(acdipole(4,ibox),sqrt(dipolex(ibox)*dipolex(ibox)+dipoley(ibox)*dipoley(ibox)+dipolez(ibox)*dipolez(ibox)),acmove)
                acdipolesq(4,ibox)=acdipolesq(1,ibox)+acdipolesq(2,ibox)+acdipolesq(3,ibox)
             end if

             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                boxlx(ibox) = cell_length(ibox,1)
                boxly(ibox) = cell_length(ibox,2)
                boxlz(ibox) = cell_length(ibox,3)
                call update_average(acboxa(ibox,1:3),cell_ang(ibox,1:3),acmove)
             end if

             call update_average(acboxl(ibox,1),boxlx(ibox),acmove)
             call update_average(acboxl(ibox,2),boxly(ibox),acmove)
             call update_average(acboxl(ibox,3),boxlz(ibox),acmove)

             call update_average(acvol(ibox),temvol,acmove)
             call update_average(acvolsq(ibox),temvol*temvol,acmove)
             call update_average(acvolume(ibox),temvol,acmove)
          end do

          if (.not.lgibbs) then
             ibox=1
             if (lnpt) then
                tmp=vbox(1,ibox)+express(ibox)*boxlx(ibox)*boxly(ibox)*boxlz(ibox)
             else
                tmp=vbox(1,ibox)
             end if
             call update_average(enthalpy,tmp,acmove)
             call update_average(enthalpy2,tmp*tmp,acmove)
          end if

          ! collect histogram data (added 8/30/99)
          if (lgrand) then
             ibox = 1
             vhist = vbox(2,ibox) + vbox(8,ibox) + vbox(11,ibox) !inter+elect+flucq
             ! KM for MPI
             ! check this if want to run grand canonical
             if (mod(acmove,ninstf).eq.0.and.myid.eq.rootid) then
                ! Formatting removed until I can figure out how to get <n>i7 to work
                ! write(io_flt, * ) acmove,(ncmt(ibox,i), i=1,nmolty),vhist
                write(io_flt,FMT='(i0,5x,'//format_n(nmolty,'i7')//',5x,g15.6)') acmove,(ncmt(ibox,i),i=1,nmolty),vhist
             end if

             if(mod(acmove,ninsth).eq.0.and.acmove.gt.nequil) then
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

             if (mod(acmove,ndumph).eq.0.and.myid.eq.rootid) then
                do i=1,nentry
                   write(io_hist, * ) (ncmt_list(i,imolty), imolty=1,nmolty), eng_list(i)
                end do
                nentry = 0

                io_ndis=get_iounit()
                do imolty=1,nmolty
                   write(file_ndis,'("n",I2.2,"dis",I1.1,A,".dat")') imolty,run_num,suffix
                   open(unit=io_ndis,access='sequential',action='write',file=file_ndis,form='formatted',iostat=jerr,status='unknown')
                   if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open ndis file '//trim(file_ndis),jerr)
                   do n=nminp(imolty),nmaxp(imolty)
                      write(io_ndis,*) n,ndist(n,imolty)
                   end do
                end do
                close(io_ndis)
             end if
          end if

          if (lstop) then
             if (acmove.ge.nstep) exit MC_cycle
          end if
       end do
!*********************************************
! ends loop over chains                     **
!*********************************************
       ! perform periodic operations
       call monper()

       time_cur = time_now()
       if (time_cur-time_prev.gt.checkpoint_interval) then
          time_prev = time_cur
          ichkpt = ichkpt+1
          if (ichkpt.gt.checkpoint_copies) ichkpt = 1
          ! write out the restart configurations to save-config file
          if (myid.eq.rootid) then
             call dump('save-config.'//integer_to_string(ichkpt))
             call write_checkpoint_main('save-stats.'//integer_to_string(ichkpt))
          end if
       end if
    end do MC_cycle
!********************************************
! ends the loop over cycles                **
!********************************************
    call cnt()

    if (lneighbor.and.myid.eq.rootid) then
       write(io_cnt,*) 'ii:',tmcc,(neigh_cnt(i),i=1,nchain)
       close(io_cnt)
    end if

    call output_safecbmc()

    ! do bin = 1,1000
    !    write(26,*) binstep*(real(bin,dp)-0.5E0_dp),profile(bin)/nstep
    ! end do

    if (myid.eq.rootid) then
       if (io_movie.ge.0) close(io_movie)
       do ibox=1,nbox
          if (io_box_movie(ibox).ge.0) close(io_box_movie(ibox))
       end do
       if (io_solute.ge.0) close(io_solute)
       if (io_cell.ge.0) close(io_cell)
       close(io_traj)

       if (lgrand) then
          close(io_flt)
          close(io_hist)
       end if

! -------------------------------------------------------------------
       ! Output MC move statistics
       write(io_output,*)
       write(io_output,*) '+++++ end of markov chain +++++'

       call output_translation_rotation_stats(io_output)

       if (pmcb .gt. 0.0E0_dp) then
          call output_cbmc_stats(io_output)
       end if

       ! write some information about volume performance ***
       if (lgibbs .or. lnpt) then
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
! -------------------------------------------------------------------
    ! Check final energies
    do ibox=1,nbox
       if ( ldielect ) then
          ! store old dipole moment
          call dipole(ibox,2)
       end if

       ! checks final value of the potential energy is consistent ***
       call sumup(ovrlap,v,ibox,.false.)
       vend(ibox) = v(1)

       ! need to check
       if (myid.eq.rootid) then
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

! -------------------------------------------------------------------
    if (myid.eq.rootid) then
       ! Write out final configurations
       write(io_output,*)
       write(io_output,"(' vstart       =',3f24.10)") (vstart(i) ,i=1,nbox)
       write(io_output,"(' vend         =',3f24.10)") (vend(i)   ,i=1,nbox)
       write(io_output,"(' vbox         =',3f24.10)") (vbox(1,i) ,i=1,nbox)
       write(io_output,*)

       ! write out the final configuration for each box, Added by Neeraj 06/26/2006 3M ***
       io_config=get_iounit()
       do ibox = 1,nbox
          write(file_config,'("box",I1.1,"config",I1.1,A,".xyz")') ibox,run_num,suffix
          open(unit=io_config,access='sequential',action='write',file=file_config,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open box config file '//trim(file_config),jerr)

          nummol = 0
          do i = 1,nchain
             if (nboxi(i).eq.ibox) then
                nummol = nummol + nunit(moltyp(i))
             end if
          end do
          write(io_config,*) nummol
          write(io_config,*)
          do i = 1,nchain
             if(nboxi(i).eq.ibox) then
                imolty = moltyp(i)
                do ii = 1,nunit(imolty)
                   ntii = ntype(imolty,ii)
                   write(io_config,'(a4,5x,3f15.4)') chemid(ntii),rxu(i,ii),ryu(i,ii),rzu(i,ii)
                end do
             end if
          end do
          close(io_config)
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

! -------------------------------------------------------------------
       ! Calculate and write out running averages ***
       do ibox=1,nbox
          if ( lpbcz ) then
             do itype = 1, nmolty
                ! number density
                acdens(ibox,itype)=1000.0E0_dp*acdens(ibox,itype)
             end do

             ! sum over all types of molecules
             temacd = sum(acdens(ibox,1:nmolty))

             ! molar volume in mL/mol
             molvol(ibox) = N_Avogadro*1E-21_dp/temacd

             ! specific density
             temspd = 0.0E0_dp
             do itype = 1, nmolty
                temspd = temspd + (acdens(ibox,itype)*masst(itype)/(N_Avogadro*1E-21_dp))
             end do
             speden(ibox) = temspd
          else
             do itype = 1, nmolty
                ! number density
                acdens(ibox,itype)=100.0E0_dp*acdens(ibox,itype)
             end do

             temacd = sum(acdens(ibox,1:nmolty))

             ! molar volume
             molvol(ibox) = 100.0E0_dp / temacd
          end if

          ! energies
          acvsq(:,ibox) = acvsq(:,ibox) - acv(:,ibox) ** 2
          acvkjmol(:,ibox) = acvkjmol(:,ibox)*R_gas/1000_dp

          ! if ( ldielect ) then
          !    flucmom(ibox) = acdipolesq(4,ibox)-acdipole(4,ibox)*acdipole(4,ibox)
          !    flucmom(ibox) = 6.9994685465110493E5_dp*flucmom(ibox)*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
          !    flucmom2(ibox) = acdipolesq(4,ibox)-acdipole(1,ibox)*acdipole(1,ibox)-acdipole(2,ibox)*acdipole(2,ibox)-acdipole(3,ibox)*acdipole(3,ibox)
          !    flucmom2(ibox) = 6.9994685465110493E5_dp*flucmom2(ibox)*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
          ! end if

          ! chemical potential
          ! This expression is suitable only for NPT ensemble
          do itype = 1, nmolty
             if(bnchem(ibox,itype).gt.0) then
                debroglie = debroglie_factor*sqrt(beta/masst(itype))
                ! determine how many steps it takes to grow molecule
                ! not counting the first inserted bead
                igrow = nugrow(itype)
                if (.not. lrigid(itype)) then
                   call schedule(igrow,itype,steps,1,0,2)
                   tmp=real(nchoi1(itype)*(nchoi(itype)**steps)*nchoih(itype),dp)*(debroglie**3)
                else
                   call schedule(igrow,itype,steps,1,0,4)
                   tmp=real(nchoi1(itype)*nchoir(itype)*(nchoi(itype)**steps)*nchoih(itype),dp)*(debroglie**3)
                end if
                acchem(ibox,itype)=(-1.0E0_dp/beta)*log((acchem(ibox,itype)/bnchem(ibox,itype))/tmp)
                itel = 2 + nEnergy + itype
                baver(itel,ibox,:)=(-1.0E0_dp/beta)*log(baver(itel,ibox,:)/tmp)
             end if
          end do
       end do

       write(io_output,"(' Averages and fluctuations',21x,4(a11,i1))") ('       Box ',i,i=1,nbox)
       write(io_output,*)
       write(io_output,"(' pressure                               [kPa] =',3f12.2)") (acpres(i),i=1,nbox)
       write(io_output,"(' pressure                  [simulation units] =',3f12.6)") ((acpres(i)*MPa2SimUnits*1E-3_dp),i=1,nbox)
       write(io_output,"(' surface tension                       [mN/m] =',3f12.4)") (acsurf(i),i=1,nbox)
       do itype = 1,nmolty
          write(io_output,"(' chem. potential of type    ',i4,'          [K] =',3f12.3)") itype,(acchem(i,itype),i=1,nbox)
       end do
       write(io_output,*)

       do i = 1,3
          write(io_output,"(' boxlength                                [A] =',3f12.3)") (acboxl(ibox,i),ibox=1,nbox)
       end do

       do ibox = 1, nbox
          if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
             do i = 1,3
                write(io_output,"(' box angle                                deg =',3f12.3)") acboxa(ibox,i)*raddeg
             end do
          end if
       end do

       do itype = 1, nmolty
          write(io_output,"(' no. of chains of type      ',i4,'              =',3f12.3)") itype,(acnbox(i,itype),i=1,nbox)
       end do
       if ( lpbcz ) then
          write(io_output,"(' molar volume                      [cm^3/mol] =',3f12.3)") (molvol(i),i=1,nbox)
          write(io_output,"(' specific density                    [g/cm^3] =',3f12.6)") (speden(i),i=1,nbox)
          do itype = 1, nmolty
             write(io_output,"(' number density of type     ',i4,' [chain/nm^3] =',3f12.5)") itype,(acdens(i,itype),i=1,nbox)
             if ( lexpand(itype) ) then
                do itype2 = 1, numcoeff(itype)
                   write(io_output,"(' number density of type     ',i4,' eetype ',i4,'  =', 2f12.5)") itype,itype2,acdens(itype,itype)*acnbox2(itype,itype,itype2)/(acnbox(itype,itype)*acmove),acnbox2(itype,itype,itype2)/(acnbox(itype,itype)*acmove)
                end do
             end if
          end do
       else
          write(io_output,"(' area per chain                     [A^2/chain] =',3f12.4)") (molvol(i),i=1,nbox)
          do itype = 1, nmolty
             write(io_output,"(' number density of type     ',i4,' [chain/nm^2] =',3f12.6)") itype,(acdens(i,itype),i=1,nbox)
          end do
       end if
       do itype = 1, nmolty
          write(io_output,"(' molfraction of type        ',i4,'              =',3f12.7)") itype,(molfra(i,itype),i=1,nbox)
       end do
       do itype = 1, nmolty
          write(io_output,"(' mean sete length of type   ',i4,'        [A^2] =',3f12.3)") itype,(asetel(i,itype),i=1,nbox)
       end do
       write(io_output,*)
       do j=1,11
          ! only 1 to 11 is the energy information
          write(io_output,"(a15,'[K per system and kJ/mol per chain] =',3(f14.2,f12.2))") vname(j),acv(j,1:nbox),acvkjmol(j,1:nbox)
       end do

       write(io_output,*)
       write(io_output,"(' fluctuation in <vtot>                          =',3f12.2)") (sqrt(acvsq(1,i)),i=1,nbox)
       write(io_output,*)

       ! Output 2nd virial coefficient data
2000   if (lvirial) then
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
             write(io_output,*) 'bvirial ',-N_Avogadro*1E-24_dp*twopi* stepvir*bvirial,' [cm^3 / mole]'

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
             write(io_output,*) 'bvirial ',N_Avogadro*1E-24_dp*bvirial,' [cm^3 / mole]'
          end do
       end if

       ! solute values
       write(io_output,*) 'type  box     vinter      vintra      vtor', '        vbend       vtail'

       do itype = 1, nmolty
          do ibox = 1, nbox
             if (solcount(ibox,itype).gt.0) then
                write(io_output,"(i5,i5,3f12.5,3f12.5,3f12.5,3f12.5,3f12.5)") itype,ibox,avsolinter(ibox,itype),avsolintra(ibox,itype),avsoltor(ibox,itype),avsolbend(ibox,itype),avsolelc(ibox,itype)
             else
                write(io_output,"(i5,i5,3f12.5,3f12.5,3f12.5,3f12.5,3f12.5)") itype,ibox,0.0,0.0,0.0,0.0,0.0
             end if
          end do
       end do

! -------------------------------------------------------------------
       ! Calculate statistical uncertainties from block averages
       if ( nblock .ge. 2 ) then
          ! global averages -
          do i = 1,nprop
             do ibox = 1,nbox
                call calculate_statistics(baver(i,ibox,1:nblock),aver(i,ibox),stdev(i,ibox),errme(i,ibox))
             end do
          end do
          do i = 1,nprop1
             do ibox = 1,nbox-1
                do jbox = ibox+1,nbox
                   call calculate_statistics(baver1(i,ibox,jbox,1:nblock),aver1(i,ibox,jbox),stdev1(i,ibox,jbox),errme1(i,ibox,jbox))
                end do
             end do
          end do

          ! write out the heat of vaporization and solubility parameters
          do ibox = 1,nbox-1
             do jbox = ibox+1,nbox
                write(io_output,"(' H_vap      [kJ/mol] btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(1,ibox,jbox),stdev1(1,ibox,jbox),errme1(1,ibox,jbox)
                write(io_output,"(' H_vap LJ  [kJ/mol] btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(2,ibox,jbox),stdev1(2,ibox,jbox),errme1(2,ibox,jbox)
                write(io_output,"(' H_vap Coul [kJ/mol] btwn box  ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(3,ibox,jbox),stdev1(3,ibox,jbox),errme1(3,ibox,jbox)
                ! write(io_output,"(' DeltaU Ext [kJ/mol] btwn box   ',i4,' and',i4, ' =',3f15.4)") ibox,jbox,acsolpar(10,ibox,jbox),stdev1(10,ibox,jbox), errme1(10,ibox,jbox)
                write(io_output,"(' pdV        [kJ/mol] btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(11,ibox,jbox),stdev1(11,ibox,jbox),errme1(11,ibox,jbox)
                write(io_output,"(' CED [cal/cc]   btwn box        ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(4,ibox,jbox),stdev1(4,ibox,jbox),errme1(4,ibox,jbox)
                write(io_output,"(' CED_LJ[cal/cc] btwn box        ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(5,ibox,jbox),stdev1(5,ibox,jbox),errme1(5,ibox,jbox)
                write(io_output,"(' CED_Coul[cal/cc] btwn box      ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(6,ibox,jbox),stdev1(6,ibox,jbox),errme1(6,ibox,jbox)
                write(io_output,"(' HSP [(cal/cc)^1/2]  btwn box   ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(7,ibox,jbox),stdev1(7,ibox,jbox),errme1(7,ibox,jbox)
                write(io_output,"(' HSP_LJ[(cal/cc)^1/2] btwn box  ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(8,ibox,jbox),stdev1(8,ibox,jbox),errme1(8,ibox,jbox)
                write(io_output,"(' HSP_Cou[(cal/cc)^1/2] btwn box ',i4,' and',i4, ' =', 3f15.4)") ibox,jbox,acsolpar(9,ibox,jbox),stdev1(9,ibox,jbox),errme1(9,ibox,jbox)
             end do
          end do

          ! specific density
          do ibox = 1, nbox
             write(io_output,"(' specific density    box ',i3, ' = ',3E12.5)") ibox,aver(1,ibox),stdev(1,ibox),errme(1,ibox)
          end do

          ! system volume
          itel = nEnergy + 4*nmolty + 4
          do ibox = 1, nbox
             write(io_output,"(' system volume       box ',i3, ' = ',3E12.5)") ibox,aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
          end do

          ! pressure
          do ibox = 1, nbox
             write(io_output,"(' pressure            box ',i3, ' = ',3g12.5)") ibox,acpres(ibox),stdev(2,ibox),errme(2,ibox)
          end do

          ! surface tension
          itel = nEnergy + 4*nmolty + 3
          do ibox = 1, nbox
             write(io_output,"(' surface tension     box ',i3, ' = ',3f12.5)") ibox,acsurf(ibox),stdev(itel,ibox),errme(itel,ibox)
          end do

          write(io_output,*)
          ! energies
          do ibox = 1, nbox
             do j=3,13
                ! only 1 to 10 is the energy information
                write(io_output,"(a15 ,' box ',i3, ' = ',3E14.5)") vname(j-2),ibox,aver(j,ibox),stdev(j,ibox),errme(j,ibox)
             end do
          end do

          write(io_output,*)

          ! Enthalpy
          do ibox = 1,nbox
             j = nEnergy + 4*nmolty + 5
             write(io_output, "(' Enthalpy Inst.[kJ/mol] for box',i3,' =',3(f12.4))") ibox,acEnthalpy(ibox),stdev(j,ibox),errme(j,ibox)
             j = nEnergy + 4*nmolty + 6
             write(io_output,"(' Enthalpy Ext. [kJ/mol] for box',i3,' =',3(f12.4))") ibox,acEnthalpy1(ibox),stdev(j,ibox),errme(j,ibox)
          end do
          write(io_output,*)

          ! residual heat capacity, in (J2/mol)
          if (.not.lgibbs) then
             tmp=(enthalpy2-enthalpy*enthalpy)/real(nchain,dp)*R_gas/(temp**2)
             if(lnpt) then
                write(io_output,*) 'Cp residual(J/Kmol) =',tmp
                write(io_output,*) ' H2=',enthalpy2
                write(io_output,*) ' H=',enthalpy
             else
                write(io_output,*) 'Cv residual(J/Kmol) =',tmp
                write(io_output,*) ' E2=',enthalpy2
                write(io_output,*) ' E=',enthalpy
             end if
          end if

          ! chemical potential
          do itype = 1, nmolty
             itel = 2 + nEnergy + itype
             do ibox = 1, nbox
                write(io_output,"(' chemical potential  itype ',i3,' box ',i3, ' = ',3f12.3)") itype,ibox,aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
             end do
          end do

          ! square end-to-end length
          do itype = 1, nmolty
             itel = 2 + nEnergy + nmolty + itype
             do ibox = 1, nbox
                write(io_output,"(' mean sete length    itype ',i3,' box ',i3, ' = ',3f12.3)") itype,ibox,aver(itel,ibox) ,stdev(itel,ibox),errme(itel,ibox)
             end do
          end do

          ! number density
          do itype = 1, nmolty
             itel = 2 + nEnergy + 2*nmolty + itype
             do ibox = 1, nbox
                if ( lpbcz ) then
                   write(io_output,"(' number density      itype ',i3,' box ',i3, ' = ',3E12.5)") itype,ibox,1.0E3_dp*aver(itel,ibox),1.0E3_dp*stdev(itel,ibox),1.0E3_dp*errme(itel,ibox)
                else
                   write(io_output,"(' number density      itype ',i3,' box ',i3, ' = ',3E12.5)") itype,ibox,1.0E2_dp*aver(itel,ibox),1.0E2_dp*stdev(itel,ibox),1.0E2_dp*errme(itel,ibox)
                end if
                if ( lexpand(itype) .and. acnbox(ibox,itype) .gt. 0.5) then
                   do itype2 = 1, numcoeff(itype)
                      molfrac = acnbox2(ibox,itype,itype2) /(acmove*acnbox(ibox,itype))
                      write(io_output,"(' number density      itype ',i3,' typ ',i3, ' = ',2E12.5)") itype,itype2,1.0E3_dp*aver(itel,ibox)*molfrac,molfrac
                   end do
                end if
             end do
          end do

          ! molfraction
          do itype = 1, nmolty
             itel = 2 + nEnergy + 3*nmolty + itype
             do ibox = 1, nbox
                write(io_output,"(' mole fraction       itype ',i3,' box ',i3, ' = ',3f12.7)") itype,ibox,aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
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
                   gconst = R_gas/(1E3_dp*beta)
                   do itype = 1,nmolty
                      itel = 2 + nEnergy + 2 * nmolty + itype
                      ostwald = aver(itel,il)/aver(itel,ig)
                      stdost  = ostwald * sqrt( (stdev(itel,il)/aver(itel,il))**2 + (stdev(itel,ig)/aver(itel,ig))**2 )
                      ! write(42,*) nunit(itype),ostwald,stdost
                      ! write(43,*) nunit(itype),-(gconst*log(ostwald)) + (eta2(ig,itype) - eta2(il,itype)) / 120.27167,gconst*stdost/ostwald
                      write(io_output,"('Ostwald Coefficient  itype ',i3,' between box ',i2, ' and ',i2,f18.6,f18.6)") itype,ig,il,ostwald,stdost
                      write(io_output,"('Free Enrgy of Transf itype ',i3,' between box ',i2, ' and ',i2,f18.6,f18.6,' kJ/mol')") itype,ig,il,-(gconst*log(ostwald))+(eta2(ig,itype)-eta2(il,itype))*R_gas*1E-3_dp,gconst*stdost/ostwald
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
                write(io_output,"(2x,i2,15(2x,e10.3))") nbl,baver(3,ibox,nbl),baver(1,ibox,nbl),baver(2,ibox,nbl),baver(3+nEnergy+4*nmolty,ibox,nbl),(baver(2+nEnergy+3*nmolty+zzz,ibox,nbl),zzz=1,nmolty)
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
       !    call analysis(2)
       ! end if

       ! ee prob
       IF(lexpee) then
          write(io_output,*)
          write(io_output,*) 'probability of each mstate in ee'
          do nnn = 1, fmstate
             write(io_output,"(i5,1x,i10)") nnn,ee_prob(nnn)
          end do
       end if

       write(io_output,*) 'Program ended at ',time_date_str()
       close(io_output)
    end if

    return
  end subroutine monola

!> \brief Read input data and initializes the positions
!>
!> reads a starting configuration from unit 7
!> calculates interaction table
  subroutine readdat(file_in)
    use const_math,only:onepi
    use const_phys,only:MPa2SimUnits,debroglie_factor
    use util_random,only:ranset
    use util_string,only:integer_to_string,real_to_string
    use util_runtime,only:err_exit
    use util_timings,only:time_date_str
    use util_files,only:get_iounit
    use util_search,only:indexOf
    use util_memory,only:reallocate
    use sim_particle,only:allocate_neighbor_list,ctrmas
    use zeolite
    use energy_kspace,only:calp,allocate_kspace
    use energy_pairwise,only:init_energy_pairwise,rzero,epsnx,type_2body
    use energy_intramolecular,only:bonds,angles,dihedrals,allocate_energy_bonded
    use energy_3body,only:readThreeBody
    use energy_4body,only:readFourBody
    use energy_sami
    use moves_simple,only:init_moves_simple,averageMaximumDisplacement
    use moves_volume,only:init_moves_volume
    use moves_cbmc,only:init_cbmc,allocate_cbmc,read_safecbmc,llplace
    use moves_ee,only:init_ee,numcoeff,sigm,epsil
    use transfer_shared,only:read_transfer
    use transfer_swap,only:init_swap
    use transfer_swatch,only:init_swatch

    character(LEN=*),intent(in)::file_in

    character(LEN=default_path_length)::file_input,file_restart,file_struct,file_run,file_movie,file_solute,file_traj
    namelist /io/ file_input,file_restart,file_struct,file_run&
     ,file_movie,file_solute,file_traj,io_output,checkpoint_interval,checkpoint_copies,use_checkpoint

    !* PARAMETER FOR IONIC SYSTEMS
    logical::lionic !< if LIONIC=.TRUE. System contains charged species, so system may not neutral

    real,allocatable::ofscale(:),ofscale2(:),qbox(:)
    integer,allocatable::ncarbon(:),inclmol(:),inclbead(:,:),inclsign(:),ainclmol(:),ainclbead(:,:),a15t(:),idummy(:),temphe(:),nures(:),k_max_l(:),k_max_m(:),k_max_n(:)
    logical,allocatable::lhere(:)

    character(LEN=default_path_length)::file_box_movie,file_cell
    integer::io_input,io_restart,jerr,seed,nijspecial,ispecial,jspecial,ndum,ij,ii,jj,i,j,k,ncres,nmtres,iensem,inpbc,nmcount,im,ibox,tcount,izz,z,zzz
    logical::lprint,lecho,lverbose,L_Ewald_Auto,lmixlb,lmixjo,lsetup,linit,lreadq,lpolar,lqqelect,lee,lxyz,lfound
    real::fqtemp,aspecd,bspecd,dum,pm,pcumu,w(3),debroglie,qtot,min_boxl
    ! variables added (3/24/05) for scaling of 1-4 interactions
    integer::nexclu,inclnum,ainclnum
    ! variables for histograms
    integer::temnc,imol,iutemp,imolty,itype,ipair
    ! Variables added (6/30/2006) for fort.4 consistency check
    integer::numvib,numbend,numtor,vib1,bend2,bend3,tor2,tor3,tor4
    integer::vibtype,bendtype,tortype
    ! real::temx(nmax,numax),temy(nmax,numax),temz(nmax,numax)
    ! KM variable added when analysis removed
    integer::nhere
! -------------------------------------------------------------------
    io_input=get_iounit()
    open(unit=io_input,access='sequential',action='read',file=file_in,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open input file '//trim(file_in),myid+1)
    end if

    file_input='fort.4'
    file_restart='fort.77'
    file_struct='input_struc.xyz'
    file_run='run1a.dat'
    file_movie='movie1a.dat'
    file_solute='fort.11'
    file_traj='fort.12'
    use_checkpoint=.false.

    read(UNIT=io_input,NML=io,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) then
       call err_exit(__FILE__,__LINE__,'reading namelist: io',jerr)
    end if

    if (myid.eq.rootid.and..not.use_checkpoint) then
       lprint=.true.
    else
       lprint=.false.
    end if

    ! rewind(io_input)
    ! read(UNIT=io_input,NML=system,iostat=jerr)
    ! if (jerr.ne.0.and.jerr.ne.-1) then
    !    call err_exit(__FILE__,__LINE__,'reading namelist: system',jerr)
    ! end if
    ! nchain=nchain+2

    ! Output unit: if 6 or 0, write to stdout/stderr; otherwise, write to a user designated file
    if (io_output.eq.5) then
       call err_exit(__FILE__,__LINE__,'unit 5 is for standard input',myid+1)
    else if(io_output.ne.6.and.io_output.ne.0.and.myid.eq.rootid) then
       io_output=get_iounit()
       open(unit=io_output,access='stream',action='write',file=file_run,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) then
          call err_exit(__FILE__,__LINE__,'cannot open output file '//trim(file_run),jerr)
       end if
    end if
    close(io_input)

! -------------------------------------------------------------------
    io_input=get_iounit()
    open(unit=io_input,access='sequential',action='read',file=file_input,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open main input file',myid+1)
    end if

    read(io_input,*)
    read(io_input,*) seed
    ! initialize random number generator
    call ranset(seed,numprocs)

! -------------------------------------------------------------------
    ! read echoing and long output flags
    read(io_input,*)
    read(io_input,*) ndum,lecho,lverbose,run_num,suffix
    read(io_input,*)
    read(io_input,*) lnpt,lgibbs,lgrand,lanes,lvirial,lmipsw,lexpee
    read(io_input,*)
    read(io_input,*) lijall,lchgall,lewald,ldielect,ltailc,lshift,ltailcZeo

    ! To add or remove helium atoms
    read(io_input,*)
    read(io_input,*) L_add,N_add,N_box2add,N_moltyp2add
    read(io_input,*)
    read(io_input,*) L_sub,N_sub,N_box2sub,N_moltyp2sub

    ! read whether to compute electrostatic interaction or not during CBMC/SWAP
    read(io_input,*)
    read(io_input,*) L_Coul_CBMC,lcutcm,ldual,lneigh
    ! read the number of unitcell replicated in each directions (a, b, and c)
    read(io_input,*)
    read(io_input,*) Num_cell_a,Num_cell_b,Num_cell_c
    ! read run information
    read(io_input,*)
    read(io_input,*) nstep, lstop, lpresim, iupdatefix
    ! read torsion, decide whether to use torsion in function form or Table
    read(io_input,*)
    read(io_input,*) lslit,lexzeo,lzgrid,lelect_field,ljoe,lsami,lmuir,lpsurf,lgraphite,lcorreg,llj,lexpsix,lmmff,lninesix,lgenlj,lgaro,lionic
    read(io_input,*)
    read(io_input,*) L_tor_table,L_spline,L_linear,L_vib_table,L_bend_table,L_vdW_table,L_elect_table

    ! KM ldielect writes to fort.27
    ! if (ldielect) then
    !    open(17,file=file_dipole,status='unknown')
    !    write(17,*) '# step  ibox   dipole_x   dipole_y   dipole_z'
    ! end if

    ! KM for MPI
    ! only processor 0 writes data
    if ( lecho.and.lprint) then
       if (lverbose) then
          write(io_output,*) 'Program started at ',time_date_str()
          write(io_output,*) 'Number of processors: ', numprocs
          write(io_output,*) thread_num,' threads per processor'
          write(io_output,*) 'Random number seed: ',seed
          write(io_output,*) 'L_Coul_CBMC:',L_Coul_CBMC
          write(io_output,*) 'Number of unit cells in a dir = ', Num_cell_a
          write(io_output,*) 'Number of unit cells in b dir = ', Num_cell_b
          write(io_output,*) 'Number of unit cells in c dir = ', Num_cell_c

          if (lstop) then
             write(io_output,*) 'number of steps:',nstep
          else
             write(io_output,*) 'number of cycles:',nstep
          end if
          write(io_output,*) 'lstep:',lstop
          write(io_output,*) 'lpresim:',lpresim
          write(io_output,*) 'iupdatefix:',iupdatefix
          write(io_output,*) 'L_tor_table:',L_tor_table
          write(io_output,*) 'L_spline:',L_spline
          write(io_output,*) 'L_linear:',L_linear
          write(io_output,*) 'L_vib_table:', L_vib_table
          write(io_output,*) 'L_bend_table:', L_bend_table
          write(io_output,*) 'L_vdW_table:', L_vdW_table
          write(io_output,*) 'L_elect_table:', L_elect_table
       else
          write(io_output,*) time_date_str(),numprocs,seed
          write(io_output,*) L_Coul_CBMC,nstep, lstop, lpresim, iupdatefix
          write(io_output,*) L_tor_table, L_spline,  L_linear
          write(io_output,*) L_vib_table, L_bend_table, L_vdW_table,  L_elect_table
       end if
    end if

    read(io_input,*)
    read(io_input,*) nbox

    nbxmax=nbox+1
    npabmax=nbxmax*(nbxmax-1)/2
    call allocate_cell()
    call allocate_sim_cell()
    call allocate_kspace()

    if ( lecho.and.lprint ) write(io_output,*) 'number of boxes in the system:',nbox

    read(io_input,*)
    read(io_input,*) (express(ibox),ibox=1,nbox)

    read(io_input,*)
    read(io_input,*) (ghost_particles(ibox),ibox=1,nbox)

    read(io_input,*)
    read(io_input,*) L_Ewald_Auto

    read(io_input,*)
    read(io_input,*) temp, fqtemp,(Elect_field(ibox),ibox=1,nbox)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'temperature:',temp,' K'
          write(io_output,*) 'external pressure:',express(1:nbox),' MPa'
          write(io_output,*) 'Ghost particles: ',ghost_particles(1:nbox)
          write(io_output,*) 'fluctuating charge temperature:',fqtemp,' K'
          write(io_output,*) 'Electric field in z direction:', Elect_field(1:nbox),' V/A'
       else
          write(io_output,*) temp, (express(ibox),ibox=1,nbox), fqtemp, Elect_field
       end if
    end if

    do ibox = 1,nbox
       express(ibox)  = express(ibox)*MPa2SimUnits
    end do

    ! read the analysis information
    read(io_input,*)
    read(io_input,*) ianalyze,nbin,lrdf,lintra,lstretch,lgvst,lbend,lete,lrhoz,bin_width
    if (lecho.and.lprint) then
       if(lverbose) then
          write(io_output,*) 'ianalyze:', ianalyze
          write(io_output,*) 'nbin', nbin
          write(io_output,*) 'lrdf:',lrdf
          write(io_output,*) 'lintra:',lintra
          write(io_output,*) 'lstretch:', lstretch
          write(io_output,*) 'lgvst:' ,lgvst
          write(io_output,*) 'lbend:' ,lbend
          write(io_output,*) 'lete:' ,lete
          write(io_output,*) 'lrhoz:', lrhoz
          write(io_output,*) 'bin_width:' ,bin_width
       else
          write(io_output,*) ianalyze,nbin,lrdf,lintra,lstretch,lgvst,lbend ,lete, lrhoz, bin_width
       end if
    end if

! -------------------------------------------------------------------
! set up constants and conversion factors ***
    lee = .false.
    qtot = 0.0E0_dp

    beta = 1.0E0_dp / temp
    fqbeta = 1.0E0_dp / fqtemp
! ------------------------------------------------------------------
    read(io_input,*)
    read(io_input,*) iprint, imv, iratio, iblock, idiele, L_movie_xyz, iheatcapacity

    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'iprint:',iprint
          write(io_output,*) 'imv:',imv
          write(io_output,*) 'iratio:',iratio
          write(io_output,*) 'iblock:',iblock
          write(io_output,*) 'idiele:',idiele
          write(io_output,*) 'L_movie_xyz',L_movie_xyz
       else
          write(io_output,*) iprint, imv, iratio, iblock, idiele,L_movie_xyz
       end if
    end if

    ! read information for histogram output (added 8/30/99)
    if (lgrand) then
       read(io_input,*)
       read(io_input,*) nequil,ninstf,ninsth, ndumph
       if ( lecho.and.lprint ) then
          if (lverbose) then
             write(io_output,*) 'nequil:',nequil
             write(io_output,*) 'ninstf:',ninstf
             write(io_output,*) 'ninsth:',ninsth
             write(io_output,*) 'ndumph:',ndumph
             write(io_output,*) 'run_num:',run_num
             write(io_output,*) 'suffix:',suffix
          else
             write(io_output,*) nequil,ninstf,ninsth,ndumph,run_num,suffix
          end if
       end if
    end if

    ! read system information
    do i = 1,nbox
       read(io_input,*)
       read(io_input,*) boxlx(i),boxly(i),boxlz(i),lsolid(i),lrect(i),kalp(i),rcut(i),rcutnn(i),numberDimensionIsIsotropic(i)
       if (i.eq.1 .and. lexzeo) then
          ! load positions of zeolite atoms
          call zeocoord(file_in,lprint)
          if (lprint) write(io_output,*) ' note zeolite determines the box size !'
       end if

       if ( lecho.and.lprint ) then
          if (lverbose) then
             write(io_output,*) 'box:',i
             write(io_output,*) '   boxlx:',boxlx(i),' A'
             write(io_output,*) '   boxly:',boxly(i),' A'
             write(io_output,*) '   boxlz:',boxlz(i),' A'
             write(io_output,*) '   lsolid:',lsolid(i)
             write(io_output,*) '   lrect:',lrect(i)
             write(io_output,*) 'neighbor list cutoff (rcutnn):',rcutnn(i), 'A'
             write(io_output,*) '   rcut:',rcut(i),'A'
             if(.not.L_Ewald_Auto) then
                write(io_output,*) '   kalp:',kalp(i)
             end if
          else
             write(io_output,*) boxlx(i),boxly(i),boxlz(i), lsolid(i),lrect(i),rcut(i),rcutnn(i)
             if (.not.L_Ewald_Auto) then
                write(io_output,*) kalp(i)
             end if
          end if
       end if
    end do

    read(io_input,*)
    read(io_input,*) nchain, nmolty

    if (lgrand) then
       nchain=nmax
       write(io_output,*)'in GCMC total number of chains set by NMAX!'
    end if

    nmax=nchain+2
    ntmax=nmolty+1
    npamax=ntmax*(ntmax-1)/2
    nprop=nEnergy+(4*ntmax)+6
    blockm=nstep/iblock
    call allocate_system()
    call allocate_neighbor_list()
    call allocate_linked_cell()
    call read_transfer(file_in,lprint)
    call init_swap()
    call init_swatch()

    allocate(io_box_movie(nbxmax),ncarbon(ntmax),idummy(ntmax),qbox(nbxmax),nures(ntmax),k_max_l(nbxmax),k_max_m(nbxmax),k_max_n(nbxmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'readdat: allocating system failed',jerr)
    end if

    ! set input arrays to zero ***
    do j=1, ntmax
       nugrow(j) = 0
       nunit(j) = 0
       do i=1, numax
          ntype(j,i) = 0
       end do
    end do

    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'number of chains:',nchain
          write(io_output,*) 'number of molecule types:',nmolty
       else
          write(io_output,*) nchain, nmolty
       end if
    end if

    read(io_input,*)
    read(io_input,*) (temtyp(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          do i = 1,nmolty
             write(io_output,*) 'number of chains of molecule type',i,':', temtyp(i)
          end do
       else
          write(io_output,*) (temtyp(i),i=1,nmolty)
       end if
    end if
    temnc = 0
    do i = 1, nmolty
       do j = 1, temtyp(i)
          temnc = temnc + 1
          moltyp(temnc) = i
       end do
    end do

    if (lgrand) then
       ! read chemical potentials (added 8/30/99 by jpotoff)
       read(io_input,*)
       read(io_input,*) (B(i),i=1,nmolty)
       if (lecho.and.lprint) then
          if (lverbose) then
             do i = 1,nmolty
                write(io_output,*) 'chemical potential for molecule type', i,':',B(i)
             end do
          else
             write(io_output,*) "B ", (B(i),i=1,nmolty)
          end if
       end if

       ! removing from here. It will be calculated once we have molecular mass
       ! to calculate debroglie wavelength.
       ! convert chemical potentials to activities
       ! do i=1,nmolty
       !    B(i) = exp(B(i)/temp)
       ! end do
    end if

    read(io_input,*)
    read(io_input,*) lmixlb, lmixjo
    if (lmixlb .and. lmixjo) then
       call err_exit(__FILE__,__LINE__,'cant use both combining rules!',myid+1)
    end if
    if ( lecho.and.lprint ) then
       if (lverbose) then
          if (lmixlb) then
             write(io_output,*) 'Lorentz-Berthelot combining rules apply'
          else
             write(io_output,*) 'Jorgensen combining rules apply'
          end if
          write(io_output,*) '   lmixlb:',lmixlb,' lmixjo:',lmixjo
       else
          write(io_output,*) lmixlb,lmixjo
       end if
    end if
    ! read special combining rule information
    read(io_input,*)
    read(io_input,*) nijspecial
    if (lecho.and.lprint) then
       if (lverbose) then
          write(io_output,*) 'number of special combining parameters:', nijspecial
       else
          write(io_output,*) nijspecial
       end if
    end if
    read(io_input,*)
    if ( nijspecial .eq. 0 ) then
       read(io_input,*)
    else
       do i = 1,nijspecial
          read(io_input,*) ispecial,jspecial,aspecd,bspecd
       end do
    end if

    read(io_input,*)
    read(io_input,*) rmin,softcut,rcutin,rbsmax,rbsmin
    if ( lecho .and.lprint) then
       if (lverbose) then
          write(io_output,*) 'minimum cutoff (rmin):',rmin,' A'
          write(io_output,*) 'softcut:',softcut
          write(io_output,*) 'CBMC inner cutoff (rcutin):',rcutin,' A'
          write(io_output,*) 'AVBMC outer cutoff (rbsmax):',rbsmax,' A'
          write(io_output,*) 'AVBMC inner cutoff (rbsmin):',rbsmin,' A'
       else
          write(io_output,*) rmin, softcut ,rcutin,rbsmax,rbsmin
       end if
    end if
    do i = 1, nbox
       if( rcut(i)/boxlx(i) .gt. 0.5E0_dp) then
          call err_exit(__FILE__,__LINE__,'rcut > 0.5*boxlx',myid+1)
       end if
    end do

    softlog = 10.0E0_dp**(-softcut)
    vol_eff = (4.0E0_dp/3.0E0_dp)*onepi*(rbsmax*rbsmax*rbsmax-rbsmin*rbsmin*rbsmin)

    ! set up the forcefield and the masses
    call init_energy_pairwise(file_in,lmixlb,lmixjo,lprint)

    allocate(lhere(nntype),temphe(nntype),stat=jerr)
    lhere=.false.
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'readdat: allocating potential failed',jerr)
    end if

    numax=0
    ! read bead potential information
    do imol = 1, nmolty
       read(io_input,*)
       read(io_input,*) nunit(imol),nugrow(imol),ncarbon(imol),nmaxcbmc(imol),iurot(imol),lelect(imol),lflucq(imol),lqtrans(imol),lexpand(imol),lavbmc1(imol),lavbmc2(imol),lavbmc3(imol),fqegp(imol)

       if (nunit(imol).gt.numax) then
          numax=nunit(imol)
          if (numax.gt.ubound(ntype,2)) then
             call reallocate(ntype,1,ntmax,1,2*numax)
             call reallocate(leaderq,1,ntmax,1,2*numax)
             call reallocate(lplace,1,ntmax,1,2*numax)
             call reallocate(lrigi,1,ntmax,1,2*numax)
             call reallocate(invib,1,ntmax,1,2*numax)
             call reallocate(itvib,1,ntmax,1,2*numax,1,6)
             call reallocate(ijvib,1,ntmax,1,2*numax,1,6)
             call reallocate(inben,1,ntmax,1,2*numax)
             call reallocate(itben,1,ntmax,1,2*numax,1,12)
             call reallocate(ijben2,1,ntmax,1,2*numax,1,12)
             call reallocate(ijben3,1,ntmax,1,2*numax,1,12)
             call reallocate(intor,1,ntmax,1,2*numax)
             call reallocate(ittor,1,ntmax,1,2*numax,1,12)
             call reallocate(ijtor2,1,ntmax,1,2*numax,1,12)
             call reallocate(ijtor3,1,ntmax,1,2*numax,1,12)
             call reallocate(ijtor4,1,ntmax,1,2*numax,1,12)
             call reallocate(irotbd,1,2*numax,1,ntmax)
             call reallocate(pmrotbd,1,2*numax,1,ntmax)
          end if
       end if

       read(io_input,*)
       read(io_input,*) maxgrow(imol),lring(imol),lrigid(imol),lrig(imol),lsetup,isolute(imol),(eta2(i,imol), i=1,nbox)

       read(io_input,*)
       read(io_input,*) lq14scale(imol),qscale(imol)

       if (lring(imol)) then
          read(io_input,*)
          read(io_input,*) iring(imol)
       else
          iring(imol) = nunit(imol)
       end if

       do i = 1, nunit(imol)
          lrigi(imol,i) = .false.
       end do

       ! irig is the site rigid sites will be grown from
       ! and frig will be the previous site (not kept rigid)

       if (lrig(imol)) then
          read(io_input,*)
          read(io_input,*) nrig(imol)

          if (nrig(imol).gt.0) then
             ! read in specific points to keep rigid in growth
             read(io_input,*)
             do i = 1, nrig(imol)
                read(io_input,*) irig(imol,i),frig(imol,i)
                lrigi(imol,irig(imol,i)) = .true.
             end do
          else
             ! we will pick irig at random in each case if nrig = 0
             read(io_input,*)
             read(io_input,*) nrigmin(imol),nrigmax(imol)

             ! nrigmin is the minimum amount of the chain to keep rigid
             ! nrigmax is the maximum

          end if
       end if

       if (lrigid(imol)) then
          read(io_input,*)
          ! number of flexible parts
          read(io_input,*) rindex(imol)
          if ( rindex(imol).gt.0) then
             do i = 1, rindex(imol)
                read(io_input,*) riutry(imol,i)
             end do
          else
             riutry(imol,1) = 1
          end if
       end if

       if ( lflucq(imol) .and. (.not. lelect(imol) ) ) then
          call err_exit(__FILE__,__LINE__,'lelect must be true if flucq is true',myid+1)
       end if
       if ( lqtrans(imol) ) then
          if (.not. lflucq(imol) ) then
             call err_exit(__FILE__,__LINE__,'lflucq must be true if interm.  CT is allowed',myid+1)
          end if
          write(io_output,*) 'Intermolecular Charge Transfer is allowed'
       end if

       lbias(imol) = .false.
       if (lavbmc1(imol) .or. lavbmc2(imol) .or. lavbmc3(imol)) then
          lbias(imol) = .true.
       end if

       ! choose only one from the three AVBMC algorithms
       if ( lavbmc1(imol) ) then
          lavbmc2(imol) = .false.
          lavbmc3(imol) = .false.
       else if ( lavbmc2(imol) ) then
          lavbmc3(imol) = .false.
       end if

       lneighbor = .false.
       if ( (lavbmc2(imol) .or. lavbmc3(imol)).and.(.not.lgaro) )  lneighbor = .true.

       if ( lecho.and.lprint ) then
          if (lverbose) then
             write(io_output,*) 'molecule type:',imol
             write(io_output,*) '   number of units:',nunit(imol)
             write(io_output,*) '   number of units for CBMC growth:', nugrow(imol)
             write(io_output,*) '   number of carbons for EH alkane:', ncarbon(imol)
             write(io_output,*) '   maximum number of units for CBMC:', nmaxcbmc(imol)
             write(io_output,*) '   iurot:',iurot(imol)
             write(io_output,*) '   lelect:',lelect(imol)
             write(io_output,*) '   lflucq:',lflucq(imol)
             write(io_output,*) '   lqtrans:',lqtrans(imol)
             write(io_output,*) '   lexpand:',lexpand(imol)
             write(io_output,*) '   lavbmc1:',lavbmc1(imol)
             write(io_output,*) '   lavbmc2:',lavbmc2(imol)
             write(io_output,*) '   lavbmc3:',lavbmc3(imol)
             write(io_output,*) '   fqegp:',fqegp(imol)
             write(io_output,*) '   lsetup:',lsetup
             write(io_output,*) '   lq14scale:',lq14scale(imol)
             write(io_output,*) '   qscale:',qscale(imol)
             do i = 1,nbox
                write(io_output,*) '   energy offset for box', i,':',eta2(i,imol),' K'
             end do
          else
             write(io_output,*) nunit(imol),nugrow(imol),ncarbon(imol) ,nmaxcbmc(imol),iurot(imol) ,lelect(imol),lflucq(imol)  ,lqtrans(imol),lexpand(imol),lavbmc1(imol),lavbmc2(imol) ,lavbmc3(imol),fqegp(imol) ,lsetup,(eta2(i,imol), i=1,nbox)
          end if
       end if
       masst(imol) = 0.0E0_dp

       if (lsetup) then
          call molsetup(io_input,imol,lprint)

          do i = 1,nunit(imol)
             lhere(ntype(imol,i)) = .true.
          end do

          goto 112
       end if

       do i = 1, nunit(imol)
          ! linear/branched chain with connectivity table -
          read(io_input,*)
          read(io_input,*) j, ntype(imol,i), leaderq(imol,i)
          ntype(imol,i)=indexOf(atoms,ntype(imol,i))
          if ( lelect(imol) .and. .not. lchgall ) then
             if ( lecho.and.lprint ) write(io_output,*) '   bead ',j,' beadtype ',atoms%list(ntype(imol,i)),chemid(ntype(imol,i)),' charge leader ',leaderq(imol,i)
             if ( leaderq(imol,i) .gt. j .and. .not. lchgall) then
                call err_exit(__FILE__,__LINE__,'group-based cut-off screwed for qq',myid+1)
             end if
          else
             if ( lecho.and.lprint ) write(io_output,*) '   bead ',j,' beadtype ',atoms%list(ntype(imol,i)),chemid(ntype(imol,i))
          end if

          IF (.not.(lmmff.or.lexpsix.or.lgaro.or.lninesix)) THEN
             if (ntype(imol,i).eq.0) then
                call err_exit(__FILE__,__LINE__,'ERROR: atom type undefined!',myid+1)
             else if (sigi(ntype(imol,i)).lt.1E-06_dp.and.epsi(ntype(imol,i)).lt.1E-06_dp.and.abs(qelect(ntype(imol,i))).lt.1E-06_dp) then
                call err_exit(__FILE__,__LINE__,'ERROR: atom type undefined!',myid+1)
             end if
          end if

          iutemp = ntype(imol,i)

          if (lpl(iutemp)) then
             lplace(imol,i) = .true.
             llplace(imol) = .true.
          else
             lplace(imol,i) = .false.
          end if

          masst(imol)=masst(imol)+mass(iutemp)
          lhere(iutemp) = .true.

          ! bond vibration -
          read(io_input,*)
          read(io_input,*) invib(imol,i)
          if ( invib(imol,i) .gt. 6 ) then
             write(io_output,*) 'imol',imol,'   i',i,'  invib',invib(imol,i)
             call err_exit(__FILE__,__LINE__,'too many vibrations',myid+1)
          end if
          do j = 1, invib(imol,i)
             read(io_input,*) ijvib(imol,i,j),itvib(imol,i,j)
             if((ijvib(imol,i,j).eq.i).or.(ijvib(imol,i,j).gt. nunit(imol))) then
                call err_exit(__FILE__,__LINE__,'check vibrations for mol type '//integer_to_string(imol)//' and bead '//integer_to_string(i),myid+1)
             end if

             itvib(imol,i,j)=indexOf(bonds,itvib(imol,i,j))
             if (itvib(imol,i,j).eq.0) then
                call err_exit(__FILE__,__LINE__,'ERROR: stretching parameters undefined!',myid+1)
             else if(brvib(itvib(imol,i,j)).lt.1E-06_dp) then
                call err_exit(__FILE__,__LINE__,'ERROR: stretching parameters undefined!',myid+1)
             end if

             if (lverbose.and.lprint) then
                write(io_output,*) '      bead',i,' bonded to bead', ijvib(imol,i,j)
                write(io_output,'(a20,i3,a13,f9.3,a5,f9.1)') '          bond type:', bonds%list(itvib(imol,i,j)),' bond length:', brvib(itvib(imol,i,j)),' k/2:', brvibk(itvib(imol,i,j))
             end if
          end do
          ! bond bending -
          read(io_input,*)
          read(io_input,*) inben(imol,i)
          ! write(io_output,*) inben(imol,i)
          if ( inben(imol,i) .gt. 12 ) then
             call err_exit(__FILE__,__LINE__,'too many bends',myid+1)
          end if
          do j = 1, inben(imol,i)
             read(io_input,*) ijben2(imol,i,j),ijben3(imol,i,j) ,itben(imol,i,j)
             if ((ijben2(imol,i,j).gt.nunit(imol)).or.( ijben3(imol,i,j).gt.nunit(imol))) then
                call err_exit(__FILE__,__LINE__,'check bending for molecule type '//integer_to_string(imol)//' bead '//integer_to_string(i),myid+1)
             end if
             if ((ijben2(imol,i,j).eq.i).or.( ijben3(imol,i,j).eq.i).or.(ijben2(imol,i,j) .eq.ijben3(imol,i,j))) then
                call err_exit(__FILE__,__LINE__,'check bending for molecule type '//integer_to_string(imol)//' bead '//integer_to_string(i),myid+1)
             end if

             itben(imol,i,j)=indexOf(angles,itben(imol,i,j))
             if (itben(imol,i,j).eq.0) then
                call err_exit(__FILE__,__LINE__,'ERROR: bending parameters undefined!',myid+1)
             else if(brben(itben(imol,i,j)).lt.1E-06_dp) then
                call err_exit(__FILE__,__LINE__,'ERROR: bending parameters undefined!',myid+1)
             end if

             if (lverbose.and.lprint) then
                write(io_output,'(1x,a10,i4,a20,a8,i4,a10,i4)') '      bead' ,i, ' bending interaction',' through',ijben2(imol,i,j ),' with bead',ijben3(imol,i,j)
                write(io_output,'(a20,i3,a13,f9.3,a5,f9.1)') '          bend type:',angles%list(itben(imol,i,j)) ,' bend angle :',brben(itben(imol,i,j))*180.0E0_dp/onepi ,' k/2:',brbenk(itben(imol,i,j))
             end if
          end do
          ! bond torsion -
          read(io_input,*)
          read(io_input,*) intor(imol,i)
          if ( intor(imol,i) .gt. 12 ) then
             call err_exit(__FILE__,__LINE__,'too many torsions',myid+1)
          end if
          do j = 1, intor(imol,i)
             read(io_input,*) ijtor2(imol,i,j),ijtor3(imol,i,j),ijtor4(imol,i,j),ittor(imol,i,j)

             if(ijtor2(imol,i,j).gt.nunit(imol).or.ijtor3(imol,i,j) .gt.nunit(imol).or.ijtor4(imol,i,j).gt.nunit(imol)) then
                call err_exit(__FILE__,__LINE__,'check torsion for molecule type '//integer_to_string(imol)//' bead '//integer_to_string(i),myid+1)
             end if

             if((ijtor2(imol,i,j).eq.i.or.ijtor3(imol,i,j).eq.i .or.ijtor4(imol,i,j).eq.i).or.(ijtor2(imol,i,j).eq. ijtor3(imol,i,j).or.ijtor2(imol,i,j).eq. (ijtor4(imol,i,j)).or.(ijtor3(imol,i,j).eq. ijtor4(imol,i,j)))) then
                call err_exit(__FILE__,__LINE__,'check torsion for molecule type '//integer_to_string(imol)//' bead '//integer_to_string(i),myid+1)
             end if

             if (lverbose.and.lprint) then
                write(io_output,'(1x,a10,i3,a22,a8,i3,a4,i3,a10,i3,a19,i4)') '      bead',i, ' torsional ','interaction through' ,ijtor2(imol,i,j),' and',ijtor3(imol,i,j) ,' with bead',ijtor4(imol,i,j),' of torsional type:',ittor(imol,i,j)
             end if

             ittor(imol,i,j)=indexOf(dihedrals,ittor(imol,i,j))
             if (ittor(imol,i,j).eq.0) then
                write(io_output,*) 'WARNING: torsion parameters undefined; set to null'
             end if
          end do
       end do

112    continue


       ! starting the self consistency check for the bond vibrations, bending, and torsions
       ! this would help in catching errors in fort.4 connectivity. Starting after continue
       ! so that if we use molsetup subroutine it will provide extra checking. (Neeraj)
       if (.not.lrigid(imol)) then
          do i = 1,nunit(imol)
             numvib=invib(imol,i)
             numbend=inben(imol,i)
             numtor=intor(imol,i)
             lfound = .false.
             do j = 1,numvib
                vib1 = ijvib(imol,i,j)
                vibtype  = itvib(imol,i,j)
                if(invib(imol,vib1).eq.0) then
                   write(io_output,*) 'Check vibration for mol. type:',imol, 'bead',vib1,'with',i
                   call err_exit(__FILE__,__LINE__,'ERROR IN FORT.4 VIBRATIONS',myid+1)
                end if
                do k =1,invib(imol,vib1)
                   if(ijvib(imol,vib1,k).eq.i) then
                      lfound = .true.
                      if(vibtype.ne.itvib(imol,vib1,k)) then
                         write(io_output,*) 'check vibration type of bead',i, 'with',vib1,'molecule type',imol,'vice versa'
                         call err_exit(__FILE__,__LINE__,'Error in fort.4 vibration specifications',myid+1)
                      end if
                   end if
                end do
                if(.not.lfound) then
                   write(io_output,*) 'Check vibration for mol. type:',imol, 'bead ',vib1,'with ',i
                   call err_exit(__FILE__,__LINE__,'Error in fort.4 vibration iformation',myid+1)
                end if
             end do
             lfound= .false.
             do j = 1,numbend
                bend2 = ijben2(imol,i,j)
                bend3 = ijben3(imol,i,j)
                bendtype = itben(imol,i,j)
                if(inben(imol,bend3).eq.0) then
                   write(io_output,*) 'Check bending for mol. type:',imol, 'bead ',bend3,'with ',i
                   call err_exit(__FILE__,__LINE__,'ERROR IN FORT.4 BENDING',myid+1)
                end if
                do k = 1,inben(imol,bend3)
                   if((ijben2(imol,bend3,k).eq.bend2).and. (ijben3(imol,bend3,k).eq.i)) then
                      lfound = .true.
                      if(itben(imol,bend3,k).ne.bendtype) then
                         write(io_output,*) 'check bending type of bead',i, 'with',bend3,'mol. typ.',imol,'and vice versa'
                         call err_exit(__FILE__,__LINE__,'Error in fort.4 bending specifications',myid+1)
                      end if
                   end if
                end do
                if(.not.lfound) then
                   write(io_output,*) 'Check bending for mol. type:',imol, 'bead ',bend3,'with ',i
                   call err_exit(__FILE__,__LINE__,'Error in fort.4 bending information',myid+1)
                end if
             end do
             lfound = .false.
             do j = 1,numtor
                tor2 = ijtor2(imol,i,j)
                tor3 = ijtor3(imol,i,j)
                tor4 = ijtor4(imol,i,j)
                tortype = ittor(imol,i,j)
                if(intor(imol,tor4).eq.0) then
                   write(io_output,*) 'Check torsion for mol. type:',imol, 'bead ',tor4,'with ',i,'and vice versa'
                   call err_exit(__FILE__,__LINE__,'ERROR IN FORT.4 TORSION',myid+1)
                end if
                do k = 1,intor(imol,tor4)
                   if((ijtor2(imol,tor4,k).eq.tor3).and.(ijtor3(imol,tor4 ,k).eq.tor2).and.(ijtor4(imol,tor4,k).eq.i)) then
                      lfound=.true.
                      if(ittor(imol,tor4,k).ne.tortype) then
                         write(io_output,*) 'check torsion type of bead',i, 'with',tor4,'mol. typ.',imol,'and vice versa'
                         call err_exit(__FILE__,__LINE__,'Error in fort.4 torsion specifications',myid+1)
                      end if
                   end if
                end do
                if(.not.lfound) then
                   write(io_output,*) 'Check torsion for mol. type:',imol, 'bead ',tor4,'with ',i
                   call err_exit(__FILE__,__LINE__,'Error in fort.4 torsion information',myid+1)
                end if
             end do
          end do
       end if

       ! Neeraj Adding molecule neutrality check
       !kea skip if lgaro
       if(.not.(lgaro.or.lionic.or.lexzeo)) then
          do i=1,nmolty
             qtot =0.0E0_dp
             do j = 1,nunit(i)
                qtot = qtot+qelect(ntype(i,j))
             end do
             if(abs(qtot).gt.1E-7_dp) then
                call err_exit(__FILE__,__LINE__,'molecule type '//integer_to_string(i)//' not neutral. check charges',myid+1)
             end if
          end do
       end if

       if ( lexpand(imol) ) then
          lee = .true.
       end if
       if ( lbias(imol) ) then
          read(io_input,*)
          read(io_input,*) pmbias(imol),(pmbsmt(ii),ii=1,nmolty),pmbias2(imol)
          if ( lecho .and.lprint) then
             if (lverbose) then
                write(io_output,*) '   AVBMC pmbias',pmbias(imol)
                do ii = 1,nmolty
                   write(io_output,*) '   AVBMC2 and 3 probability for', ' molecule',' type',ii,':',pmbsmt(ii)
                end do
                write(io_output,*) '   AVBMC3 pmbias2:',pmbias2(imol)
             else
                write(io_output,*) '   AVBMC bias for cluster formation and' ,' destruction'
                write(io_output,*) pmbias(imol),(pmbsmt(ii),ii=1,nmolty) ,pmbias2(imol)
             end if
          end if
          if (rbsmax .lt. rbsmin) then
             call err_exit(__FILE__,__LINE__,'rbsmax should be greater than rbsmin',myid+1)
          end if
       end if
       !kea 6/4/09 -- added for multiple rotation centers
       ! To assign multiple rotation centers, set iurot(imol) < 0
       ! Add line after molecule specification, avbmc parameters
       ! First, number of rotation centers
       ! Second, identity of centers (0=COM,integer::> 0 = bead number)
       ! Third, give probability to rotate around different centers
       if(iurot(imol).lt.0) then
          read(io_input,*)
          read(io_input,*) nrotbd(imol),irotbd(1:nrotbd(imol),imol),pmrotbd(1:nrotbd(imol),imol)
          if( lecho.and.lprint ) then
             if ( lverbose ) then
                write(io_output,*) ' Multiple rotation centers',nrotbd(imol)
                do ii=1,nrotbd(imol)
                   write(io_output,*) 'Rotation center',ii,':', irotbd(ii,imol),'    probability to select:', pmrotbd(ii,imol)
                end do
             end if
          end if
       end if
    end do

    call allocate_molecule()
    call allocate_energy_bonded()
    call init_moves_simple()
    call init_moves_volume()
    call init_cbmc()
    call init_ee()

    allocate(inclmol(ntmax*numax*numax),inclbead(ntmax*numax*numax,2),inclsign(ntmax*numax*numax),ofscale(ntmax*numax*numax),ofscale2(ntmax*numax*numax),ainclmol(ntmax*numax*numax),ainclbead(ntmax*numax*numax,2),a15t(ntmax*numax*numax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'readdat: allocating molecule failed',jerr)
    end if

    ! check whether there is a polarizable molecule
    lpolar = .false.
    lqqelect = .false.

    do imol = 1, nmolty
       if (lflucq(imol)) lpolar = .true.
       if (lelect(imol)) lqqelect = .true.
    end do

    if ( .not. lqqelect ) then
       if ( lewald .or. lchgall ) then
          if (lprint) write(io_output,*) 'No charges in the system -> lewald is now turned off'
       end if
    end if

    if ( .not. lpolar ) then
       if ( lanes  ) then
          call err_exit(__FILE__,__LINE__,'lanes should be false for nonpolarizable systems!',myid+1)
       end if
       if ( lfepsi ) then
          call err_exit(__FILE__,__LINE__,'lfepsi should be false for nonpolarizable systems!',myid+1)
       end if
    end if

    ! This B(i) goes in the acceptance rules

    if (lgrand) then
       do i=1,nmolty
          debroglie = debroglie_factor*sqrt(beta/masst(i))
          B(i) = exp(B(i)/temp)/(debroglie*debroglie*debroglie)
       end do
    end if

    ! read linkcell information
    read(io_input,*)
    read(io_input,*) licell,rintramax,boxlink

    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'licell:',licell
          write(io_output,*) 'rintramax:',rintramax,' A'
          write(io_output,*) 'boxlink:',boxlink
       else
          write(io_output,*) licell,rintramax,boxlink
       end if
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
    IF (boxlink .LE. nbxmax)THEN
       if (lsolid(boxlink).and.(.not.lrect(boxlink))) then
          call err_exit(__FILE__,__LINE__,'Linkcell not implemented for nonrectangular boxes',myid+1)
       end if
    end if

    ! read the atomic displacements
    read(io_input,*)
    read(io_input,*) Armtrax, Armtray, Armtraz

    if(lecho.and.lprint) then
       if(lverbose) then
          write(io_output,'(a41,f8.4,f8.4,f8.4)') 'initial maximum displacements for atoms:',Armtrax, Armtray , Armtraz
       else
          write(io_output,*) Armtrax, Armtray, Armtraz
       end if
    end if

    ! read displacement information
    read(io_input,*)
    read(io_input,*) rmtrax(1,1),rmtray(1,1),rmtraz(1,1)
    do im = 1,nbox
       do imol = 1,nmolty
          rmtrax(imol,im) = rmtrax(1,1)
          rmtray(imol,im) = rmtray(1,1)
          rmtraz(imol,im) = rmtraz(1,1)
       end do
    end do
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,'(a41,f8.4,f8.4,f8.4)') ' initial maximum x, y and z displacement:',rmtrax(1,1) ,rmtray(1,1),rmtraz(1,1)
       else
          write(io_output,*) rmtrax(1,1), rmtray(1,1), rmtraz(1,1)
       end if
    end if

    read(io_input,*)
    read(io_input,*) rmrotx(1,1),rmroty(1,1),rmrotz(1,1)
    do im = 1,nbox
       do imol = 1,nmolty
          rmrotx(imol,im) = rmrotx(1,1)
          rmroty(imol,im) = rmroty(1,1)
          rmrotz(imol,im) = rmrotz(1,1)
       end do
    end do
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,'(a41,f8.4,f8.4,f8.4)') ' initial maximum x, y and z rotation:    ',rmrotx(1,1) ,rmroty(1,1),rmrotz(1,1)
       else
          write(io_output,*) rmrotx(1,1), rmroty(1,1), rmrotz(1,1)
       end if
    end if
    read(io_input,*)
    read(io_input,*) tatra,tarot
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'target translational acceptance ratio:',tatra
          write(io_output,*) 'target rotational acceptance ratio:',tarot
       else
          write(io_output,*) tatra, tarot
       end if
    end if

    ! read initial setup information
    read(io_input,*)
    read(io_input,*) linit,lreadq
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'linit:',linit
          write(io_output,*) 'lreadq:',lreadq
       else
          write(io_output,*) linit, lreadq
       end if
    end if
    read(io_input,*)
    read(io_input,*) (lbranch(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          do i = 1,nmolty
             write(io_output,*) 'lbranch for molecule type',i,':',lbranch(i)
          end do
       else
          write(io_output,*) (lbranch(i),i=1,nmolty)
       end if
    end if
    do i = 1, nbox
       read(io_input,*)
       read(io_input,*) (ininch(j,i),j=1,nmolty)
       if ( lecho.and.lprint ) then
          if (lverbose) then
             write(io_output,*) 'box:',i
             do j = 1,nmolty
                write(io_output,*) '   initial number of chains of type', j,':',ininch(j,i)
             end do
          else
             write(io_output,*) 'box:',i,(ininch(j,i),j=1,nmolty)
          end if
       end if
       read(io_input,*)
       read(io_input,*) inix(i),iniy(i),iniz(i),inirot(i),inimix(i),zshift(i),dshift(i),nchoiq(i)
       if ( lecho.and.lprint ) then
          if (lverbose) then
             write(io_output,'(a36,a18,3i5)') '    initial number of chains in x, y' ,' and z directions:',inix(i),iniy(i),iniz(i)
             write(io_output,*) '   initial rotational displacement:', inirot(i)
             write(io_output,*) '   inimix:',inimix(i)
             write(io_output,*) '   zshift:',zshift(i)
             write(io_output,*) '   dshift:',dshift(i)
             write(io_output,*) '   nchoiq:',nchoiq(i)
          else
             write(io_output,*) inix(i),iniy(i),iniz(i),inirot(i), inimix(i),zshift(i),dshift(i),nchoiq(i)
          end if
       end if
    end do

    ! read ensemble specific information
    read(io_input,*)
    read(io_input,*) rmvol(1), tavol, iratv, iratp, rmflcq(1,1), taflcq
    do izz = 1,nmolty
       do zzz = 1,nbox
          rmflcq(izz,zzz) = rmflcq(1,1)
       end do
    end do
    if ( lgibbs ) then
       do zzz = 2, nbox
          rmvol(zzz) = rmvol(1)
       end do
    end if
    if ( lecho.and.lprint ) then
       if (lverbose) then
          do izz = 1,nbox
             write(io_output,*) 'initial maximum volume displacement ' ,'(rmvol) in box',izz,':',rmvol(izz)
          end do
          write(io_output,*) 'target volume acceptance ratio (tavol):',tavol
          write(io_output,*) 'iratv:',iratv
          write(io_output,*) 'iratp:',iratp
          do izz = 1,nmolty
             do zzz = 1,nbox
                write(io_output,*) 'initial maximum fluct. charge', ' displ. for chain type',izz,' box', zzz,':',rmflcq(izz,zzz)
             end do
          end do
          write(io_output,*) 'target fluctuating charge acceptance ratio', ' (taflcq):',taflcq
       else
          write(io_output,*) rmvol, tavol, iratv, iratp
          write(io_output,*) 'rmflcq', ((rmflcq(izz,zzz),zzz=1,nbox),izz=1,nmolty),taflcq
       end if
    end if

    read(io_input,*)
    read(io_input,*) pmvol,(pmvlmt(j),j=1,nbox)
    if ( .not. lfold )  call err_exit(__FILE__,__LINE__,'volume move only correct with folded coordinates',myid+1)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmvol:',pmvol
          do j = 1,nbox
             write(io_output,*) '   pmvlmt for box',j,':',pmvlmt(j)
          end do
       else
          write(io_output,*) 'pmvol',pmvol,(pmvlmt(j),j=1,nbox)
       end if
    end if
    read(io_input,*)
    read(io_input,*) nvolb,(pmvolb(j),j=1,nvolb)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'nvolb:',nvolb
          do j = 1,nvolb
             write(io_output,*) '   pmvolb:',pmvolb(j)
          end do
       else
          write(io_output,*) '   nvolb',nvolb,(pmvolb(j),j=1,nvolb)
       end if
    end if
    read(io_input,*)
    do j = 1,nvolb
       read(io_input,*) box5(j),box6(j)

       if ((lsolid(box5(j)) .and. .not. lrect(box5(j))) .and. (lsolid(box6(j)) .and. .not. lrect(box6(j)))) then
          call err_exit(__FILE__,__LINE__,'can not perform volume move between two non-rectangular boxes',myid+1)
       end if

       if ( lecho.and.lprint ) then
          if (lverbose) then
             write(io_output,*) '   box pair for volume move number',j,':', box5(j),box6(j)
          else
             write(io_output,*) box5(j),box6(j)
          end if
       end if
    end do

    lxyz = .false.
    do j = 1,nbox
       if (lsolid(j) .and. .not. lxyz) then
          lxyz = .true.
          read(io_input,*)
          read(io_input,*) pmvolx,pmvoly
          if (lecho.and.lprint) then
             if (lverbose) then
                write(io_output,*) 'pmvolx:',pmvolx
                write(io_output,*) 'pmvoly:',pmvoly
             else
                write(io_output,*) pmvolx,pmvoly
             end if
          end if
       end if
    end do

    ! read swatch information
    read(io_input,*)
    read(io_input,*) pmswat,nswaty
    if ( nswaty .gt. npamax ) then
       call err_exit(__FILE__,__LINE__,'nswaty gt npamax',myid+1)
    end if

    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmswat:',pmswat
          write(io_output,*) '   number of swatch pairs (nswaty):',nswaty
       else
          write(io_output,*) 'pmswat, nswaty',pmswat,nswaty
       end if
    end if

    if (nswaty .gt. 0) then
       read(io_input,*)
       read(io_input,*) ((nswatb(i,j),j=1,2),i=1,nswaty)
       read(io_input,*)
       read(io_input,*) (pmsatc(i),i=1,nswaty)

       ! safety checks on swatch
       do i = 1,nswaty
          if ( nswatb(i,1) .eq. nswatb(i,2) ) then
             write(io_output,*) 'nswaty ',i,' has identical moltyp'
             call err_exit(__FILE__,__LINE__,'cannot swatch identical moltyp',myid+1)
          end if
       end do

       if( lecho.and.lprint ) then
          if (lverbose) then
             do i = 1,nswaty
                write(io_output,*) '   swatch molecule type pairs:', (nswatb(i,j),j=1,2)
             end do
             do i = 1,nswaty
                write(io_output,*) '   probability of each swatch pair:', pmsatc(i)
             end do
          else
             do i=1,nswaty
                write(io_output,*) 'nswatb', (nswatb(i,j),j=1,2)
                write(io_output,*) 'pmsatc',pmsatc(i)
             end do
          end if
       end if

       do i=1,nswaty
          ! number of beads that remain in the same position
          read(io_input,*)
          read(io_input,*) nsampos(i),(ncut(i,j),j=1,2)
          if (lecho.and.lprint) then
             if (lverbose) then
                write(io_output,*) '   nsampos:',nsampos(i)
                write(io_output,*) '   ncut:',(ncut(i,j),j=1,2)
             else
                write(io_output,*) 'nsampos',nsampos(i),' ncut', (ncut(i,j),j=1,2)
             end if
          end if
          ! bead number
          read(io_input,*)
          do j = 1,nsampos(i)
             read(io_input,*) (splist(i,j,k),k=1,2)
             if (lecho.and.lprint) then
                if (lverbose) then
                   write(io_output,*) '   splist:', (splist(i,j,k),k=1,2)
                else
                   write(io_output,*) 'splist',(splist(i,j,k),k=1,2)
                end if
             end if
          end do

          read(io_input,*)
          read(io_input,*) (( gswatc(i,j,k), k=1,2*ncut(i,j) ), j=1,2 )
          ! if (lecho.and.lprint) then
          !    if (lverbose) then
          !       do izz = 1,ncut(i,j)
          !          write(io_output,*) '   grow from and prev for ncut',izz,':',(gswatc(i,j,izz),j=1,2)
          !       end do
          !    else
          !       write(io_output,*) 'gswatc',(( gswatc(i,j,k),k=1,2*ncut(i,j) ), j=1,2 )
          !    end if
          ! end if

          read(io_input,*)
          read(io_input,*) nswtcb(i), (pmswtcb(i,ipair), ipair=1,nswtcb(i))
          read(io_input,*)
          do ipair = 1,nswtcb(i)
             read(io_input,*) box3(i,ipair),box4(i,ipair)
             if (lecho.and.lprint) then
                if (lverbose) then
                   write(io_output,*) '   box pair:', box3(i,ipair),box4(i,ipair)
                else
                   write(io_output,*) box3(i,ipair),box4(i,ipair)
                end if
             end if
          end do
       end do
    else
       ! skip past all of the swatch info
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
       read(io_input,*)
    end if

    ! read swap info
    read(io_input,*)
    read(io_input,*) pmswap, (pmswmt(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmswap:',pmswap
          do i = 1,nmolty
             write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   swap probability for molecule type ',' ',i ,' (pmswmt):',pmswmt(i)
          end do
       else
          write(io_output,*) 'pmswap',pmswap,(pmswmt(i),i=1,nmolty)
       end if
    end if
    do i = 1, nmolty
       read(io_input,*)
       read(io_input,*) nswapb(i), (pmswapb(i,ipair),ipair=1,nswapb(i))
       if ( lecho.and.lprint ) then
          if (lverbose) then
             write(io_output,*) '   number of swap moves for molecule type', i,':',nswapb(i)
             do ipair = 1,nswapb(i)
                write(io_output,*) '      pmswapb:',pmswapb(i,ipair)
             end do
          else
             write(io_output,*) nswapb(i),  (pmswapb(i,ipair),ipair=1,nswapb(i))
          end if
       end if
       read(io_input,*)
       do ipair = 1, nswapb(i)
          read(io_input,*) box1(i,ipair), box2(i,ipair)
          if ( lecho.and.lprint) then
             if (lverbose) then
                write(io_output,*) '      box pair:', box1(i,ipair), box2(i,ipair)
             else
                write(io_output,*) 'box pair:', box1(i,ipair), box2(i,ipair)
             end if
          end if
       end do
    end do

    ! read cbmc info
    read(io_input,*)
    read(io_input,*) pmcb, (pmcbmt(i),i=1,nmolty)
    read(io_input,*)
    read(io_input,*) (pmall(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmcb:',pmcb
          do i = 1,nmolty
             write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   CBMC probability for molecule type ',' ',i ,' (pmcbmt):',pmcbmt(i)
          end do
          do i = 1,nmolty
             write(io_output,*) '   pmall for molecule type',i,':',pmall(i)
          end do
       else
          write(io_output,*) 'pmcb',pmcb,(pmcbmt(i),i=1,nmolty), 'pmall',(pmall(i),i=1,nmolty)
       end if
    end if
    read(io_input,*)
    read(io_input,*) (pmfix(i),i=1,nmolty)
    if (lecho.and.lprint) then
       if (lverbose) then
          do i = 1,nmolty
             write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   probability for SAFE-CBMC for ','molecule typ e',i ,'  (pmfix):',pmfix(i)
          end do
       else
          write(io_output,*) (pmfix(i),i=1,nmolty)
       end if
    end if
    do i = 1, nmolty
       if (lring(i).and.pmfix(i).lt.0.999.and. .not.lrig(i)) then
          call err_exit(__FILE__,__LINE__,'a ring can only be used with safe-cbmc',myid+1)
       end if
    end do
    ! read fluctuating charge info
    read(io_input,*)
    read(io_input,*) pmflcq, (pmfqmt(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmflcq:',pmflcq
          do i = 1,nmolty
             write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   flcq probability for molecule type ',' ',i ,' (pmfqmt):',pmfqmt(i)
          end do
       else
          write(io_output,*) 'pmflcq',pmflcq,(pmfqmt(i),i=1,nmolty)
       end if
    end if
    ! read expanded-coefficient move info
    read(io_input,*)
    read(io_input,*) pmexpc, (pmeemt(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmexpc:',pmexpc
          do i = 1,nmolty
             write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   expanded ens. prob. for molecule ','type ',i ,' (pmeemt):',pmeemt(i)
          end do
       else
          write(io_output,*) 'pmexpc',pmexpc,(pmeemt(i),i=1,nmolty)
       end if
    end if
    ! read new expanded ensemble info
    read(io_input,*)
    read(io_input,*) pmexpc1
    ! read atom translation probability
    read(io_input,*)
    read(io_input,*) pm_atom_tra
    ! read translation info
    read(io_input,*)
    read(io_input,*) pmtra,(pmtrmt(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmtra:',pmtra
          do i = 1,nmolty
             write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   translation probability for molecule','   typ e',i ,' (pmtrmt):',pmtrmt(i)
          end do
       else
          write(io_output,*) 'pmtra',pmtra,(pmtrmt(i),i=1,nmolty)
       end if
    end if

    ! read rotation info
    read(io_input,*)
    read(io_input,*) (pmromt(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'pmrot:',1.0E0_dp
          do i = 1,nmolty
             write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   rotational probability for molecule','    typ e',i ,' (pmromt):',pmromt(i)
          end do
       else
          write(io_output,*) (pmromt(i),i=1,nmolty)
       end if
    end if

    ! writeout probabilities, in percent
    if (lverbose.and.lprint) then
       write(io_output,*)
       write(io_output,*) 'percentage move probabilities:'
       write(io_output,'(1x,a19,f8.2,a2)') 'volume move       :', 100.0E0_dp*pmvol,' %'
       pcumu = pmvol
       if (pmswat .gt. pmvol) then
          pm = pmswat - pcumu
          pcumu = pcumu + pm
       else
          pm = 0.0E0_dp
       end if
       write(io_output,'(1x,a19,f8.2,a2)') 'swatch move       :', 100.0E0_dp*pm,' %'
       if (pmswap .gt. pmswat) then
          pm = pmswap - pcumu
          pcumu = pcumu + pm
       else
          pm = 0.0E0_dp
       end if
       write(io_output,'(1x,a19,f8.2,a2)') 'swap move         :', 100.0E0_dp*pm,' %'
       if (pmcb .gt. pmswap) then
          pm = pmcb - pcumu
          pcumu = pcumu + pm
       else
          pm = 0.0E0_dp
       end if
       write(io_output,'(1x,a19,f8.2,a2)') 'CBMC move         :', 100.0E0_dp*pm,' %'
       if (pmflcq .gt. pmcb) then
          pm = pmflcq - pcumu
          pcumu = pcumu + pm
       else
          pm = 0.0E0_dp
       end if
       write(io_output,'(1x,a19,f8.2,a2)') 'fluct charge move :', 100.0E0_dp*pm,' %'
       if (pmexpc .gt. pmflcq) then
          pm = pmexpc - pcumu
          pcumu = pcumu + pm
       else
          pm = 0.0E0_dp
       end if
       write(io_output,'(1x,a19,f8.2,a2)') 'expanded ens move :', 100.0E0_dp*pm,' %'
       if (pmtra .gt. pmexpc) then
          pm = pmtra - pcumu
          pcumu = pcumu + pm
       else
          pm = 0.0E0_dp
       end if
       write(io_output,'(1x,a19,f8.2,a2)') 'translation move  :', 100.0E0_dp*pm,' %'
       pm = 1.0E0_dp - pmtra
       write(io_output,'(1x,a19,f8.2,a2)') 'rotation move     :', 100.0E0_dp*pm,' %'
       write(io_output,*)
       write(io_output,*) 'Fraction of atom translations move', pm_atom_tra
    end if

    ! read growth details
    read(io_input,*)
    read(io_input,*) (nchoi1(i),i=1,nmolty),(nchoi(i),i=1,nmolty) ,(nchoir(i),i=1,nmolty) ,(nchoih(i),i=1,nmolty),(nchtor(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*)
          write(io_output,*) 'molecule type :',(i,'  ',i=1,nmolty)
          write(io_output,*) '     nchoi1   :',(nchoi1(i),' ',i=1,nmolty)
          write(io_output,*) '     nchoi    :',(nchoi(i),' ',i=1,nmolty)
          write(io_output,*) '     nchoir   :',(nchoir(i),' ',i=1,nmolty)
          write(io_output,*) '     nchoih   :',(nchoih(i),' ',i=1,nmolty)
          write(io_output,*) '     nchtor   :',(nchtor(i),i=1,nmolty)
       else
          write(io_output,*) 'nchoi1',(nchoi1(i),i=1,nmolty)
          write(io_output,*) 'nchoi',(nchoi(i),i=1,nmolty)
          write(io_output,*) 'nchoir',(nchoir(i),i=1,nmolty)
          write(io_output,*) 'nchoih',(nchoih(i),i=1,nmolty)
          write(io_output,*) 'nchtor',(nchtor(i),i=1,nmolty)
       end if
    end if

    nchmax=0
    nchtor_max=0
    do imol=1,nmolty
       if ( nchoi1(imol) .gt. nchmax ) then
          nchmax=nchoi1(imol)
       end if
       if (nchoi(imol) .gt. nchmax ) then
          nchmax=nchoi(imol)
       end if
       if (nchoir(imol) .gt. nchmax ) then
          nchmax=nchoir(imol)
       end if
       if ( nchoih(imol) .ne. 1 .and. nunit(imol) .eq. nugrow(imol) )  then
          call err_exit(__FILE__,__LINE__,'nchoih must be 1 (one) if nunit = nugrow',myid+1)
       end if
       if (nchtor(imol) .gt. nchtor_max ) then
          nchtor_max=nchtor(imol)
       end if
    end do

    call read_safecbmc()

    ! read information for CBMC bond angle growth
    read(io_input,*)
    read(io_input,*) (nchbna(i),i=1,nmolty),(nchbnb(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          do i = 1,nmolty
             write(io_output,*) 'nchbna and nchbnb for molecule type',i,':', nchbna(i),nchbnb(i)
          end do
       else
          write(io_output,*) 'nchbna ',(nchbna(i),i=1,nmolty)
          write(io_output,*) 'nchbnb ',(nchbnb(i),i=1,nmolty)
       end if
    end if

    read(io_input,*)
    read(io_input,*) (icbdir(i),i=1,nmolty),(icbsta(i),i=1,nmolty)
    if ( lecho.and.lprint ) then
       if (lverbose) then
          do i = 1,nmolty
             write(io_output,*) 'icbdir for molecule type',i,':',icbdir(i)
          end do
          do i = 1,nmolty
             write(io_output,*) 'icbsta for molecule type',i,':',icbsta(i)
          end do
       else
          write(io_output,*) 'icbdir',(icbdir(i),i=1,nmolty)
          write(io_output,*) 'icbsta',(icbsta(i),i=1,nmolty)
       end if
    end if

    nchbn_max=0
    do i=1,nmolty
       if ( nchbna(i) .gt. nchbn_max ) then
          nchbn_max=nchbna(i)
       end if
       if ( nchbnb(i) .gt. nchbn_max ) then
          nchbn_max=nchbnb(i)
       end if
       if ( abs(icbsta(i)) .gt. nunit(i) ) then
          call err_exit(__FILE__,__LINE__,'icbsta > nunit for molecule '//integer_to_string(i),myid+1)
       end if
    end do

    call allocate_cbmc()

    ! read exclusion table for intermolecular interactions
    read(io_input,*)
    read(io_input,*) nexclu
    do i = 1, nmolty
       do ii = 1, numax
          do j = 1, nmolty
             do jj = 1, numax
                lexclu(i,ii,j,jj) = .false.
             end do
          end do
       end do
    end do
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*)'nexclu:',nexclu
       else
          write(io_output,*) nexclu
       end if
    end if
    if (nexclu .ne. 0) then
       do ndum = 1, nexclu
          read(io_input,*) i, ii, j, jj
          lexclu(i,ii,j,jj) = .true.
          lexclu(j,jj,i,ii) = .true.
          if (lverbose.and.lprint) then
             write(io_output,*) 'excluding interactions between bead',ii, ' on chain',i,' and bead',jj,' on chain',j
          end if
       end do
    else
       read(io_input,*)
    end if

    ! read inclusion list
    read(io_input,*)
    read(io_input,*) inclnum
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'inclnum:',inclnum
       else
          write(io_output,*) inclnum
       end if
    end if
    if (inclnum .ne. 0) then
       do ndum = 1, inclnum
          read(io_input,*) inclmol(ndum),inclbead(ndum,1),inclbead(ndum,2),inclsign(ndum),ofscale(ndum),ofscale2(ndum)
          if (lecho.and.lprint ) then
             if (lverbose) then
                if (inclsign(ndum) .eq. 1) then
                   write(io_output,*) 'including intramolecular ' ,'interactions for chain type',inclmol(ndum), ' between beads',inclbead(ndum,1),' and', inclbead(ndum,2)
                else
                   write(io_output,*) 'excluding intramolecular ' ,'interactions for chain type',inclmol(ndum), ' between beads',inclbead(ndum,1),' and', inclbead(ndum,2)

                end if
             else
                write(io_output,*) inclmol(ndum),inclbead(ndum,1) ,inclbead(ndum,2),inclsign(ndum),ofscale(ndum),ofscale2(ndum)
             end if
          end if
       end do
    else
       read(io_input,*)
    end if

    ! read a15 inclusion list
    read(io_input,*)
    read(io_input,*) ainclnum
    if ( lecho.and.lprint ) then
       if (lverbose) then
          write(io_output,*) 'ainclnum:',ainclnum
       else
          write(io_output,*) ainclnum
       end if
    end if
    if (ainclnum .ne. 0) then
       do ndum = 1, ainclnum
          read(io_input,*) ainclmol(ndum),ainclbead(ndum,1),ainclbead(ndum,2) ,a15t(ndum)
          if (lecho.and.lprint) then
             if (lverbose) then
                write(io_output,*) 'repulsive 1-5 OH interaction for', ' chain type',ainclmol(ndum),' between beads', ainclbead(ndum,1),' and',ainclbead(ndum,2), ' of type:',a15t(ndum)
             else
                write(io_output,*) ainclmol(ndum),ainclbead(ndum,1) ,ainclbead(ndum,2),a15t(ndum)
             end if
          end if
       end do
    else
       read(io_input,*)
    end if

    ! set up the inclusion table
    call inclus(inclnum,inclmol,inclbead,inclsign,ncarbon,ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2,lprint)

    ! read in information on the chemical potential checker
    read(io_input,*)
    read(io_input,*) lucall
    read(io_input,*)
    read(io_input,*) (ucheck(jj),jj=1,nmolty)
    if (lecho.and.lprint) then
       if (lverbose) then
          write(io_output,*) 'lucall:',lucall
          do jj = 1,nmolty
             write(io_output,*) '   ucheck for molecule type', jj,':',ucheck(jj)
          end do
       else
          write(io_output,*) lucall,(ucheck(jj),jj=1,nmolty)
       end if
    end if

    ! read information for virial coefficient calculation
    read(io_input,*)
    read(io_input,*) nvirial,starvir,stepvir
    if (lecho.and.lprint) then
       if (lverbose) then
          write(io_output,*) 'nvirial:',nvirial
          write(io_output,*) 'starvir:',starvir
          write(io_output,*) 'stepvir:',stepvir
       else
          write(io_output,*) nvirial,starvir,stepvir
       end if
    end if

    if (lvirial) then
       if ( nvirial .gt. maxvir ) then
          call err_exit(__FILE__,__LINE__,'nvirial .gt. maxvir',myid+1)
       end if

       read(io_input,*)
       read(io_input,*) ntemp,(virtemp(jj),jj=1,ntemp)
       if (lecho.and.lprint) then
          if (lverbose) then
             write(io_output,*) 'ntemp:',ntemp
             write(io_output,*) 'calculation of virial coefficient ', 'at the following temperatures:', (virtemp(jj),jj=1,ntemp)
          else
             write(io_output,*) ntemp, 'Calculation of virial coefficient ', 'at the following temperatures:', (virtemp(jj),jj=1,ntemp)
          end if
       end if
    end if

    ! JLR 11-11-09
    ! reading in extra variables for RPLC simulations
    read(io_input,*)
    read(io_input,*) (lideal(i),i=1,nbox)
    if (lecho.and.lprint)  then
       write(io_output,*) 'lideal: ', (lideal(i),i=1,nbox)
    end if
    do i = 1,nbox
       if (lideal(i) .and. lexpee) then
          call err_exit(__FILE__,__LINE__,'Cannot have lideal and lexpee both true (if you want this you will have change code)',myid+1)
       end if
    end do
    read(io_input,*)
    read(io_input,*) (ltwice(i),i=1,nbox)
    read(io_input,*)
    read(io_input,*) (lrplc(i),i=1,nmolty)
    if (lecho.and.lprint) then
       write(io_output,*) 'ltwice: ', (ltwice(i),i=1,nbox)
       write(io_output,*) 'lrplc: ', (lrplc(i),i=1,nmolty)
    end if
    ! END JLR 11-11-09

    ! JLR 11-11-09
    ! skip this part if analysis will not be done (if ianalyze .gt. nstep)
    ! KM 01/10 remove analysis
    ! if (ianalyze.lt.nstep) then
    !    nhere = 0
    !    do izz=1,nntype
    !       if ( lhere(izz) ) then
    !          nhere = nhere + 1
    !          temphe(nhere) = izz
    !          beadtyp(nhere)=izz
    !       end if
    !    end do
    !    do izz = 1,nhere
    !       atemp = temphe(izz)
    !       decode(atemp) = izz
    !    end do
    ! end if
    ! END JLR 11-11-09
    close(io_input)
! -------------------------------------------------------------------
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
    !    call err_exit(__FILE__,__LINE__,'lchgall is true and lewald is false. Not checked for accuracy!',myid+1)
    ! end if

    ! check information of .INC-files
    if (lprint) then
       write(io_output,*)
       write(io_output,*) '***** program   =  THE MAGIC BLACK BOX'
       iensem = 0
       if ( lgrand ) then
          iensem = iensem + 1
          write(io_output,*) 'grand-canonical ensemble'
       end if
       if ( lvirial ) then
          write(io_output,*) 'Computing Second Virial Coefficient'
          if ( nchain .ne. 2) then
             call err_exit(__FILE__,__LINE__,'nchain must equal 2',myid+1)
          end if
       end if
       if ( lgibbs ) then
          iensem = iensem + 1
          if ( lnpt ) then
             write(io_output,*) 'NPT Gibbs ensemble'
          else
             write(io_output,*) 'NVT Gibbs ensemble'
          end if
       end if
       if ( iensem .eq. 0 ) then
          if ( lnpt ) then
             write(io_output,*) 'Isobaric-isothermal ensemble'
             iensem = iensem + 1
          else
             write(io_output,*) 'Canonical ensemble'
             iensem = iensem + 1
          end if
       end if
       if ( iensem .gt. 1 ) then
          call err_exit(__FILE__,__LINE__,'INCONSISTENT ENSEMBLE SPECIFICATION',myid+1)
       end if

       inpbc = 0
       if ( lpbc ) then
          write(io_output,*) 'using periodic boundaries'
          if ( lpbcx ) then
             inpbc = inpbc + 1
             write(io_output,*) 'in x-direction'
          end if
          if ( lpbcy ) then
             inpbc = inpbc + 1
             write(io_output,*) 'in y-direction'
          end if
          if ( lpbcz ) then
             inpbc = inpbc + 1
             write(io_output,*) 'in z-direction'
          end if
          if ( inpbc .eq. 0 ) then
             call err_exit(__FILE__,__LINE__,'INCONSISTENT PBC SPECIFICATION',myid+1)
          end if
          write(io_output,*) inpbc,'-dimensional periodic box'
       else
          write(io_output,*) 'cluster mode (no pbc)'
          if ( lgibbs .or. lgrand ) then
             call err_exit(__FILE__,__LINE__,'INCONSISTENT SPECIFICATION OF LPBC  AND ENSEMBLE',myid+1)
          end if
       end if

       if ( lfold ) then
          write(io_output,*) 'particle coordinates are folded into central box'
          if ( .not. lpbc ) then
             call err_exit(__FILE__,__LINE__,'INCONSISTENT SPECIFICATION OF LPBC AND LFOLD',myid+1)
          end if
       end if

       if ( lijall ) then
          write(io_output,*) 'all i-j-interactions are considered', '  (no potential cut-off)'
       end if

       if ( lchgall ) then
          write(io_output,*) 'all the inter- and intramolecular Coulombic', ' interactions are considered (no group-based cutoff)'
       end if
       if ( lcutcm ) then
          write(io_output,*) 'additional (COM) cutoff on computed rcmu'
          if ( lijall ) then
             call err_exit(__FILE__,__LINE__,'cannot have lijall with lcutcm',myid+1)
          end if
          ! if ( lchgall ) call err_exit(__FILE__,__LINE__,'cannot have lchgall with lcutcm',myid+1)
       end if
       if ( ldual ) then
          write(io_output,*) 'Dual Cutoff Configurational-bias Monte Carlo'
       end if

       write(io_output,*)  'CBMC simultaneously grows all beads conected to the same bead'
       write(io_output,*) 'with bond angles generated from Gaussian', ' distribution'

       if ( ljoe ) then
          write(io_output,*) 'external 12-3 potential for SAM (Joe Hautman)'
       end if
       if (lslit) then
          write(io_output,*) '10-4-3 graphite slit pore potential'
       end if
       if ( lsami ) then
          write(io_output,*) 'external potential for Langmuir films (Sami)'
          write(io_output,*) 'WARNING: LJ potential defined in SUSAMI'
          write(io_output,*) 'WARNING: sets potential cut-off to 2.5sigma'
          write(io_output,*) 'WARNING: has build-in tail corrections'
          if ( ltailc ) then
             call err_exit(__FILE__,__LINE__,'INCONSISTENT SPECIFICATION OF LTAILC AND LSAMI',myid+1)
          end if
       end if

       if ( lexzeo ) then
          write(io_output,*) 'external potential for zeolites'
       end if

       write(io_output,*) 'Program will call Explct.f if needed'

       if ( lexpsix ) then
          write(io_output,*) 'Exponential-6 potential for bead-bead'
       else if ( lmmff ) then
          write(io_output,*) 'Buffered 14-7 potential for bead-bead'
       else if (lninesix) then
          write(io_output,*) '9-6 potential for bead-bead'
       else if (lgenlj) then
          write(io_output,*) 'Generalized Lennard Jones potential for bead-bead'
       else if (lgaro) then
          write(io_output,*) 'Feuston-Garofalini potential for bead-bead'
       else
          write(io_output,*) 'Lennard-Jones potential for bead-bead'
       end if

       if ( ltailc ) then
          write(io_output,*) 'with additional tail corrections'
       end if

       if ( lshift) then
          write(io_output,*) 'using a shifted potential'
       end if

       if (lshift.and.ltailc) then
          call err_exit(__FILE__,__LINE__,'lshift.and.ltailc!',myid+1)
          return
       end if

       write(io_output,*) 'Coulombic inter- and intramolecular interactions'
       if ( lewald ) write(io_output,*) 'Ewald-sum will be used to calculate Coulombic interactions'
       write(io_output,*)
       write(io_output,*) 'MOLECULAR MASS:', (masst(i),i=1,nmolty)
    end if

! -------------------------------------------------------------------

    if (lgrand .and. .not.(lslit)) then
       ! volume ideal gas box is set arbitry large!
       boxlx(2)=1000*boxlx(1)
       boxly(2)=1000*boxly(1)
       boxlz(2)=1000*boxlz(1)
    end if
    if (linit) then
       do ibox = 1,nbox
          if (lsolid(ibox).and..not.lrect(ibox) .and..not.(ibox.eq.1.and.lexzeo)) then
             call err_exit(__FILE__,__LINE__,'Cannot initialize non-rectangular system',myid+1)
          end if
       end do
    end if

    if ( linit ) then
       call initia(file_struct)
       nnstep = 0
    else
       if (use_checkpoint) file_restart='save-config'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM needed for long input records, like
       ! nboxi and moltyp for a lot of molecules
       io_restart=get_iounit()
       open(unit=io_restart,access='sequential',action='read',file=file_restart,form='formatted',iostat=jerr,recl=4096,status='old')
       if (jerr.ne.0) then
          call err_exit(__FILE__,__LINE__,'cannot open main restart file',myid+1)
       end if
       read(io_restart,*) nnstep
       read(io_restart,*) Armtrax, Armtray, Armtraz
       do im = 1,nbox
          do imol = 1,nmolty
             read(io_restart,*) rmtrax(imol,im), rmtray(imol,im), rmtraz(imol,im)
             read(io_restart,*) rmrotx(imol,im), rmroty(imol,im), rmrotz(imol,im)
          end do
          call averageMaximumDisplacement(im,imol)
       end do

       if (lprint) then
          write(io_output,*)  'new maximum displacements read from restart-file'
          do im = 1,nbox
             write(io_output,*)'box      #',im
             do imol = 1,nmolty
                write(io_output,*) 'molecule type',imol
                write(io_output,"(' max trans. displacement:        ',3f10.6)") rmtrax(imol,im), rmtray(imol,im), rmtraz(imol,im)
                write(io_output,"(' max rot. displacement:          ',3f10.6)") rmrotx(imol,im), rmroty(imol,im), rmrotz(imol,im)
             end do
          end do
       end if

       do im=1,nbox
          read(io_restart,*) (rmflcq(i,im),i=1,nmolty)
       end do
       if ( lecho.and.lprint) then
          do im=1,nbox
             write(io_output,*) 'maximum fluc q displacements: Box #',im
             write(io_output,*) (rmflcq(i,im),i=1,nmolty)
          end do
       end if
       ! changed so fort.77 the same for all ensembles
       ! 06/08/09 KM
       read(io_restart,*) (rmvol(ibox), ibox = 1,nbox)
       if (lprint) then
          write(io_output,"(' max volume displacement:        ',3E12.4)") (rmvol(ibox), ibox = 1,nbox)
          write(io_output,*)
          write(io_output,*)
       end if
       do ibox = 1,nbox
          if (lsolid(ibox) .and. .not. lrect(ibox)) then
             read(io_restart,*) (rmhmat(ibox,j),j=1,9)
             if (lprint) write(io_output,*) (rmhmat(ibox,j),j=1,9)
          end if
       end do

       if (lprint) then
          write(io_output,*)
          write(io_output,*) 'new box size read from restart-file'
          write(io_output,*)
       end if

       do ibox = 1,nbox
          if (ibox.eq.1.and.lexzeo) then
             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                read(io_restart,*) dum,dum,dum,dum,dum,dum,dum,dum,dum
             else
                read(io_restart,*) dum,dum,dum
             end if
          else if (lsolid(ibox) .and. .not. lrect(ibox)) then

             read(io_restart,*) (hmat(ibox,j),j=1,9)
             if (lprint) then
                write(io_output,*)
                write(io_output,*) 'HMAT COORDINATES FOR BOX',ibox
                write(io_output,*) (hmat(ibox,j),j=1,3)
                write(io_output,*) (hmat(ibox,j),j=4,6)
                write(io_output,*) (hmat(ibox,j),j=7,9)
                write(io_output,*)
             end if

             call matops(ibox)

             if (lprint) then
                write(io_output,*)
                write(io_output,"(' width of  box ',i1,':',3f12.6)") ibox ,min_width(ibox,1),min_width(ibox,2) ,min_width(ibox,3)
             end if

             w(1) = min_width(ibox,1)
             w(2) = min_width(ibox,2)
             w(3) = min_width(ibox,3)

             if (rcut(ibox)/w(1) .gt. 0.5E0_dp .or.  rcut(ibox)/w(2) .gt. 0.5E0_dp .or.  rcut(ibox)/w(3) .gt. 0.5E0_dp) then
                call err_exit(__FILE__,__LINE__,'rcut > half cell width',myid+1)
             end if

             if (lprint) then
                write(io_output,*)
                write(io_output,*) 'ibox:  ', ibox
                write(io_output,'("cell length |a|:",2x,f12.3)') cell_length(ibox,1)
                write(io_output,'("cell length |b|:",2x,f12.3)')  cell_length(ibox,2)
                write(io_output,'("cell length |c|:",2x,f12.3)')  cell_length(ibox,3)

                write(io_output,*)
                write(io_output,'("cell angle alpha:",2x,f12.3)')  cell_ang(ibox,1)*180.0E0_dp/onepi
                write(io_output,'("cell angle beta: ",2x,f12.3)') cell_ang(ibox,2)*180.0E0_dp/onepi
                write(io_output,'("cell angle gamma:",2x,f12.3)') cell_ang(ibox,3)*180.0E0_dp/onepi

                ! write(io_output,"(' angle of  box ',i1,'  :  ','   alpha:   ',f12.6,'
                !     &  beta: ', f12.6, '    gamma:   ',f12.6)") ibox,cell_ang(ibox,1)*180.0E0_dp/onepi,
                !     &                cell_ang(ibox,2)*180.0E0_dp/onepi,cell_ang(ibox,3)
                !     &                *180/onepi
             end if
          else
             read(io_restart,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
             if (lprint) then
                write(io_output,*)
                write(io_output,*)
                write(io_output,"(' dimension box ',i1,'  :','  a:   ',f12.6,'   b:   ',f12.6,'   c   :  ' ,f12.6)") ibox,boxlx(ibox),boxly(ibox),boxlz(ibox)
             end if
             do i = 1, nbox
                if( (rcut(i)/boxlx(i) .gt. 0.5E0_dp).or. (rcut(i)/boxly(i) .gt. 0.5E0_dp).or. (rcut(i)/boxlz(i) .gt. 0.5E0_dp)) then
                   call err_exit(__FILE__,__LINE__,'rcut > 0.5*boxlx',myid+1)
                end if
             end do
          end if
       end do

       if (lprint) then
          write(io_output,*)
          write(io_output,*) 'Finished writing simulation box related info'
       end if

       read(io_restart,*) ncres
       read(io_restart,*) nmtres
       ! check that number of particles in fort.4 & fort.77 agree ---
       if ( ncres .ne. nchain .or. nmtres .ne. nmolty ) then
          write(io_output,*) 'nchain',nchain,'ncres',ncres
          write(io_output,*) 'nmolty',nmolty,'nmtres',nmtres
          call err_exit(__FILE__,__LINE__,'conflicting information in restart and control files',myid+1)
       end if
       read(io_restart,*) (nures(i),i=1,nmtres)

       ! do i = 1, nmtres
       !    if ( nures(i) .ne. nunit(i) ) then
       !       write(io_output,*) 'unit',i,'nunit',nunit(i),'nures',nures(i)
       !       call err_exit(__FILE__,__LINE__,'conflicting information in restart and control files',myid+1)
       !    end if
       ! end do
       ! write(io_output,*) 'ncres',ncres,'   nmtres',nmtres

       read(io_restart,*) (moltyp(i),i=1,ncres)
       read(io_restart,*) (nboxi(i),i=1,ncres)
       if ( lee ) then
          do i = 1, nmtres
             if ( lexpand(i) ) read(io_restart,*) eetype(i)
          end do
          do i = 1, nmtres
             if ( lexpand(i) ) read(io_restart,*) rmexpc(i)
          end do
       end if

       ! write(io_output,*) 'start reading coordinates'
       ! check that particles are in correct boxes ---
       ! obtain ncmt values
       do ibox = 1,nbox
          nchbox(ibox) = 0
       end do
       do i = 1, nmolty
          do ibox = 1,nbox
             ncmt(ibox,i) = 0
          end do
          if ( lexpand(i) ) then
             !!??? problem in expand ensemble
             do j = 1, numcoeff(i)
                do ibox = 1,2
                   ncmt2(ibox,i,j) = 0
                end do
             end do
          end if
       end do

       do i = 1, ncres
          if( nboxi(i) .le. nbox ) then
             ibox = nboxi(i)
             if ( ibox .ne. 1 .and. .not. lgibbs  .and. .not.lgrand) then
                call err_exit(__FILE__,__LINE__,'Particle found outside BOX 1',myid+1)
             end if
             nchbox(ibox) = nchbox(ibox) + 1
             imolty = moltyp(i)
             ncmt(ibox,imolty) = ncmt(ibox,imolty) + 1
             if ( lexpand(imolty) ) then
                if ( ibox .gt. 2 ) then
                   call err_exit(__FILE__,__LINE__,'put in box 1 and 2 for such  molecules',myid+1)
                end if
                itype = eetype(imolty)
                ncmt2(ibox,imolty,itype) =  ncmt2(ibox,imolty,itype) + 1
                do j = 1,nunit(imolty)
                   sigma_f(imolty,j) = sigm(imolty,j,itype)
                   epsilon_f(imolty,j) = epsil(imolty,j,itype)
                end do
             end if
          else
             write(io_output,*) 'i:',i,'nboxi(i)',nboxi(i)
             call err_exit(__FILE__,__LINE__,'Particle found in ill-defined box',myid+1)
          end if

       end do

       ! check that number of particles of each type is consistent
       do i = 1, nmolty
          tcount = 0
          do ibox = 1, nbox
             tcount = tcount + ncmt(ibox,i)
          end do
          if ( tcount .ne. temtyp(i) ) then
             write(io_output,*) 'type',i
             write(io_output,*) 'ncmt',(ncmt(ibox,i), ibox = 1,nbox)
             write(io_output,*) 'temtyp', temtyp(i)
             call err_exit(__FILE__,__LINE__,'Particle type number inconsistency',myid+1)
          end if
       end do

       ! write(io_output,*) 'particles found in correct box with correct type'

       do i = 1,nbxmax
          qbox(i) = 0.0E0_dp
       end do

       do i = 1, nchain
          ! write(io_output,*) 'reading coord of chain i',i
          imolty = moltyp(i)
          do j = 1, nunit(imolty)
             read(io_restart,*) rxu(i,j), ryu(i,j), rzu(i,j),qqu(i,j)
             if (.not. lreadq) then
                qqu(i,j) = qelect(ntype(imolty,j))
             end if
             qbox(nboxi(i)) = qbox(nboxi(i)) + qqu(i,j)
          end do
       end do

       close(io_restart)

       do i = 1,nbxmax
          if ( abs(qbox(i)) .gt. 1E-6_dp ) then
             if (i.eq.1.and.lexzeo) cycle
             call err_exit(__FILE__,__LINE__,'box '//integer_to_string(i)//' has a net charge of '//real_to_string(qbox(i)),myid+1)
          end if
       end do
    end if

    if (L_Ewald_Auto) then
       if ( lchgall ) then
          if (lewald) then
             do ibox = 1,nbox
                if (lsolid(ibox).and.(.not.lrect(ibox))) then
                   min_boxl = min(min_width(ibox,1),min_width(ibox,2), min_width(ibox,3))
                else
                   min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
                end if
                kalp(ibox) = 6.40E0_dp
                calp(ibox) = kalp(ibox)/min_boxl
             end do
          else
             call err_exit(__FILE__,__LINE__,'lewald should be true when lchgall is true',myid+1)
          end if
       else
          do ibox = 1, nbox
             if (lsolid(ibox).and.(.not.lrect(ibox))) then
                min_boxl = min(min_width(ibox,1),min_width(ibox,2), min_width(ibox,3))
             else
                min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
             end if
             ! rcut(ibox) = 0.4E0_dp*min_boxl
             calp(ibox) = 3.2E0_dp/rcut(ibox)
          end do
       end if
    else
       if ( lchgall ) then
          ! if real space term are summed over all possible pairs in the box
          ! kalp(1) & kalp(2) are fixed while calp(1) & calp(2) change
          ! according to the boxlength, kalp(1) should have a value greater than
          ! 5.0
          do ibox = 1, nbox
             if (lsolid(ibox).and.(.not.lrect(ibox))) then
                min_boxl = min(min_width(ibox,1),min_width(ibox,2), min_width(ibox,3))
             else
                min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
             end if
             calp(ibox) = kalp(ibox)/min_boxl
             if ( lewald ) then
                if ( kalp(ibox) .lt. 5.6E0_dp ) then
                   call err_exit(__FILE__,__LINE__,'Warning, kalp is too small',myid+1)
                end if
             else
                !kea
                if(.not. lgaro) then
                   call err_exit(__FILE__,__LINE__,'lewald should be true when lchgall is true',myid+1)
                end if
             end if
          end do
       else
          ! if not lchgall, calp(1) & calp(2) are fixed
          do ibox = 1, nbox
             calp(ibox) = kalp(ibox)
             if ( lewald ) then
                if (calp(ibox)*rcut(ibox).lt.3.2E0_dp.and.lprint) then
                   write(io_output,*) 'Warning, kalp too small in box',ibox
                   write(io_output,*) ibox,calp(ibox),rcut(ibox)
                   ! JLR 11-24-09
                   ! you may want a smaller kalp, e.g. when comparing to previous work
                   ! This does not need to be an error
                   ! call err_exit(__FILE__,__LINE__,'kalp is too small, set to 3.2/rcutchg',myid+1)
                   write(io_output,*) 'kalp is too small, set to 3.2/rcutchg'
                   ! END JLR 11-24-09
                end if
             end if
          end do
       end if
    end if

    if (L_Ewald_Auto) then
       do ibox = 1,nbox
          if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
             k_max_l(ibox) = aint(boxlx(ibox)*calp(ibox))+1
             k_max_m(ibox) = aint(boxly(ibox)*calp(ibox))+1
             k_max_n(ibox) = aint(boxlz(ibox)*calp(ibox))+1
          else
             k_max_l(ibox) = aint(hmat(ibox,1)*calp(ibox))+2
             k_max_m(ibox) = aint(hmat(ibox,5)*calp(ibox))+2
             k_max_n(ibox) = aint(hmat(ibox,9)*calp(ibox))+2
          end if
       end do
    end if

    if (L_Ewald_Auto.and.lprint) then
       write(io_output,*)
       write(io_output,*) '****Ewald Parameters*****'
       write(io_output,*) 'ibox   calp(ibox)  kmaxl(ibox)   kmaxm(ibox)', '   kmaxn(ibox)   rcut(ibox)'
       do ibox = 1,nbox
          write(io_output,'(i4,5x,f12.6,3i12,12x,f12.4)') ibox, calp(ibox),  k_max_l(ibox),k_max_m(ibox), k_max_n(ibox),rcut(ibox)
       end do
       write(io_output,*)
    end if

    ! book keeping arrays
    do ibox=1,nbox
       do imolty=1,nmolty
          idummy(imolty)=0
       end do
       do i = 1, nchain
          if ( nboxi(i) .eq. ibox ) then
             imolty=moltyp(i)
             idummy(imolty) = idummy(imolty)+1
             parbox(idummy(imolty),ibox,imolty)=i
             ! pparbox(i,ibox) = nmcount
          end if
       end do
       nmcount = 0
       do imolty=1,nmolty
          nmcount = nmcount + idummy(imolty)
       end do
       if ( nmcount .ne. nchbox(ibox) ) then
          write(io_output,*) 'nmcount: ',nmcount,'; nchbox: ',nchbox
          call err_exit(__FILE__,__LINE__,'Readdat: nmcount ne nchbox',myid+1)
       end if
    end do
    ! set idummy counter to 0
    do imolty=1,nmolty
       idummy(imolty)=0
    end do
    ! set up parall
    do i =1,nchain
       imolty = moltyp(i)
       idummy(imolty) = idummy(imolty) + 1
       parall(imolty,idummy(imolty)) = i
    end do

! -------------------------------------------------------------------
    ! read/produce initial/starting configuration
    ! zeolite external potential
    if ( lexzeo ) call suzeo(lprint)
    call readThreeBody(file_in)
    call readFourBody(file_in)
! ----------------------------------------------------------------

    ! reordering of numbers for charmm
    ! do i = 1, nchain
    !    imolty = moltyp(i)
    !    temtyp(i) = imolty
    !    do j = 1, nunit(imolty)
    !       temx(i,j) = rxu(i,j)
    !       temy(i,j) = ryu(i,j)
    !       temz(i,j) = rzu(i,j)
    !    end do
    ! end do
    ! innew = 0
    ! do it = 1, nmolty
    !    do i = 1, nchain
    !       imolty = temtyp(i)
    !       if ( imolty .eq. it ) then
    !          innew = innew + 1
    !          moltyp(innew) = it
    !          do j = 1, nunit(imolty)
    !             rxu(innew,j) = temx(i,j)
    !             ryu(innew,j) = temy(i,j)
    !             rzu(innew,j) = temz(i,j)
    !          end do
    !       end if
    !    end do
    !    write(io_output,*) 'it =',it,'   innew =',innew
    ! end do

! -------------------------------------------------------------------

    ! set the centers of mass if LFOLD = .TRUE.

    if ( lfold ) then
       do ibox = 1, nbox
          call ctrmas(.true.,ibox,0,6)
       end do
    end if

    ! check that rintramax is really valid
    if (licell) then
       do i = 1,nchain
          if (2.0E0_dp*rcmu(i) .gt. rintramax) then
             call err_exit(__FILE__,__LINE__,'rintramax for the linkcell list too small',myid+1)
          end if
       end do
    end if

    if (myid.eq.rootid) then
       ! write out movie-header
       !*** The very first frame (nnstep=0) is no longer written out because it's usually useless
       if (imv.le.nstep) then
          io_movie=get_iounit()
          open(unit=io_movie,access='stream',action='write',file=file_movie,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open movie file '//trim(file_movie),jerr)
          if (lprint) then
             nhere = 0
             do izz=1,nntype
                if ( lhere(izz) ) then
                   nhere = nhere + 1
                   temphe(nhere) = izz
                end if
             end do
             write(io_movie,*) nstep/imv,nchain,nmolty,nbox,nhere
             write(io_movie,*) (rcut(ibox),ibox=1,nbox)
             write(io_movie,*) (temphe(izz),izz=1,nhere)

             do imolty = 1,nmolty
                write(io_movie,*) nunit(imolty)
                ! output bond connectivity information
                do ii=1,nunit(imolty)
                   write(io_movie,*) invib(imolty,ii),(ijvib(imolty,ii,z),z=1,invib(imolty,ii))
                end do

                ! output torsional connectivity information
                do j = 1,nunit(imolty)
                   write(io_movie,*) intor(imolty,j),(ijtor2(imolty,j,ii),ijtor3(imolty,j,ii),ijtor4(imolty,j,ii),ii=1,intor(imolty,j))
                end do
             end do
          end if
       else
          io_movie=-1
       end if

       ! open xyz movie file for individual simulation box
       if (L_movie_xyz) then
          do ibox = 1,nbox
             write(file_box_movie,'("box",I1.1,"movie",I1.1,A,".xyz")') ibox ,run_num,suffix
             io_box_movie(ibox)=get_iounit()
             open(unit=io_box_movie(ibox),access='stream',action='write',file=file_box_movie,form='formatted',iostat=jerr,status='unknown')
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open box movie file '//trim(file_box_movie),jerr)
          end do
       else
          io_box_movie=-1
       end if

       ! write out isolute movie header
       if (ANY(isolute(1:nmolty).le.nstep)) then
          io_solute=get_iounit()
          open(unit=io_solute,access='stream',action='write',file=file_solute,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open solute movie file '//trim(file_solute),jerr)
          if (lprint) then
             write(io_solute,*) nmolty
             do imol = 1, nmolty
                write(io_solute,*) imol,nunit(imol),(nstep/isolute(imol))*temtyp(imol)
             end do
          end if
       else
          io_solute=-1
       end if

       ! cell parameters using crystallographic convention
       if (ANY(lsolid.and..not.lrect)) then
          write(file_cell,'("cell_param",I1.1,A,".dat")') run_num,suffix
          io_cell=get_iounit()
          open(unit=io_cell,access='stream',action='write',file=file_cell,form='formatted',iostat=jerr,status='unknown')
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open cell file '//trim(file_cell),jerr)
       else
          io_cell=-1
       end if

       ! set up info at beginning of fort.12 for analysis
       io_traj=get_iounit()
       open(unit=io_traj,access='stream',action='write',file=file_traj,form='formatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open traj file '//trim(file_traj),jerr)
       if (lprint) write(io_traj,*) nstep,nmolty,(masst(i),i=1,nmolty)
    end if
! -------------------------------------------------------------------

    ! write out connectivity and bonded interactions
    if (lverbose.and.lprint) then
       do imol = 1,nmolty
          write(io_output,*) 'molecule type',imol
          if (nunit(imol) .gt. 1) then
             write(io_output,*) '   i   j   type_i type_j   bond length', '        k/2'
          end if
          do i = 1,nunit(imol)
             do j = 1, invib(imol,i)
                write(io_output,'(i5,i4,i7,i7,f13.4,f14.1)') i,ijvib(imol,i ,j),atoms%list(ntype(imol,i)),atoms%list(ntype(imol,ijvib(imol,i,j))),brvib(itvib(imol,i,j)),brvibk(itvib(imol,i,j))
             end do
          end do

          if (nunit(imol) .gt. 2) then
             write(io_output,*)
             write(io_output,*) '   i   j   k   type_i type_j type_k', '     angle      k/2'
          end if
          do i = 1,nunit(imol)
             do j = 1,inben(imol,i)
                write(io_output,'(i5,i4,i4,i7,i7,i7,f12.2,f12.1)') i ,ijben2(imol,i,j),ijben3(imol,i,j),atoms%list(ntype(imol,i))&
                 ,atoms%list(ntype(imol,ijben2(imol,i,j))),atoms%list(ntype(imol,ijben3(imol,i,j)))&
                 ,brben(itben(imol,i,j))*180.0E0_dp/onepi,brbenk(itben(imol,i,j))
             end do
          end do

          if (nunit(imol) .gt. 3) then
             write(io_output,*)
             write(io_output,*) '   i   j   k   l    type_i type_j type_k', ' type_l     torsion type'
          end if

          do i = 1,nunit(imol)
             do j = 1, intor(imol,i)
                write(io_output,'(i5,i4,i4,i4,i8,i7,i7,i7,i14)') i,ijtor2(imol,i,j),ijtor3(imol,i,j),ijtor4(imol,i,j)&
                 ,atoms%list(ntype(imol,i)),atoms%list(ntype(imol,ijtor2(imol,i,j))),atoms%list(ntype(imol,ijtor3(imol,i,j)))&
                 ,atoms%list(ntype(imol,ijtor4(imol,i,j))),dihedrals%list(ittor(imol,i,j))
             end do
          end do

       end do
    end if

    ! write out non-bonded interaction table
    if (lprint) then
       write(io_output,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
       ! added a flag for lgaro...need to add for all potentials, probably
       if ((.not.lexpsix).and.(.not.lmmff).and.(.not.lgaro)) then
          if (lninesix) then
             write(io_output,*) '     i   j    r_0,ij     epsij         q0(i)          q0(j)'
          else
             write(io_output,*) '     i   j    sig2ij     epsij         q0(i)          q0(j)'
          end if
          do i = 1, nntype
             do j = 1,nntype
                if (lhere(i).and.lhere(j)) then
                   ij=type_2body(i,j)
                   if (lninesix) then
                      write(io_output,'(3x,2i4,2f10.5,2f15.6)') atoms%list(i),atoms%list(j) ,rzero(ij),epsnx(ij),qelect(i),qelect(j)
                   else
                      write(io_output,'(3x,2i4,2f10.5,2f15.6)') atoms%list(i),atoms%list(j) ,sqrt(sig2ij(ij)),epsij(ij),qelect(i),qelect(j)
                   end if
                end if
             end do
          end do
       end if
    end if

    ! write input data to unit 6 for control ***
    if (lprint) then
       write(io_output,*)
       write(io_output,*) 'number of mc cycles:            ', nstep
       write(io_output,*) 'number of chains:               ', nchain
       write(io_output,*)
       write(io_output,*) 'temperature:                    ', temp
       write(io_output,*)
       write(io_output,*) 'ex-pressure:                    ', (express(ibox),ibox=1,nbox)
       write(io_output,*)
    end if

    deallocate(ncarbon,idummy,qbox,nures,k_max_l,k_max_m,k_max_n,lhere,temphe,inclmol,inclbead,inclsign,ofscale,ofscale2,ainclmol,ainclbead,a15t,stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'readdat: deallocation failed',jerr)
    end if

    return
  end subroutine readdat

!> \brief Perform periodic updates and block averages
  subroutine monper()
    use const_math,only:raddeg
    use const_phys,only:N_Avogadro,R_gas,MPa2SimUnits
    use util_math,only:update_average,store_block_average
    use util_string,only:format_n
    use energy_intramolecular,only:U_bonded
    use energy_pairwise,only:energy,coru
    use moves_simple,only:update_translation_rotation_max_displacement
    use moves_volume,only:update_volume_max_displacement
    use moves_cbmc,only:opt_safecbmc
    use transfer_shared,only:opt_bias,lopt_bias,freq_opt_bias
    use transfer_swap,only:acchem,bnchem
    use prop_pressure,only:pressure

    logical::lfq,ovrlap
    integer::im,i,ibox,jbox,Temp_nmol,m,mm,imolty,nummol,ii,ntii,intg,ilunit,k,j,itype,itel,zzz,jmolty
    real::ratflcq,press1,surf,temvol,Temp_Mol_Vol,Temp_Energy,Heat_vapor_T,Heat_vapor_LJ,Heat_vapor_COUL,CED_T,CED_LJ,CED_COUL,HSP_T,HSP_LJ,HSP_COUL,pdV,setx,sety,setz,setel,v(nEnergy),vol,rho,temmass,dpr,dpp

#ifdef __DEBUG__
    write(io_output,*) 'begin MONPER in ',myid
#endif
! -------------------------------------------------------------------
    ! Optimize and output MC move parameters
    if (ANY(lopt_bias).and.mod(nnn,freq_opt_bias).eq.0) then
       call opt_bias()
    end if

    if (mod(nnn,iratio).eq.0) then
       call update_translation_rotation_max_displacement(io_output)

       ! adjust maximum charge displacement for fluc Q
       lfq = .false.
       do im = 1,nbox
          do i = 1,nmolty
             if ( bnflcq(i,im) .gt. 0.5E0_dp ) then
                lfq = .true.
                ratflcq = bsflcq(i,im)/(bnflcq(i,im)*taflcq)
                if ( ratflcq .lt. 0.1E0_dp ) then
                   rmflcq(i,im) = rmflcq(i,im) * 0.1E0_dp
                else
                   rmflcq(i,im) = rmflcq(i,im) * ratflcq
                end if
             end if
             ! accumulate flcq info for final output
             bsflcq2(i,im) = bsflcq2(i,im) + bsflcq(i,im)
             bnflcq2(i,im) = bnflcq2(i,im) + bnflcq(i,im)
             ! rezero flcq
             bsflcq(i,im) = 0.0E0_dp
             bnflcq(i,im) = 0.0E0_dp
          end do
       end do
       if ( lfq.and.myid.eq.rootid ) then
          ! write out information about fluctuating charge success
          write(io_output,*) 'Box:   rmflcq for moltyps'
          do im =1,nbox
             write(io_output,*) im,(rmflcq(i,im),i=1,nmolty)
          end do
       end if
    end if

    if (mod(nnn,iupdatefix).eq.0) then
       call opt_safecbmc()
    end if

    if ((lgibbs.or.lnpt).and.(mod(nnn,iratv).eq.0).and.(pmvol.gt.0.0E0_dp)) then
       call update_volume_max_displacement(io_output)
    end if
! -------------------------------------------------------------------
    ! Calculate pressure and other related properties
    if (mod(nnn,iratp).eq.0) then
       ! calculate pressure ***
       acnp = acnp + 1
       do ibox = 1, nbox
          call pressure( press1, surf, ibox )
          pres(ibox) = press1
          call update_average(acpres(ibox),press1,acnp)
          call update_average(acsurf(ibox),surf,acnp)

          if (lsolid(ibox) .and. .not. lrect(ibox)) then
             temvol = cell_vol(ibox)
          else
             if ( lpbcz ) then
                temvol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
             else
                temvol = boxlx(ibox)*boxly(ibox)
             end if
          end if

          ! Enthalpy calculation
          Temp_nmol = sum(ncmt(ibox,1:nmolty))
          ! molar volume in m3/mol, energies in kJ/mol
          Temp_Mol_Vol = temvol/Temp_nmol*N_Avogadro*1E-30_dp
          Temp_Energy  = vbox(1,ibox)/Temp_nmol*R_gas/1000_dp
          call update_average(acEnthalpy(ibox),Temp_Energy+pres(ibox)*Temp_Mol_Vol,acnp)
          call update_average(acEnthalpy1(ibox),Temp_Energy+(express(ibox)*1E3_dp/MPa2SimUnits)*Temp_Mol_Vol,acnp)
       end do

       ! cannot calculate a heat of vaporization for only one box,
       ! and some compilers choke because Heat_vapor_T will not be
       ! defined if nbox == 1
       if (lgibbs) then
          do ibox = 1,nbox-1
             do jbox = ibox+1,nbox
                ! WRITE(io_output,*) 'ieouwfe ',ibox,jbox
                call calcsolpar(pres,Heat_vapor_T,Heat_vapor_LJ,Heat_vapor_COUL,pdV,CED_T,CED_LJ,CED_COUL,HSP_T,HSP_LJ,HSP_COUL,ibox,jbox)

                ! Heat of vaporization
                call update_average(acsolpar(1,ibox,jbox),Heat_vapor_T,acnp)
                call update_average(acsolpar(2,ibox,jbox),Heat_vapor_LJ,acnp)
                call update_average(acsolpar(3,ibox,jbox),Heat_vapor_COUL,acnp)
                call update_average(acsolpar(4,ibox,jbox),CED_T,acnp)
                call update_average(acsolpar(5,ibox,jbox),CED_LJ,acnp)
                call update_average(acsolpar(6,ibox,jbox),CED_COUL,acnp)
                call update_average(acsolpar(7,ibox,jbox),HSP_T,acnp)
                call update_average(acsolpar(8,ibox,jbox),HSP_LJ,acnp)
                call update_average(acsolpar(9,ibox,jbox),HSP_COUL,acnp)
                !call update_average(acsolpar(10,ibox,jbox),DeltaU_Ext,acnp)
                call update_average(acsolpar(11,ibox,jbox),pdV,acnp)
             end do
          end do
       end if
    end if
! -------------------------------------------------------------------
    ! Print out summary current simulation status
    if ((mod(nnn,iprint).eq.0).and.(myid.eq.rootid)) then
       ! write out runtime information ***
       write(io_output,FMT='(i6,i8,e12.4,f10.3,f12.1,'//format_n(nmolty,"i5")//')') nnn,tmcc,vbox(1,1),boxlx(1),pres(1),(ncmt(1,imolty),imolty=1,nmolty)
       if ( lgibbs ) then
          do ibox = 2, nbox
             write(io_output,FMT='(14x,e12.4,f10.3,f12.1,'//format_n(nmolty,"i5")//')') vbox(1,ibox),boxlx(ibox),pres(ibox) ,(ncmt(ibox,imolty),imolty=1,nmolty)
          end do
       end if
    end if

    if ((lgibbs.or.lnpt).and.(.not.lvirial).and.(myid.eq.rootid)) then
       do ibox = 1,nbox
          if ( lpbcz ) then
             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                write(io_cell,'(i8,6f12.4)') tmcc,cell_length(ibox,1)/Num_cell_a,cell_length(ibox,2)/Num_cell_b,cell_length(ibox,3)/Num_cell_c,cell_ang(ibox,1)*raddeg,cell_ang(ibox,2)*raddeg,cell_ang(ibox,3)*raddeg
                write(io_traj,FMT='(7E13.5,'//format_n(nmolty,"i5")//')') hmat(ibox,1),hmat(ibox,4),hmat(ibox,5),hmat(ibox,7),hmat(ibox,8),hmat(ibox,9),vbox(1,ibox),(ncmt(ibox,itype),itype=1,nmolty)
             else
                write(io_traj,'(4E13.5,'//format_n(nmolty,"i5")//')') boxlx(ibox),boxly(ibox),boxlz(ibox),vbox(1,ibox),(ncmt(ibox,itype),itype=1,nmolty)
             end if
          else
             write(io_traj,'(2E12.5,'//format_n(nmolty,"i4")//')') boxlx(ibox)*boxly(ibox),vbox(1,ibox),(ncmt(ibox,itype),itype=1,nmolty)
          end if
       end do
    end if
! -------------------------------------------------------------------
    ! Write out movie frames
    if (mod(nnn,imv).eq.0) then
       if ( lvirial ) then
          call virial(binvir,binvir2)
       else if (myid.eq.rootid) then
          ! write out the movie configurations ***
          write(io_movie,*) nnn
          do ibox = 1, nbox
             write(io_movie,*) (ncmt(ibox,zzz),zzz=1,nmolty)

             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                write(io_movie,*) (hmat(ibox,zzz),zzz=1,9)
             else
                write(io_movie,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
             end if
          end do
          do m = 1, nchain
             imolty = moltyp(m)
             write(io_movie,'(4(1x,i5),3(1x,f16.6))') m,imolty,nunit(imolty),nboxi(m),xcm(m),ycm(m),zcm(m)
             do mm = 1, nunit(imolty)
                write(io_movie,'(4(1x,f14.6),i5)') rxu(m,mm),ryu(m,mm),rzu(m,mm),qqu(m,mm),ntype(imolty,mm)
             end do
          end do

          ! KM for MPI
          if (L_movie_xyz) then
             do ibox = 1,nbox
                nummol = 0
                do i = 1,nchain
                   if (nboxi(i).eq.ibox) then
                      nummol = nummol + nunit(moltyp(i))
                   end if
                end do
                write(io_box_movie(ibox),*) nummol
                write(io_box_movie(ibox),*)
                do i = 1,nchain
                   if(nboxi(i).eq.ibox) then
                      imolty = moltyp(i)
                      do ii = 1,nunit(imolty)
                         ntii = ntype(imolty,ii)
                         write(io_box_movie(ibox),'(a4,5x,3f15.4)') chemid(ntii), rxu(i,ii), ryu(i,ii), rzu(i,ii)
                      end do
                   end if
                end do
             end do
          end if
       end if
    end if
! -------------------------------------------------------------------
    if (lucall) then
       call err_exit(__FILE__,__LINE__,'not recently checked for accuracy',myid+1)
       ! do j = 1,nmolty
       !    if ( ucheck(j) .gt. 0 ) then
       !       call chempt(bsswap,j,ucheck(j))
       !    end if
       ! end do
    end if

    ! calculate the integrand of thermosynamic integration
    if (lmipsw.and.(mod(nnn,iratipsw).eq.0)) then
       acipsw = acipsw+1
       call deriv(1)
       call update_average(acdvdl,dvdl,acipsw)
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
       mnbox(ibox,imolty) = mnbox(ibox,imolty) + 1
       call update_average(asetel(ibox,imolty),setel,mnbox(ibox,imolty))
    end do

    do imolty = 1, nmolty
       ! KM for MPI
       ! only do anything for lsolute if monper is not called from readdat (nnn.ne.0)
       ! calculate energy and write out movie for lsolute
       if ((nnn.ne.0).and.mod(nnn,isolute(imolty)).eq.0) then
          do ibox = 1, nbox
             do k = 1, ncmt(ibox,imolty)
                i = parbox(k,ibox,imolty)
                ! set coords for energy and write out conformations
                if (myid.eq.rootid) write(io_solute,*) imolty,ibox,nunit(imolty)
                do j = 1, nunit(imolty)
                   rxuion(j,1) = rxu(i,j)
                   ryuion(j,1) = ryu(i,j)
                   rzuion(j,1) = rzu(i,j)
                   if (myid.eq.rootid) write(io_solute,*) ntype(imolty,j),rxuion(j,1),ryuion(j,1),rzuion(j,1),qqu(j,1)
                end do

                call energy(i,imolty,v,1,ibox,1,nunit(imolty),.true.,ovrlap,.false.,.false.,.false.,.false.)
                if (ovrlap) write(io_output,*)  '*** DISASTER, OVERLAP IN MONPER'

                if (ltailc) then
                   ! tail corrections
                   if (lsolid(ibox).and..not.lrect(ibox)) then
                      vol = cell_vol(ibox)
                   else
                      vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
                   end if
                   v(3) = 0.0E0_dp
                   do jmolty = 1, nmolty
                      rho = ncmt(ibox,jmolty) / vol
                      v(3) = v(3) + ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
                   end do
                end if

                call U_bonded(i,imolty,v(5),v(6),v(7))

                solcount(ibox,imolty) = solcount(ibox,imolty) + 1
                call update_average(avsolinter(ibox,imolty),v(2)/2.0_dp+v(3),solcount(ibox,imolty))
                call update_average(avsolintra(ibox,imolty),v(4),solcount(ibox,imolty))
                call update_average(avsoltor(ibox,imolty),v(7),solcount(ibox,imolty))
                call update_average(avsolbend(ibox,imolty),v(6),solcount(ibox,imolty))
                call update_average(avsolelc(ibox,imolty),v(8)+v(14),solcount(ibox,imolty))
             end do
          end do
       end if
    end do
! -------------------------------------------------------------------
    ! calculation of block averages
    if (mod(nnn,iblock).eq.0)  then
       nblock = nblock + 1
       do ibox=1,nbox
          ! specific density
          temmass = 0.0E0_dp
          if ( lpbcz ) then
             do itype = 1, nmolty
                temmass = temmass + masst(itype)*acdens(ibox,itype)
             end do
             dpr = temmass*1E24_dp/N_Avogadro
          else
             do itype = 1, nmolty
                temmass = temmass + acdens(ibox,itype)
             end do
             dpr = temmass
          end if
          call store_block_average(baver(1,ibox,nblock),dpr,acmove,bccold(1,ibox),nccold(1,ibox))

          ! pressure
          call store_block_average(baver(2,ibox,nblock),acpres(ibox),acnp,bccold(2,ibox),nccold(2,ibox))

          ! surface tension
          itel = nEnergy + 4*nmolty + 3
          call store_block_average(baver(itel,ibox,nblock),acsurf(ibox),acnp,bccold(itel,ibox),nccold(itel,ibox))

          ! box volume
          itel = nEnergy + 4*nmolty + 4
          call store_block_average(baver(itel,ibox,nblock),acvolume(ibox),acmove,bccold(itel,ibox),nccold(itel,ibox))

          ! energies
          do j=1,nEnergy
             itel = 2 + j
             call store_block_average(baver(itel,ibox,nblock),acv(j,ibox),acmove,bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! chemical potential
          do itype = 1, nmolty
             itel = 2 + nEnergy + itype
             dpp = acchem(ibox,itype)/bnchem(ibox,itype)
             call store_block_average(baver(itel,ibox,nblock),dpp,bnchem(ibox,itype),bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! square end-to-end length
          do itype = 1, nmolty
             itel = 2 + nEnergy + nmolty + itype
             call store_block_average(baver(itel,ibox,nblock),asetel(ibox,itype),mnbox(ibox,itype),bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! number density
          do itype = 1, nmolty
             itel = 2 + nEnergy + 2*nmolty + itype
             call store_block_average(baver(itel,ibox,nblock),acdens(ibox,itype),acmove,bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! mol fraction
          do itype = 1, nmolty
             itel = 2 + nEnergy + 3*nmolty + itype
             call store_block_average(baver(itel,ibox,nblock),molfra(ibox,itype),acmove,bccold(itel,ibox),nccold(itel,ibox))
          end do

          ! Enthalpy
          itel = nEnergy + 4*nmolty + 5
          call store_block_average(baver(itel,ibox,nblock),acEnthalpy(ibox),acnp,bccold(itel,ibox),nccold(itel,ibox))
          itel = nEnergy + 4*nmolty + 6
          call store_block_average(baver(itel,ibox,nblock),acEnthalpy1(ibox),acnp,bccold(itel,ibox),nccold(itel,ibox))
       end do
       do ibox = 1,nbox-1
          do jbox = ibox+1,nbox
             do i = 1,nprop1
                call store_block_average(baver1(i,ibox,jbox,nblock),acsolpar(i,ibox,jbox),acnp,bccold1(i,ibox,jbox),nccold1(i,ibox,jbox))
             end do
          end do
       end do
    end if

    if (ldielect.and.(mod(nnn,idiele).eq.0).and.myid.eq.rootid) then
       ! use fort.27 to calculate dielectric constant
       !---------not currently used---------!
       ! If you really want this quantity comment should be taken out**
       ! dielect = 6.9994685465110493E5_dp*beta/(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
       ! write(14,*) nnn,dielect*acdipolesq(4,ibox)/acmove
       ! write(15,*) nnn,dielect*(acdipolesq(4,ibox)/acmove-(acdipole(1,ibox)/acmove)**2-(acdipole(2,ibox)/acmove)**2 - (acdipole(3,ibox)/acmove)**2)
       ! write(16,*) nnn,acdipole(1,ibox)/acmove,acdipole(2,ibox)/acmove,acdipole(3,ibox)/acmove
       ! write(25,*) nnn+nnstep,vbox(1,1)
       !---------not currently used---------!

       ! output the fluctuation information
       write(14,*) nnn,acvol(ibox)
       write(15,*) nnn,acvolsq(ibox)
       write(16,*) nnn,acv(1,ibox)
       write(17,*) nnn,acvsq(1,ibox)
       write(18,*) nnn,acvol(ibox)*acv(1,ibox)
       write(19,*) nnn,acv(14,ibox)-acvol(ibox)*acv(1,ibox)
       do ibox = 1,nbox
          write(27,*) dipolex(ibox),dipoley(ibox),dipolez(ibox)
       end do
    end if

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

    ! Output running averages of residual heat capacity, in (J2/mol)
    if (mod(nnn,iheatcapacity).eq.0.and..not.lgibbs.and.myid.eq.rootid) then
       write(56,'(I12,F18.6,F18.2,F18.6)') nnn,(enthalpy2-enthalpy*enthalpy)/real(nchain,dp)*R_gas/(temp**2),enthalpy2,enthalpy
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end MONPER in ',myid
#endif
    return
  end subroutine monper

  subroutine read_checkpoint_main(file_chkpt)
    use util_mp,only:mp_bcast
    use moves_simple,only:read_checkpoint_simple
    use moves_cbmc,only:read_checkpoint_cbmc
    use moves_volume,only:read_checkpoint_volume
    use transfer_swap,only:read_checkpoint_swap
    use transfer_swatch,only:read_checkpoint_swatch
    use transfer_shared,only:read_checkpoint_transfer_shared
    use moves_ee,only:read_checkpoint_ee
    character(LEN=*),intent(in)::file_chkpt
    integer::io_chkpt,jerr,i,pos_output,pos_movie,pos_box_movie(nbox),pos_solute,pos_cell,pos_traj

    if (myid.eq.rootid) then
       io_chkpt=get_iounit()
       open(unit=io_chkpt,access='stream',action='read',file=file_chkpt,form='unformatted',iostat=jerr,status='unknown')
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot read checkpoint file '//trim(file_chkpt),jerr)

       read(io_chkpt) pos_output,pos_movie,pos_box_movie,pos_solute,pos_cell,pos_traj,nnstep,nnn,acmove,vbox,vstart,acdens,molfra,acnbox,acnbox2,acv,acvsq,acvkjmol,acdipole,acdipolesq,acboxa,acboxl,acvol,acvolsq,acvolume,enthalpy,enthalpy2,acnp,acpres,acsurf,acEnthalpy,acEnthalpy1,solcount,acsolpar,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc,mnbox,asetel,acipsw,acdvdl,nblock,nccold1,bccold1,baver1,nccold,bccold,baver
    end if

    call read_checkpoint_simple(io_chkpt)
    call read_checkpoint_cbmc(io_chkpt)
    call read_checkpoint_volume(io_chkpt)
    call read_checkpoint_swap(io_chkpt)
    call read_checkpoint_swatch(io_chkpt)
    call read_checkpoint_transfer_shared(io_chkpt)

    ! fluctuating charges
    if (myid.eq.rootid) read(io_chkpt) bnflcq,bsflcq,bnflcq2,bsflcq2

    call read_checkpoint_ee(io_chkpt)

    call mp_bcast(nnstep,1,rootid,groupid)
    call mp_bcast(nnn,1,rootid,groupid)
    call mp_bcast(acmove,1,rootid,groupid)
    call mp_bcast(vbox,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(vstart,nbxmax,rootid,groupid)
    call mp_bcast(acdens,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(molfra,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acnbox,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acnbox2,nbxmax*ntmax*20,rootid,groupid)
    call mp_bcast(acv,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(acvsq,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(acvkjmol,nEnergy*nbxmax,rootid,groupid)
    call mp_bcast(acdipole,4*nbxmax,rootid,groupid)
    call mp_bcast(acdipolesq,4*nbxmax,rootid,groupid)
    call mp_bcast(acboxa,nbxmax*3,rootid,groupid)
    call mp_bcast(acboxl,nbxmax*3,rootid,groupid)
    call mp_bcast(acvol,nbxmax,rootid,groupid)
    call mp_bcast(acvolsq,nbxmax,rootid,groupid)
    call mp_bcast(acvolume,nbxmax,rootid,groupid)
    call mp_bcast(enthalpy,1,rootid,groupid)
    call mp_bcast(enthalpy2,1,rootid,groupid)
    call mp_bcast(acnp,1,rootid,groupid)
    call mp_bcast(acpres,nbxmax,rootid,groupid)
    call mp_bcast(acsurf,nbxmax,rootid,groupid)
    call mp_bcast(acEnthalpy,nbxmax,rootid,groupid)
    call mp_bcast(acEnthalpy1,nbxmax,rootid,groupid)
    call mp_bcast(solcount,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acsolpar,nprop1*nbxmax*nbxmax,rootid,groupid)
    call mp_bcast(avsolinter,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsolintra,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsolbend,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsoltor,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(avsolelc,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(mnbox,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(asetel,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acipsw,1,rootid,groupid)
    call mp_bcast(acdvdl,1,rootid,groupid)
    call mp_bcast(nblock,1,rootid,groupid)
    call mp_bcast(nccold1,nprop1*nbxmax*nbxmax,rootid,groupid)
    call mp_bcast(bccold1,nprop1*nbxmax*nbxmax,rootid,groupid)
    call mp_bcast(baver1,nprop1*nbxmax*nbxmax*blockm,rootid,groupid)
    call mp_bcast(nccold,nprop*nbxmax,rootid,groupid)
    call mp_bcast(bccold,nprop*nbxmax,rootid,groupid)
    call mp_bcast(baver,nprop*nbxmax*blockm,rootid,groupid)
    call mp_bcast(bnflcq,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bsflcq,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bnflcq2,ntmax*nbxmax,rootid,groupid)
    call mp_bcast(bsflcq2,ntmax*nbxmax,rootid,groupid)

    if (myid.eq.rootid) then
       close(io_chkpt)
       write(UNIT=io_output,FMT="()",ADVANCE='NO',POS=pos_output)
       if (io_movie.ge.0) write(UNIT=io_movie,FMT="()",ADVANCE='NO',POS=pos_movie)
       do i=1,nbox
          if (io_box_movie(i).ge.0) write(UNIT=io_box_movie(i),FMT="()",ADVANCE='NO',POS=pos_box_movie(i))
       end do
       if (io_solute.ge.0) write(UNIT=io_solute,FMT="()",ADVANCE='NO',POS=pos_solute)
       if (io_cell.ge.0) write(UNIT=io_cell,FMT="()",ADVANCE='NO',POS=pos_cell)
       write(UNIT=io_traj,FMT="()",ADVANCE='NO',POS=pos_traj)
    end if
  end subroutine read_checkpoint_main

  subroutine write_checkpoint_main(file_chkpt)
    use util_files,only:flush_force
    use moves_simple,only:write_checkpoint_simple
    use moves_cbmc,only:write_checkpoint_cbmc
    use moves_volume,only:write_checkpoint_volume
    use transfer_swap,only:write_checkpoint_swap
    use transfer_swatch,only:write_checkpoint_swatch
    use transfer_shared,only:write_checkpoint_transfer_shared
    use moves_ee,only:write_checkpoint_ee
    character(LEN=*),intent(in)::file_chkpt
    integer::io_chkpt,jerr,i,pos_output,pos_movie,pos_box_movie(nbox),pos_solute,pos_cell,pos_traj

    pos_output=flush_force(io_output)
    pos_movie=flush_force(io_movie)
    do i=1,nbox
       pos_box_movie(i)=flush_force(io_box_movie(i))
    end do
    pos_solute=flush_force(io_solute)
    pos_cell=flush_force(io_cell)
    pos_traj=flush_force(io_traj)

    io_chkpt=get_iounit()
    open(unit=io_chkpt,access='stream',action='write',file=file_chkpt,form='unformatted',iostat=jerr,status='unknown')
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'cannot open checkpoint file '//trim(file_chkpt),jerr)

    write(io_chkpt) pos_output,pos_movie,pos_box_movie,pos_solute,pos_cell,pos_traj,nnstep,nnn,acmove,vbox,vstart,acdens,molfra,acnbox,acnbox2,acv,acvsq,acvkjmol,acdipole,acdipolesq,acboxa,acboxl,acvol,acvolsq,acvolume,enthalpy,enthalpy2,acnp,acpres,acsurf,acEnthalpy,acEnthalpy1,solcount,acsolpar,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc,mnbox,asetel,acipsw,acdvdl,nblock,nccold1,bccold1,baver1,nccold,bccold,baver

    call write_checkpoint_simple(io_chkpt)
    call write_checkpoint_cbmc(io_chkpt)
    call write_checkpoint_volume(io_chkpt)
    call write_checkpoint_swap(io_chkpt)
    call write_checkpoint_swatch(io_chkpt)
    call write_checkpoint_transfer_shared(io_chkpt)

    ! fluctuating charges
    write(io_chkpt) bnflcq,bsflcq,bnflcq2,bsflcq2

    call write_checkpoint_ee(io_chkpt)

    close(io_chkpt)
  end subroutine write_checkpoint_main
END MODULE topmon_main
