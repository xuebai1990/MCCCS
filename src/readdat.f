subroutine readdat(file_in,lucall,ucheck,nvirial,starvir,stepvir)
  use const_math,only:onepi
  use const_phys,only:MPa2SimUnits
  use util_random,only:ranset
  use util_runtime,only:err_exit
  use util_timings,only:time_date_str
  use util_files,only:get_iounit
  use util_search,only:initiateTable,addToTable
  use sim_system
  use sim_cell
  use zeolite
  use energy_kspace,only:calp
  use energy_intramolecular,only:read_tor_table,read_vib_table,read_bend_table
  use energy_pairwise,only:suijtab,read_vdW_table,read_elect_table,rzero,epsnx
  use energy_3body,only:buildTripletTable
  use energy_4body,only:buildQuadrupletTable
  use energy_sami
  use moves_cbmc,only:llplace
  use moves_ee,only:read_expand,numcoeff,sigm,epsil
  use transfer_shared,only:read_transfer
  implicit none

      character(LEN=*),intent(in)::file_in
      character(LEN=default_path_length)::fileout
      integer::io_input,io_restart,jerr
      integer::seed
      logical::lecho,lverbose
      logical::L_Ewald_Auto
      real::fqtemp
      logical::lmixlb,lmixjo
      integer::nijspecial,ispecial,jspecial
      real::aspecd,bspecd
      integer::ncarbon(ntmax)
      logical::lsetup,linit,lreadq
! -- variables added (3/24/05) for scaling of 1-4 interactions
      integer::nexclu,inclnum,inclmol(ntmax*numax*numax),inclbead(ntmax*numax*numax,2),inclsign(ntmax*numax*numax),ainclnum
      real::ofscale(ntmax*numax*numax),ofscale2(ntmax*numax*numax)
      integer::ainclmol(ntmax*numax*numax),ainclbead(ntmax*numax*numax,2),a15t(ntmax*numax*numax)
      logical::lucall
      integer::ucheck(ntmax)
      integer::nvirial
      real::starvir,stepvir

      real::rcnnsq,umatch,dum ,pm,pcumu,qbox,w(3)
      integer::ndum,ij,ji,ii,jj
      logical::lsolute(ntmax),lhere(nntype)

! -- variables for histograms
      integer::temnc, imol, iutemp, imolty, itype,ipair ,bdum,bin
      integer::idummy(ntmax)

      integer::i,j,k,ncres, nmtres, iensem, inpbc, nmcount
      integer::im,nures,ibox,tcount,nnframe

      integer::izz,temphe,z,zzz
      integer::k_max_l,k_max_m,k_max_n

      real::debroglie,qtot,min_boxl

      logical::lpolar,lqqelect,lee,lratfix
      logical::lprint,lxyz,lfound

      dimension lratfix(ntmax)
      dimension qbox(nbxmax)

! -- Variables added (6/30/2006) for fort.4 consistency check
      integer::numvib,numbend,numtor,vib1,bend2,bend3,tor2,tor3,tor4
      integer::vibtype,bendtype,tortype

!      real::temx,temy,temz

      dimension nures(ntmax)
      dimension temphe(nntype)
!      dimension temx(nmax,numax),temy(nmax,numax),temz(nmax,numax)
      dimension k_max_l(nbxmax),k_max_m(nbxmax),k_max_n(nbxmax)

! KEA torsion variables
      integer::mmm,ttor
! KM tabulated potential variables
      integer::tvib, tbend
! KM variable added when analysis removed
      integer::nhere

! -- reads input data and initializes the positions
!
! --------------------------------------------------------------------

! *** set input arrays to zero ***
      do j=1, ntmax
         lratfix(j) = .false.
         nugrow(j) = 0
         nunit(j) = 0
         do i=1, numax
            ntype(j,i) = 0
         end do
      end do
      do i = 1,nntype**2
         lspecial(i) = .false.
      end do
      do i = 1, nntype
         lhere(i) = .false.
      end do

      lee = .false.
      qtot = 0.0d0
      ldie = .false.

! -------------------------------------------------------------------
      io_input=get_iounit()
      open(unit=io_input,access='sequential',action='read',file=file_in,form='formatted',iostat=jerr,status='old')
      if (jerr.ne.0) then
         call err_exit('cannot open input file')
      end if
      call read_system(io_input)
      call read_transfer(io_input)
      close(io_input)

! -------------------------------------------------------------------
      io_input=get_iounit()
      open(unit=io_input,access='sequential',action='read',file=file_input,form='formatted',iostat=jerr,status='old')
      if (jerr.ne.0) then
         call err_exit('cannot open main input file')
      end if

      read(io_input,*)
      read(io_input,*) seed
! --- initialize random number generator
      call ranset(seed)

! -------------------------------------------------------------------
! *** Output unit (if 2, write to runXX.dat file; if 6, write to stdout/screen; outherwise,
! *** user designate a file to which to write) KEA 6/3/09 (defined in control.inc)
! *** read echoing and long output flags
      read(io_input,*)
      read(io_input,*) ndum,lecho,lverbose,run_num,suffix
      read(io_input,*)
      read(io_input,*) lnpt,lgibbs,lgrand,lanes,lvirial,lmipsw,lexpee
      read(io_input,*)
      read(io_input,*) lijall,lchgall,lewald,ldielect,ltailc,lshift,ltailcZeo

! *** To add or remove helium atoms
      read(io_input,*)
      read(io_input,*) L_add,N_add,N_box2add,N_moltyp2add
      read(io_input,*)
      read(io_input,*) L_sub,N_sub,N_box2sub,N_moltyp2sub

! - read whether to compute electrostatic interaction or not during CBMC/SWAP
      read(io_input,*)
      read(io_input,*) L_Coul_CBMC,lcutcm,ldual,lneigh
!-- read the number of unitcell replicated in each directions (a, b, and c)
      read(io_input,*)
      read(io_input,*) Num_cell_a,Num_cell_b,Num_cell_c
! - read run information
      read(io_input,*)
      read(io_input,*) nstep, lstop, lpresim, iupdatefix
! - read torsion, decide whether to use torsion in function form or Table
      read(io_input,*)
      read(io_input,*) lslit,lexzeo,lzgrid,lelect_field,ljoe,lsami,lmuir,lpsurf,lgraphite,lcorreg,llj,lexpsix,lmmff,lninesix,lgenlj,lgaro,lionic
      read(io_input,*)
      read(io_input,*) L_tor_table,L_spline,L_linear,L_vib_table,L_bend_table,L_vdW_table,L_elect_table

! KM ldielect writes to fort.27
!      if (ldielect) then
!        open(17,file=file_dipole,status='unknown')
!        write(17,*) '# step  ibox   dipole_x   dipole_y   dipole_z'
!      end if

! KM for MPI
! only processor 0 writes data
      if ( lecho.and.myid.eq.0) then
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
      if ( lecho.and.myid.eq.0 ) write(io_output,*) 'number of boxes in the system:',nbox

      read(io_input,*)
      read(io_input,*) (express(ibox),ibox=1,nbox)

      read(io_input,*)
      read(io_input,*) (ghost_particles(ibox),ibox=1,nbox)

      read(io_input,*)
      read(io_input,*) L_Ewald_Auto

      read(io_input,*)
      read(io_input,*) temp, fqtemp,(Elect_field(ibox),ibox=1,nbox)
      if ( lecho.and.myid.eq.0 ) then
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

!$$   read the analysis information
      read(io_input,*)
      read(io_input,*)ianalyze,nbin,lrdf,lintra,lstretch,lgvst,lbend,lete, lrhoz,bin_width
      if (lecho.and.myid.eq.0) then
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

! *** set up constants and conversion factors ***
      beta = 1.0d0 / temp
      fqbeta = 1.0d0 / fqtemp

! ------------------------------------------------------------------
      read(io_input,*)
      read(io_input,*) iprint, imv, iratio, iblock, idiele, L_movie_xyz, iheatcapacity

      if (L_movie_xyz) then
        do ibox = 1,nbox
           write(fileout,'("box",I1.1,"movie",I1.1,A,".xyz")') ibox ,run_num,suffix
           if (myid.eq.0) then
              open (210+ibox,FILE=fileout,status='unknown')
           end if
        end do
      end if

      if ( lecho.and.myid.eq.0 ) then
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
! - read information for histogram output (added 8/30/99)
      if (lgrand) then
         read(io_input,*)
         read(io_input,*) nequil,ninstf, ninsth, ndumph
         if ( lecho.and.myid.eq.0 ) then
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

! KM for MPI
! jobs stop in monola so that all processors die
      if (dint(dble(nstep)/dble(iblock)) .gt. 100) then
         write(io_output,*) 'too many blocks'
         ldie = .true.
         return
      end if

! - read system information
      do i = 1,nbox
         read(io_input,*)
         read(io_input,*) boxlx(i),boxly(i),boxlz(i),lsolid(i),lrect(i), kalp(i),rcut(i),rcutnn(i),numberDimensionIsIsotropic(i)
         if (i.eq.1 .and. lexzeo) then
! === load positions of zeolite atoms
            call zeocoord(file_in,lhere)
            if (myid.eq.0) write(io_output,*) ' note zeolite determines the box size !'
         end if

         if ( lecho.and.myid.eq.0 ) then
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
      if ( nmolty .gt. ntmax ) then
         write(io_output,*) 'nmolty gt ntmax'
         ldie = .true.
         return
      end if
      if ( nchain .gt. nmax-2 ) then
         write(io_output,*) 'nchain gt nmax-2'
         ldie = .true.
         return
      end if
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'number of chains:',nchain
            write(io_output,*) 'number of molecule types:',nmolty
         else
            write(io_output,*) nchain, nmolty
         end if
      end if

      call initiateTable(atoms,nmolty)

      read(io_input,*)
      read(io_input,*) (temtyp(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
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
! - read chemical potentials (added 8/30/99 by jpotoff)
         read(io_input,*)
         read(io_input,*) (B(i),i=1,nmolty)
         if (lecho.and.myid.eq.0) then
            if (lverbose) then
               do i = 1,nmolty
                  write(io_output,*) 'chemical potential for molecule type', i,':',B(i)
               end do
            else
               write(io_output,*) "B ", (B(i),i=1,nmolty)
            end if
         end if

! --- removing from here. It will be calculated once we have molecular mass
! --- to calculate debroglie wavelength.
! - convert chemical potentials to activities
!         do i=1,nmolty
!            B(i) = exp(B(i)/temp)
!         end do
      end if

      if (lgrand) then
         nchain=nmax
         write(io_output,*)'in GCMC total number of chains set by NMAX!'
      end if

      read(io_input,*)
      read(io_input,*) lmixlb, lmixjo
      if (lmixlb .and. lmixjo) then
         write(io_output,*) 'cant use both combining rules!'
         ldie = .true.
         return
      end if
      if ( lecho.and.myid.eq.0 ) then
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
! --- read special combining rule information
      read(io_input,*)
      read(io_input,*) nijspecial
      if (lecho.and.myid.eq.0) then
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
            ij = (ispecial-1)*nntype + jspecial
            ji = (jspecial-1)*nntype + ispecial
            aspecial(ij) = aspecd
            bspecial(ij) = bspecd
            aspecial(ji) = aspecd
            bspecial(ji) = bspecd
            lspecial(ij) = .true.
            lspecial(ji) = .true.
            if (lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(io_output,*) 'special parameter number',i
                  write(io_output,*) '   ispecial:',ispecial
                  write(io_output,*) '   jspecial:',jspecial
                  write(io_output,*) '   aspecd:',aspecd
                  write(io_output,*) '   bspecd:',bspecd
               else
                  write(io_output,*) ispecial,jspecial,aspecd,bspecd
               end if
            end if
         end do
      end if

      read(io_input,*)
      read(io_input,*) rmin, softcut,rcutin, rbsmax,rbsmin
      if ( lecho .and.myid.eq.0) then
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
         if( rcut(i)/boxlx(i) .gt. 0.5d0) then
            write(io_output,*) 'rcut > 0.5*boxlx'
            ldie = .true.
            return
         end if
      end do

      softlog = 10.0d0**(-softcut)
      vol_eff = (4.0d0/3.0d0)*onepi*(rbsmax*rbsmax*rbsmax-rbsmin*rbsmin*rbsmin)

! - set up the forcefield and the masses
      call suijtab(file_in,lmixlb,lmixjo)

! - read bead potential information
      do imol = 1, nmolty
         read(io_input,*)
         read(io_input,*) nunit(imol),nugrow(imol),ncarbon(imol),nmaxcbmc(imol),iurot(imol),lelect(imol),lflucq(imol),lqtrans(imol),lexpand(imol),lavbmc1(imol),lavbmc2(imol),lavbmc3(imol) ,fqegp(imol)
         read(io_input,*)
         read(io_input,*) maxgrow(imol),lring(imol),lrigid(imol) ,lrig(imol),lsetup,isolute(imol),(eta2(i,imol), i=1,nbox)

         read(io_input,*)
         read(io_input,*) lq14scale(imol),qscale(imol)

         if (isolute(imol).lt.nstep) then
            lsolute(imol) = .true.
         else
            lsolute(imol) = .false.
         end if

         if (lring(imol)) then
            read(io_input,*)
            read(io_input,*) iring(imol)
         else
            iring(imol) = nunit(imol)
         end if

         do i = 1, nunit(imol)
            lrigi(imol,i) = .false.
         end do

!     *** irig is the site rigid sites will be grown from
!     *** and frig will be the previous site (not kept rigid)

         if (lrig(imol)) then
            read(io_input,*)
            read(io_input,*) nrig(imol)

            if (nrig(imol).gt.0) then
!     --- read in specific points to keep rigid in growth
               read(io_input,*)
               do i = 1, nrig(imol)
                  read(io_input,*) irig(imol,i),frig(imol,i)
                  lrigi(imol,irig(imol,i)) = .true.
               end do
            else
!     --- we will pick irig at random in each case if nrig = 0
               read(io_input,*)
               read(io_input,*) nrigmin(imol),nrigmax(imol)

!     --- nrigmin is the minimum amount of the chain to keep rigid
!     --- nrigmax is the maximum

            end if
         end if

         if (lrigid(imol)) then
            read(io_input,*)
!     - number of flexible parts
            read(io_input,*) rindex(imol)
            if ( rindex(imol).gt.0) then
               do i = 1, rindex(imol)
                  read(io_input,*) riutry(imol,i)
               end do
            else
               riutry(imol,1) = 1
            end if
         end if


         if ( nunit(imol) .gt. numax ) then
            write(io_output,*) 'nunit gt numax'
            ldie = .true.
            return
         end if
         if ( lflucq(imol) .and. (.not. lelect(imol) ) ) then
            write(io_output,*) 'lelect must be true if flucq is true'
            ldie = .true.
            return
         end if
         if ( lqtrans(imol) ) then
            if (.not. lflucq(imol) ) then
               write(io_output,*) 'lflucq must be true if interm.  CT is allowed'
               ldie = .true.
               return
            end if
            write(io_output,*) 'Intermolecular Charge Transfer is allowed'
         end if

         lbias(imol) = .false.
         if (lavbmc1(imol) .or. lavbmc2(imol) .or. lavbmc3(imol)) then
            lbias(imol) = .true.
         end if

         ! *** choose only one from the three AVBMC algorithms
         if ( lavbmc1(imol) ) then
            lavbmc2(imol) = .false.
            lavbmc3(imol) = .false.
         else if ( lavbmc2(imol) ) then
            lavbmc3(imol) = .false.
         end if

         lneighbor = .false.
         if ( (lavbmc2(imol) .or. lavbmc3(imol)).and.(.not.lgaro) )  lneighbor = .true.

         if ( lecho.and.myid.eq.0 ) then
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
         masst(imol) = 0.0d0

         if (lsetup) then
            call molsetup(io_input,imol)

            do i = 1,nunit(imol)
               lhere(ntype(imol,i)) = .true.
            end do

            goto 112
         end if

         do i = 1, nunit(imol)
! - linear/branched chain with connectivity table -
            read(io_input,*)
            if ( lelect(imol) .and. .not. lchgall ) then
               read(io_input,*) j, ntype(imol,i), leaderq(imol,i)
               ndum=addToTable(atoms,ntype(imol,i),expand=.true.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
! for other potentials, sigi isn't defined, so I added this if/else statement
! maybe should have a switch for every potential flag as a check?
               IF(llj)THEN
               if(sigi(ntype(imol,i)).lt.1d-06.and.epsi(ntype(imol,i)) .lt.1d-06.and.abs(qelect(ntype(imol,i))).lt.1d-06) then
                  write(io_output,*)
                  write(io_output,*) '****PROBLEM IN SUIJTAB****'
                  write(io_output,*) 'check if the beadtyp',ntype(imol,i), 'is defined'
                 ldie = .true.
                 return
               end if
! KM this isn't necessary
!               ELSE
!                  if (myid.eq.0)  write(io_output,*)
!     &                 'Confirm that your parameters are defined.'
               end if

               if ( lecho.and.myid.eq.0 )  write(io_output,*) '   bead ',j,' beadtype ',ntype(imol,i) ,' charge leader ',leaderq(imol,i)
               if ( leaderq(imol,i) .gt. j .and. .not. lchgall) then
                  write(io_output,*) 'group-based cut-off screwed for qq'
                  ldie = .true.
                  return
               end if
            else
               read(io_input,*) j, ntype(imol,i)
               if ( lecho.and.myid.eq.0 ) then
                  if (lverbose) then
                     write(io_output,*) '   bead ',j,' beadtype ', ntype(imol,i),chname(ntype(imol,i))
                  else
                     write(io_output,*) '   bead ',j,' beadtype ', ntype(imol,i)
                  end if
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

! - bond vibration -
            read(io_input,*)
            read(io_input,*) invib(imol,i)
            if ( invib(imol,i) .gt. 6 ) then
               write(io_output,*) 'imol',imol,'   i',i,'  invib',invib(imol,i)
               write(io_output,*) 'too many vibrations'
               ldie = .true.
               return
            end if
            do j = 1, invib(imol,i)
               read(io_input,*) ijvib(imol,i,j),itvib(imol,i,j)
                if(brvib(itvib(imol,i,j)).lt.1d-06) then
                  write(io_output,*)
                  write(io_output,*) '****PROBLEM IN SUVIBE****'
                  write(io_output,*) 'check if the vib type ',itvib(imol,i,j), 'is defined'
                 ldie = .true.
                 return
               end if

               if((ijvib(imol,i,j).eq.i).or.(ijvib(imol,i,j).gt. nunit(imol))) then
                 write(io_output,*) 'check vibrations for mol type',imol, 'and bead',i
                 ldie = .true.
                 return
               end if
               if (lverbose.and.myid.eq.0) then
!                  write(io_output,*) '      bead',i,' bonded to bead',
!     &                 ijvib(imol,i,j),' with bond type:',
!     &                 itvib(imol,i,j)
                  write(io_output,*) '      bead',i,' bonded to bead', ijvib(imol,i,j)

                  write(io_output,'(a20,i3,a13,f9.3,a5,f9.1)') '          bond type:', itvib(imol,i,j),' bond length:', brvib(itvib(imol,i,j)),' k/2:', brvibk(itvib(imol,i,j))
               end if
            end do
! - bond bending -
            read(io_input,*)
            read(io_input,*) inben(imol,i)
!            write(io_output,*) inben(imol,i)
            if ( inben(imol,i) .gt. 12 ) then
               write(io_output,*) 'too many bends'
               ldie = .true.
               return
            end if
            do j = 1, inben(imol,i)
               read(io_input,*) ijben2(imol,i,j),ijben3(imol,i,j) ,itben(imol,i,j)
               if(brben(itben(imol,i,j)).lt.1d-06) then
                 write(io_output,*)
                 write(io_output,*) '****PROBLEM IN SUVIBE****'
                 write(io_output,*) 'check if thE bend type',itben(imol,i,j), 'is defined'
                 ldie = .true.
                 return
               end if


               if ((ijben2(imol,i,j).gt.nunit(imol)).or.( ijben3(imol,i,j).gt.nunit(imol))) then
                   write(io_output,*) 'check bending for the mol type',imol, 'bead',i
                   ldie = .true.
                   return
               end if
               if ((ijben2(imol,i,j).eq.i).or.( ijben3(imol,i,j).eq.i).or.(ijben2(imol,i,j) .eq.ijben3(imol,i,j))) then
                   write(io_output,*) 'check bending for the mol type',imol, 'bead',i
                   ldie = .true.
                   return
               end if


               if (lverbose.and.myid.eq.0) then
!                  write(io_output,*) '      bead',i, ' bending interaction',
!     &                 ' through',ijben2(imol,i,j),' with bead',
!     &                 ijben3(imol,i,j),' of bend type:',itben(imol,i,j)
                  write(io_output,'(1x,a10,i4,a20,a8,i4,a10,i4)') '      bead' ,i, ' bending interaction',' through',ijben2(imol,i,j ),' with bead',ijben3(imol,i,j)
                  write(io_output,'(a20,i3,a13,f9.3,a5,f9.1)') '          bend type:',itben(imol,i,j) ,' bend angle :',brben(itben(imol,i,j))*180.0d0/onepi ,' k/2:',brbenk(itben(imol,i,j))
               end if
            end do
! - bond torsion -
            read(io_input,*)
            read(io_input,*) intor(imol,i)
            if ( intor(imol,i) .gt. 12 ) then
               write(io_output,*) 'too many torsions'
               ldie = .true.
               return
            end if
            do j = 1, intor(imol,i)
               read(io_input,*) ijtor2(imol,i,j),ijtor3(imol,i,j), ijtor4(imol,i,j),ittor(imol,i,j)

               if(ijtor2(imol,i,j).gt.nunit(imol).or.ijtor3(imol,i,j) .gt.nunit(imol).or.ijtor4(imol,i,j).gt.nunit(imol)) then
                   write(io_output,*) 'check torsion for the mol type',imol, 'bead',i
                   ldie = .true.
                   return
               end if

               if((ijtor2(imol,i,j).eq.i.or.ijtor3(imol,i,j).eq.i .or.ijtor4(imol,i,j).eq.i).or.(ijtor2(imol,i,j).eq. ijtor3(imol,i,j).or.ijtor2(imol,i,j).eq. (ijtor4(imol,i,j)).or.(ijtor3(imol,i,j).eq. ijtor4(imol,i,j)))) then
                   write(io_output,*) 'check torsion for the mol type',imol, 'bead',i
                   ldie = .true.
                   return
               end if

               if (lverbose.and.myid.eq.0) then
                  write(io_output,'(1x,a10,i3,a22,a8,i3,a4,i3,a10,i3,a19,i4)') '      bead',i, ' torsional ','interaction through' ,ijtor2(imol,i,j),' and',ijtor3(imol,i,j) ,' with bead',ijtor4(imol,i,j),' of torsional type:' ,ittor(imol,i,j)
               end if
            end do
         end do

 112     continue


! -- starting the self consistency check for the bond vibrations, bending, and torsions
! -- this would help in catching errors in fort.4 connectivity. Starting after continue
! -- so that if we use molsetup subroutine it will provide extra checking. (Neeraj)
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
                     write(io_output,*) 'ERROR IN FORT.4 VIBRATIONS'
                     write(io_output,*) 'Check vibration for mol. type:',imol, 'bead',vib1,'with',i
                     ldie = .true.
                     return
                  end if
                  do k =1,invib(imol,vib1)
                     if(ijvib(imol,vib1,k).eq.i) then
                        lfound = .true.
                        if(vibtype.ne.itvib(imol,vib1,k)) then
                           write(io_output,*) 'Error in fort.4 vibration', ' specifications'
                           write(io_output,*) 'check vibration type of bead',i, 'with',vib1,'molecule type',imol,'vice versa'
                           ldie = .true.
                           return
                        end if
                     end if
                  end do
                  if(.not.lfound) then
                     write(io_output,*) 'Error in fort.4 vibration iformation'
                     write(io_output,*) 'Check vibration for mol. type:',imol, 'bead ',vib1,'with ',i
                     ldie = .true.
                     return
                  end if
               end do
               lfound= .false.
               do j = 1,numbend
                  bend2 = ijben2(imol,i,j)
                  bend3 = ijben3(imol,i,j)
                  !            if ( i .eq. 17) then
                  !              write(io_output,*) bend2, bend3, numbend
                  !            end if
                  bendtype = itben(imol,i,j)
                  if(inben(imol,bend3).eq.0) then
                     write(io_output,*) 'ERROR IN FORT.4 BENDING'
                     write(io_output,*) 'Check bending for mol. type:',imol, 'bead ',bend3,'with ',i
                     ldie = .true.
                     return
                  end if
                  do k = 1,inben(imol,bend3)
                     if((ijben2(imol,bend3,k).eq.bend2).and. (ijben3(imol,bend3,k).eq.i)) then
                        lfound = .true.
                        if(itben(imol,bend3,k).ne.bendtype) then
                           write(io_output,*) 'Error in fort.4 bending', ' specifications'
                           write(io_output,*) 'check bending type of bead',i, 'with',bend3,'mol. typ.',imol,'and vice versa'
                           ldie = .true.
                           return
                        end if
                     end if
                  end do
                  if(.not.lfound) then
                     write(io_output,*) 'Error in fort.4 bending iformation'
                     write(io_output,*) 'Check bending for mol. type:',imol, 'bead ',bend3,'with ',i
                     ldie = .true.
                     return
                  end if
               end do
               lfound = .false.
               do j = 1,numtor
                  tor2 = ijtor2(imol,i,j)
                  tor3 = ijtor3(imol,i,j)
                  tor4 = ijtor4(imol,i,j)
                  tortype = ittor(imol,i,j)
                  if(intor(imol,tor4).eq.0) then
                     write(io_output,*) 'ERROR IN FORT.4 TORSION'
                     write(io_output,*) 'Check torsion for mol. type:',imol, 'bead ',tor4,'with ',i,'and vice versa'
                     ldie = .true.
                     return
                  end if
                  do k = 1,intor(imol,tor4)
                     if((ijtor2(imol,tor4,k).eq.tor3).and.(ijtor3(imol,tor4 ,k).eq.tor2).and.(ijtor4(imol,tor4,k).eq.i)) then
                        lfound=.true.
                        if(ittor(imol,tor4,k).ne.tortype) then
                           write(io_output,*) 'Error in fort.4 torsion', ' specifications'
                           write(io_output,*) 'check torsion type of bead',i, 'with',tor4,'mol. typ.',imol,'and vice versa'
                           ldie = .true.
                           return
                        end if
                     end if
                  end do
                  if(.not.lfound) then
                     write(io_output,*) 'Error in fort.4 torsion iformation'
                     write(io_output,*) 'Check torsion for mol. type:',imol, 'bead ',tor4,'with ',i
                     ldie = .true.
                     return
                  end if
               end do
            end do
         end if

! -- Neeraj Adding molecule neutrality check
!kea skip if lgaro
         if(.not.(lgaro .or.lionic)) then
            do i=1,nmolty
               qtot =0.0d0
               do j = 1,nunit(i)
                  qtot = qtot+qelect(ntype(i,j))
               end do
!            if(dabs(qtot).gt.1d-7) then
!               write(io_output,*)'molecule type',i,'not neutral check charges'
!               ldie = .true.
!               return
!            end if
            end do
         end if


         if ( lexpand(imol) ) then
            lee = .true.
         end if
         if ( lbias(imol) ) then
            read(io_input,*)
            read(io_input,*) pmbias(imol),(pmbsmt(ii),ii=1,nmolty),pmbias2(imol)
            if ( lecho .and.myid.eq.0) then
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
               write(io_output,*)'rbsmax should be greater than rbsmin'
               ldie = .true.
               return
            end if
         end if
!kea 6/4/09 -- added for multiple rotation centers
! -- To assign multiple rotation centers, set iurot(imol) < 0
! -- Add line after molecule specification, avbmc parameters
! -- First, number of rotation centers
! -- Second, identity of centers (0=COM,integer::> 0 = bead number)
! -- Third, give probability to rotate around different centers
         if(iurot(imol).lt.0) then
            read(io_input,*)
            read(io_input,*) nrotbd(imol),irotbd(1:nrotbd(imol),imol), pmrotbd(1:nrotbd(imol),imol)
            if( lecho.and.myid.eq.0 ) then
               if ( lverbose ) then
                  write(io_output,*) ' Multiple rotation centers',nrotbd(imol)
                  do ii=1,nrotbd(imol)
                     write(io_output,*) 'Rotation center',ii,':', irotbd(ii,imol),'    probability to select:', pmrotbd(ii,imol)
                  end do
               end if
            end if
         end if
      end do

      call read_expand()

      if (L_tor_table) then
         call read_tor_table(io_output)
      end if

      if (L_vib_table) then
         call read_vib_table(io_output)
      end if

      if (L_bend_table) then
         call read_bend_table(io_output)
      end if

      if (L_vdW_table) then
         call read_vdW_table(io_output)
      end if

      if (L_elect_table) then
         call read_elect_table(io_output)
      end if

! -- check whether there is a polarizable molecule

      lpolar = .false.
      lqqelect = .false.

      do imol = 1, nmolty
         if (lflucq(imol)) lpolar = .true.
         if (lelect(imol)) lqqelect = .true.
      end do

      if ( .not. lqqelect ) then
         if ( lewald .or. lchgall ) then
            write(io_output,*) 'no charges in the system and turn off lewald'
            ldie = .true.
            return
         end if
      end if

      if ( .not. lpolar ) then
         if ( lanes  ) then
            write(io_output,*) 'lanes should be false for nonpolarizable systems!'
            ldie = .true.
            return
         end if
         if ( lfepsi ) then
            write(io_output,*) 'lfepsi should be false for nonpolarizable systems!'
            ldie = .true.
            return
         end if
      end if

! This B(i) goes in the acceptance rules

      if (lgrand) then
        do i=1,nmolty
            debroglie = 17.458d0/( dsqrt(masst(i)/beta ))
            B(i) = exp(B(i)/temp)/(debroglie*debroglie*debroglie)
         end do
      end if

! - read linkcell information
      read(io_input,*)
      read(io_input,*) licell,rintramax,boxlink

      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'licell:',licell
            write(io_output,*) 'rintramax:',rintramax,' A'
            write(io_output,*) 'boxlink:',boxlink
         else
            write(io_output,*) licell,rintramax,boxlink
         end if
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
      IF(boxlink .LE. nbxmax)THEN
         if (lsolid(boxlink).and.(.not.lrect(boxlink))) then
             write(io_output,*)  'Linkcell not implemented for nonrectangular boxes'
             ldie = .true.
             return
         end if
      end if

! -- read the atomic displacements
      read(io_input,*)
      read(io_input,*) Armtrax, Armtray, Armtraz

      if(lecho.and.myid.eq.0) then
         if(lverbose) then
            write(io_output,'(a41,f8.4,f8.4,f8.4)') 'initial maximum displacements for atoms:',Armtrax, Armtray , Armtraz
         else
           write(io_output,*) Armtrax, Armtray, Armtraz
         end if
       end if

! - read displacement information
      read(io_input,*)
      read(io_input,*) rmtrax(1,1),rmtray(1,1),rmtraz(1,1)
      do im = 1,nbox
         do imol = 1,nmolty
            rmtrax(imol,im) = rmtrax(1,1)
            rmtray(imol,im) = rmtray(1,1)
            rmtraz(imol,im) = rmtraz(1,1)
         end do
      end do
      if ( lecho.and.myid.eq.0 ) then
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
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,'(a41,f8.4,f8.4,f8.4)') ' initial maximum x, y and z rotation:    ',rmrotx(1,1) ,rmroty(1,1),rmrotz(1,1)
         else
            write(io_output,*) rmrotx(1,1), rmroty(1,1), rmrotz(1,1)
         end if
      end if
      read(io_input,*)
      read(io_input,*) tatra,tarot
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'target translational acceptance ratio:',tatra
            write(io_output,*) 'target rotational acceptance ratio:',tarot
         else
            write(io_output,*) tatra, tarot
         end if
      end if

! - read initial setup information
      read(io_input,*)
      read(io_input,*) linit,lreadq
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'linit:',linit
            write(io_output,*) 'lreadq:',lreadq
         else
            write(io_output,*) linit, lreadq
         end if
      end if
      read(io_input,*)
      read(io_input,*) (lbranch(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
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
         if ( lecho.and.myid.eq.0 ) then
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
         read(io_input,*) inix(i),iniy(i),iniz(i),inirot(i),inimix(i), zshift(i),dshift(i),nchoiq(i)
         if ( lecho.and.myid.eq.0 ) then
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

! - read ensemble specific information
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
      if ( lecho.and.myid.eq.0 ) then
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
      if ( lecho.and.myid.eq.0 ) then
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
      if ( lecho.and.myid.eq.0 ) then
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
         if ( lecho.and.myid.eq.0 ) then
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
            if (lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(io_output,*) 'pmvolx:',pmvolx
                  write(io_output,*) 'pmvoly:',pmvoly
               else
                  write(io_output,*) pmvolx,pmvoly
               end if
            end if
         end if
      end do

!     --- read swatch information
      read(io_input,*)
      read(io_input,*) pmswat,nswaty
      if ( nswaty .gt. npamax ) then
         write(io_output,*) 'nswaty gt npamax'
         ldie = .true.
         return
      end if

      if ( lecho.and.myid.eq.0 ) then
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

!        --- safety checks on swatch
         do i = 1,nswaty
            if ( nswatb(i,1) .eq. nswatb(i,2) ) then
               write(io_output,*) 'nswaty ',i,' has identical moltyp'
               write(io_output,*) 'cannot swatch identical moltyp'
               ldie = .true.
               return
            end if
         end do

         if( lecho.and.myid.eq.0 ) then
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
!           --- number of beads that remain in the same position
            read(io_input,*)
            read(io_input,*) nsampos(i),(ncut(i,j),j=1,2)
            if (lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(io_output,*) '   nsampos:',nsampos(i)
                  write(io_output,*) '   ncut:',(ncut(i,j),j=1,2)
               else
                  write(io_output,*) 'nsampos',nsampos(i),' ncut', (ncut(i,j),j=1,2)
               end if
            end if
!           --- bead number
            read(io_input,*)
            do j = 1,nsampos(i)
               read(io_input,*) (splist(i,j,k),k=1,2)
               if (lecho.and.myid.eq.0) then
                  if (lverbose) then
                     write(io_output,*) '   splist:', (splist(i,j,k),k=1,2)
                  else
                     write(io_output,*) 'splist',(splist(i,j,k),k=1,2)
                  end if
               end if
            end do

            read(io_input,*)
            read(io_input,*) (( gswatc(i,j,k), k=1,2*ncut(i,j) ), j=1,2 )
!!            if (lecho.and.myid.eq.0) then
!!               if (lverbose) then
!!                  do izz = 1,ncut(i,j)
!!                     write(io_output,*) '   grow from and prev for ncut',izz,':',
!!     &                    (gswatc(i,j,izz),j=1,2)
!!                  end do
!!               else
!!                  write(io_output,*) 'gswatc',(( gswatc(i,j,k),
!!     &                 k=1,2*ncut(i,j) ), j=1,2 )
!!               end if
!!            end if

            read(io_input,*)
            read(io_input,*) nswtcb(i), (pmswtcb(i,ipair), ipair=1,nswtcb(i))
            read(io_input,*)
            do ipair = 1,nswtcb(i)
               read(io_input,*) box3(i,ipair),box4(i,ipair)
               if (lecho.and.myid.eq.0) then
                  if (lverbose) then
                     write(io_output,*) '   box pair:', box3(i,ipair),box4(i,ipair)
                  else
                     write(io_output,*) box3(i,ipair),box4(i,ipair)
                  end if
               end if
            end do
         end do

      else
!        --- skip past all of the swatch info
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

! --- read swap info
      read(io_input,*)
      read(io_input,*) pmswap, (pmswmt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
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
         if ( lecho.and.myid.eq.0 ) then
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
            if ( lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(io_output,*) '      box pair:', box1(i,ipair), box2(i,ipair)
               else
                  write(io_output,*) 'box pair:', box1(i,ipair), box2(i,ipair)
               end if
            end if
         end do
      end do

! --- read cbmc info
      read(io_input,*)
      read(io_input,*) pmcb, (pmcbmt(i),i=1,nmolty)
      read(io_input,*)
      read(io_input,*) (pmall(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
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
      if (lecho.and.myid.eq.0) then
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
            write(io_output,*) 'a ring can only be used with safe-cbmc'
            ldie = .true.
            return
         end if
      end do
! --- read fluctuating charge info
      read(io_input,*)
      read(io_input,*) pmflcq, (pmfqmt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'pmflcq:',pmflcq
            do i = 1,nmolty
               write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   flcq probability for molecule type ',' ',i ,' (pmfqmt):',pmfqmt(i)
            end do
         else
            write(io_output,*) 'pmflcq',pmflcq,(pmfqmt(i),i=1,nmolty)
         end if
      end if
! --- read expanded-coefficient move info
      read(io_input,*)
      read(io_input,*) pmexpc, (pmeemt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
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
! -- read atom translation probability
      read(io_input,*)
      read(io_input,*) pm_atom_tra
! --- read translation info
      read(io_input,*)
      read(io_input,*) pmtra,(pmtrmt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'pmtra:',pmtra
            do i = 1,nmolty
               write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   translation probability for molecule','   typ e',i ,' (pmtrmt):',pmtrmt(i)
            end do
         else
            write(io_output,*) 'pmtra',pmtra,(pmtrmt(i),i=1,nmolty)
         end if
      end if

! --- read rotation info
      read(io_input,*)
      read(io_input,*) (pmromt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'pmrot:',1.0d0
            do i = 1,nmolty
               write(io_output,'(1x,a41,a5,i4,a10,f8.4)') '   rotational probability for molecule','    typ e',i ,' (pmromt):',pmromt(i)
            end do
         else
            write(io_output,*) (pmromt(i),i=1,nmolty)
         end if
      end if

! *** writeout probabilities, in percent
      if (lverbose.and.myid.eq.0) then
         write(io_output,*)
         write(io_output,*) 'percentage move probabilities:'
         write(io_output,'(1x,a19,f8.2,a2)') 'volume move       :', 100.0d0*pmvol,' %'
         pcumu = pmvol
         if (pmswat .gt. pmvol) then
            pm = pmswat - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         end if
         write(io_output,'(1x,a19,f8.2,a2)') 'swatch move       :', 100.0d0*pm,' %'
         if (pmswap .gt. pmswat) then
            pm = pmswap - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         end if
         write(io_output,'(1x,a19,f8.2,a2)') 'swap move         :', 100.0d0*pm,' %'
         if (pmcb .gt. pmswap) then
            pm = pmcb - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         end if
         write(io_output,'(1x,a19,f8.2,a2)') 'CBMC move         :', 100.0d0*pm,' %'
         if (pmflcq .gt. pmcb) then
            pm = pmflcq - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         end if
         write(io_output,'(1x,a19,f8.2,a2)') 'fluct charge move :', 100.0d0*pm,' %'
         if (pmexpc .gt. pmflcq) then
            pm = pmexpc - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         end if
         write(io_output,'(1x,a19,f8.2,a2)') 'expanded ens move :', 100.0d0*pm,' %'
         if (pmtra .gt. pmexpc) then
            pm = pmtra - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         end if
         write(io_output,'(1x,a19,f8.2,a2)') 'translation move  :', 100.0d0*pm,' %'
         pm = 1.0d0 - pmtra
         write(io_output,'(1x,a19,f8.2,a2)') 'rotation move     :', 100.0d0*pm,' %'
         write(io_output,*)
         write(io_output,*) 'Fraction of atom translations move', pm_atom_tra
      end if

! --- read growth details
      read(io_input,*)
      read(io_input,*) (nchoi1(i),i=1,nmolty),(nchoi(i),i=1,nmolty) ,(nchoir(i),i=1,nmolty) ,(nchoih(i),i=1,nmolty),(nchtor(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
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

      do imol = 1, nmolty
         if (pmfix(imol).gt.0) then
            read(23,*) counttot
!     --- read in from fort.23 the bead-bead distribution
            do i = 1, iring(imol)
               do j = 1, iring(imol)
                  if (i.eq.j) goto 110
                  do bin = 1, maxbin
                     read(23,*) bdum,probf(i,j,bin)
                     hist(i,j,bin) = 0
                  end do
 110              continue
               end do
            end do
         end if
      end do

!     --- read information for CBMC bond angle growth
      read(io_input,*)
      read(io_input,*) (nchbna(i),i=1,nmolty),(nchbnb(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
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
      if ( lecho.and.myid.eq.0 ) then
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

!     ---- Error checking
      do i=1,nmolty
         if ( nchoi1(i) .gt. nchmax ) then
            write(io_output,*) 'nchoi1 gt nchmax'
            ldie = .true.
            return
         end if
         if (nchoi(i) .gt. nchmax ) then
            write(io_output,*) 'nchoi gt nchmax'
            ldie = .true.
            return
         end if
         if ( nchoih(i) .ne. 1 .and. nunit(i) .eq. nugrow(i) )  then
            write(io_output,*) ' nchoih must be one if nunit = nugrow'
            ldie = .true.
            return
         end if
         if ( nchtor(i) .gt. nchtor_max ) then
            write(io_output,*) 'nchtor gt nchtor_max'
            ldie = .true.
            return
         end if
         if ( nchbna(i) .gt. nchbn_max ) then
            write(io_output,*) 'nchbna gt nchbn_max'
            ldie = .true.
            return
         end if
         if ( nchbnb(i) .gt. nchbn_max ) then
            write(io_output,*) 'nchbnb gt nchbn_max'
            ldie = .true.
            return
         end if
         if ( icbsta(i) .gt. numax ) then
            write(io_output,*) 'icbsta gt numax'
            ldie = .true.
            return
         end if
      end do

! --- read exclusion table for intermolecular interactions
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
      if ( lecho.and.myid.eq.0 ) then
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
            if (lverbose.and.myid.eq.0) then
               write(io_output,*) 'excluding interactions between bead',ii, ' on chain',i,' and bead',jj,' on chain',j
            end if
         end do
      else
         read(io_input,*)
      end if

! --- read inclusion list
      read(io_input,*)
      read(io_input,*) inclnum
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'inclnum:',inclnum
         else
            write(io_output,*) inclnum
         end if
      end if
      if (inclnum .ne. 0) then
         do ndum = 1, inclnum
            read(io_input,*) inclmol(ndum),inclbead(ndum,1),inclbead(ndum,2),inclsign(ndum),ofscale(ndum),ofscale2(ndum)
            if (lecho.and.myid.eq.0 ) then
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

! --- read a15 inclusion list
      read(io_input,*)
      read(io_input,*) ainclnum
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(io_output,*) 'ainclnum:',ainclnum
         else
            write(io_output,*) ainclnum
         end if
      end if
      if (ainclnum .ne. 0) then
         do ndum = 1, ainclnum
            read(io_input,*) ainclmol(ndum),ainclbead(ndum,1),ainclbead(ndum,2) ,a15t(ndum)
            if (lecho.and.myid.eq.0) then
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

! - set up the inclusion table
      call inclus(inclnum,inclmol,inclbead,inclsign,ncarbon,ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2)

! -  read in information on the chemical potential checker
      read(io_input,*)
      read(io_input,*) lucall
      read(io_input,*)
      read(io_input,*) (ucheck(jj),jj=1,nmolty)
      if (lecho.and.myid.eq.0) then
         if (lverbose) then
            write(io_output,*) 'lucall:',lucall
            do jj = 1,nmolty
               write(io_output,*) '   ucheck for molecule type', jj,':',ucheck(jj)
            end do
         else
            write(io_output,*) lucall,(ucheck(jj),jj=1,nmolty)
         end if
      end if

! -   read information for virial coefficient calculation
      read(io_input,*)
      read(io_input,*) nvirial,starvir,stepvir
      if (lecho.and.myid.eq.0) then
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
            write(io_output,*) 'nvirial .gt. maxvir'
            ldie = .true.
            return
         end if

         read(io_input,*)
         read(io_input,*) ntemp,(virtemp(jj),jj=1,ntemp)
         if (lecho.and.myid.eq.0) then
            if (lverbose) then
               write(io_output,*) 'ntemp:',ntemp
               write(io_output,*) 'calculation of virial coefficient ', 'at the following temperatures:', (virtemp(jj),jj=1,ntemp)
            else
               write(io_output,*) ntemp, 'Calculation of virial coefficient ', 'at the following temperatures:', (virtemp(jj),jj=1,ntemp)
            end if
         end if
      end if

! --- JLR 11-11-09
! --- reading in extra variables for RPLC simulations
! --- if lideal=.true. then intermolecular interactions are not computed
! --- if ltwice=.true. then mimage is applied twice
! --- if lrplc=.true. there are some special rules in CBMC for how to grow chains
      read(io_input,*)
      read(io_input,*) (lideal(i),i=1,nbox)
      if (lecho.and.myid.eq.0)  then
         write(io_output,*) 'lideal: ', (lideal(i),i=1,nbox)
      end if
      do i = 1,nbox
         if (lideal(i) .and. lexpee) then
            write(io_output,*) 'cannot have lideal and lexpee both true'
            write(io_output,*) 'If you want this you will have change code'
            ldie = .true.
            return
         end if
      end do
      read(io_input,*)
      read(io_input,*) (ltwice(i),i=1,nbox)
      read(io_input,*)
      read(io_input,*) (lrplc(i),i=1,nmolty)
      if (lecho.and.myid.eq.0) then
         write(io_output,*) 'ltwice: ', (ltwice(i),i=1,nbox)
         write(io_output,*) 'lrplc: ', (lrplc(i),i=1,nmolty)
      end if
! --- END JLR 11-11-09

!--- JLR 11-11-09
!--- skip this part if analysis will not be done (if ianalyze .gt. nstep)
! KM 01/10 remove analysis
!      if (ianalyze.lt.nstep) then
!         nhere = 0
!         do izz=1,nntype
!            if ( lhere(izz) ) then
!               nhere = nhere + 1
!               temphe(nhere) = izz
!               beadtyp(nhere)=izz
!            end if
!         end do
!         do izz = 1,nhere
!            atemp = temphe(izz)
!            decode(atemp) = izz
!         end do
!      end if
! --- END JLR 11-11-09
      close(io_input)

! -------------------------------------------------------------------

! *** restart saver
      if ( nstep .gt. 100 ) then
         irsave = nstep / 50
      else
         irsave = nstep + 1
      end if

! -------------------------------------------------------------------

! * check information of .INC-files
      if (myid.eq.0) then
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
               write(io_output,*) 'nchain must equal 2'
               ldie = .true.
               return
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
            write(io_output,*) 'INCONSISTENT ENSEMBLE SPECIFICATION'
            ldie = .true.
            return
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
               write(io_output,*) 'INCONSISTENT PBC SPECIFICATION'
               ldie = .true.
               return
            end if
            write(io_output,*) inpbc,'-dimensional periodic box'
         else
            write(io_output,*) 'cluster mode (no pbc)'
            if ( lgibbs .or. lgrand ) then
               write(io_output,*) 'INCONSISTENT SPECIFICATION OF LPBC  AND ENSEMBLE'
               ldie = .true.
               return
            end if
         end if

         if ( lfold ) then
            write(io_output,*) 'particle coordinates are folded into central box'
            if ( .not. lpbc ) then
               write(io_output,*) 'INCONSISTENT SPECIFICATION OF LPBC AND  LFOLD'
               ldie = .true.
               return
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
               write(io_output,*) 'cannot have lijall with lcutcm'
               ldie = .true.
               return
            end if
!         if ( lchgall ) call err_exit('cannot have lchgall with lcutcm')
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
               write(io_output,*) 'INCONSISTENT SPECIFICATION OF LTAILC AND  LSAMI'
               ldie = .true.
               return
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
            if (lshift) then
               write(io_output,*) ' lshift.and.ltailc!'
               ldie = .true.
               return
            end if
         end if

         if ( lshift) then
            write(io_output,*) 'using a shifted potential'
            if (ltailc) then
               write(io_output,*) ' lshift.and.ltailc!'
               ldie = .true.
               return
            end if
         end if

         write(io_output,*) 'Coulombic inter- and intramolecular interactions'
         if ( lewald ) write(io_output,*) 'Ewald-sum will be used to calculate Coulombic interactions'
         write(io_output,*)
         write(io_output,*) 'MOLECULAR MASS:', (masst(i),i=1,nmolty)
      end if

! -------------------------------------------------------------------

      if (lgrand .and. .not.(lslit)) then
! ---     volume ideal gas box is set arbitry large!
         boxlx(2)=1000*boxlx(1)
         boxly(2)=1000*boxly(1)
         boxlz(2)=1000*boxlz(1)
      end if
      if (linit) then
         do ibox = 1,nbox
            if (lsolid(ibox).and..not.lrect(ibox) .and..not.(ibox.eq.1.and.lexzeo)) then
               write(io_output,*) 'Cannot initialize non-rectangular system'
               ldie = .true.
               return
            end if
         end do
      end if

      if ( linit ) then
         call initia
         nnstep = 0
      else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM needed for long input records, like
! nboxi and moltyp for a lot of molecules
         io_restart=get_iounit()
         open(unit=io_restart,access='sequential',action='read',file=file_restart,form='formatted',iostat=jerr,recl=4096,status='old')
         if (jerr.ne.0) then
            call err_exit('cannot open main restart file')
         end if
         read(io_restart,*) nnstep
         read(io_restart,*) Armtrax, Armtray, Armtraz
         do im = 1,nbox
            do imol = 1,nmolty
               read(io_restart,*) rmtrax(imol,im), rmtray(imol,im) , rmtraz(imol,im)
               read(io_restart,*) rmrotx(imol,im), rmroty(imol,im) , rmrotz(imol,im)
            end do
         end do

         call averageMaximumDisplacement()

         if (myid.eq.0) then
            write(io_output,*)  'new maximum displacements read from restart-file'
            do im = 1,nbox
               if (myid.eq.0) write(io_output,*)'box      #',im
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
         if ( lecho.and.myid.eq.0) then
            do im=1,nbox
               write(io_output,*) 'maximum fluc q displacements: Box #',im
               write(io_output,*) (rmflcq(i,im),i=1,nmolty)
            end do
         end if
! -- changed so fort.77 the same for all ensembles
! -- 06/08/09 KM
!         if ( lgibbs .or. lgrand .or. lnpt ) then
         read(io_restart,*) (rmvol(ibox), ibox = 1,nbox)
         if (myid.eq.0) then
            write(io_output,"(' max volume displacement:        ',3e12.4)") (rmvol(ibox), ibox = 1,nbox)
            write(io_output,*)
            write(io_output,*)
         end if
         do ibox = 1,nbox
            if (lsolid(ibox) .and. .not. lrect(ibox)) then
               read(io_restart,*) (rmhmat(ibox,j),j=1,9)
               if (myid.eq.0) write(io_output,*) (rmhmat(ibox,j),j=1,9)
            end if
         end do

            if (myid.eq.0) then
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
                  if (myid.eq.0) then
                     write(io_output,*)
                     write(io_output,*) 'HMAT COORDINATES FOR BOX',ibox
                     write(io_output,*) (hmat(ibox,j),j=1,3)
                     write(io_output,*) (hmat(ibox,j),j=4,6)
                     write(io_output,*) (hmat(ibox,j),j=7,9)
                     write(io_output,*)
                  end if

                  call matops(ibox)

                  if (myid.eq.0) then
                     write(io_output,*)
                     write(io_output,"(' width of  box ',i1,':',3f12.6)") ibox ,min_width(ibox,1),min_width(ibox,2) ,min_width(ibox,3)
                  end if

                  w(1) = min_width(ibox,1)
                  w(2) = min_width(ibox,2)
                  w(3) = min_width(ibox,3)

                  if (rcut(ibox)/w(1) .gt. 0.5d0 .or.  rcut(ibox)/w(2) .gt. 0.5d0 .or.  rcut(ibox)/w(3) .gt. 0.5d0) then
                     write(io_output,*) 'rcut > half cell width'
                     ldie = .true.
                     return
                  end if

                  if (myid.eq.0) then
                     write(io_output,*)
                     write(io_output,*) 'ibox:  ', ibox
                     write(io_output,'("cell length |a|:",2x,f12.3)') cell_length(ibox,1)
                     write(io_output,'("cell length |b|:",2x,f12.3)')  cell_length(ibox,2)
                     write(io_output,'("cell length |c|:",2x,f12.3)')  cell_length(ibox,3)

                     write(io_output,*)
                     write(io_output,'("cell angle alpha:",2x,f12.3)')  cell_ang(ibox,1)*180.0d0/onepi
                     write(io_output,'("cell angle beta: ",2x,f12.3)') cell_ang(ibox,2)*180.0d0/onepi
                     write(io_output,'("cell angle gamma:",2x,f12.3)') cell_ang(ibox,3)*180.0d0/onepi

!                     write(io_output,"(' angle of  box ',i1,'  :  ','   alpha:   ',f12.6,'
!     &  beta: ', f12.6, '    gamma:   ',f12.6)") ibox,cell_ang(ibox,1)*180.0d0/onepi,
!     &                cell_ang(ibox,2)*180.0d0/onepi,cell_ang(ibox,3)
!     &                *180/onepi
                  end if
               else
                  read(io_restart,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
                  if (myid.eq.0) then
                     write(io_output,*)
                     write(io_output,*)
                     write(io_output,"(' dimension box ',i1,'  :','  a:   ',f12.6,'   b:   ',f12.6,'   c   :  ' ,f12.6)") ibox,boxlx(ibox),boxly(ibox),boxlz(ibox)
                  end if
                  do i = 1, nbox
                     if( (rcut(i)/boxlx(i) .gt. 0.5d0).or. (rcut(i)/boxly(i) .gt. 0.5d0).or. (rcut(i)/boxlz(i) .gt. 0.5d0)) then
                        write(io_output,*) 'rcut > 0.5*boxlx'
                        ldie = .true.
                        return
                     end if
                  end do
               end if
            end do

!         end if ! end if ( lgibbs .or. lgrand .or. lnpt )
         if (myid.eq.0) then
            write(io_output,*)
            write(io_output,*) 'Finished writing simulation box related info'
         end if

         read(io_restart,*) ncres
         read(io_restart,*) nmtres
! --- check that number of particles in fort.4 & fort.77 agree ---
         if ( ncres .ne. nchain .or. nmtres .ne. nmolty ) then
            write(io_output,*) 'conflicting information in restart and control files'
            write(io_output,*) 'nchain',nchain,'ncres',ncres
            write(io_output,*) 'nmolty',nmolty,'nmtres',nmtres
            ldie = .true.
            return
         end if
         read(io_restart,*) (nures(i),i=1,nmtres)

!         do i = 1, nmtres
!            if ( nures(i) .ne. nunit(i) ) then
!               write(io_output,*)
!     +           'conflicting information in restart and control files'
!               write(io_output,*) 'unit',i,'nunit',nunit(i),'nures',nures(i)
!               call err_exit('')
!            end if
!         end do
!         write(io_output,*) 'ncres',ncres,'   nmtres',nmtres

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

!         write(io_output,*) 'start reading coordinates'
! --- check that particles are in correct boxes ---
! --- obtain ncmt values
         do ibox = 1,nbox

            nchbox(ibox) = 0

         end do
         do i = 1, nmolty
            do ibox = 1,nbox
               ncmt(ibox,i) = 0
            end do
            if ( lexpand(i) ) then

! ??? problem in expand ensemble

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
                  write(io_output,*) 'Particle found outside BOX 1'
                  ldie = .true.
                  return
               end if
               nchbox(ibox) = nchbox(ibox) + 1
               imolty = moltyp(i)
               ncmt(ibox,imolty) = ncmt(ibox,imolty) + 1
               if ( lexpand(imolty) ) then
                  if ( ibox .gt. 2 ) then
                     write(io_output,*) 'put in box 1 and 2 for such  molecules'
                     ldie = .true.
                     return
                  end if
                  itype = eetype(imolty)
                  ncmt2(ibox,imolty,itype) =  ncmt2(ibox,imolty,itype) + 1
                  do j = 1,nunit(imolty)
                     sigma(imolty,j) = sigm(imolty,j,itype)
                     epsilon(imolty,j) = epsil(imolty,j,itype)
                  end do
               end if
            else
               write(io_output,*) 'i:',i,'nboxi(i)',nboxi(i)
               write(io_output,*) 'Particle found in ill-defined box'
               ldie = .true.
               return
            end if

         end do

! --- check that number of particles of each type is consistent
         do i = 1, nmolty
            tcount = 0
            do ibox = 1, nbox
               tcount = tcount + ncmt(ibox,i)
            end do
            if ( tcount .ne. temtyp(i) ) then
               write(io_output,*) 'Particle type number inconsistency'
               write(io_output,*) 'type',i
               write(io_output,*) 'ncmt',(ncmt(ibox,i), ibox = 1,nbox)
               write(io_output,*) 'temtyp', temtyp(i)
               ldie = .true.
               return
            end if
         end do

!         write(io_output,*) 'particles found in correct box with correct type'

         do i = 1,nbxmax
            qbox(i) = 0.0d0
         end do

         do i = 1, nchain
!            write(io_output,*) 'reading coord of chain i',i
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
            if ( dabs(qbox(i)) .gt. 1d-6 ) then
               if (i.eq.1.and.lexzeo) cycle
               write(io_output,*) 'box',i,' has a net charge of',qbox(i)
               ldie = .true.
               return
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
                  kalp(ibox) = 6.40d0
                  calp(ibox) = kalp(ibox)/min_boxl
               end do
           else
              write(io_output,*) 'lewald should be true when lchgall is true'
              ldie = .true.
              return
           end if
        else
           do ibox = 1, nbox
              if (lsolid(ibox).and.(.not.lrect(ibox))) then
                  min_boxl = min(min_width(ibox,1),min_width(ibox,2), min_width(ibox,3))
              else
                  min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
              end if
!              rcut(ibox) = 0.4d0*min_boxl
              calp(ibox) = 3.2d0/rcut(ibox)
           end do
        end if
      else
        if ( lchgall ) then
! *** if real::space term are summed over all possible pairs in the box
! *** kalp(1) & kalp(2) are fixed while calp(1) & calp(2) change
! *** according to the boxlength, kalp(1) should have a value greater than
! *** 5.0
          do ibox = 1, nbox
             if (lsolid(ibox).and.(.not.lrect(ibox))) then
                min_boxl = min(min_width(ibox,1),min_width(ibox,2), min_width(ibox,3))
             else
                min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
             end if
             calp(ibox) = kalp(ibox)/min_boxl
             if ( lewald ) then
                if ( kalp(ibox) .lt. 5.6d0 ) then
                   write(io_output,*) 'Warning, kalp is too small'
                   ldie = .true.
                   return
                end if
             else
!kea
                if(.not. lgaro) then
                   write(io_output,*) 'lewald should be true when lchgall  is true'
                   ldie = .true.
                   return
                end if
             end if
          end do
        else
! *** if not lchgall, calp(1) & calp(2) are fixed
         do ibox = 1, nbox
            calp(ibox) = kalp(ibox)
            if ( lewald ) then
               if (calp(ibox)*rcut(ibox).lt.3.2d0.and.myid.eq.0) then
                  write(io_output,*) 'Warning, kalp too small in box',ibox
                  write(io_output,*) ibox,calp(ibox),rcut(ibox)
!cc --- JLR 11-24-09
!cc --- you may want a smaller kalp, e.g. when comparing to previous work
!cc --- This does not need to be a call err_exit('')
!                  call err_exit('kalp is too small, set to 3.2/rcutchg')
                  write(io_output,*) 'kalp is too small, set to 3.2/rcutchg'
!cc --- END JLR 11-24-09
               end if
            end if
         end do
        end if
      end if

      if (L_Ewald_Auto) then
        do ibox = 1,nbox
           if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
              k_max_l(ibox) = dint(boxlx(ibox)*calp(ibox))+1
              k_max_m(ibox) = dint(boxly(ibox)*calp(ibox))+1
              k_max_n(ibox) = dint(boxlz(ibox)*calp(ibox))+1
           else
              k_max_l(ibox) = dint(hmat(ibox,1)*calp(ibox))+2
              k_max_m(ibox) = dint(hmat(ibox,5)*calp(ibox))+2
              k_max_n(ibox) = dint(hmat(ibox,9)*calp(ibox))+2
           end if
        end do
      end if

      if (L_Ewald_Auto.and.myid.eq.0) then
         write(io_output,*)
         write(io_output,*) '****Ewald Parameters*****'
         write(io_output,*) 'ibox   calp(ibox)  kmaxl(ibox)   kmaxm(ibox)', '   kmaxn(ibox)   rcut(ibox)'
         do ibox = 1,nbox
            write(io_output,'(i4,5x,f12.6,3i12,12x,f12.4)') ibox, calp(ibox),  k_max_l(ibox),k_max_m(ibox), k_max_n(ibox),rcut(ibox)
         end do
         write(io_output,*)
      end if

!--- book keeping arrays
      do ibox=1,nbox
         do imolty=1,nmolty
            idummy(imolty)=0
         end do
         do i = 1, nchain
            if ( nboxi(i) .eq. ibox ) then
               imolty=moltyp(i)
               idummy(imolty) = idummy(imolty)+1
               parbox(idummy(imolty),ibox,imolty)=i
!               pparbox(i,ibox) = nmcount
            end if
         end do
         nmcount = 0
         do imolty=1,nmolty
            nmcount = nmcount + idummy(imolty)
         end do
         if ( nmcount .ne. nchbox(ibox) ) then
            write(io_output,*) 'Readdat: nmcount ne nchbox', nmcount, nchbox
            ldie = .true.
            return
         end if
      end do
!     set idummy counter to 0
      do imolty=1,nmolty
         idummy(imolty)=0
      end do
!     set up parall
      do i =1,nchain
         imolty = moltyp(i)
         idummy(imolty) = idummy(imolty) + 1
         parall(imolty,idummy(imolty)) = i
      end do

! -------------------------------------------------------------------
! *** read/produce initial/starting configuration ***
! *** zeolite external potential
      if ( lexzeo ) call suzeo(lhere)
      call buildTripletTable(file_in)
      call buildQuadrupletTable(file_in)
! ----------------------------------------------------------------

! - reordering of numbers for charmm
!      do i = 1, nchain
!         imolty = moltyp(i)
!         temtyp(i) = imolty
!         do j = 1, nunit(imolty)
!            temx(i,j) = rxu(i,j)
!            temy(i,j) = ryu(i,j)
!            temz(i,j) = rzu(i,j)
!         end do
!      end do
!      innew = 0
!      do it = 1, nmolty
!         do i = 1, nchain
!            imolty = temtyp(i)
!            if ( imolty .eq. it ) then
!               innew = innew + 1
!               moltyp(innew) = it
!               do j = 1, nunit(imolty)
!                  rxu(innew,j) = temx(i,j)
!                  ryu(innew,j) = temy(i,j)
!                  rzu(innew,j) = temz(i,j)
!               end do
!            end if
!         end do
!         write(io_output,*) 'it =',it,'   innew =',innew
!      end do

! -------------------------------------------------------------------

! --- set the centers of mass if LFOLD = .TRUE.

      if ( lfold ) then
         do ibox = 1, nbox
            call ctrmas(.true.,ibox,0,6)
         end do
      end if

! * check that rintramax is really valid
      if (licell) then
         do i = 1,nchain
            if (2.0d0*rcmu(i) .gt. rintramax) then
               write(io_output,*) 'rintramax for the linkcell list too small'
               ldie = .true.
               return
            end if
         end do
      end if

! * calculate number of frames *
      nnframe = nstep / imv

! *** write out movie-header ***
      if (myid.eq.0) then
         open(10, file=file_movie, status='unknown')
         if ( nnframe .ne. 0 ) then
            nhere = 0
            do izz=1,nntype
               if ( lhere(izz) ) then
                  nhere = nhere + 1
                  temphe(nhere) = izz
               end if
            end do
            write(10,*) nnframe,nchain,nmolty,nbox,nhere
            write(10,*) (rcut(ibox),ibox=1,nbox)
            write(10,*) (temphe(izz),izz=1,nhere)

            do imolty = 1,nmolty
               write(10,*) nunit(imolty)
!     output bond connectivity information
               do ii=1,nunit(imolty)
                  write(10,*) invib(imolty,ii) ,(ijvib(imolty,ii,z),z=1,invib(imolty,ii))
               end do

!     output torsional connectivity information
               do j = 1,nunit(imolty)
                  write(10,*) intor(imolty,j),(ijtor2(imolty,j,ii), ijtor3(imolty,j,ii),ijtor4(imolty,j,ii),ii=1,intor(imolty,j))
               end do
            end do

         end if
      end if

!     --- write out isolute movie header
      lprint = .false.
      do imol = 1,nmolty
         if (lsolute(imol)) then
            lprint = .true.
         end if
      end do

      if (lprint.and.myid.eq.0) then
         write(11,*) nmolty
         do imol = 1, nmolty
            write(11,*) imol,nunit(imol),(nstep / isolute(imol)) * temtyp(imol)
         end do
      end if


! *** write out initial configuration for first movie frame ***
      if (nnstep .eq. 0) then
         dum = 1.0d0
! * fixed by adding nbox, why the hell didn't this cause errors before?
! * KM fixed by adding dum for acsolpar
            call monper(dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,nbox,nnstep,dum,.false.,.false.,.false. ,.true.,.false.,.false.,lratfix,lsolute,dum,0.0d0,0.0d0)
      end if

! *** calculate constants for lmuir external potential ***
      if ( lmuir ) then
         sigpri = 0.715d0 * dsqrt( 3.8d0 * 3.93d0 )
         c9ch2 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*47.0d0) ) * sigpri**9
         c3ch2 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*47.0d0) ) * sigpri**3
         c9ch3 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*114.0d0) ) * sigpri**9
         c3ch3 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*114.0d0) ) * sigpri**3
         zprmin = ( 3.0d0**(1/6.0d0) ) * sigpri
         v2prmin = c9ch2 / zprmin**9 - c3ch2 / zprmin**3
         v3prmin = c9ch3 / zprmin**9 - c3ch3 / zprmin**3
         betac2 = beta1 - v2prmin
         betac3 = beta1 - v3prmin
         write(io_output,*) 'external potential for Langmuir monolayers used'
         write(io_output,*) 'zprmin',zprmin
         write(io_output,*) 'v2prmin',v2prmin,'v3prmin',v3prmin
      end if

! -------------------------------------------------------------------

! * write out connectivity and bonded interactions
      if (lverbose.and.myid.eq.0) then
         do imol = 1,nmolty
            write(io_output,*) 'molecule type',imol
            if (nunit(imol) .gt. 1) then
               write(io_output,*) '   i   j   type_i type_j   bond length', '        k/2'
            end if
            do i = 1,nunit(imol)
               do j = 1, invib(imol,i)
                  write(io_output,'(i5,i4,i7,i7,f13.4,f14.1)') i,ijvib(imol,i ,j),ntype(imol,i),ntype(imol,ijvib(imol,i,j)) ,brvib(itvib(imol,i,j)),brvibk(itvib(imol,i,j))
               end do
            end do

            if (nunit(imol) .gt. 2) then
               write(io_output,*)
               write(io_output,*) '   i   j   k   type_i type_j type_k', '     angle      k/2'
            end if
            do i = 1,nunit(imol)
               do j = 1,inben(imol,i)
                  write(io_output,'(i5,i4,i4,i7,i7,i7,f12.2,f12.1)') i ,ijben2(imol,i,j),ijben3(imol,i,j),ntype(imol,i) ,ntype(imol,ijben2(imol,i,j)),ntype(imol,ijben3(imol ,i,j)),brben(itben(imol,i,j))*180.0d0/onepi ,brbenk(itben(imol,i,j))
               end do
            end do

            if (nunit(imol) .gt. 3) then
               write(io_output,*)
               write(io_output,*) '   i   j   k   l    type_i type_j type_k', ' type_l     torsion type'
            end if

            do i = 1,nunit(imol)
               do j = 1, intor(imol,i)
                  write(io_output,'(i5,i4,i4,i4,i8,i7,i7,i7,i14)') i ,ijtor2(imol,i,j),ijtor3(imol,i,j),ijtor4(imol,i,j) ,ntype(imol,i),ntype(imol,ijtor2(imol,i,j)) ,ntype(imol,ijtor3(imol,i,j)),ntype(imol,ijtor4(imol ,i,j)),ittor(imol,i,j)
               end do
            end do

         end do
      end if

! * write out non-bonded interaction table
      if (myid.eq.0) then
         write(io_output,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
! added a flag for lgaro...need to add for all potentials, probably
         if ((.not.lexpsix).and.(.not.lmmff).and.(.not.lgaro)) then
            if (lninesix) then
               write(io_output,*) '     i   j    r_0,ij     epsij         q0(i)          q0(j)'
            else if (lgenlj) then
               write(io_output,*) '     i   j    sig2ij     epsij         q0(i)          q0(j)'
            else
               write(io_output,*)  '     i   j    sig2ij     epsij         q0(i)          q0(j)'
            end if
            do i = 1, nntype
               do j = 1,nntype
                  if (lhere(i).and.lhere(j)) then
                     if (lninesix) then
                        ij = (i-1)*nxatom + j
                        write(io_output,'(3x,2i4,2f10.5,2f15.6)') i,j ,rzero(ij),epsnx(ij),qelect(i),qelect(j)
                     else if (lgenlj) then
                        ij = (i-1)*nntype + j
                        write(io_output,'(3x,2i4,2f10.5,2f15.6)') i,j ,dsqrt(sig2ij(ij)) , epsij(ij),qelect(i),qelect(j)
                     else
                        ij = (i-1)*nntype + j
                        write(io_output,'(3x,2i4,2f10.5,2f15.6)') i,j ,dsqrt(sig2ij(ij)) , epsij(ij),qelect(i),qelect(j)
                     end if
                  end if
               end do
            end do
         end if
      end if

      if ( lneigh ) then
! *** If using neighbour list make sure the rcut & rcutnn is the same
! *** for all the boxes
         do ibox = 2,nbox
            if ((dabs(rcut(1)-rcut(ibox)).gt.1.0d-10).and. (dabs(rcutnn(1)-rcut(ibox)).gt.1.0d-10)) then
               write(io_output,*) 'Keep rcut and rcutnn for all the  boxes same'
               ldie = .true.
               return
            end if
         end do

! *** calculate squares of nn-radius and max. update displacement ***
         rcnnsq = rcutnn(1) * rcutnn(1)
         upnn = ( rcutnn(1) - rcut(1) ) / 3.0d0
         upnnsq = upnn * upnn

! *** calculate max. angular displacement that doesn't violate upnn ***
! *** calculate max. all-trans chain length ( umatch ) ***
         umatch = 0.0d0
         do j = 1, nunit(1) - 1
            umatch = umatch + brvib(1)
         end do
         upnndg = asin( upnn / umatch )

        do im=1,2
           do imol=1,nmolty
              if ( rmtrax(imol,im) .gt. upnn ) then
                 write(io_output,*) ' rmtrax greater than upnn',im,imol
                 rmtrax(imol,im) = upnn
              end if
              if ( rmtray(imol,im) .gt. upnn ) then
                 write(io_output,*) ' rmtray greater than upnn',im,imol
                 rmtray(imol,im) = upnn
              end if
              if ( rmtraz(imol,im) .gt. upnn ) then
                 write(io_output,*) ' rmtraz greater than upnn',im,imol
                 rmtraz(imol,im) = upnn
              end if

              if ( rmrotx(imol,im) .gt. upnndg ) then
                 write(io_output,*) ' rmrotx greater than upnndg',im,imol
                 rmrotx(imol,im) = upnndg
              end if
              if ( rmroty(imol,im) .gt. upnndg ) then
                 write(io_output,*) ' rmroty greater than upnndg',im,imol
                 rmroty(imol,im) = upnndg
              end if
              if ( rmrotz(imol,im) .gt. upnndg ) then
                 write(io_output,*) ' rmrotz greater than upnndg',im,imol
                 rmrotz(imol,im) = upnndg
              end if
           end do
        end do
      end if

! *** write input data to unit 6 for control ***
      if (myid.eq.0) then
         write(io_output,*)
         write(io_output,*) 'number of mc cycles:            ', nstep
         write(io_output,*) 'number of chains:               ', nchain
         write(io_output,*)
         write(io_output,*) 'temperature:                    ', temp
         write(io_output,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
         write(io_output,*) 'ex-pressure:                    ',  (express(ibox),ibox=1,nbox)
         write(io_output,*)
      end if


! -------------------------------------------------------------------

      do i = 1,nbox
         if ( rcut(i) .ge. rcutnn(i) .and. lneigh ) then
            write(io_output,*) ' rcut greater equal rcutnn for box',i
            ldie = .true.
            return
         end if
      end do

      return
      end
