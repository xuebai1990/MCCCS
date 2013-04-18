MODULE energy_pairwise
  use var_type,only:dp,default_string_length
  use const_math,only:onepi,twopi,sqrtpi,degrad
  use const_phys,only:qqfact
  use util_math,only:erfunc
  use util_runtime,only:err_exit
  use util_string,only:uppercase
  use util_files,only:readLine
  use util_mp,only:mp_sum,mp_lor,mp_allgather
  use sim_system
  use sim_cell
  use energy_kspace,only:calp,sself,correct
  use energy_intramolecular,only:lininter_bend
  use energy_external,only:U_ext
  use energy_sami,only:ljsami,ljmuir
  use energy_garofalini
  use energy_3body,only:hasThreeBody
  use energy_4body,only:hasFourBody
  implicit none
  private
  save
  public::sumup,energy,boltz,coru,read_ff,init_ff,exsix,ljpsur,type_2body

  real,parameter::a15(2)=(/4.0E7_dp,7.5E7_dp/) !< 1-5 correction term for unprotected hydrogen-oxygen interaction; 1 for ether oxygens, 2 for alcohol oxygens
  !< OLD VALUES: a15(2)=/17.0**6,16.0**6/)
  integer,allocatable::atom_type(:),nonbond_type(:)

  integer,allocatable::vdWsplits(:,:),electsplits(:,:)
  real,allocatable::rvdW(:,:,:),tabvdW(:,:,:),relect(:,:,:),tabelect(:,:,:)
  integer::ntabvdW,ntabelect

! EXPSIX.INC
  integer,parameter::natom=3 !< for exsix
  real,public::aexsix(natom),bexsix(natom),cexsix(natom),sexsix(natom),consu(natom),consp(natom)

! NSIX.INC
! ***************************************************
! stores interaction parameters for 9-6 potential
! ***************************************************
  integer,parameter::nxatom=5 !< for ninesix
  real,public::rzero(nxatom*nxatom),epsnx(nxatom*nxatom),shiftnsix(nxatom*nxatom)

! ***********************************************************
! parameters for Generalized Lennard Jones Potential
  real,public::n0,n1
! repulsive part
  parameter(n0=12.0E0_dp)
! attractive part
  parameter(n1=6.0E0_dp)
! Ref:  J. Chem. Phys. 120, 4994 (2004)
! ***********************************************************

! MERCK.INC
  integer,parameter::natomtyp=3
  real,public::epsimmff(natomtyp),sigimmff(natomtyp),smmff(natomtyp),ammff(natomtyp),nmmff(natomtyp),gmmff(natomtyp),sigisq(natomtyp),alphammff(natomtyp),coru_cons(natomtyp),corp_cons(natomtyp)

contains
!*****************************************************************
!> \brief Calculates the total potential energy for a configuration.
!>
!> \param ovrlap logical, true for substantial atom overlap
!> \param v* energies
!> \param ibox box number
!> \param lvol true if called from volume.f, no output of summary infomation
!******************************************************************
  subroutine sumup(ovrlap,v,ibox,lvol)
    use energy_kspace,only:recipsum
    use energy_intramolecular,only:U_bonded
    use energy_3body,only:U3System
    use energy_4body,only:U4System

    real(kind=dp)::v(nEnergy),vrecipsum,vwell,my_velect
    logical::ovrlap,lvol
    logical::lexplt,lqimol,lqjmol,lcoulo(numax,numax),lij2,liji,lqchgi
    integer::i,imolty,ii,j,jmolty,jj,ntii,ntjj,ntij,iunit,ibox,nmcount,ntj,k,mmm
    real::rcutsq,rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij,rijsq,rho,rij,rbcut,calpi,qqii
    ! real::vtemp
    real::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,vol
    ! Neeraj & RP for MPI
    real::sum_vvib,sum_vbend,sum_vtg
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start SUMUP in ',myid,' for box ', ibox
#endif

    if ( lpbc ) call setpbc(ibox)

    rbcut = rcut(ibox)
    rcutsq = rbcut*rbcut
    calpi=calp(ibox)
    rminsq = rmin * rmin

    ! KM for MPI
    my_velect = 0.0E0_dp

    ovrlap = .false.
    v = 0.0E0_dp
    !kea - 3body garofalini term
    v3garo = 0.0E0_dp
    vwell = 0.0E0_dp
    ! check the molecule count ***
    nmcount = 0
    do i = 1, nchain
       if ( nboxi(i) .eq. ibox ) then
          nmcount=nmcount+1
       end if
       neigh_cnt(i) = 0
    end do
    if ( nmcount .ne. nchbox(ibox) ) then
       call err_exit(__FILE__,__LINE__,'SUMUP: nmcount ne nchbox'//integer_to_string(nmcount)//integer_to_string(nchbox(ibox)),myid+1)
    end if

! ###############################################################
! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************

    ! loop over all chains i
    if (.not.(lgrand.and.ibox.eq.2) .and. .not.lideal(ibox)) then
       ! MPI
       do i = 1, nchain - 1
          ! check if i is in relevant box ###
          if ( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             lqimol = lelect(imolty)

             if ( lcutcm .and. lvol ) then
                xcmi = xcm(i)
                ycmi = ycm(i)
                zcmi = zcm(i)
                rcmi = rcmu(i)
             else
                lij2 = .true.
             end if

             if ( nugrow(imolty) .eq. nunit(imolty) ) then
                lexplt = .false.
             else
                lexplt = .true.
             end if

             ! loop over all chains j with j>i
             molecule2: do j = i + 1 + myid, nchain, numprocs
                ! check for simulation box ###
                if ( nboxi(j) .eq. ibox ) then
                   jmolty = moltyp(j)
                   lqjmol = lelect(jmolty)

                   if (lcutcm .and. lvol ) then
                      ! check if ctrmas within rcmsq
                      rxuij = xcmi - xcm(j)
                      ryuij = ycmi - ycm(j)
                      rzuij = zcmi - zcm(j)
                      ! minimum image the ctrmas pair separations
                      if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)

                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rij  = sqrt(rijsq)
                      rcm = rbcut + rcmi + rcmu(j)
                      rcmsq = rcm*rcm

                      ! if ( lneighbor .and. rcmsq .lt. rbsmax**2 .and. rcmsq .gt. rbsmin**2 ) then
                      ! neigh_cnt(i,jmolty)=neigh_cnt(i,jmolty)+1
                      ! neighbor(neigh_cnt(i,jmolty),i,jmolty)=j
                      ! neigh_cnt(j,imolty)=neigh_cnt(j,imolty)+1
                      ! neighbor(neigh_cnt(j,imolty),j,imolty)=i
                      ! end if

                      if ( rijsq .gt. rcmsq ) then
                         if ( lqimol .and. lqjmol .and. lchgall ) then
                            lij2 = .false.
                         else
                            cycle molecule2
                         end if
                      else
                         lij2 = .true.
                      end if
                   end if

                   do ii = 1,nunit(imolty)
                      ntii = ntype(imolty,ii)
                      liji = lij(ntii)
                      lqchgi = lqchg(ntii)
                      rxui = rxu(i,ii)
                      ryui = ryu(i,ii)
                      rzui = rzu(i,ii)

                      ! loop over all beads jj of chain j
                      bead2: do jj = 1, nunit(jmolty)
                         ! check exclusion table
                         if (lexclu(imolty,ii,jmolty,jj)) cycle bead2

                         ntjj = ntype(jmolty,jj)
                         if ( lij2 ) then
                            if ((.not.(liji.and.lij(ntjj))).and.(.not.(lqchgi.and.lqchg(ntjj)))) cycle bead2
                         else
                            if (.not.(lqchgi.and.lqchg(ntjj))) cycle bead2
                         end if

                         ntij=type_2body(ntii,ntjj)

                         if (lexpee) rminsq=rminee(ntij)*rminee(ntij)

                         rxuij = rxui - rxu(j,jj)
                         ryuij = ryui - ryu(j,jj)
                         rzuij = rzui - rzu(j,jj)
                         ! minimum image the pair separations ***
                         if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)

                         rijsq=(rxuij*rxuij)+(ryuij*ryuij)+(rzuij*rzuij)
                         rij=sqrt(rijsq)

                         if (rijsq.lt.rminsq .and. .not.(lexpand(imolty).or.lexpand(jmolty))) then
                            if ( .not. lvol .and.myid.eq.rootid) then
                               write(io_output,*) 'overlap inter'
                               write(io_output,*)'rijsq rminsq',rijsq,rminsq
                               write(io_output,*) 'i ii', i, ii
                               write(io_output,*) 'i-pos', rxui,ryui,rzui
                               write(io_output,*) 'j jj', j, jj
                               write(io_output,*) 'j-pos',  rxu(j,jj),ryu(j,jj),rzu(j,jj)
                            end if
                            ovrlap = .true.
                            ! RP added for MPI to compensate ovrlap
                            ! return
                            goto 199
                         else if (rijsq.lt.rcutsq .or. lijall) then
                            v(2)=v(2)+U2(rij,rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                         end if

                         v(8)=v(8)+Q2(rij,rijsq,rcutsq,i,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo)

!cc  KM for MPI
!cc  all processors need to know neighbor information
!cc  lneighbor and lgaro will not work in parallel
!cc  calculation of neighbors assumes everything is sequential
                         if ( lneighbor .and. ii .eq. 1 .and.  jj .eq. 1 .and. rijsq .lt. rbsmax**2  .and. rijsq .gt. rbsmin**2 ) then
                            ! neigh_cnt(i,jmolty)=neigh_cnt(i,jmolty)+1
                            ! neighbor(neigh_cnt(i,jmolty),i,jmolty)=j
                            ! neigh_cnt(j,imolty)=neigh_cnt(j,imolty)+1
                            ! neighbor(neigh_cnt(j,imolty),j,imolty)=i
                            neigh_cnt(i)=neigh_cnt(i)+1
                            neighbor(neigh_cnt(i),i)=j
                            neigh_cnt(j)=neigh_cnt(j)+1
                            neighbor(neigh_cnt(j),j)=i
                         else if(lgaro) then
                            if((ntij.eq.4.and.rijsq.lt.grijsq(2,1)) .or.(ntij.eq.6.and.rijsq.lt. grijsq(3,1))) then
                               write(64,*) 'neighbor',i,' (', neigh_cnt(i)+1,')',j,' (', neigh_cnt(j)+1,')'
                               neigh_cnt(i)=neigh_cnt(i)+1
                               neighbor(neigh_cnt(i),i)=j
                               neigh_cnt(j)=neigh_cnt(j)+1
                               neighbor(neigh_cnt(j),j)=i
                               ndij(neigh_cnt(i),i) = rij
                               ndij(neigh_cnt(j),j) =  ndij(neigh_cnt(i),i)
                               nxij(neigh_cnt(i),i) = rxuij
                               nyij(neigh_cnt(i),i) = ryuij
                               nzij(neigh_cnt(i),i) = rzuij
                               nxij(neigh_cnt(j),j) = -rxuij
                               nyij(neigh_cnt(j),j) = -ryuij
                               nzij(neigh_cnt(j),j) = -rzuij
                            end if
                         end if
                      end do bead2
                   end do !do ii = 1,nunit(imolty)
                end if
             end do molecule2
          end if
       end do !do i = 1, nchain - 1

! Returning from ovrlap--------------
199    continue
! KM don't check overlap until allreduce is finished
! if(ovrlap .eq. .true.)then
! write(io_output,*)'521: in sumup ovrlap=',ovrlap,'myid=',myid
! end if
! -----------------------------------------
       call mp_lor(ovrlap,1,groupid)
       if(ovrlap)then
          ! write(io_output,*)'530: in sumup ovrlap=',ovrlap,'myid=',myid
          return
       end if

       call mp_sum(v(2),1,groupid)
       call mp_sum(v(8),1,groupid)

       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro .and..not.L_vdW_table) then
          v(2) = 4.0E0_dp * v(2)
       end if

! KEA garofalini 3 body potential
       if (lgaro) then
          call triad
          call vthreebody(v3garo)
       end if

       if (hasThreeBody) v(2)=v(2)+U3System(ibox)
       if (hasFourBody) v(2)=v(2)+U4System(ibox)

       if (ltailc) then
          ! add tail corrections for the Lennard-Jones energy
          if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
             vol = cell_vol(ibox)
          else
             vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
          end if
          do jmolty = 1, nmolty
             rho = ncmt(ibox,jmolty) / vol
             do imolty = 1, nmolty
                v(3) = v(3) +  ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
             end do
          end do
          v(2) = v(2) + v(3)
       end if
    end if

!$$$c      write(io_output,*)
!$$$c      write(io_output,*) '+++++++'
!$$$c      vtemp = v(8)
!$$$c      write(io_output,*) 'direct space part:',v(8)*qqfact

    if ( ldielect ) then
       call dipole(ibox,0)
    end if

    if (lewald.and..not.lideal(ibox)) then
       call recipsum(ibox,vrecipsum)
       ! update self terms and correction terms
       sself = 0.0E0_dp
       correct = 0.0E0_dp
       ! combine to reduce numerical error
       ! vsc = 0.0E0_dp

       ! MPI
       do imolty=1,nmolty
          do nmcount=myid+1,ncmt(ibox,imolty),numprocs
             i=parbox(nmcount,ibox,imolty)
             do ii = 1,nunit(imolty)
                sself = sself + qqu(i,ii)*qqu(i,ii)
                ! 1.772.. is the square root of pi
                ! vsc = vsc - qqu(i,ii)*qqu(i,ii)*calpi/1.772453851E0_dp
                do jj = ii+1,nunit(imolty)
                   rxuij = rxu(i,ii) - rxu(i,jj)
                   ryuij = ryu(i,ii) - ryu(i,jj)
                   rzuij = rzu(i,ii) - rzu(i,jj)
                   ! JLR 11-17-09  need call to mimage for intrachain
                   if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
                   ! END JLR 11-17-09
                   rij = sqrt(rxuij*rxuij + ryuij*ryuij + rzuij*rzuij)

                   ! correct should only be calculated if ii and jj should NOT interact,
                   ! so only calculating it if lqinclu is false
                   ! this part is 1,2 and 1,3
                   if (.not. lqinclu(imolty,ii,jj)) then
                      correct = correct + qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                      ! vsc = vsc + qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                   else
                      correct=correct+(1.0E0_dp - qscale2(imolty,ii,jj))*qqu(i,ii)*qqu(i,jj)* (erfunc(calpi*rij)-1.0E0_dp)/rij
                      ! vsc = vsc + (1.0E0_dp - qscale2(imolty,ii,jj))*qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                   end if
                end do
             end do
          end do
       end do

       call mp_sum(correct,1,groupid)
       call mp_sum(sself,1,groupid)

! vdipole = (dipolex*dipolex+dipoley*dipoley+dipolez*dipolez)*(2.0E0_dp*onepi)/(3.0E0_dp*boxlx(ibox)**3.0E0_dp)
! write(io_output,*) dipolex,dipoley,dipolez
       sself = -sself*calpi/sqrtpi
       v(8) = v(8) + sself + correct + vrecipsum
! v(8) = v(8) + vsc + vrecipsum
    end if

!$$$c at this point velect contains all intermolecular charge interactions,
!$$$c plus the ewald self term and intramolecular corrections
!$$$
! write(io_output,*)
! write(io_output,*) '== After Inter === velect is:',v(8)*qqfact
!$$$
!$$$       vtemp = v(8)

! ################################################################

! have to recalculate ewald terms if volume changes
    if ( .not. lvol .or. (lvol .and. lewald) ) then
! *******************************
! INTRACHAIN INTERACTIONS ***
! *******************************
       ! write(io_output,*) 'starting intrachain'
       ! loop over all chains i
       ! MPI
       do imolty=1,nmolty
          do nmcount=myid+1,ncmt(ibox,imolty),numprocs
             i=parbox(nmcount,ibox,imolty)
             do ii = 1, nunit(imolty)-1
! write(io_output,*) 'ntype(imolty,ii),ii',ntype(imolty,ii),ii
                ntii = ntype(imolty,ii)
                rxui = rxu(i,ii)
                ryui = ryu(i,ii)
                rzui = rzu(i,ii)
                do jj = ii+1, nunit(imolty)
                   if (linclu(imolty,ii,jj).or.lqinclu(imolty,ii,jj)) then
                      ntjj = ntype(imolty,jj)

                      ntij=type_2body(ntii,ntjj)

                      if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

                      rxuij = rxu(i,ii) - rxu(i,jj)
                      ryuij = ryu(i,ii) - ryu(i,jj)
                      rzuij = rzu(i,ii) - rzu(i,jj)
                      ! JLR 11-17-09  need call to mimage for intrachain
                      if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
                      ! END JLR 11-17-09
                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rij  = sqrt(rijsq)
                      if (linclu(imolty,ii,jj)) then

                         if (rijsq.lt.rminsq.and..not.lexpand(imolty)) then
                            if ( .not. lvol ) then
                               write(io_output,*) 'overlap intra'
                               write(io_output,*) 'rijsq rminsq', rijsq,  rminsq
                               write(io_output,*) 'i ii', i, ii
                               write(io_output,*) 'i-pos', rxui,ryui,rzui
                               write(io_output,*) 'jj', jj
                               write(io_output,*) 'j-pos' ,rxu(i,jj),ryu(i,jj),rzu(i,jj)
                            end if
                            ovrlap = .true.
                            ! RP added for MPI to compensate ovrlap
                            ! return
                            goto 299
                            ! -------------------------------
                         else if (rijsq.lt.rcutsq.or.lijall) then
                            !> \bug skip intra if it is bending 1-3 and using a table?
                            if (L_bend_table) then
                               do mmm=1,inben(imolty,ii)
                                  if (ijben3(imolty,ii,mmm).eq.jj)then
                                     v(4) = v(4) + lininter_bend(rij,itben(imolty,ii,mmm))
                                     goto 94
                                  end if
                               end do
                            end if

                            v(4)=v(4)+U2(rij,rijsq,i,imolty,ii,ntii,i,imolty,jj,ntjj,ntij)
                         end if
                      end if !if (linclu(imolty,ii,jj))

94                    if (lqinclu(imolty,ii,jj)) then
                         ! calculate intramolecular charge interaction
                         my_velect=my_velect+qscale2(imolty,ii,jj)*Q2(rij,rijsq,rcutsq,i,imolty,ii,ntii,lqchg(ntii),i,imolty,jj,ntjj,calpi,lcoulo)
                      end if
                   end if
                end do
             end do
          end do
       end do

! RP added for MPI----- Returning from ovrlap--------------
299    continue

       call mp_lor(ovrlap,1,groupid)
       if(ovrlap)then
          ! write(io_output,*)'941: in sumup ovrlap=',ovrlap,'myid=',myid
          return
       end if
! -----------------------------------------
       call mp_sum(v(4),1,groupid)
       call mp_sum(my_velect,1,groupid)
       v(8) = v(8) + my_velect

!$$$       vtemp = v(8) - vtemp
!$$$
!$$$c       write(io_output,*) '== Intra Velect ===',vtemp*qqfact
! write(io_output,*) '== After Intra  === velect is:',v(8)*qqfact
! write(io_output,*) 'vintra ', v(4)
! write(io_output,*) 'vinter ', v(2)
! write(io_output,*) 'test', v(4)

       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro .and..not.L_vdW_table) then
          v(4) = 4.0E0_dp * v(4)
       end if

! ################################################################

! *************************************
! INTRACHAIN FLUCQ INTERACTIONS ***
! *************************************
!c RP added for MPI
!> \bug removed at some point?
       do i = 1,nchain
! do i = myid+1,nchain,numprocs
! calculate intramolecular flucq energy for chain i
          if( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             if ( lelect(imolty) ) then

                if ( lflucq(imolty) ) then
                   iunit = nunit(imolty)
                   do ii = 1, iunit

                      ntii = ntype(imolty,ii)
                      qqii = qqu(i,ii)
                      do jj = ii, iunit

                         if ( ii .eq. jj) then
                            v(11) = v(11) + xiq(ntii)*qqii + jayself(ntii)*qqii*qqii
                         else
                            ntjj = ntype(imolty,jj)
                            ntij = type_2body(ntii,ntjj)

                            v(11) = v(11)  + jayq(ntij)*qqii*qqu(i,jj)
                         end if
                      end do
                   end do
                   v(11) = v(11) - fqegp(imolty)
                else
                   v(11) = 0.0E0_dp
                end if
             end if
          end if
       end do
! ------------------------------------------

! **************************************************
! CALCULATION OF VIB. + BEND. + TORS. ENERGY ***
! **************************************************

       ! NOTE here virtual coordinates can be used!!!
       ! MPI
       do imolty=1,nmolty
          do nmcount=myid+1,ncmt(ibox,imolty),numprocs
             i=parbox(nmcount,ibox,imolty)
             call U_bonded(i,imolty,sum_vvib,sum_vbend,sum_vtg)
             v(5)=v(5)+sum_vvib
             v(6)=v(6)+sum_vbend
             v(7)=v(7)+sum_vtg
          end do
       end do
       call mp_sum(v(5),1,groupid)
       call mp_sum(v(6),1,groupid)
       call mp_sum(v(7),1,groupid)
    end if !if ( .not. lvol .or. (lvol .and. lewald) )

! ################################################################

! ***************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

! for adsorption isotherms, don't calculate energy w/surface
! in box 2
    if ((lelect_field) .or. ((ibox .eq. 1) .and. (ljoe .or. lsami .or. lmuir .or.  lexzeo .or. lgraphite .or. lslit))) then
       ! MPI
       do imolty=1,nmolty
          do nmcount=myid+1,ncmt(ibox,imolty),numprocs
             i=parbox(nmcount,ibox,imolty)
             do j = 1, nunit(imolty)
                ntj = ntype(imolty,j)
                v(9)=v(9)+U_ext(ibox,i,j,ntj)
             end do
          end do
       end do

       call mp_sum(v(9),1,groupid)
    end if

! --------------------------------------------------------------------
! calculation of additional gaussian potential needed in thermodynamic
! integration in stages b and c
! --------------------------------------------------------------------
    if (lmipsw) then
       ! MPI
       do imolty=1,nmolty
          !> \bug most likely only works when nbox = 1
       do nmcount=myid+1,ncmt(ibox,imolty),numprocs
          i=parbox(nmcount,ibox,imolty)
          if (lwell(imolty)) then
             rxui = xcm(i)
             ryui = ycm(i)
             rzui = zcm(i)
             do j = 1, nwell(imolty)*nunit(imolty)
                k = j-int(j/nunit(imolty))*nunit(imolty)
                if (k.eq.0) k = nunit(imolty)
                rxuij = rxui-rxwell(j,imolty)
                ryuij = ryui-rywell(j,imolty)
                rzuij = rzui-rzwell(j,imolty)
                call mimage(rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                rcm = rcut(ibox)+rcmu(i)
                rcmsq = rcm*rcm
                if (rijsq.lt.rcmsq) then
                   do ii = 1, nunit(imolty)
                      if (awell(ii,k,imolty).lt.1.0E-6_dp) cycle
                      rxui = rxu(i,ii)
                      ryui = ryu(i,ii)
                      rzui = rzu(i,ii)
                      rxuij = rxui-rxwell(j,imolty)
                      ryuij = ryui-rywell(j,imolty)
                      rzuij = rzui-rzwell(j,imolty)
                      call mimage(rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                      vwell = vwell-awell(ii,k,imolty)*exp(-bwell*rijsq)
                   end do
                end if
             end do
          end if
       end do
       end do

       call mp_sum(vwell,1,groupid)
    end if

! ----------------------------------------------------------------------------

! write(io_output,*) 'self,corr:',(v(8)-vrecipsum)*qqfact
! write(io_output,*) 'vsc, new self cor:',vsc*qqfact
! write(io_output,*) 'recip space part :',vrecipsum*qqfact
! write(io_output,*) 'sc and recip:',(vsc+vrecipsum)*qqfact

    if (.not.L_elect_table) then
       v(8) = v(8)*qqfact
    end if

    v(1) = v(2) + v(4) + v(9) + v(8) + v(11) + v3garo

! write(io_output,*) 'v in sumup',v(1)

    vipsw = v(1)
    vwellipsw = vwell

    if (lstagea) then
       v(1) = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*v(1)
    else if (lstageb) then
       v(1) = etais*v(1)+lambdais*vwell
    else if (lstagec) then
       v(1) = (etais+(1.0E0_dp-etais)*lambdais)*v(1)+(1.0E0_dp-lambdais)*vwell
    end if

    v(1) = v(1) + v(5) + v(6) + v(7)

    if ( .not. lvol.and.myid.eq.rootid ) then
       write(io_output,*)
       write(io_output,*) 'sumup control'
       write(io_output,*) 'number of chains', nchbox(ibox)
       do i = 1, nmolty
          write(io_output,*) 'number of chains of type',i,ncmt(ibox,i)
       end do
       write(io_output,*) 'inter lj energy ', v(2)
       write(io_output,*) 'intra lj energy ', v(4)
       if (ltailc) write(io_output,*) 'Tail correction ', v(3)
       write(io_output,*) 'bond vibration  ', v(5)
       write(io_output,*) 'bond bending    ', v(6)
       write(io_output,*) 'torsional       ', v(7)
       write(io_output,*) 'external        ', v(9)
       write(io_output,*) 'coulombic energy', v(8)
       ! write(io_output,*) 'exact energy    ', 1.74756*1.67*831.441/3.292796
       write(io_output,*) 'fluc Q energy   ', v(11)
       write(io_output,*) 'well energy     ', vwellipsw
       if(lgaro) write(io_output,*) '3-body garo     ', v3garo
       write(io_output,*) 'total energy    ', v(1)
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end SUMUP in ',myid
#endif

    return
  end subroutine sumup

!*****************************************************************
!> \brief Calculates the total potential energy for a configuration.
!>
!> \param i calculate the energies associated with chain i
!> \param imolty molecule type of chain i
!> \param v* energies
!> \param flagon flag for old(flagon=1)/new(flagon=2) configurations
!> \param ibox box number of chain i
!> \param istart, iuend calculate energies from bead istart to bead iuend for chain i
!> \param lljii whether to include intramolecular LJ interactions
!> \param ovrlap atom overlap
!> \param ltors whether to calculate torsional energy
!> \param lcharge_table whether need to set up charge interaction table; true if called from CBMC
!> \param lfavor
!*****************************************************************
  subroutine energy(i,imolty,v,flagon,ibox,istart,iuend,lljii,ovrlap,ltors,lcharge_table,lfavor,lAtom_traxyz)
    use sim_particle,only:lnn,ctrmas
    use energy_intramolecular,only:U_torsion
    use energy_3body,only:U3MolSys
    use energy_4body,only:U4MolSys
    use energy_garofalini,only:triad_en

    logical::lqimol,lqjmol,lexplt,lcoulo(numax,numax),lfavor,lij2,liji,lqchgi
    logical::lljii,ovrlap,ltors,lcharge_table,lAtom_traxyz

    integer::growii,growjj,k,jcell(nmax),nmole
    integer::i,ibox,istart,iuend,ii,ntii,flagon,jjj,iii,mmm,j,jj,ntjj,ntij,ntj,imolty,jmolty,jjend
    integer::nchp2

    real::v(nEnergy),rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq,ryuij,rzuij,rij,rijsq,rbcut,calpi
    real::vwell
    real::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq

    ! KEA
    integer::neigh_j,neighj(maxneigh)
    real::ndijj(maxneigh),nxijj(maxneigh),nyijj(maxneigh),nzijj(maxneigh)
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start ENERGY in ',myid,' for molecule ',i,' in box ',ibox
#endif

    if ( lpbc ) call setpbc(ibox)

    rbcut = rcut(ibox)
    rcutsq = rbcut * rbcut
    calpi = calp(ibox)
    if (ldual) rcinsq = rcutin*rcutin
    rminsq = rmin * rmin

    ovrlap = .false.
    v = 0.0E0_dp
    sself  = 0.0E0_dp
    correct = 0.0E0_dp
    !kea
    v3garo = 0.0E0_dp

    if ( istart .eq. 1 .and. flagon .eq. 2) then
       neigh_icnt = 0
    else if(lgaro.and.flagon.eq.1) then
       neigh_j = 0
    end if

! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************
    lqimol = lelect(imolty)

    if (nugrow(imolty) .eq. nunit(imolty)) then
       lexplt = .false.
    else
       lexplt = .true.
       growii = nugrow(imolty)
    end if

    nchp2 = nchain + 2
    do ii = 1,nunit(imolty)
       rxu(nchp2,ii) = rxuion(ii,flagon)
       ryu(nchp2,ii) = ryuion(ii,flagon)
       rzu(nchp2,ii) = rzuion(ii,flagon)
       qqu(nchp2,ii) = qquion(ii,flagon)
    end do
    nboxi(nchp2) = ibox
    moltyp(nchp2) = imolty

    if ( lcutcm .or. lfavor ) then
       ! calculate the center of mass of chain i and give it a dummy #
       call ctrmas(.false.,ibox,nchp2,9)
       xcmi = xcm(nchp2)
       ycmi = ycm(nchp2)
       zcmi = zcm(nchp2)
       rcmi = rcmu(nchp2)
       ! write(io_output,*) 'rcmi:',rcmi
    else
       lij2 = .true.
    end if

    if (licell.and.(ibox.eq.boxlink)) then
       ii=1
       rxui = rxuion(ii,flagon)
       ryui = ryuion(ii,flagon)
       rzui = rzuion(ii,flagon)
       call get_cell_neighbors(rxui,ryui,rzui,ibox,jcell,nmole)
    else
       nmole = nchain
    end if

! loop over all chains except i - not for grand can. with ibox=2 !
! JLR 11-24-09  also don't loop if box is ideal gas
    if (.not.(lgrand.and.(ibox.eq.2)) .and. .not.(lideal(ibox))) then
! END JLR 11-24-09
! RP added for MPI
       do k = myid+1,nmole,numprocs
! do k = 1, nmole
          if (licell.and.(ibox.eq.boxlink)) then
             j = jcell(k)
          else
             j = k
          end if

          jmolty = moltyp(j)
          lqjmol = lelect(jmolty)
          growjj = nugrow(jmolty)

! check for simulation box ###
          if ( ( ibox .eq. nboxi(j) ) .and. (i .ne. j )) then
             if ( lneigh ) then
                if ( .not. lnn(j,i) ) cycle
             end if

             if (lcutcm .or. lfavor) then
                ! check if ctrmas within rcmsq
                rxuij = xcmi - xcm(j)
                ryuij = ycmi - ycm(j)
                rzuij = zcmi - zcm(j)
                ! minimum image the ctrmas pair separations ***
                if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                rij  = sqrt(rijsq)
                rcm = rbcut + rcmi + rcmu(j)
                rcmsq = rcm*rcm
                ! write(io_output,*) rcm,rcmi,rcmu(j)
                if ( lfavor ) then
                   favor(j) = (rminsq/rijsq)**2*5.0E0_dp
                   favor2(j) = rminsq/rijsq
                end if
                if ( rijsq .gt. rcmsq .and. lcutcm) then
                   if ( lqimol .and. lqjmol .and. lchgall ) then
                      lij2 = .false.
                   else
                      cycle
                   end if
                else
                   lij2 = .true.
                end if
             end if

             if ( lcharge_table .and. (.not. lchgall) ) then
                ! called from CBMC and must set up charge-interaction table ---
                do ii = 1,nugrow(imolty)
                   do jj = 1,nugrow(jmolty)
                      iii = leaderq(imolty,ii)
                      jjj = leaderq(jmolty,jj)
                      if ( iii .eq. ii .and. jjj .eq. jj ) then
                         rxuij = rxuion(ii,flagon) - rxu(j,jj)
                         ryuij = ryuion(ii,flagon) - ryu(j,jj)
                         rzuij = rzuion(ii,flagon) - rzu(j,jj)
                         if ( lpbc )  call mimage(rxuij,ryuij,rzuij,ibox)
                         rijsq = rxuij*rxuij + ryuij*ryuij  + rzuij*rzuij
                         if ((rijsq .lt. rcutsq) .or. lijall) then
                            lcoulo(ii,jj) = .true.
                         else
                            lcoulo(ii,jj) = .false.
                         end if
                      end if
                   end do
                end do
             end if

             ! loop over all beads ii of chain i
             do ii = istart, iuend
                ntii = ntype(imolty,ii)
                liji = lij(ntii)
                lqchgi = lqchg(ntii)
                rxui = rxuion(ii,flagon)
                ryui = ryuion(ii,flagon)
                rzui = rzuion(ii,flagon)

                ! loop over all beads jj of chain j
                do jj = 1, nunit(jmolty)
                   ! check exclusion table
                   if ( lexclu(imolty,ii,jmolty,jj) ) cycle

                   ntjj = ntype(jmolty,jj)
                   if ( lij2 ) then
                      if ((.not.(liji.and.lij(ntjj))).and.(.not.(lqchgi.and.lqchg(ntjj))))  cycle
                   else
                      if (.not.(lqchgi.and.lqchg(ntjj))) cycle
                   end if

                   ntij=type_2body(ntii,ntjj)

                   if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

                   rxuij = rxui - rxu(j,jj)
                   ryuij = ryui - ryu(j,jj)
                   rzuij = rzui - rzu(j,jj)
                   ! minimum image the pair separations ***
                   if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                   rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   rij  = sqrt(rijsq)
                   if (rijsq.lt.rminsq .and. .not.(lexpand(imolty).or.lexpand(jmolty))) then
                      ovrlap = .true.
                      ! write(io_output,*) 'inter ovrlap:',i,j, myid
                      ! write(io_output,*) 'i xyz',rxui,ryui,rzui
                      ! write(io_output,*) 'j xyz',rxu(j,jj),ryu(j,jj),rzu(j,jj)
                      ! write(io_output,*) 'ii:',ii,'jj:',jj
                      ! write(io_output,*) 'distance', sqrt(rijsq)
                      ! RP added for MPI
                      ! return
                      goto 99
                   else if ((rijsq .lt. rcutsq) .or. lijall) then
                      v(2)=v(2)+U2(rij,rijsq,nchp2,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                   end if

                   v(8)=v(8)+Q2(rij,rijsq,rcutsq,nchp2,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo)

                   ! KM lneighbor and lgaro does not work in parallel
                   if ( lneighbor .and. ii .eq. 1 .and.  jj .eq. 1 .and. flagon .eq. 2 .and. rijsq .lt. rbsmax**2  .and. rijsq .gt. rbsmin**2) then
                      ! neigh_icnt1(jmolty)=neigh_icnt1(jmolty)+1
                      ! neighi1(neigh_icnt1(jmolty),jmolty)=j
                      neigh_icnt=neigh_icnt+1
                      neighi(neigh_icnt)=j
                   else if(lgaro) then
                      if((ntij.eq.4.and.rijsq.lt.grijsq(2,1)).or. (ntij.eq.6.and.rijsq.lt.grijsq(3,1))) then
                         if(flagon.eq.2) then
                            neigh_icnt=neigh_icnt+1
                            neighi(neigh_icnt)=j
                            ndiji(neigh_icnt) = rij
                            nxiji(neigh_icnt) = rxuij
                            nyiji(neigh_icnt) = ryuij
                            nziji(neigh_icnt) = rzuij
                         else if(flagon.eq.1) then
                            neigh_j = neigh_j+1
                            neighj(neigh_j) = j
                            ndijj(neigh_j) = rij
                            nxijj(neigh_j) = rxuij
                            nyijj(neigh_j) = ryuij
                            nzijj(neigh_j) = rzuij
                         end if
                      end if
                   end if
                end do
             end do
          end if
       end do
    end if

! RP added for MPI
! Returning from ovrlap--------------
99  continue

    call mp_lor(ovrlap,1,groupid)
    if(ovrlap)then
       ! write(io_output,*)'630 in energy ovrlap=',ovrlap,'myid=',myid
       return
    end if

    call mp_sum(v(2),1,groupid)
    call mp_sum(v(8),1,groupid)
! -----------------------------------------------

    if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro .and..not.L_vdW_table) then
       v(2) = 4.0E0_dp * v(2)
    end if

!kea - garo: add three body loop for intermolecular interactions
    if (lgaro.and..not.lideal(ibox)) then
       if(flagon.eq.2) then
          call triad_en(i,v3garo,neigh_icnt,neighi,ndiji,nxiji,nyiji,nziji,.true.)
       else if(flagon.eq.1) then
          call triad_en(i,v3garo,neigh_j,neighj,ndijj,nxijj,nyijj,nzijj,.false.)
       end if
    end if

    if (hasThreeBody) v(2)=v(2)+U3MolSys(i,istart,iuend,flagon)
    if (hasFourBody) v(2)=v(2)+U4MolSys(i,istart,iuend,flagon)

! ################################################################

! the intramolecular van der waals and ewald terms have to be calculated
! for the explicit atom placement models
! *******************************
! INTRACHAIN INTERACTIONS ***
! *******************************

! JLR 11-19-09 commenting this out, alway do mimage for intrachain
! for expanded ensemble
! lmim = .false.
! mlen2 = rcmu(nchp2)*2E0_dp
! if ( mlen2>boxlx(ibox) .or. mlen2>boxly(ibox) .or. mlen2>boxlz(ibox)) lmim = .true.
! END JLR 11-19-09

! calculate intramolecular energy correction for chain i
    do ii = istart, iuend
       ntii = ntype(imolty,ii)
       rxui = rxuion(ii,flagon)
       ryui = ryuion(ii,flagon)
       rzui = rzuion(ii,flagon)
       if (lAtom_traxyz) then
          jjend=nunit(imolty)
       else
          jjend=ii-1
       end if
       do jj = 1,jjend
          if (jj.eq.ii) cycle
          ntjj = ntype(imolty,jj)

          ntij=type_2body(ntii,ntjj)

          if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

          rxuij = rxuion(ii,flagon) - rxuion(jj,flagon)
          ryuij = ryuion(ii,flagon) - ryuion(jj,flagon)
          rzuij = rzuion(ii,flagon) - rzuion(jj,flagon)
          ! JLR 11-19-09 always do mimage for intrachain
          ! if (lpbc .and. lmim) call mimage( rxuij,ryuij,rzuij,ibox)
          if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
          ! END JLR 11-19-09
          rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
          rij = sqrt(rijsq)

          if (lqinclu(imolty,ii,jj) ) then
             ! calculation of intramolecular electrostatics
             v(8)=v(8)+qscale2(imolty,ii,jj)*Q2(rij,rijsq,rcutsq,nchp2,imolty,ii,ntii,lqchg(ntii),nchp2,imolty,jj,ntjj,calpi,lcoulo)
          end if

          ! calculation of other non-bonded interactions
          if ( linclu(imolty,ii,jj) ) then
             if (lljii) then
                if (rijsq.lt.rminsq .and. .not.lexpand(imolty)) then
                   ovrlap = .true.
                   ! write(io_output,*) 'intra ovrlap:',ii,jj
                   return
                else if (rijsq.lt.rcutsq .or. lijall) then
                   !> \bug skip intra if it is bending 1-3 and using a table?
                   if (L_bend_table) then
                      do mmm=1,inben(imolty,ii)
                         if (ijben3(imolty,ii,mmm).eq.jj) then
                            v(4) = v(4) + lininter_bend(rij,itben(imolty,ii,mmm))
                            goto 96
                         end if
                      end do
                   end if

                   v(4)=v(4)+U2(rij,rijsq,nchp2,imolty,ii,ntii,nchp2,imolty,jj,ntjj,ntij)
                end if
             end if
          end if

96        if (lewald.and..not.lideal(ibox)) then
             ! compute the ewald intramolecular (self and correction) terms for
             ! the interactions of the placed atoms with themselves, and with the
             ! rest of their own molecule, if there's no interaction
             ! these are 1,2 and 1,3
             if (.not. lqinclu(imolty,ii,jj)) then
                correct=correct+qquion(ii,flagon)*qquion(jj,flagon)*(erfunc(calpi*rij)-1.0E0_dp)/rij
                ! 1,4 interaction which we scale by qscale
             else
                correct=correct+(1.0E0_dp-qscale2(imolty,ii,jj))*qquion(ii,flagon)*qquion(jj,flagon)*(erfunc(calpi*rij)-1.0E0_dp)/rij
             end if
          end if
       end do
       if (lewald.and..not.lideal(ibox)) then
          sself = sself + qquion(ii,flagon)*qquion(ii,flagon)
       end if
    end do
    if (lewald.and..not.lideal(ibox)) then
       sself = -sself * calpi/sqrtpi
       v(14) = sself + correct
    end if
    if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj  .and. .not. lninesix .and..not.L_vdW_table) then
       v(4) = 4.0E0_dp * v(4)
    end if

! ################################################################

! ***************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

    if ((lelect_field.and.lqimol) .or. ((ibox .eq. 1) .and. (ljoe .or. lsami .or. lmuir .or. lexzeo .or. lgraphite .or. lslit))) then
       do j = istart,iuend
          ntj = ntype(imolty,j)
          v(9)=v(9)+U_ext(ibox,nchp2,j,ntj)
       end do
    end if

! *********************************************************************
! calculation of torsion energy for explicit atom methyl groups ****
! *********************************************************************
    if ( ltors ) then
       v(7)=U_torsion(nchp2,imolty,nugrow(imolty)+1,.true.)
    end if

! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
    vwell = 0.0E0_dp
    if (lwell(imolty).and.lmipsw) then
       rxui = xcmi
       ryui = ycmi
       rzui = zcmi
       do j = 1, nwell(imolty)*nunit(imolty)
          k = j - int(j/nunit(imolty))*nunit(imolty)
          if (k.eq.0) k = nunit(imolty)
          rxuij = rxui-rxwell(j,imolty)
          ryuij = ryui-rywell(j,imolty)
          rzuij = rzui-rzwell(j,imolty)
          call mimage(rxuij,ryuij,rzuij,ibox)
          rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
          rcm = rcut(ibox)+rcmi
          rcmsq = rcm*rcm
          if (rijsq.lt.rcmsq) then
             do ii = 1, nunit(imolty)
                if (awell(ii,k,imolty).lt.1.0E-6_dp) cycle
                rxui = rxuion(ii,flagon)
                ryui = ryuion(ii,flagon)
                rzui = rzuion(ii,flagon)
                rxuij = rxui-rxwell(j,imolty)
                ryuij = ryui-rywell(j,imolty)
                rzuij = rzui-rzwell(j,imolty)
                call mimage(rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                vwell = vwell-awell(ii,k,imolty)*exp(-bwell*rijsq)
             end do
          end if
       end do
    end if

! ----------------------------------------------------------------------------

    if (.not.L_elect_table) then
       v(8) = v(8)*qqfact
       v(14) = v(14)*qqfact
    end if

! note that vintra is only computed when the flag lljii is true
    v(1) = v(2) + v(9) + v(4) + v(8) + v(14) + v3garo
! write(io_output,*) 'vinter:',v(2),'vext:',v(9),'vintra:',v(4),'velect',v(8),'vewald:',v(14),'v'

    if (flagon.eq.1) then
       vipswo = v(1)
       vwellipswo = vwell
    else
       vipswn = v(1)
       vwellipswn = vwell
    end if

    if (lmipsw) then
       if (lstagea) then
          v(1) = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*v(1)
       else if (lstageb) then
          v(1) = etais*v(1)+lambdais*vwell
       else if (lstagec) then
          v(1) = (etais+(1.0E0_dp-etais)*lambdais)*v(1)+(1.0E0_dp-lambdais)*vwell
       end if
    end if

#ifdef __DEBUG__
    ! write(io_output,*) 'v :', v(1)
    write(io_output,*) 'end ENERGY in ',myid
#endif
    return
  end subroutine energy

!*****************************************************************
!> \brief Calculates the potential energy and the boltzmann factor
!>       for ichoi trial positions.
!>
!> \param lnew true for new configurations
!> \param lfirst true for insertion of the first bead in swap moves
!> \param ovrlap logical variable, true for walk termination
!> \param i calculates Boltzmann weights for the newly-grown beads in chain i
!> \param icharge usually identical to i
!> \param imolty molecule type of chain i
!> \param ibox box number of chain i
!> \param ichoi number of trial positions
!> \param iufrom the bead from which the new beads were grown
!> \param ntogrow number of new beads that have been grown
!> \param glist the list of new beads that have been grown; ntogrow entries
!> \param maxlen maximum possible distance of the newly-grown beads from iufrom
!*****************************************************************
  subroutine boltz(lnew,lfirst,ovrlap,i,icharge,imolty,ibox,ichoi,iufrom,ntogrow,glist,maxlen)
    use sim_particle,only:lnn
    use util_mp,only:mp_set_displs

    logical::lnew,ovrlap,lcmno(nmax),lfirst
    logical::lqimol,lqjmol,liji,lqchgi
    integer::ichoi,growjj,igrow,count,glist(numax),icharge,cnt,jcell(nmax)
    integer::i,imolty,ibox,ntogrow,itrial,ntii,j,jj,ntjj,ntij,iu,jmolty,iufrom,ii,k,nmole
    ! integer::NRtype
    real::rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij,rij,rijsq,maxlen,rcm,rcmsq,corr,rcutmax
    real::v(nEnergy),vwell,rcutsq,rcinsq
    integer::mmm

    ! RP added for MPI
    integer::rcounts(numprocs),displs(numprocs),my_start,my_end,blocksize,my_itrial
    real::my_vtry(nchmax),my_vtrintra(nchmax),my_vtrext(nchmax),my_vtrinter(nchmax),my_vtrelect(nchmax),my_vtrewald(nchmax),my_bfac(nchmax),my_vipswot(nchmax),my_vwellipswot(nchmax),my_vipswnt(nchmax),my_vwellipswnt(nchmax)
    logical::my_lovr(nchmax)
! ------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start BOLTZ in ',myid
#endif

    if ( lpbc ) call setpbc(ibox)

    ! determine the potential cutoffs
    rcutsq = rcut(ibox)*rcut(ibox)

    ! KM initialize variables
    do j=1,ichoi
       my_lovr(j) = .false.
       lovr(j) = .false.
       my_vtry(j) = 0.0E0_dp
       my_vtrintra(j) = 0.0E0_dp
       my_vtrelect(j) = 0.0E0_dp
       my_vtrext(j) = 0.0E0_dp
       my_vtrinter(j) = 0.0E0_dp
       my_vtrewald(j) = 0.0E0_dp
       my_bfac(j) = 0.0E0_dp
       my_vipswot(j) = 0.0E0_dp
       my_vwellipswot(j) = 0.0E0_dp
       my_vipswnt(j) = 0.0E0_dp
       my_vwellipswnt(j) = 0.0E0_dp
    end do

    if ( ldual ) then
       ! use rcutin for both types of interactions (except intra)
       rcinsq = rcutin*rcutin
    else
       ! compute the cutoffs squared for each interaction
       rcinsq =  rcutsq
       if ( lcutcm ) then
          ! not needed when ldual is true since will use rcutin then
          rcutmax = rcut(ibox)
       end if
    end if

    ! compute minimum cutoff squared
    rminsq = rmin * rmin

    lqimol = lelect(imolty)
    igrow = nugrow(imolty)

    if ( lcutcm .and. (.not. lfirst) ) then
       ! check previous bead (iufrom) COM for each molecule in the box
       if ( lnew ) then
          ! for trial chain ###
          rxui  = rxnew(iufrom)
          ryui  = rynew(iufrom)
          rzui  = rznew(iufrom)
       else
          ! for old chain ###
          rxui  = rxu(i,iufrom)
          ryui  = ryu(i,iufrom)
          rzui  = rzu(i,iufrom)
       end if

!> \todo to be MPI parallelized
       do j = 1,nchain
          lcmno(j) = .false.
          if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
             rxuij = rxui-xcm(j)
             ryuij = ryui-ycm(j)
             rzuij = rzui-zcm(j)
             ! minimum image the pseudo-ctrmas pair separation
             if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)

             rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
             rij  = sqrt(rijsq)

             if ( ldual ) then
                rcm = rcutin + rcmu(j) + maxlen
                rcmsq = rcm*rcm
             else
                rcm = rcutmax + rcmu(j) + maxlen
                rcmsq = rcm*rcm
             end if
             if (rijsq .gt. rcmsq ) lcmno(j) = .true.
          end if
       end do
    end if

    ! RP added for MPI
    blocksize = ichoi/numprocs
    rcounts = blocksize
    blocksize = ichoi - blocksize * numprocs
    if (blocksize.gt.0) rcounts(1:blocksize) = rcounts(1:blocksize) + 1
    call mp_set_displs(rcounts,displs,blocksize,numprocs)
    my_start = displs(myid+1) + 1
    my_end = my_start + rcounts(myid+1) - 1
    !if (ldebug) write(myid+100,*)'boltz: my_start=',my_start,'; my_end=',my_end,'; ichoi=',ichoi,'; rcounts=',rcounts,'; displs=',displs

    my_itrial = 0
    ! do itrial = 1, ichoi
    do itrial = my_start,my_end
       my_itrial = my_itrial + 1
       ! lovr(itrial) = .false.
       my_lovr(my_itrial) = .false.

       v(2) = 0.0E0_dp
       v(4) = 0.0E0_dp
       v(9) = 0.0E0_dp
       v(8) = 0.0E0_dp
       v(14) = 0.0E0_dp

       ! Only if L_Coul_CBMC is true, then compute electrostatic interactions/corrections
       if(L_Coul_CBMC.and.lewald.and..not.lideal(ibox)) then
          do count = 1,ntogrow
             ii = glist(count)
             ! This part does not change for fixed charge moves, but is
             ! used in the swap rosenbluth weight. - ewald self term
             v(14) = v(14) - qqu(icharge,ii)*qqu(icharge,ii)*calp(ibox)/sqrtpi
          end do
       end if

       ! no intramolecular interactions if this is the first bead
       if ( .not. lfirst ) then
! *****************************************
! INTRACHAIN BEAD-BEAD INTERACTIONS ***
! *****************************************
          ! cycle through molecule and check bead by bead
          do iu = 1, igrow
             ! see if iu exists in the new chain yet
             if (.not. lexist(iu)) cycle
             ntjj = ntype(imolty,iu)

             ! loop over all the grown beads
             do count = 1,ntogrow
                ii = glist(count)
                ! see if iu has nonbonded intramolecular interaction with ii
                if (linclu(imolty,ii,iu).or.(L_Coul_CBMC.and.(lqinclu(imolty,ii,iu).or.lewald))) then
                   ! assign bead type for ii,iu, and the cross term
                   ntii = ntype(imolty,ii)
                   ntij=type_2body(ntii,ntjj)

                   if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
                   ! determine distances
                   if ( lnew ) then
                      ! use new trial chain coordinates
                      rxuij  = rxnew(iu) - rxp(count,itrial)
                      ryuij  = rynew(iu) - ryp(count,itrial)
                      rzuij  = rznew(iu) - rzp(count,itrial)
                   else
                      ! use old chain coordinates
                      rxuij  = rxu(i,iu) - rxp(count,itrial)
                      ryuij  = ryu(i,iu) - ryp(count,itrial)
                      rzuij  = rzu(i,iu) - rzp(count,itrial)
                   end if
                   if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
                   rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   rij = sqrt(rijsq)
                end if
                if ( linclu(imolty,ii,iu) .or. lqinclu(imolty,ii,iu)) then
                   if ( linclu(imolty,ii,iu) ) then
                      if (rijsq.lt.rminsq.and..not.lexpand(imolty)) then
                         ! RP added for MPI
                         my_lovr(my_itrial) = .true.
                         ! write(io_output,*) 'intra overlap'
                         goto 19
                      else if (rijsq.lt.rcutsq .or. lijall) then
                         !> \bug skip intra if it is bending 1-3 and using a table?
                         if (L_bend_table) then
                            do mmm=1,inben(imolty,ii)
                               if (ijben3(imolty,ii,mmm).eq.iu) then
                                  v(4) = v(4) + lininter_bend(rij,itben(imolty,ii,mmm))
                                  goto 96
                               end if
                            end do
                         end if

                         v(4)=v(4)+U2(rij,rijsq,i,imolty,ii,ntii,i,imolty,iu,ntjj,ntij)
                      end if
                   end if

                   ! intramolecular charge interaction
                   ! compute velect (coulomb and ewald)
96                 if (L_Coul_CBMC.and.lqinclu(imolty,ii,iu).and.lqchg(ntii).and.lqchg(ntjj).and.rijsq.lt.rcinsq) then
                      ! boltz.f has problem to compute the electrostatic interactions
                      ! in a group-based way because the leader q might not be grown at
                      ! present, so it calculates electrostatic interaction not based on
                      ! group but on its own distance in SC, but should be corrected
                      ! later by calling energy subroutine.
                      if (L_elect_table) then
                         v(8) = v(8) + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(icharge,iu)*lininter_elect(rij,ntii,ntjj)
                      else if (lewald) then
                         ! compute real space term of vewald
                         v(8) = v(8) + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(icharge,iu)*erfunc(calp(ibox)*rij)/ rij
                         ! ewald sum correction term
                         corr = (1.0E0_dp - qscale2(imolty,ii,iu))*qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0E0_dp) /rij
                         v(14) = v(14) + corr
                      else
                         v(8) = v(8) + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(i,iu)/rij
                      end if
                   end if
                   ! end charge calculation
                   ! will only add correction if lqinclu is false.
                else if (L_Coul_CBMC.and.lewald) then
                   ! ewald sum correction term
                   corr = qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0E0_dp) /rij
                   v(14) = v(14) + corr
                end if
             end do
          end do
          !> \bug double Ewald correction?
          if (L_Coul_CBMC.and.lewald.and.ntogrow.gt.1.and..not.lideal(ibox)) then
             ! ewald sum correction term for interactions of the
             ! growing beads with each other
             ! this is 1-3, so don't need to consult lqinclu
             ! should change this since it corrects for all currently grown beads in this
             ! step, which could at somepoint be further than 1-3 apart!!! (say for rigrot...)
             do cnt = 1,ntogrow-1
                iu = glist(cnt)
                do count = cnt+1,ntogrow
                   ii = glist(count)
                   ! determine distances - use trial chain coordinates
                   rxuij  = rxp(cnt,itrial) - rxp(count,itrial)
                   ryuij  = ryp(cnt,itrial) - ryp(count,itrial)
                   rzuij  = rzp(cnt,itrial) - rzp(count,itrial)
                   !> \bug no mimage convertion?
                   rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   rij   = sqrt(rijsq)
                   ! ewald sum correction term
                   corr = qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0E0_dp)/rij
                   v(14) = v(14) + corr
                end do
             end do
          end if

          if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj .and. .not. lninesix .and..not.L_vdW_table.and..not.L_bend_table) then
             v(4) = 4.0E0_dp * v(4)
          end if
       end if

! grand-canonical: if ibox = 2 (ideal gas box) only intra-chain
! JLR 11-24-09 don't compute if lideal
       if (.not.(lgrand.and.ibox.eq.2).and..not.lideal(ibox)) then
! END JLR 11-24-09
          if (licell.and.(ibox.eq.boxlink)) then
! we only use count = 1, the rest should be taken care of
! with rintramax
             count = 1
             rxui = rxp(count,itrial)
             ryui = ryp(count,itrial)
             rzui = rzp(count,itrial)
             call get_cell_neighbors(rxui,ryui,rzui,ibox,jcell,nmole)
          else
             nmole = nchain
          end if

! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************
! loop over all chains except i
          do_nmole:do k = 1, nmole
             if (licell.and.(ibox.eq.boxlink)) then
                j = jcell(k)
             else
                j = k
             end if

             ! check for simulation box
             if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
                ! check neighbor list and COM table calculated above
                if ((lneigh.and.(.not.lnew).and.(.not.lnn(j,i))).or.((.not.lfirst).and.(lcmno(j).and.lcutcm))) cycle do_nmole

                jmolty = moltyp(j)
                lqjmol = lelect(jmolty)
                growjj = nugrow(jmolty)

                ! loop over all beads of molecule i grown this step
108             do count = 1,ntogrow
                   ! assign bead type for ii
                   ii = glist(count)
                   ntii = ntype(imolty,ii)
                   liji = lij(ntii)
                   lqchgi = lqchg(ntii)

                   ! assign positions to r*ui
                   rxui = rxp(count,itrial)
                   ryui = ryp(count,itrial)
                   rzui = rzp(count,itrial)

                   if ( lfirst .and. lcutcm ) then
                      ! check if ctrmas within rcmsq
                      rxuij = rxui-xcm(j)
                      ryuij = ryui-ycm(j)
                      rzuij = rzui-zcm(j)
                      ! minimum image the ctrmas pair separations
                      if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rij   = sqrt(rijsq)
                      ! determine cutoff
                      if ( ldual ) then
                         ! must be lfirst so no previous bead
                         rcm = rcutin + rcmu(j)
                      else
                         ! standard lcutcm cutoff
                         rcm = rcutmax + rcmu(j)
                      end if
                      ! check if interaction distance is greater than cutoff
                      rcmsq = rcm*rcm
                      if ( rijsq .gt. rcmsq ) cycle do_nmole
                   end if

                   ! loop over all beads jj of chain j
                   do jj = 1,nunit(jmolty)
                      ! check exclusion table
                      if ( lexclu(imolty,ii,jmolty,jj) ) cycle
                      ! start iswatch add-on ***
                      !> \todo is there a way to pull this out of the loops?
                      if (liswatch.and.j.eq.other.and.(.not.liswinc(jj,jmolty))) then
                         ! write(io_output,*) 'iSwatch-skipping:',jj
                         cycle
                      end if
                      ! end iswatch add-on ***

                      ntjj = ntype(jmolty,jj)
                      if ((.not.(liji.and.lij(ntjj))).and.(.not.(lqchgi.and.lqchg(ntjj))).and..not.L_vdW_table) cycle

                      ntij=type_2body(ntii,ntjj)

                      if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

                      rxuij = rxui - rxu(j,jj)
                      ryuij = ryui - ryu(j,jj)
                      rzuij = rzui - rzu(j,jj)
                      ! minimum image the pair separations ***
                      if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rij   = sqrt(rijsq)
                      ! compute vinter (eg. lennard-jones)
                      if (rijsq.lt.rminsq.and..not.(lexpand(imolty).or.lexpand(jmolty))) then
                         my_lovr(my_itrial) = .true.
                         ! write(io_output,*) 'j:',j,jj
                         ! write(io_output,*) 'rjsq:',rijsq,rminsq
                         goto 19
                      else if (rijsq.lt.rcinsq.or.lijall) then
                         v(2)=v(2)+U2(rij,rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                      end if

                      ! compute velect (coulomb and ewald)
                      if (L_Coul_CBMC.and.lqchg(ntii).and.lqchg(ntjj).and.rijsq.lt.rcinsq) then
                         ! boltz.f has problem to compute the electrostatic interactions
                         ! in a group-based way because the leader q might not be grown at
                         ! present, so it calculates electrostatic interaction not based on
                         ! group but on its own distance in SC, but should be corrected
                         ! later by calling energy subroutine.
                         if (L_elect_table) then
                            v(8) = v(8) + qqu(icharge,ii)*qqu(j,jj)*lininter_elect(rij,ntii,ntjj)
                         else if (lewald) then
                            ! compute real space term of velect
                            v(8) = v(8) + qqu(icharge,ii)*qqu(j,jj)*erfunc(calp(ibox)*rij)/rij
                         else
                            ! compute all electrostatic interactions
                            v(8) = v(8) + qqu(icharge,ii)*qqu(j,jj)/ rij
                         end if
                      end if
                   end do
                end do
             end if
          end do do_nmole

          if (.not.lsami.and..not.lexpsix.and..not.lmmff.and..not.lgenlj.and..not.lninesix.and..not.L_vdW_table) then
             v(2) = 4.0E0_dp * v(2)
          end if
       end if
! ################################################################

! **************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

! not for grand can. with ibox=2 !
! required for histogram reweighting to work for monolayer
! phase diagrams.
! not used for adsorption isotherms
       if (ibox .eq. 1) then
          if ((lelect_field.and.lqimol).or.((ibox.eq.1).and.(ljoe .or. lsami .or. lmuir .or. lexzeo .or. lgraphite .or. lslit))) then
             do count = 1,ntogrow
                ! assign bead type for ii
                ii = glist(count)
                ntii = ntype(imolty,ii)
                rxu(nchain+2,ii) = rxp(count,itrial)
                ryu(nchain+2,ii) = ryp(count,itrial)
                rzu(nchain+2,ii) = rzp(count,itrial)
                v(9)=v(9)+U_ext(ibox,nchain+2,ii,ntii)
             end do
          end if
       end if

! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
       vwell = 0.0E0_dp
       if (lwell(imolty).and.lmipsw) then
          rxui = xcm(i)
          ryui = ycm(i)
          rzui = zcm(i)
          do j = 1, nwell(imolty)*nunit(imolty)
             k = j - int(j/nunit(imolty))*nunit(imolty)
             if (k.eq.0) k = nunit(imolty)
             rxuij = rxui-rxwell(j,imolty)
             ryuij = ryui-rywell(j,imolty)
             rzuij = rzui-rzwell(j,imolty)
             call mimage(rxuij,ryuij,rzuij,ibox)
             rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
             rcm = rcut(ibox)+rcmu(i)+maxlen
             rcmsq = rcm*rcm
             if (rijsq.lt.rcmsq) then
                do count = 1, ntogrow
                   ii = glist(count)
                   if (awell(ii,k,imolty).lt.1.0E-6_dp) cycle
                   rxui = rxp(count,itrial)
                   ryui = ryp(count,itrial)
                   rzui = rzp(count,itrial)
                   rxuij = rxui-rxwell(j,imolty)
                   ryuij = ryui-rywell(j,imolty)
                   rzuij = rzui-rzwell(j,imolty)
                   call mimage(rxuij,ryuij,rzuij,ibox)
                   rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                   vwell = vwell-awell(ii,k,imolty)*exp(-bwell*rijsq)
                end do
             end if
          end do
       end if
! ----------------------------------------------------------------------------

! *********************************************
! CALCULATION OF TOTAL POTENTIAL ENERGY ***
! *********************************************
! write(23,*)
! write(23,*) 'wirting out total energy'
! write(23,*) 'Lnew', lnew
! if (NRtype.eq.1) then
! write(23,*) 'called from swap'
! else
! write(23,*) 'called from rigrot'
! end if

19     if ( my_lovr(my_itrial) ) then
          my_bfac(my_itrial) = 0.0E0_dp
       else
          if (.not.L_elect_table) then
             v(8) = v(8)*qqfact
             v(14) = v(14)*qqfact
          end if
          v(1) = v(2)+v(4)+v(9)+v(8)+v(14)
          if (.not.lnew) then
             my_vipswot(my_itrial) = v(1)
             my_vwellipswot(my_itrial) = vwell
          else
             my_vipswnt(my_itrial) = v(1)
             my_vwellipswnt(my_itrial) = vwell
          end if

          if (lstagea) then
             v(1) = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*v(1)
          else if (lstageb) then
             v(1) = etais*v(1)+lambdais*vwell
          else if (lstagec) then
             v(1) = (etais+(1.0E0_dp-etais)*lambdais)*v(1)+ (1.0E0_dp-lambdais)*vwell
          end if

          my_vtry(my_itrial) = v(1)
          my_vtrintra(my_itrial) = v(4)
          my_vtrext(my_itrial)   = v(9)
          my_vtrinter(my_itrial) = v(2)
          my_vtrelect(my_itrial) = v(8)
          my_vtrewald(my_itrial) = v(14)
          ! write(23,*) 'itrial' ,itrial
          ! write(23,*) vtr(1,itrial), vtr(4,itrial), vtr(9,itrial),
          !   &       vtr(2,itrial),vtr(8,itrial), vtr(14,itrial)

          if ((my_vtry(my_itrial)*beta).gt.(2.3E0_dp*softcut))then
             ! write(io_output,*) 'caught by softcut',vtr(1,itrial)*beta
             my_lovr(my_itrial) = .true.
             my_bfac(my_itrial) = 0.0E0_dp
          else if((my_vtry(my_itrial)*beta).lt.-2.303E0_dp*308)then
             ! write(io_output,*) '### warning: weight too big out of range'
             my_lovr(my_itrial) = .true.
             my_bfac(my_itrial) = 0.0E0_dp
          else
             my_bfac(my_itrial) = exp ( -(my_vtry(my_itrial)*beta) )
          end if
       end if
    end do

    call mp_allgather(my_vtry,vtr(1,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrintra,vtr(4,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrext,vtr(9,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrinter,vtr(2,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrelect,vtr(8,:),rcounts,displs,groupid)
    call mp_allgather(my_vtrewald,vtr(14,:),rcounts,displs,groupid)
    call mp_allgather(my_bfac,bfac,rcounts,displs,groupid)
    call mp_allgather(my_vipswot,vipswot,rcounts,displs,groupid)
    call mp_allgather(my_vwellipswot,vwellipswot,rcounts,displs,groupid)
    call mp_allgather(my_vipswnt,vipswnt,rcounts,displs,groupid)
    call mp_allgather(my_vwellipswnt,vwellipswnt,rcounts,displs,groupid)
    call mp_allgather(my_lovr,lovr,rcounts,displs,groupid)
    ovrlap = .true.
    if (ANY(.not.lovr(1:ichoi))) ovrlap=.false.
! ----------------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'end BOLTZ in ',myid
#endif
    return
  end subroutine boltz

! **********************************
!> \brief tail-corrections in energy ***
! **********************************
  function coru(imolty,jmolty,rho,ibox)

      real::coru,rci3,rho,epsilon2,sigma2
      real::rci1
      integer::imolty,jmolty,ii,jj, ntii, ntjj, ntij ,ibox

      coru = 0.0E0_dp
      do ii = 1, nunit(imolty)
         ntii = ntype(imolty,ii)
         do jj = 1, nunit(jmolty)
            ntjj = ntype(jmolty,jj)
            ntij = type_2body(ntii,ntjj)
            if (lexpsix) then
               coru = coru + rho*consu(ntij)
            else if (lmmff) then
               coru = coru + rho * epsimmff(ntij) * coru_cons(ntij) * sigimmff(ntij)**3.0E0_dp*twopi
            else if (lninesix) then
               coru = coru + 8.0E0_dp*onepi*rho*epsnx(ntij)* rzero(ntij)**3*(rzero(ntij)/rcut(ibox))**3* ((rzero(ntij)/rcut(ibox))**3/3.0E0_dp - 1.0E0_dp)
            else if (lgenlj) then
               rci3 = sig2ij(ntij)**(3.0E0_dp/2.0E0_dp) / rcut(ibox)**3
               rci1 = rci3 **(1.0E0_dp/3.0E0_dp)

               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigma_f(jmolty,jj))**2/4.0E0_dp
                  epsilon2 = sqrt(epsilon_f(imolty,ii) *epsilon_f(jmolty,jj))
               else if ( lexpand(imolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigi(ntjj))**2/4.0E0_dp
                  epsilon2 = sqrt(epsilon_f(imolty,ii)*epsi(ntjj))
               else if ( lexpand(jmolty) ) then
                  sigma2 = (sigma_f(jmolty,jj)+sigi(ntii))**2/4.0E0_dp
                  epsilon2 = sqrt(epsilon_f(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               end if
               coru = coru + 2.0E0_dp * onepi * epsilon2 * sigma2 ** (1.50E0_dp) * rho&
                * (  (( (2.0E0_dp**(4.0E0_dp*n1/n0))/(2.0E0_dp*n1-3.0E0_dp)) * rci1 **(2.0E0_dp*n1-3.0E0_dp) )&
                - ( (2.0E0_dp**((2.0E0_dp*n1/n0)+1.0E0_dp))/(n1-3.0E0_dp)) * rci1 **(n1-3.0E0_dp) )

            else
               rci3 = sig2ij(ntij)**(3.0E0_dp/2.0E0_dp) / rcut(ibox)**3
               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigma_f(jmolty,jj))**2/4.0E0_dp
                  epsilon2 = sqrt(epsilon_f(imolty,ii) *epsilon_f(jmolty,jj))
               else if ( lexpand(imolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigi(ntjj))**2/4.0E0_dp
                  epsilon2 = sqrt(epsilon_f(imolty,ii)*epsi(ntjj))
               else if ( lexpand(jmolty) ) then
                  sigma2 = (sigma_f(jmolty,jj)+sigi(ntii))**2/4.0E0_dp
                  epsilon2 = sqrt(epsilon_f(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               end if
               coru = coru +  8.0E0_dp * onepi * epsilon2 *  sigma2**(1.5E0_dp) *rho *  (rci3 * rci3 * rci3 / 9.0E0_dp - rci3 / 3.0E0_dp)
            end if
         end do
      end do
      return
  end function coru

!> \brief Read in force field parameters
!>
!> \todo Need to reorganize pairwise interactions so that more than one functional type can be used
  subroutine read_ff(io_ff,lmixlb,lmixjo)
    use util_search,only:initiateTable,addToTable
    use util_memory,only:reallocate
    use energy_intramolecular,only:read_energy_bonded

    integer,INTENT(IN)::io_ff
    logical,intent(in)::lmixlb,lmixjo
    integer,parameter::initial_size=20
    character(LEN=default_string_length)::line_in
    integer::jerr,i,j,ij,it,nmix
    real::sigmaTmp,epsilonTmp

    call initiateTable(atoms,initial_size)

    call init_tabulated_potential_pair()

    !> Looking for section ATOMS
    nntype=0
    rewind(io_ff)
    CYCLE_READ_ATOMS:DO
       call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) exit cycle_read_atoms

       if (UPPERCASE(line_in(1:5)).eq.'ATOMS') then
          allocate(atom_type(1:initial_size),sigi(1:initial_size),epsi(1:initial_size),qelect(1:initial_size),mass(1:initial_size),lij(1:initial_size),lqchg(1:initial_size),chemid(1:initial_size),stat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_pairwise: atoms allocation failed',jerr)
          sigi = 0.0E0_dp
          epsi = 0.0E0_dp
          qelect = 0.0E0_dp
          mass = 0.0E0_dp
          lij = .true.
          lqchg = .false.

          do
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section ATOMS',jerr)
             if (UPPERCASE(line_in(1:9)).eq.'END ATOMS') exit
             nntype=nntype+1
             read(line_in,*) i
             i=addToTable(atoms,i,expand=.true.)
             if (i.gt.ubound(atom_type,1)) then
                call reallocate(atom_type,1,2*ubound(atom_type,1))
                call reallocate(sigi,1,2*ubound(sigi,1))
                call reallocate(epsi,1,2*ubound(epsi,1))
                call reallocate(qelect,1,2*ubound(qelect,1))
                call reallocate(mass,1,2*ubound(mass,1))
                call reallocate(lij,1,2*ubound(lij,1))
                call reallocate(lqchg,1,2*ubound(lqchg,1))
                call reallocate(chemid,1,2*ubound(chemid,1))
             end if
             read(line_in,*) j,atom_type(i),sigi(i),epsi(i),qelect(i),mass(i),chemid(i)
             !read(line_in,'(I,4F,2A)') i,sigi(i),epsi(i),qelect(i),mass(i),chemid(i),chname(i)
             if (qelect(i).ne.0) then
                lqchg(i)=.true.
             else
                lqchg(i)=.false.
             end if
             if (sigi(i).eq.0 .and. epsi(i).eq.0) then
                lij(i)=.false.
             else
                lij(i)=.true.
             end if
          end do
          exit cycle_read_atoms
       end if
    END DO CYCLE_READ_ATOMS

    nmix=nntype*nntype

    allocate(lpl(1:nntype),xiq(1:nntype),jayself(1:nntype),nonbond_type(1:nmix),sig2ij(1:nmix),epsij(1:nmix),rminee(1:nmix),ecut(1:nmix),jayq(1:nmix),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_pairwise: nonbond allocation failed',jerr)

    xiq = 0.0E0_dp
    lpl = .false.
    sig2ij=0.0E0_dp
    epsij=0.0E0_dp
    jayq=0.0E0_dp

    ! Computation of un-like interactions
    ! convert input data to program units
    ! calculate square sigmas and epsilons for lj-energy subroutines
    do i = 1, nntype
       do j = 1, nntype
          ij = (i-1)*nntype + j
          if (lmixlb) then
             ! Lorentz-Berthelot rules --- sig_ij = 0.5 [ sig_i + sig_j ]
             sig2ij(ij) =(0.5E0_dp*(sigi(i)+sigi(j)))**2
             if (sigi(i).eq.0.0E0_dp.or.sigi(j).eq.0.0E0_dp) sig2ij(ij) = 0.0E0_dp
          else if (lmixjo) then
             ! Jorgensen mixing rules --- sig_ij = [ sig_i * sig_j ]^(1/2)
             sig2ij(ij) = sigi(i) * sigi(j)
          end if
          epsij(ij) = sqrt(epsi(i)*epsi(j))
       end do
    end do

    !> Looking for section NONBOND
    REWIND(io_ff)
    CYCLE_READ_NONBOND:DO
       call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) exit cycle_read_nonbond

       if (UPPERCASE(line_in(1:7)).eq.'NONBOND') then
          nmix=0
          do
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section NONBOND',jerr)
             if (UPPERCASE(line_in(1:11)).eq.'END NONBOND') exit
             nmix=nmix+1
             read(line_in,*) i,j,it,sigmaTmp,epsilonTmp
             ij=(i-1)*nntype+j
             nonbond_type(ij)=it
             nonbond_type((j-1)*nntype+i)=it
             sig2ij(ij)=sigmaTmp*sigmaTmp
             sig2ij((j-1)*nntype+i)=sig2ij(ij)
             epsij(ij)=epsilonTmp
             epsij((j-1)*nntype+i)=epsilonTmp
          end do
          exit cycle_read_nonbond
       end if
    END DO CYCLE_READ_NONBOND

    !> set up the strectching and bending constants
    call read_energy_bonded(io_ff)

    close(io_ff)
  end subroutine read_ff

  subroutine init_ff(io_input,lprint)
    use energy_external,only:init_energy_external
    use energy_sami,only:susami
    INTEGER,INTENT(IN)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::i,j,ij
    real::sr2,sr6,rzeronx(nxatom),epsilonnx(nxatom),rcheck,rs1,rs7,sr7

    do i = 1, nntype
       do j = 1, nntype
          ij = (i-1)*nntype + j
          if (lshift) then
             sr2 = sig2ij(ij)/(rcut(1)*rcut(1))
             sr6 = sr2 * sr2 * sr2
             ecut(ij)= sr6*(sr6-1.0E0_dp)*epsij(ij)
          end if
       end do
    end do

    !> read external potentials
    call init_energy_external(io_input,lprint)

    if ( lsami ) then
       call susami()
       rcheck = 2.5E0_dp * 3.527E0_dp
       if ( rcut(1) .ne. rcheck ) then
          if (lprint) write(io_output,*) 'WARNING ### rcut set to 2.5sigma for SAMI'
          rcut(1) = rcheck
       end if
    else if (lgaro) then
! KEA adding garofalini silica/water potential
! see J. Phys. Chem. 94 5351 (1990)
! form:
! U(2) = A exp(-rij/rho)+[zi zj erfc (rij/beta)]/rij + a/[1+exp(b/(rij-c))]
! U(3) = h3(rij rik thetajik) + h3(rjk rji thetakji) + h3(rki rkj thetaijk)
! h3(rij rik thetajik) = lambda exp[(gamma/(rij-rij*))+(gamma/(rik-rik*))]*
!             [cos(theta)-cos(theta*)]**2     for rij<rij* and rik<rik*
! 0 otherwise
       do i=1,6
          galpha(i) = 0.0E0_dp
          grho(i) = 0.0E0_dp
          gbeta(i) = 0.0E0_dp
          ecut(i) = 0.0E0_dp
          do j=1,3
             ga(i,j) = 0.0E0_dp
             gb(i,j) = 0.0E0_dp
             gc(i,j) = 0.0E0_dp
          end do
       end do
       do i=1,4
          glambda(i) = 0.0E0_dp
          do j=1,2
             grij(i,j) = 0.0E0_dp
             grijsq(i,j) = 0.0E0_dp
             ggamma(i,j) = 0.0E0_dp
          end do
          gtheta(i) = 0.0E0_dp
       end do

! Parameters (galpha,grho,gbeta,ga,gb,gc; lambda,grij,ggamma,gtheta)
! Si-Si
       galpha(1) = 13597175.7E0_dp
       grho(1) = 0.29E0_dp
       gbeta(1) = 2.29E0_dp
       lqchg(1) = .true.
       qelect(1) = 4.0E0_dp
       mass(1) = 28.09E0_dp
       ecut(1) = garofalini(rcut(1)*rcut(1),1,qelect(1),qelect(1),1,1)
       chemid(1) = 'Si '

! O-O
       galpha(2) = 5251972.5E0_dp
       grho(2) = 0.29E0_dp
       gbeta(2) = 2.34E0_dp
       lqchg(2) = .true.
       qelect(2) = -2.0E0_dp
       mass(2) = 16.00E0_dp
       ecut(2) = garofalini(rcut(1)*rcut(1),2,qelect(2),qelect(2) ,2,2)
       chemid(2) = 'O  '

! H-H
       galpha(3) = 246299.4E0_dp
       grho(3) = 0.35E0_dp
       gbeta(3) = 2.1E0_dp
       ga(3,1) = -38243.8E0_dp
       gb(3,1) = 6.0E0_dp
       gc(3,1) = 1.51E0_dp
       ga(3,2) = 2515.9E0_dp
       gb(3,2) = 2.0E0_dp
       gc(3,2) = 2.42E0_dp
       lqchg(3) = .true.
       qelect(3) = 1.0E0_dp
       mass(3) = 1.0078E0_dp
       ecut(3) = garofalini(rcut(1)*rcut(1),3,qelect(3),qelect(3) ,3,3)
       chemid(3) = 'H  '

! Si-O
       galpha(4) =  21457024.2E0_dp
       grho(4) = 0.29E0_dp
       gbeta(4) = 2.34E0_dp
       ecut(4) = garofalini(rcut(1)*rcut(1),4,qelect(1),qelect(2) ,1,2)

! Si-H
       galpha(5) = 499842.9E0_dp
       grho(5) = 0.29E0_dp
       gbeta(5) = 2.31E0_dp
       ga(5,1) = -33715.5E0_dp
       gb(5,1) = 6.0E0_dp
       gc(5,1) = 2.2E0_dp
       ecut(5) = garofalini(rcut(1)*rcut(1),5,qelect(1),qelect(3) ,1,3)

! 0-H
       galpha(6) = 2886049.4E0_dp
       grho(6) = 0.29E0_dp
       gbeta(6) = 2.26E0_dp
       ga(6,1) = -15096.7E0_dp
       gb(6,1) = 15.0E0_dp
       gc(6,1) = 1.05E0_dp
       ga(6,2) = 55353.6E0_dp
       gb(6,2) = 3.2E0_dp
       gc(6,2) = 1.50E0_dp
       ga(6,3) = -6038.7E0_dp
       gb(6,3) = 5.0E0_dp
       gc(6,3) = 2.0E0_dp
       ecut(6) = garofalini(rcut(1)*rcut(1),6,qelect(2),qelect(3) ,2,3)

! Si-O-Si
       glambda(1) = 21732.3E0_dp
       ggamma(1,1) = 2.0E0_dp
       grij(1,1) = 2.6E0_dp
       grijsq(1,1) = grij(1,1)*grij(1,1)
       gtheta(1) = cos(109.5E0_dp*degrad)

! O-Si-O
       glambda(2) = 1376379.0E0_dp
       ggamma(2,1) = 2.8E0_dp
       grij(2,1) = 3.0E0_dp
       grijsq(2,1) = grij(2,1)*grij(2,1)
       gtheta(2) = cos(109.5E0_dp*degrad)

! H-O-H
       glambda(3) = 2535435.0E0_dp
       ggamma(3,1) = 1.3E0_dp
       grij(3,1) = 1.6E0_dp
       grijsq(3,1) = grij(3,1)*grij(3,1)
       gtheta(3) = cos(104.5E0_dp*degrad)

! Si-O-H
       glambda(4) = 362205.0E0_dp
       ggamma(4,1) = 2.0E0_dp
       ggamma(4,2) = 1.2E0_dp
       grij(4,1) = grij(1,1)
       grij(4,2) = 1.5E0_dp
       grijsq(4,1) = grij(4,1)*grij(4,1)
       grijsq(4,2) = grij(4,2)*grij(4,2)
       gtheta(4) = cos(109.5E0_dp*degrad)
    else if(lexpsix) then
! Explicit atom carbon Williams Exp-6 potential
! J.Chem.Phys 47 11 4680 (1967) paramter set IV
! note that the combination of C--H has to be the average of C
! and H the way it is set up in the energy subroutines
! natom is set in expsix.inc
! U(r) = A*r^(-6) + B*exp[C*r]
       do i = 1,natom
          aexsix(i) = 0.0E0_dp
          bexsix(i) = 0.0E0_dp
          cexsix(i) = 0.0E0_dp
          qelect(i) = 0.0E0_dp
          xiq(i) = 0.0E0_dp
          lqchg(i) = .false.
          lij(i) = .true.
       end do

! O--- nonbonded interaction
! BKS [O]
       qelect(1) = -1.2
       lqchg(1) = .true.
       aexsix(1) = -175.0000E0_dp*1.602177E-19_dp/1.3806503E-23_dp
       bexsix(1) = 1388.7730E0_dp*1.602177E-19_dp/1.3806503E-23_dp
       cexsix(1) = -2.76000E0_dp
       mass(1) = 15.9994

! O---Si nonbonded interaction
       aexsix(2) = -133.5381E0_dp*1.602177E-19_dp/1.3806503E-23_dp
       bexsix(2) = 18003.7572E0_dp*1.602177E-19_dp/1.3806503E-23_dp
       cexsix(2) = -4.87318E0_dp
       lqchg(2) = .true.

! Si---Si nonbonded interaction
! BKS [Si]
       qelect(3) = 2.4
       lqchg(3) = .true.
       aexsix(3) = 0.
       bexsix(3) = 0.
       cexsix(3) = 0.
       mass(3) = 28.0899

       if (lshift) then
          do i=1,natom
             sexsix(i) = aexsix(i)*(rcut(1)**(-6)) + bexsix(i)*Exp(cexsix(i)*rcut(1))
          end do
       else
          do i=1,natom
             consp(i) = (2.0E0_dp/3.0E0_dp)*onepi*(2.0E0_dp*aexsix(i)/(rcut(1) *rcut(1)*rcut(1))&
              +bexsix(i)*exp(cexsix(i)*rcut(1)) *(-6.0E0_dp/(cexsix(i)*cexsix(i)*cexsix(i))+6.0E0_dp *rcut(1)&
              /(cexsix(i)*cexsix(i))-3.0E0_dp*rcut(1)* rcut(1)/ cexsix(i)+rcut(1)*rcut(1)*rcut(1)))
             consu(i) = 2.0E0_dp*onepi*(aexsix(i)/(3.0E0_dp*rcut(1)*rcut(1)* rcut(1)) +(-rcut(1)*rcut(1)&
              +2.0E0_dp*rcut(1)/cexsix(i)-2.0E0_dp/ (cexsix(i)* cexsix(i)))*bexsix(i)*exp(cexsix(i)*rcut(1))/ cexsix(i))
          end do
! write(11,*) 'consp(i)',consp
! write(11,*) 'consu(i)',consu
       end if
    else if (lmmff) then
! Merk Molecular Force Field (MMFF)
! J. Am. Chem. Soc. 1992 Vol.114 P7827-7843 (Thomas A. Halgren)
! natomtyp is set in mmff.inc
! U(r) = epsi*(1.07/(rs+0.07))^7 * (1.12/(rs^7+0.12)-2)
! rs = r / sigimmff

! C---C nonbonded interaction ***
       alphammff(1) = 1.050E0_dp
       nmmff(1) = 2.490E0_dp
       ammff(1) = 3.890E0_dp
       gmmff(1) = 1.282E0_dp
       sigimmff(1) = ammff(1)*sqrt(sqrt(alphammff(1)))
       epsimmff(1) = 45582.6E0_dp*gmmff(1)*gmmff(1)*sqrt(nmmff(1)) /(ammff(1)**6.0E0_dp)
! sigimmff(1) = 1.126763255691509E0_dp
! epsimmff(1) = 0.9994354715851470E0_dp
       sigisq(1) = sigimmff(1)*sigimmff(1)
       mass(1) = 12.011E0_dp
! mass(1) = 1.0E0_dp

! H---H nonbonded interaction ***
       alphammff(3) = 0.250E0_dp
       nmmff(3) = 0.800E0_dp
       ammff(3) = 4.200E0_dp
       gmmff(3) = 1.209E0_dp
       sigimmff(3) = ammff(3)*sqrt(sqrt(alphammff(3)))
       epsimmff(3) = 45582.6E0_dp*gmmff(3)*gmmff(3)*sqrt(nmmff(3)) /(ammff(3)**6.0E0_dp)
       sigisq(3) = sigimmff(3)*sigimmff(3)
       mass(3) = 1.0078E0_dp

! C---H nonbonded interaction by using cominbination rule ***
       sigimmff(2) = 0.5E0_dp*(sigimmff(1)+sigimmff(3))*(1.0E0_dp +  0.2E0_dp*(1.0E0_dp-exp(-12.E0_dp*(((sigimmff(1)-sigimmff(3))/ (sigimmff(3)+sigimmff(1)))**2.0E0_dp))))
       epsimmff(2) = 91165.1E0_dp*gmmff(1)*gmmff(3)*alphammff(1)* alphammff(3)/((sqrt(alphammff(1)/nmmff(1))+ sqrt(alphammff(3)/nmmff(3)))*(sigimmff(2)**6.0E0_dp))
       sigisq(2) = sigimmff(2)*sigimmff(2)

       if (lshift) then
          do i=1,natomtyp
             rs1 = rcut(1)/sigimmff(i)
             rs7 = rs1**7.0E0_dp
             sr7 = (1.07E0_dp/(rs1+0.07E0_dp))**7.0E0_dp
             smmff(i) = epsimmff(i)*sr7*(1.12E0_dp/(rs7+0.12)-2)
          end do
       else
          do i = 1,natomtyp
             smmff(i) = 0.0E0_dp
          end do
          coru_cons(1) = -2.4837937263569310E-02_dp
          ! coru_cons(1) = -5.8244592746534724E-03_dp
          coru_cons(2) = -1.7583010189381791E-02_dp
          coru_cons(3) = -8.3770412792126582E-03_dp
          corp_cons(1) = 0.1696349613545569E0_dp
          ! corp_cons(1) = 4.0098456560842058E-02_dp
          corp_cons(2) = 0.1203650950025348E0_dp
          corp_cons(3) = 5.7576802340310304E-02_dp
       end if
    else if (lninesix) then
! special potential for all-atom formic acid from llnl 4/6/04 jms
! sp2 carbon site H-[C](=O)-OH
       rzeronx(1) = 4.1834161E0_dp
       epsilonnx(1) = 45.224E0_dp
       qelect(1) = 0.44469E0_dp
       lqchg(1) = .true.
       mass(1) = 12.011E0_dp
       chemid(1) = 'C'

! hydroxyl oxygen site H-C(=O)-[O]H
       rzeronx(2) = 3.5694293136E0_dp
       epsilonnx(2) = 47.151E0_dp
       qelect(2) = -0.55296E0_dp
       lqchg(2) = .true.
       mass(2) = 15.9996E0_dp
       chemid(2) = 'OH'

! carbonyl oxygen site H-C[=O]-OH
       rzeronx(3) = 3.0014635172E0_dp
       epsilonnx(3) = 146.008E0_dp
       qelect(3) = -0.43236E0_dp
       lqchg(3) = .true.
       mass(3) = 15.9996E0_dp
       chemid(3) = 'O='

! hydrogen site [H]-C(=O)-OH
       rzeronx(4) = 0.8979696387E0_dp
       epsilonnx(4) = 2.4054E0_dp
       qelect(4) = 0.10732E0_dp
       lqchg(4) = .true.
       mass(4) = 1.00794E0_dp
       chemid(4) = 'HC'

! acidic hydrogen site H-C(=O)-O[H]
       rzeronx(5) = 1.115727276E0_dp
       epsilonnx(5) = 12.027E0_dp
       qelect(5) = 0.43331E0_dp
       lqchg(5) = .true.
       mass(5) = 1.00794E0_dp
       chemid(5) = 'HO'

! calculate all site-site parameters via Lorentz-Berthelot rules
       do i = 1,nxatom
          do j = 1,nxatom
             ij = (i-1)*nxatom + j
             rzero(ij) = 0.5E0_dp*(rzeronx(i) + rzeronx(j))
             epsnx(ij) = sqrt(epsilonnx(i)*epsilonnx(j))
             if (lshift) then
                shiftnsix(ij) = 4.0E0_dp*epsnx(ij)*(rzero(ij) /rcut(1))**6 * (2.0E0_dp*(rzero(ij)/rcut(1))**3 - 3.0E0_dp)
             end if
          end do
       end do
    end if

! ! --- TraPPE-UA? Methane [CH4] sp3 charged with polarizability
! sigi(28) = 3.73E0_dp
! epsi(28) = 148.0E0_dp
! ! is this correct?
! mass(28) = 16.043E0_dp
! qelect(28) = -0.572E0_dp
! lqchg(28) = .true.
! jayself(28) = 0.5E0_dp*117403E0_dp
! xiq(28) = 9449.3E0_dp
! chname(28) = 'Tr C CH4 chg pol '
! chemid(28)  = 'C  '

! ! --- Methane hydrogen charged with polarizibility
! sigi(29) = 0.0E0_dp
! epsi(29) = 0.0E0_dp
! mass(29) = 1.0078E0_dp
! qelect(29) = 0.143E0_dp
! lqchg(29) = .true.
! jayself(29) = 0.5E0_dp*177700E0_dp
! xiq(29) = 0.0E0_dp
! lij(29) = .false.
! chname(29) = 'Tr H CH4 chg pol '
! chemid(29)  = 'H  '

! ! --- SPC-FQ oxygen [O]   S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(109) = 3.176
! epsi(109) = 148.0E0_dp
! mass(109) = 15.999E0_dp
! qelect(109) = -0.672123708
! lqchg(109) = .true.
! xiq(109) = 36899.0E0_dp
! jayself(109) = (0.5E0_dp)*(503.2E0_dp)*(367.0E0_dp)
! chname(109) = 'SPC-FQ O water   '
! chemid(109)  = '0  '

! ! --- SPC-FQ hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(110) = 0.0E0_dp
! epsi(110) = 0.0E0_dp
! mass(110) = 1.0079E0_dp
! qelect(110) = 0.336061854
! lij(110) = .false.
! lqchg(110) = .true.
! xiq(110) = 0.0E0_dp
! jayself(110) = (0.5E0_dp)*(503.2E0_dp)*(392.2E0_dp)
! chname(110) = 'SPC-FQ H water   '
! chemid(110)  = 'H  '

! ! --- TIP4P-FQ Oxygen [O] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(111) = 3.159E0_dp
! epsi(111) = 144.1E0_dp
! !      epsi(111) = 105.0E0_dp
! mass(111) = 15.999E0_dp
! chname(111) = 'TIP4P-FQ O water '
! chemid(111)  = 'O  '

! ! --- TIP4P-FQ Hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(112) = 0.0E0_dp
! epsi(112) = 0.0E0_dp
! mass(112) = 1.0079E0_dp
! qelect(112) = 0.35E0_dp
! lij(112) = .false.
! lqchg(112) = .true.
! xiq(112) = 0.0E0_dp
! jayself(112) = (0.5E0_dp)*(503.2E0_dp)*(353.0E0_dp)
! chname(112) = 'TIP4P-FQ H water '
! chemid(112)  = 'H  '

! ! --- TIP4P-FQ Charge [Q] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(113) = 0.0E0_dp
! epsi(113) = 0.0E0_dp
! mass(113) = 0.0E0_dp
! qelect(113) = -0.70E0_dp
! lij(113) = .false.
! lqchg(113) = .true.
! xiq(113) = 34464.0E0_dp
! jayself(113) = (0.5E0_dp)*(503.2E0_dp)*(371.6E0_dp)
! chname(113) = 'TIP4P-FQ M water '
! chemid(113)  = 'M  '

! ! --- TraPPE carbon dioxide carbon in [C]O2-fq (jpotoff 2/15/00)
! sigi(131) = 2.80E0_dp
! epsi(131) = 28.5E0_dp
! mass(131) = 12.011E0_dp
! qelect(131) = 0.6512E0_dp
! lqchg(131) = .true.
! !      xiq(131) = (503.2E0_dp)*123.2E0_dp
! xiq(131) = 0.0E0_dp
! jayself(131) = (0.5E0_dp)*(503.2E0_dp)*(233.5E0_dp)
! chname(131) = 'Tr-FQ C in CO2   '
! chemid(131)  = 'C  '

! ! --- TraPPE carbon dioxide oxygen in C[O]2-fq (jpotoff 2/15/00)
! sigi(132) = 3.06E0_dp
! epsi(132) = 80.5E0_dp
! mass(132) = 15.999E0_dp
! qelect(132) = -0.3256E0_dp
! lqchg(132) = .true.
! !      xiq(132) = (503.2E0_dp)*201.56E0_dp
! xiq(132) = 39430.75E0_dp
! jayself(132) = (0.5E0_dp)*(503.2E0_dp)*(308.17E0_dp)
! chname(132) = 'Tr-FQ O in CO2   '
! chemid(132)  = 'O  '

! ! - CO2-FQ Carbon-Oxygen cross term (JCO)
! i = 131
! j = 132
! djay = (503.2E0_dp)*(133.905E0_dp)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ij) = djay
! jayq(ji) = djay

! ! - CO2-FQ Oxygen-Oxygen cross term (JOO)
! i = 132
! j = 132
! djay = (503.2E0_dp)*(1.09E0_dp)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- SPC-FQ water Oxygen-Hydrogen cross term
! i = 109
! j = 110
! djay = (503.2E0_dp)*(276.0E0_dp)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ij) = djay
! jayq(ji) = djay

! ! --- SPC-FQ water Hydrogen-Hydrogen cross term
! i = 110
! j = 110
! djay = (503.2E0_dp)*(196.0E0_dp)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- TIP4P water Charge-Hydrogen cross term
! i = 112
! j = 113
! djay = (503.2E0_dp)*(286.4E0_dp)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ji) = djay
! jayq(ij) = djay

! ! --- TIP4P water Hydrogen-Hydrogen cross term
! i = 112
! j = 112
! djay = (503.2E0_dp)*(203.6E0_dp)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- Methane C-H cross term
! i = 28
! j = 29
! ij = (i-1)*nntype + j
! ji = (j-i)*nntype + i
! jayq(ji) = 114855.0E0_dp
! jayq(ij) = 114855.0E0_dp

! ! --- Methane H-H cross term
! i = 29
! j = 29
! ij = (i-1)*nntype + j
! jayq(ij) = 112537.0E0_dp
  end subroutine init_ff

  function lininter_vdW(r,typi,typj) result(tabulated_vdW)
    use util_math,only:polint
    use util_search,only:indexOf,LOCATE
    real::tabulated_vdW
    real,intent(in)::r
    integer,intent(in)::typi,typj
    integer::low,high

    low=locate(rvdW(:,typi,typj),vdWsplits(typi,typj),r,2)
    high=low+1
    if (rvdW(low,typi,typj).gt.r.or.rvdW(high,typi,typj).lt.r) then
       write(io_output,*) 'problem in lininter_vdW!'
       write(io_output,*) 'r', r, ' typi', typi, ' typj ', typj
       write(io_output,*) 'low ', low, rvdW(low, typi, typj)
       write(io_output,*) 'high ', high, rvdW(high, typi, typj)
       write(io_output,*)
    end if
    call polint(rvdW(low:high,typi,typj),tabvdW(low:high,typi,typj),2,r,tabulated_vdW)
    return
  end function lininter_vdW

  function lininter_elect(r,typi,typj) result(tabulated_elect)
    use util_math,only:polint
    use util_search,only:indexOf,LOCATE
    real::tabulated_elect
    real,intent(in)::r
    integer,intent(in)::typi,typj
    integer::low,high

    low=locate(relect(:,typi,typj),electsplits(typi,typj),r,2)
    high=low+1
    if (relect(low,typi,typj).gt.r.or.relect(high,typi,typj).lt.r) then
       write(io_output,*) 'problem in lininter_elect!'
       write(io_output,*) 'r', r, ' typi', typi, ' typj ', typj
       write(io_output,*) 'low ', low, relect(low, typi, typj)
       write(io_output,*) 'high ', high, relect(high, typi, typj)
       write(io_output,*)
    end if
    call polint(relect(low:high,typi,typj),tabelect(low:high,typi,typj),2,r,tabulated_elect)
    return
  end function lininter_elect

!> \brief calculates the energy using the exp-6 potential
!> parameters defined in suijtab.f
!> \author M.G. Martin
  function exsix(rijsq,ntij)
      real::rijsq,rij,exsix
      integer::ntij

      rij=sqrt(rijsq)
      exsix = aexsix(ntij)/(rijsq*rijsq*rijsq) + bexsix(ntij)*exp(cexsix(ntij)*rij)
      if (lshift) exsix = exsix-sexsix(ntij)
      return
  end function exsix

!> \brief calculates the energy of the 9-6 potential
!> parameters defined in suijtab.f
!> \author JMS
  function ninesix(rijsq,ntij)
      real::rijsq,rij,ror,ninesix
      integer::ntij

      rij=sqrt(rijsq)
      ror = rzero(ntij)/rij
      ninesix = 4.0E0_dp*epsnx(ntij)*ror**6*(2.0E0_dp*ror**3 - 3.0E0_dp)
      if (lshift) then
         ninesix = ninesix - shiftnsix(ntij)
      end if

      return
  end function ninesix

!> \brief calculates the energy of the Generalized Lennard Jones
!> potential
!>
!> parameters defined  poten.inc \n
!> Also you should make sure that the rcut you choose
!> is .gt. sigma*(2**(2/n0)) (because currently tail corrections
!> are evaluated with this assumption which is rarely incorrect.)
  function genlj (rijsq,sr2,epsilon2)
      real::rijsq,rij,srij,sr2,epsilon2,genlj

      rij=sqrt(rijsq)
      srij=sqrt(sr2)

      if (rij.le.rij*srij*2.0_dp**(2.0_dp/n0)) then
         genlj = 4.0E0_dp*epsilon2*(srij**n0-srij**(n0/2.0E0_dp))
      else
         genlj = epsilon2*(2.0E0_dp**(4.0E0_dp*n1/n0)*srij**(2.0_dp*n1)-2.0_dp**(2.0_dp*n1/n0+1.0E0_dp)*srij**n1)
      end if

! In reduced units
! xij=(1.0E0_dp/sqrt(sr2))
!
! if ( (xij) .le. 2.0E0_dp**(2.0E0_dp/n0) ) then
! genlj = 4.0E0_dp * (((1/xij)**n0)-((1/xij)**(n0/2.0E0_dp)))
! else
! genlj =(( (2.0E0_dp**((4.0E0_dp*n1/n0)))*(1/xij)**(2.0E0_dp*n1))-
!     & (2.0E0_dp**((2.0E0_dp*(n1/n0))+1.0E0_dp)*((1/xij)**(n1))))
! end if
! v(2)=v(2)+e

! not usually used in MC
! if (lshift) then
! genlj = genlj - shiftgenlj()
! end if

      return
  end function genlj

!> \brief calculates energy for a polymeric surfactant bead.
  function ljpsur ( rijsq, ntij )

      real::ljpsur, rijsq, sr, sr6
      integer::ntij

! --------------------------------------------------------------------
! AT PRESENT: all sigma = 1.0
! epsilon   = 1.0
! --------------------------------------------------------------------

      sr = 1.0E0_dp / rijsq
      sr6 = sr**3

      if ( ntij .eq. 1 ) then
! nonpolar-nonpolar interaction ( LJ interaction )
         ljpsur = sr6*sr6 - sr6
      else
! polar-polar or polar-nonpolar interaction ( repulsive LJ interaction )
         if ( rijsq .le. 1.259921E0_dp ) then
            ljpsur = sr6*sr6 - sr6 + 0.25E0_dp
         else
            ljpsur = 0.0E0_dp
         end if
      end if

      return
  end function ljpsur

!> \brief calculates the energy using the Buffered 14-7 potential
!> parameters defined in suijtab.f
!> \author Bin Chen
  function mmff(rijsq,ntij)
      real::rijsq,mmff,rs2,rs1,sr1,sr7,rs7
      integer::ntij

        rs2 = rijsq / (sigisq(ntij))
        rs1 = sqrt(rs2)
        rs7 = rs1*rs2*rs2*rs2
        sr1 = 1.07E0_dp/(rs1+0.07E0_dp)
        sr7 = sr1**7.0E0_dp
        mmff = sr7*(1.12E0_dp/(rs7+0.12E0_dp) - 2.0E0_dp) *epsimmff(ntij)

      if (lshift) mmff = mmff-smmff(ntij)
      return
  end function mmff

!DEC$ ATTRIBUTES FORCEINLINE :: type_2body
  function type_2body(ntii,ntjj)
    integer::type_2body
    integer,intent(in)::ntii,ntjj

    if (lexpsix.or.lmmff) then
       type_2body = (ntii+ntjj)/2
    else if (lninesix) then
       type_2body = (ntii-1)*nxatom + ntjj
! KEA garofalini
    else if (lgaro) then
       if (ntii.eq.ntjj) then
          type_2body = ntii
       else
          type_2body = ntii+ntjj+1
       end if
    else
       type_2body = (ntii-1)*nntype + ntjj
    end if

  end function type_2body

!DEC$ ATTRIBUTES FORCEINLINE :: U2
  function U2(rij,rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
    real::U2
    real,intent(in)::rij,rijsq
    integer,intent(in)::i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij

    real::epsilon2,sigma2,sr2,sr6,qave

    U2=0.0E0_dp
    if (L_vdW_table) then
       U2=lininter_vdW(rij,ntii,ntjj)
    else if (lsami) then
       U2=ljsami(rijsq,ntij)
    else if (lexpsix) then
       U2=exsix(rijsq,ntij)
    else if (lmmff) then
       U2=mmff(rijsq,ntij)
    else if (lninesix) then
       U2=ninesix(rijsq,ntij)
    else if (lgenlj) then
       sr2 = sig2ij(ntij) / rijsq
       epsilon2=epsij(ntij)
       U2=genlj(rijsq,sr2,epsilon2)
    else if ( lmuir ) then
       U2=ljmuir(rijsq,ntij)
    else if ( lpsurf ) then
       U2=ljpsur(rijsq,ntij)
    else if (lgaro) then
       ! KEA garofalini potential; do NOT work with CBMC/boltz
       U2=garofalini(rijsq,ntij,qqu(i,ii),qqu(j,jj),i,j)
       if(lshift) then
          U2=U2-ecut(ntij)
       end if
    else if ((lij(ntii).and.lij(ntjj)).or.(i.eq.j)) then
       if (lexpand(imolty).and.lexpand(jmolty)) then
          sigma2=(sigma_f(imolty,ii)+sigma_f(jmolty,jj))/2.0E0_dp
          sr2 = sigma2*sigma2/rijsq
          epsilon2=sqrt(epsilon_f(imolty,ii)*epsilon_f(jmolty,jj))
       else if (lexpand(imolty)) then
          sigma2=(sigma_f(imolty,ii)+ sigi(ntjj))/2.0E0_dp
          sr2 = sigma2*sigma2/rijsq
          epsilon2=sqrt(epsilon_f(imolty,ii)*epsi(ntjj))
       else if (lexpand(jmolty)) then
          sigma2=(sigma_f(jmolty,jj)+ sigi(ntii))/2.0E0_dp
          sr2 = sigma2*sigma2/rijsq
          epsilon2=sqrt(epsi(ntii)*epsilon_f(jmolty,jj))
       else
          sr2 = sig2ij(ntij) / rijsq
          epsilon2=epsij(ntij)
       end if

       if (lfepsi.and.i.ne.j) then ! NOT for intramolecular interactions
          sr6 = rijsq*rijsq*rijsq
          if ((.not.lqchg(ntii)).and.(.not.lqchg(ntjj))) then
             if (nunit(imolty).eq.4) then
                !> \bug TIP-4P structure (temperary use?)
                qave=(qqu(i,4)+qqu(j,4))/2.0E0_dp
             else
                qave=(qqu(i,4)+qqu(i,5)+qqu(j,4)+qqu(j,5))*0.85E0_dp
             end if
          else
             qave=(qqu(i,ii)+qqu(j,jj))/2.0E0_dp
          end if
          U2=((aslope*(qave-a0)*(qave-a0)+ashift)/sr6-(bslope*(qave- b0)*(qave-b0)+bshift))/sr6*epsilon2
       else
          sr6 = sr2 * sr2 * sr2
          U2=sr6*(sr6-1.0E0_dp)*epsilon2
          if (lshift) then
             U2=U2-ecut(ntij)
          end if
          if (i.eq.j) then
             ! intramolecular interactions
             U2=U2*ljscale(imolty,ii,jj)
             if (lainclu(imolty,ii,jj)) then
                ! OH 1-5 interaction
               U2=U2+0.25E0_dp*a15(a15type(imolty,ii,jj))/((rijsq**2)*(rijsq**2)*(rijsq**2))
            end if
          end if
       end if
    end if
  end function U2

!> \brief Direct-space charge interactions
!DEC$ ATTRIBUTES FORCEINLINE :: Q2
  function Q2(rij,rijsq,rcutsq,i,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo)
    real::Q2
    real,intent(in)::rij,rijsq,rcutsq,calpi
    integer,intent(in)::i,imolty,ii,ntii,j,jmolty,jj,ntjj
    logical,intent(in)::lqchgi
    logical,intent(inout)::lcoulo(numax,numax)

    integer::iii,jjj

    Q2=0.0E0_dp
    !kea - skip for garofalini; included in vinter
    if(.not.lgaro.and.lqchgi.and.lqchg(ntjj)) then
       if (.not.lewald) then
          if (.not.lchgall) then
             ! All-Atom charges (charge-group look-up table)
             iii = leaderq(imolty,ii)
             jjj = leaderq(jmolty,jj)
             if (iii.eq.ii .and. jjj.eq.jj)then
                ! set up the charge-interaction table
                if ( rijsq .lt. rcutsq ) then
                   lcoulo(iii,jjj) = .true.
                else
                   lcoulo(iii,jjj) = .false.
                end if
             end if
             ! set up table for neighboring groups- make sure they interact when
             ! leaderqs are only 2 bonds apart. For intramolecular charge interactions
             if (i.eq.j.and..not.lqinclu(imolty,iii,jjj)) then
                lcoulo(iii,jjj)  = .true.
             end if
          end if
          if (lchgall.or.lcoulo(iii,jjj) ) then
             if (L_elect_table) then
                Q2=qqu(i,ii)*qqu(j,jj)*lininter_elect(rij,ntii,ntjj)
             else
                Q2=qqu(i,ii)*qqu(j,jj)/rij
             end if
          end if
       else if (lchgall.or.rijsq.lt.rcutsq) then
          Q2=qqu(i,ii)*qqu(j,jj)*erfunc(calpi*rij)/rij
       end if
    end if
  end function Q2

!> \brief Read in tabulated potential for nonbonded pair interactions (vdW and elect) and set up linear interpolation
  subroutine read_tabulated_potential_pair(file_tab,ntab,r,tab,splits,lists)
    use util_search,only:addToTable
    use util_memory,only:reallocate
    use util_files,only:get_iounit
    character(len=*),intent(in)::file_tab
    integer,intent(out)::ntab
    integer,allocatable,intent(inout)::splits(:,:)
    real,allocatable,intent(inout)::r(:,:,:),tab(:,:,:)
    type(LookupTable),intent(inout)::lists

    integer::io_tab,mmm,ii,jj,i,jerr

    io_tab=get_iounit()
    open(unit=io_tab,access='sequential',action='read',file=file_tab,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open tabulated potential file: '//file_tab,myid+1)
    end if

    read(io_tab,*) ntab
    do mmm=1,ntab
       ! ii and jj are bead types
       read(io_tab,*) ii, jj
       ii=addToTable(lists,ii,expand=.true.)
       jj=addToTable(lists,jj,expand=.true.)
       if (ii.gt.size(splits,1).or.jj.gt.size(splits,1)) then
          call reallocate(splits,1,2*size(splits,1),1,2*size(splits,2))
          call reallocate(r,1,size(r,1),1,2*size(r,2),1,2*size(r,3))
          call reallocate(tab,1,size(tab,1),1,2*size(tab,2),1,2*size(tab,3))
       end if
       i=1
       do
          if (i.gt.size(r,1)) then
             call reallocate(r,1,2*size(r,1),1,size(r,2),1,size(r,3))
             call reallocate(tab,1,2*size(tab,1),1,size(tab,2),1,size(tab,3))
          end if
          read(io_tab,*,end=17) r(i,ii,jj),  tab(i,ii,jj)
          if (r(i,ii,jj).eq.1000) exit
          ! write(io_tab+10,*) i,r(i,ii,jj),tab(i,ii,jj)
          i=i+1
       end do
17     splits(ii,jj)=i-1
    end do
    close(io_tab)
  end subroutine read_tabulated_potential_pair

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Calculates nonbonding van der Waals and electrostatic potential
!>  using linear interpolation between two points \n
!>  requires file: fort.43 -- vdW, fort.44 -- Q \n
!>    fort.43: number of tabulated potentials, potential number,
!>  number of points per angstrom, tabulated potential
!>  (repeat last three parts for each additional potential) \n
!>  for unlike interactions, list 1-2 and 2-1 \n
!>  separate potentials with 1000 \n
!>  make sure potential does not go up to infinity! \n
!>  bead type numbers should be defined in suijtab, but it doesn't
!>  matter what they are (the parameters aren't used) \n
!>  KM 12/03/08 \n
!>    fort.44: number of tabulated potentials, potential number,
!>  number of points per angstrom, tabulated potential
!>  (repeat last three parts for each additional potential) \n
!>  for unlike interactions, list 1-2 and 2-1 \n
!>  separate potentials with 1000 \n
!>  KM 04/23/09
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine init_tabulated_potential_pair()
    use sim_system,only:L_vdW_table,L_elect_table
    integer,parameter::initial_size=10,grid_size=1500
    integer::jerr

    if (L_vdW_table) then
       allocate(vdWsplits(1:initial_size,1:initial_size),rvdW(1:grid_size,1:initial_size,1:initial_size),tabvdW(1:grid_size,1:initial_size,1:initial_size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_pair: allocation failed for vdW_table',myid+1)
       call read_tabulated_potential_pair('fort.43',ntabvdW,rvdW,tabvdW,vdWsplits,atoms)
    end if

    if (L_elect_table) then
       allocate(electsplits(1:initial_size,1:initial_size),relect(1:grid_size,1:initial_size,1:initial_size),tabelect(1:grid_size,1:initial_size,1:initial_size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_pair: allocation failed for elect_table',myid+1)
       call read_tabulated_potential_pair('fort.44',ntabelect,relect,tabelect,electsplits,atoms)
    end if
  end subroutine init_tabulated_potential_pair
end MODULE energy_pairwise
