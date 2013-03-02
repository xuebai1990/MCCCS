MODULE energy_pairwise
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
  public::sumup,energy,boltz,coru,init_energy_pairwise,exsix,ljpsur,type_2body

  real,parameter::a15(2)=(/4.0d7,7.5d7/) !< 1-5 correction term for unprotected hydrogen-oxygen interaction; 1 for ether oxygens, 2 for alcohol oxygens
  !OLD VALUES: a15(2)=/17.0**6,16.0**6/)
  integer,allocatable::atom_type(:),nonbond_type(:)

  integer,allocatable::vdWsplits(:,:),electsplits(:,:)
  real,allocatable::rvdW(:,:,:),tabvdW(:,:,:),relect(:,:,:),tabelect(:,:,:)
  integer::ntabvdW,ntabelect

! EXPSIX.INC
  integer,parameter::natom=3 !< for exsix
  real,public::aexsix,bexsix,cexsix,sexsix,consp,consu
  dimension aexsix(natom),bexsix(natom),cexsix(natom),sexsix(natom),consu(natom),consp(natom)

! NSIX.INC
! ***************************************************
! stores interaction parameters for 9-6 potential *
! ***************************************************
  integer,parameter::nxatom=5 !< for ninesix
  real,public::rzero(nxatom*nxatom),epsnx(nxatom*nxatom),shiftnsix(nxatom*nxatom)

! ***********************************************************
! parameters for Generalized Lennard Jones Potential  ***
  real,public::n0,n1
! repulsive part
  parameter(n0=12.0d0)
! attractive part
  parameter(n1=6.0d0)
! Ref:  J. Chem. Phys. 120, 4994 (2004)         ****
! ***********************************************************

! MERCK.INC
  integer,parameter::natomtyp=3
  real,public::epsimmff(natomtyp),sigimmff(natomtyp),smmff(natomtyp),ammff(natomtyp),nmmff(natomtyp),gmmff(natomtyp),sigisq(natomtyp),alphammff(natomtyp),coru_cons(natomtyp),corp_cons(natomtyp)

contains
!*****************************************************************
!> \brief Calculates the total potential energy for a configuration.
!
!> ovrlap: logical, true for substantial atom overlap
!> v*: energies
!> ibox: box number
!> lvol: true if called from volume.f, no output of summary infomation
!******************************************************************
  subroutine sumup(ovrlap,v,vinter,vtail,vintra,vvib,vbend,vtg,vext,velect,vflucq,ibox,lvol)
    use energy_kspace,only:recipsum
    use energy_intramolecular,only:U_bonded
    use energy_3body,only:U3System
    use energy_4body,only:U4System

    logical::ovrlap,lvol
    logical::lexplt,lqimol,lqjmol,lcoulo(numax,numax),lij2,liji,lqchgi
    integer::i,imolty,ii,j,jmolty,jj,ntii,ntjj,ntij,iunit,ibox,nmcount,ntj,k,mmm
    real::v,vinter,vintra,vtail,vvib,vbend,vtg,vext,velect,vflucq,qqii
    real::rcutsq,rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij,rijsq,rho,rij,vrecipsum,rbcut,vwell,calpi
    ! real::vtemp
    real::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq
    real::vintera,velecta,vol
    ! Neeraj & RP for MPI
    real::sum_vvib,sum_vbend,sum_vtg,my_velect
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start SUMUP in ',myid,' for box ', ibox
#endif

    if ( lpbc ) call setpbc(ibox)

    rbcut = rcut(ibox)
    rcutsq = rbcut*rbcut
    calpi=calp(ibox)
    rminsq = rmin * rmin

    vintera = 0.0d0
    velecta = 0.0d0

    ! KM for MPI
    my_velect = 0.0d0

    ovrlap = .false.
    v = 0.0d0
    vinter = 0.0d0
    vintra = 0.0d0
    vtail = 0.0d0
    vtg = 0.0d0
    vbend = 0.0d0
    vvib = 0.0d0
    vext = 0.0d0
    velect = 0.0d0
    vflucq = 0.0d0
    !kea - 3body garofalini term
    v3garo = 0.0d0
    vwell = 0.0d0
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
! not if lgrand and ibox =2
    if (.not.(lgrand.and.ibox.eq.2) .and. .not.lideal(ibox)) then
! RP added for MPI
! do i = 1, nchain - 1
       do i = myid+1,nchain-1,numprocs
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
             molecule2: do j = i + 1, nchain
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
                      rij  = dsqrt(rijsq)
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
                            if ((.not.(liji .and. lij(ntjj))).and.(.not.(lqchgi.and.lqchg(ntjj)))) cycle bead2
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
                         rij=dsqrt(rijsq)

                         ! if ( i .eq. 12 .and. ii .eq. 6 .and. j .eq. 95 .and. jj .eq. 1 ) then
                         ! write(io_output,*) 'CONTROL CONTROL CONTROL'
                         ! write(io_output,*) 'box',ibox,nboxi(i),nboxi(j)
                         ! write(io_output,*) 'i xyz',rxui,ryui,rzui
                         ! write(io_output,*) 'j xyz',rxu(j,jj),ryu(j,jj),rzu(j,jj)
                         ! write(io_output,*) 'r*uij',rxuij,ryuij,rzuij
                         ! write(io_output,*) 'dist2',rijsq
                         ! write(io_output,*) 'distance', dsqrt(rijsq)
                         ! end if

                         if (rijsq.lt.rminsq .and. .not.(lexpand(imolty).or.lexpand(jmolty))) then
                            if ( .not. lvol .and.myid.eq.0) then
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
                            vinter=vinter+U2(rij,rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                         end if

                         velect=velect+Q2(rij,rijsq,rcutsq,i,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo)

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

       call mp_sum(vinter,1,groupid)
       call mp_sum(velect,1,groupid)

       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro .and..not.L_vdW_table) then
          vinter = 4.0d0 * vinter
       end if

! KEA garofalini 3 body potential
       if (lgaro) then
          call triad
          call vthreebody(v3garo)
       end if

       if (hasThreeBody) vinter=vinter+U3System(ibox)
       if (hasFourBody) vinter=vinter+U4System(ibox)

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
                vtail = vtail +  ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
             end do
          end do
          vinter = vinter + vtail
       end if
    end if

!$$$c      write(io_output,*)
!$$$c      write(io_output,*) '+++++++'
!$$$c      vtemp = velect
!$$$c      write(io_output,*) 'direct space part:',velect*qqfact

    if ( ldielect ) then
       call dipole(ibox,0)
    end if

    if (lewald.and..not.lideal(ibox)) then
       call recipsum(ibox,vrecipsum)
       ! update self terms and correction terms
       sself = 0.0d0
       correct = 0.0d0
       ! combine to reduce numerical error
       ! vsc = 0.0d0

       ! RP added for MPI
       ! do i = 1,nchain
       do i = myid+1,nchain,numprocs
          if (nboxi(i) .eq. ibox) then
             imolty = moltyp(i)
             do ii = 1,nunit(imolty)
                sself = sself + qqu(i,ii)*qqu(i,ii)
                ! 1.772.. is the square root of pi
                ! vsc = vsc - qqu(i,ii)*qqu(i,ii)*calpi/1.772453851d0
                do jj = ii+1,nunit(imolty)
                   rxuij = rxu(i,ii) - rxu(i,jj)
                   ryuij = ryu(i,ii) - ryu(i,jj)
                   rzuij = rzu(i,ii) - rzu(i,jj)
                   ! JLR 11-17-09  need call to mimage for intrachain
                   if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
                   ! END JLR 11-17-09
                   rij = dsqrt(rxuij*rxuij + ryuij*ryuij + rzuij*rzuij)

                   ! correct should only be calculated if ii and jj should NOT interact,
                   ! so only calculating it if lqinclu is false
                   ! this part is 1,2 and 1,3
                   if (.not. lqinclu(imolty,ii,jj)) then
                      correct = correct + qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0d0)/rij
                      ! vsc = vsc + qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0d0)/rij
                   else
                      correct=correct+(1.0d0 - qscale2(imolty,ii,jj))*qqu(i,ii)*qqu(i,jj)* (erfunc(calpi*rij)-1.0d0)/rij
                      ! vsc = vsc + (1.0d0 - qscale2(imolty,ii,jj))*qqu(i,ii)*qqu(i,jj)*(erfunc(calpi*rij)-1.0d0)/rij
                   end if
                end do
             end do
          end if
       end do

       call mp_sum(correct,1,groupid)
       call mp_sum(sself,1,groupid)

! vdipole = (dipolex*dipolex+dipoley*dipoley+dipolez*dipolez)*(2.0d0*onepi)/(3.0d0*boxlx(ibox)**3.0d0)
! write(io_output,*) dipolex,dipoley,dipolez
       sself = -sself*calpi/sqrtpi
       velect = velect + sself + correct + vrecipsum
! velect = velect + vsc + vrecipsum
    end if

!$$$c at this point velect contains all intermolecular charge interactions,
!$$$c plus the ewald self term and intramolecular corrections
!$$$
! write(io_output,*)
! write(io_output,*) '== After Inter === velect is:',velect*qqfact
!$$$
!$$$       vtemp = velect

! ################################################################

! have to recalculate ewald terms if volume changes
    if ( .not. lvol .or. (lvol .and. lewald) ) then
! *******************************
! INTRACHAIN INTERACTIONS ***
! *******************************
! write(io_output,*) 'starting intrachain'
! loop over all chains i
! RP added for MPI
! do i = 1, nchain
       do i = myid+1,nchain,numprocs
! check if i is in relevant box ###
! write(io_output,*) 'nboxi(i),i,ibox',nboxi(i),i,ibox
          if ( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
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
                      rij  = dsqrt(rijsq)
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
                            !!?? skip intra if it is bending 1-3 and using a table??
                            if (L_bend_table) then
                               do mmm=1,inben(imolty,ii)
                                  if (ijben3(imolty,ii,mmm).eq.jj)then
                                     vintra = vintra + lininter_bend(rij,itben(imolty,ii,mmm))
                                     goto 94
                                  end if
                               end do
                            end if

                            vintra=vintra+U2(rij,rijsq,i,imolty,ii,ntii,i,imolty,jj,ntjj,ntij)
                         end if
                      end if !if (linclu(imolty,ii,jj))

94                    if (lqinclu(imolty,ii,jj)) then
                         ! calculate intramolecular charge interaction
                         my_velect=my_velect+qscale2(imolty,ii,jj)*Q2(rij,rijsq,rcutsq,i,imolty,ii,ntii,lqchg(ntii),i,imolty,jj,ntjj,calpi,lcoulo)
                      end if
                   end if
                end do
             end do
          end if
       end do

! RP added for MPI----- Returning from ovrlap--------------
299    continue

       call mp_lor(ovrlap,1,groupid)
       if(ovrlap)then
          ! write(io_output,*)'941: in sumup ovrlap=',ovrlap,'myid=',myid
          return
       end if
! -----------------------------------------
       call mp_sum(vintra,1,groupid)
       call mp_sum(my_velect,1,groupid)
       velect = velect + my_velect

!$$$       vtemp = velect - vtemp
!$$$
!$$$c       write(io_output,*) '== Intra Velect ===',vtemp*qqfact
! write(io_output,*) '== After Intra  === velect is:',velect*qqfact
! write(io_output,*) 'vintra ', vintra
! write(io_output,*) 'vinter ', vinter
! write(io_output,*) 'test', vintra

       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro .and..not.L_vdW_table) then
          vintra = 4.0d0 * vintra
       end if

! ################################################################

! *************************************
! INTRACHAIN FLUCQ INTERACTIONS ***
! *************************************
!c RP added for MPI
! removed at some point?
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
                            vflucq = vflucq + xiq(ntii)*qqii + jayself(ntii)*qqii*qqii
                         else
                            ntjj = ntype(imolty,jj)
                            ntij = type_2body(ntii,ntjj)

                            vflucq = vflucq  + jayq(ntij)*qqii*qqu(i,jj)
                         end if
                      end do
                   end do
                   vflucq = vflucq - fqegp(imolty)
                else
                   vflucq = 0.0d0
                end if
             end if
          end if
       end do
! ------------------------------------------

! **************************************************
! CALCULATION OF VIB. + BEND. + TORS. ENERGY ***
! **************************************************

! NOTE here virtual coordinates can be used!!!
! RP added for MPI
! do i = 1, nchain
       do i = myid+1,nchain,numprocs
          imolty = moltyp(i)
          ! check if i is in relevant box ###
          if ( nboxi(i) .eq. ibox ) then
             call U_bonded(i,imolty,sum_vvib,sum_vbend,sum_vtg)
             vvib=vvib+sum_vvib
             vbend=vbend+sum_vbend
             vtg=vtg+sum_vtg
          end if
       end do
       call mp_sum(vvib,1,groupid)
       call mp_sum(vbend,1,groupid)
       call mp_sum(vtg,1,groupid)
    end if !if ( .not. lvol .or. (lvol .and. lewald) )

! ################################################################

! ***************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

! for adsorption isotherms, don't calculate energy w/surface
! in box 2
    if ((lelect_field) .or. ((ibox .eq. 1) .and. (ljoe .or. lsami .or. lmuir .or.  lexzeo .or. lgraphite .or. lslit))) then
! RP added for MPI
! do i = 1, nchain
       do i = myid+1, nchain,numprocs
          ! check if i is in relevant box ###
          if ( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             do j = 1, nunit(imolty)
                ntj = ntype(imolty,j)
                vext=vext+U_ext(ibox,i,j,ntj)
             end do
          end if
       end do

       call mp_sum(vext,1,groupid)
    end if

! --------------------------------------------------------------------
! calculation of additional gaussian potential needed in thermodynamic
! integration in stages b and c
! --------------------------------------------------------------------
    if (lmipsw) then
! RP added for MPI
! do i = 1, nchain
       do i = myid+1,nchain,numprocs
          imolty = moltyp(i)
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
                      if (awell(ii,k,imolty).lt.1.0d-6) cycle
                      rxui = rxu(i,ii)
                      ryui = ryu(i,ii)
                      rzui = rzu(i,ii)
                      rxuij = rxui-rxwell(j,imolty)
                      ryuij = ryui-rywell(j,imolty)
                      rzuij = rzui-rzwell(j,imolty)
                      call mimage(rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                      vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
                   end do
                end if
             end do
          end if
       end do

       call mp_sum(vwell,1,groupid)
    end if

! ----------------------------------------------------------------------------

! write(io_output,*) 'self,corr:',(velect-vrecipsum)*qqfact
! write(io_output,*) 'vsc, new self cor:',vsc*qqfact
! write(io_output,*) 'recip space part :',vrecipsum*qqfact
! write(io_output,*) 'sc and recip:',(vsc+vrecipsum)*qqfact

    if (.not.L_elect_table) then
       velect = velect*qqfact
    end if

    v = vinter + vintra + vext + velect + vflucq + v3garo

! write(io_output,*) 'v in sumup',v

    vipsw = v
    vwellipsw = vwell

    if (lstagea) then
       v = (1.0d0-lambdais*(1.0d0-etais))*v
    else if (lstageb) then
       v = etais*v+lambdais*vwell
    else if (lstagec) then
       v = (etais+(1.0d0-etais)*lambdais)*v+(1.0d0-lambdais)*vwell
    end if

    v = v + vvib + vbend + vtg

    if ( .not. lvol.and.myid.eq.0 ) then
       write(io_output,*)
       write(io_output,*) 'sumup control'
       write(io_output,*) 'number of chains', nmcount
       do i = 1, nmolty
          write(io_output,*) 'number of chains of type',i,ncmt(ibox,i)
       end do
       write(io_output,*) 'inter lj energy ', vinter
       write(io_output,*) 'intra lj energy ', vintra
       if (ltailc) write(io_output,*) 'Tail correction ', vtail
       write(io_output,*) 'bond vibration  ', vvib
       write(io_output,*) 'bond bending    ', vbend
       write(io_output,*) 'torsional       ', vtg
       write(io_output,*) 'external        ', vext
       write(io_output,*) 'coulombic energy', velect
       ! write(io_output,*) 'exact energy    ', 1.74756*1.67*831.441/3.292796
       write(io_output,*) 'fluc Q energy   ', vflucq
       write(io_output,*) 'well energy     ', vwellipsw
       if(lgaro) write(io_output,*) '3-body garo     ', v3garo
       write(io_output,*) 'total energy    ', v
    end if

#ifdef __DEBUG__
    write(io_output,*) 'end SUMUP in ',myid
#endif

    return
  end subroutine sumup

!*****************************************************************
!> \brief Calculates the total potential energy for a configuration.
!>
!> i: calculate the energies associated with chain i
!> imolty: molecule type of chain i
!> v*: energies
!> flagon: flag for old(flagon=1)/new(flagon=2) configurations
!> ibox: box number of chain i
!> istart, iuend: calculate energies from bead istart to bead iuend for chain i
!> lljii: whether to include intramolecular LJ interactions
!> ovrlap: atom overlap
!> ltors: whether to calculate torsional energy
!> lcharge_table: whether need to set up charge interaction table; true if called from CBMC
!> lfavor:
!*****************************************************************
  subroutine energy(i,imolty,v,vintra,vinter,vext,velect,vewald,flagon,ibox,istart,iuend,lljii,ovrlap,ltors,vtors,lcharge_table,lfavor,lAtom_traxyz)
    use sim_particle,only:lnn
    use energy_intramolecular,only:U_torsion
    use energy_3body,only:U3MolSys
    use energy_4body,only:U4MolSys
    use energy_garofalini,only:triad_en

    logical::lqimol,lqjmol,lexplt,lcoulo(numax,numax),lfavor,lij2,liji,lqchgi
    logical::lljii,ovrlap,ltors,lcharge_table,lAtom_traxyz

    integer::growii,growjj,k,jcell(nmax),nmole
    integer::i,ibox,istart,iuend,ii,ntii,flagon,jjj,iii,mmm,j,jj,ntjj,ntij,ntj,imolty,jmolty,jjend
    integer::nchp2

    real::v,vintra,vinter,vext,rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq,ryuij,rzuij,rij,rijsq,vtors,velect,vewald,rbcut,calpi
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
    v = 0.0d0
    vinter = 0.0d0
    vintra = 0.0d0
    vext = 0.0d0
    velect = 0.0d0
    vewald = 0.0d0
    vtors = 0.0d0
    sself  = 0.0d0
    correct = 0.0d0
    !kea
    v3garo = 0.0d0

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
                rij  = dsqrt(rijsq)
                rcm = rbcut + rcmi + rcmu(j)
                rcmsq = rcm*rcm
                ! write(io_output,*) rcm,rcmi,rcmu(j)
                if ( lfavor ) then
                   favor(j) = (rminsq/rijsq)**2*5.0d0
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
                   rij  = dsqrt(rijsq)
                   if (rijsq.lt.rminsq .and. .not.(lexpand(imolty).or.lexpand(jmolty))) then
                      ovrlap = .true.
                      ! write(io_output,*) 'inter ovrlap:',i,j, myid
                      ! write(io_output,*) 'i xyz',rxui,ryui,rzui
                      ! write(io_output,*) 'j xyz',rxu(j,jj),ryu(j,jj),rzu(j,jj)
                      ! write(io_output,*) 'ii:',ii,'jj:',jj
                      ! write(io_output,*) 'distance', dsqrt(rijsq)
                      ! RP added for MPI
                      ! return
                      goto 99
                   else if ((rijsq .lt. rcutsq) .or. lijall) then
                      vinter=vinter+U2(rij,rijsq,nchp2,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                   end if

                   velect=velect+Q2(rij,rijsq,rcutsq,nchp2,imolty,ii,ntii,lqchgi,j,jmolty,jj,ntjj,calpi,lcoulo)

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

    call mp_sum(vinter,1,groupid)
    call mp_sum(velect,1,groupid)
! -----------------------------------------------

    if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro .and..not.L_vdW_table) then
       vinter = 4.0d0 * vinter
    end if

!kea - garo: add three body loop for intermolecular interactions
    if (lgaro.and..not.lideal(ibox)) then
       if(flagon.eq.2) then
          call triad_en(i,v3garo,neigh_icnt,neighi,ndiji,nxiji,nyiji,nziji,.true.)
       else if(flagon.eq.1) then
          call triad_en(i,v3garo,neigh_j,neighj,ndijj,nxijj,nyijj,nzijj,.false.)
       end if
    end if

    if (hasThreeBody) vinter=vinter+U3MolSys(i,istart,iuend,flagon)
    if (hasFourBody) vinter=vinter+U4MolSys(i,istart,iuend,flagon)

! ################################################################

! the intramolecular van der waals and ewald terms have to be calculated
! for the explicit atom placement models
! *******************************
! INTRACHAIN INTERACTIONS ***
! *******************************

! JLR 11-19-09 commenting this out, alway do mimage for intrachain
! for expanded ensemble
! lmim = .false.
! mlen2 = rcmu(nchp2)*2d0
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
          rij = dsqrt(rijsq)

          if (lqinclu(imolty,ii,jj) ) then
             ! calculation of intramolecular electrostatics
             velect=velect+qscale2(imolty,ii,jj)*Q2(rij,rijsq,rcutsq,nchp2,imolty,ii,ntii,lqchg(ntii),nchp2,imolty,jj,ntjj,calpi,lcoulo)
          end if

          ! calculation of other non-bonded interactions
          if ( linclu(imolty,ii,jj) ) then
             if (lljii) then
                if (rijsq.lt.rminsq .and. .not.lexpand(imolty)) then
                   ovrlap = .true.
                   ! write(io_output,*) 'intra ovrlap:',ii,jj
                   return
                else if (rijsq.lt.rcutsq .or. lijall) then
                   !!?? skip intra if it is bending 1-3 and using a table??
                   if (L_bend_table) then
                      do mmm=1,inben(imolty,ii)
                         if (ijben3(imolty,ii,mmm).eq.jj) then
                            vintra = vintra + lininter_bend(rij,itben(imolty,ii,mmm))
                            goto 96
                         end if
                      end do
                   end if

                   vintra=vintra+U2(rij,rijsq,nchp2,imolty,ii,ntii,nchp2,imolty,jj,ntjj,ntij)
                end if
             end if
          end if

96        if (lewald.and..not.lideal(ibox)) then
             ! compute the ewald intramolecular (self and correction) terms for
             ! the interactions of the placed atoms with themselves, and with the
             ! rest of their own molecule, if there's no interaction
             ! these are 1,2 and 1,3
             if (.not. lqinclu(imolty,ii,jj)) then
                correct=correct+qquion(ii,flagon)*qquion(jj,flagon)*(erfunc(calpi*rij)-1.0d0)/rij
                ! 1,4 interaction which we scale by qscale
             else
                correct=correct+(1.0d0-qscale2(imolty,ii,jj))*qquion(ii,flagon)*qquion(jj,flagon)*(erfunc(calpi*rij)-1.0d0)/rij
             end if
          end if
       end do
       if (lewald.and..not.lideal(ibox)) then
          sself = sself + qquion(ii,flagon)*qquion(ii,flagon)
       end if
    end do
    if (lewald.and..not.lideal(ibox)) then
       sself = -sself * calpi/sqrtpi
       vewald = sself + correct
    end if
    if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj  .and. .not. lninesix .and..not.L_vdW_table) then
       vintra = 4.0d0 * vintra
    end if

! ################################################################

! ***************************************************************
! CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************

    if ((lelect_field.and.lqimol) .or. ((ibox .eq. 1) .and. (ljoe .or. lsami .or. lmuir .or. lexzeo .or. lgraphite .or. lslit))) then
       do j = istart,iuend
          ntj = ntype(imolty,j)
          vext=vext+U_ext(ibox,nchp2,j,ntj)
       end do
    end if

! *********************************************************************
! calculation of torsion energy for explicit atom methyl groups ****
! *********************************************************************
    if ( ltors ) then
       vtors=U_torsion(nchp2,imolty,nugrow(imolty)+1,.true.)
    end if

! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
    vwell = 0.0d0
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
                if (awell(ii,k,imolty).lt.1.0d-6) cycle
                rxui = rxuion(ii,flagon)
                ryui = ryuion(ii,flagon)
                rzui = rzuion(ii,flagon)
                rxuij = rxui-rxwell(j,imolty)
                ryuij = ryui-rywell(j,imolty)
                rzuij = rzui-rzwell(j,imolty)
                call mimage(rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
             end do
          end if
       end do
    end if

! ----------------------------------------------------------------------------

    if (.not.L_elect_table) then
       velect = velect*qqfact
       vewald = vewald*qqfact
    end if

! note that vintra is only computed when the flag lljii is true
    v = vinter + vext + vintra + velect + vewald + v3garo
! write(io_output,*) 'vinter:',vinter,'vext:',vext,'vintra:',vintra,'velect',velect,'vewald:',vewald,'v'

    if (flagon.eq.1) then
       vipswo = v
       vwellipswo = vwell
    else
       vipswn = v
       vwellipswn = vwell
    end if

    if (lmipsw) then
       if (lstagea) then
          v = (1.0d0-lambdais*(1.0d0-etais))*v
       else if (lstageb) then
          v = etais*v+lambdais*vwell
       else if (lstagec) then
          v = (etais+(1.0d0-etais)*lambdais)*v+(1.0d0-lambdais)*vwell
       end if
    end if

#ifdef __DEBUG__
    ! write(io_output,*) 'v :', v
    write(io_output,*) 'end ENERGY in ',myid
#endif
    return
  end subroutine energy

!*****************************************************************
!> \brief Calculates the potential energy and the boltzmann factor
!>       for ichoi trial positions.
!>
!> lnew: true for new configurations
!> lfirst: true for insertion of the first bead in swap moves
!> ovrlap: logical variable, true for walk termination
!> i: calculates Boltzmann weights for the newly-grown beads in chain i
!> icharge: usually identical to i
!> imolty: molecule type of chain i
!> ibox: box number of chain i
!> ichoi: number of trial positions
!> iufrom: the bead from which the new beads were grown
!> ntogrow: number of new beads that have been grown
!> glist: the list of new beads that have been grown; ntogrow entries
!> maxlen: maximum possible distance of the newly-grown beads from iufrom
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
    real::vinter,vintra,vext,velect,vewald,vwell,v,rcutsq,rcinsq
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
       my_vtry(j) = 0.0d0
       my_vtrintra(j) = 0.0d0
       my_vtrelect(j) = 0.0d0
       my_vtrext(j) = 0.0d0
       my_vtrinter(j) = 0.0d0
       my_vtrewald(j) = 0.0d0
       my_bfac(j) = 0.0d0
       my_vipswot(j) = 0.0d0
       my_vwellipswot(j) = 0.0d0
       my_vipswnt(j) = 0.0d0
       my_vwellipswnt(j) = 0.0d0
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

!!?? to be MPI parallelized
       do j = 1,nchain
          lcmno(j) = .false.
          if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
             rxuij = rxui-xcm(j)
             ryuij = ryui-ycm(j)
             rzuij = rzui-zcm(j)
             ! minimum image the pseudo-ctrmas pair separation
             if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)

             rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
             rij  = dsqrt(rijsq)

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

       vinter = 0.0d0
       vintra = 0.0d0
       vext = 0.0d0
       velect = 0.0d0
       vewald = 0.0d0

       ! Only if L_Coul_CBMC is true, then compute electrostatic interactions/corrections
       if(L_Coul_CBMC.and.lewald.and..not.lideal(ibox)) then
          do count = 1,ntogrow
             ii = glist(count)
             ! This part does not change for fixed charge moves, but is
             ! used in the swap rosenbluth weight. - ewald self term
             vewald = vewald - qqu(icharge,ii)*qqu(icharge,ii)*calp(ibox)/sqrtpi
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
                   rij = dsqrt(rijsq)
                end if
                if ( linclu(imolty,ii,iu) .or. lqinclu(imolty,ii,iu)) then
                   if ( linclu(imolty,ii,iu) ) then
                      if (rijsq.lt.rminsq.and..not.lexpand(imolty)) then
                         ! RP added for MPI
                         my_lovr(my_itrial) = .true.
                         ! write(io_output,*) 'intra overlap'
                         goto 19
                      else if (rijsq.lt.rcutsq .or. lijall) then
                         !!?? skip intra if it is bending 1-3 and using a table??
                         if (L_bend_table) then
                            do mmm=1,inben(imolty,ii)
                               if (ijben3(imolty,ii,mmm).eq.iu) then
                                  vintra = vintra + lininter_bend(rij,itben(imolty,ii,mmm))
                                  goto 96
                               end if
                            end do
                         end if

                         vintra=vintra+U2(rij,rijsq,i,imolty,ii,ntii,i,imolty,iu,ntjj,ntij)
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
                         velect = velect + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(icharge,iu)*lininter_elect(rij,ntii,ntjj)
                      else if (lewald) then
                         ! compute real space term of vewald
                         velect = velect + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(icharge,iu)*erfunc(calp(ibox)*rij)/ rij
                         ! ewald sum correction term
                         corr = (1.0d0 - qscale2(imolty,ii,iu))*qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0d0) /rij
                         vewald = vewald + corr
                      else
                         velect = velect + qscale2(imolty,ii,iu)*qqu(icharge,ii)*qqu(i,iu)/rij
                      end if
                   end if
                   ! end charge calculation
                   ! will only add correction if lqinclu is false.
                else if (L_Coul_CBMC.and.lewald) then
                   ! ewald sum correction term
                   corr = qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0d0) /rij
                   vewald = vewald + corr
                end if
             end do
          end do
          !!?? double Ewald correction??
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
                   !!?? no mimage convertion??
                   rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   rij   = dsqrt(rijsq)
                   ! ewald sum correction term
                   corr = qqu(icharge,ii)*qqu(icharge,iu)*(erfunc(calp(ibox)*rij)-1.0d0)/rij
                   vewald = vewald + corr
                end do
             end do
          end if

          if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff  .and. .not. lgenlj .and. .not. lninesix .and..not.L_vdW_table.and..not.L_bend_table) then
             vintra = 4.0d0 * vintra
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
                      rij   = dsqrt(rijsq)
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
                      ! is there a way to pull this out of the loops? ***
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
                      rij   = dsqrt(rijsq)
                      ! compute vinter (eg. lennard-jones)
                      if (rijsq.lt.rminsq.and..not.(lexpand(imolty).or.lexpand(jmolty))) then
                         my_lovr(my_itrial) = .true.
                         ! write(io_output,*) 'j:',j,jj
                         ! write(io_output,*) 'rjsq:',rijsq,rminsq
                         goto 19
                      else if (rijsq.lt.rcinsq.or.lijall) then
                         vinter=vinter+U2(rij,rijsq,i,imolty,ii,ntii,j,jmolty,jj,ntjj,ntij)
                      end if

                      ! compute velect (coulomb and ewald)
                      if (L_Coul_CBMC.and.lqchg(ntii).and.lqchg(ntjj).and.rijsq.lt.rcinsq) then
                         ! boltz.f has problem to compute the electrostatic interactions
                         ! in a group-based way because the leader q might not be grown at
                         ! present, so it calculates electrostatic interaction not based on
                         ! group but on its own distance in SC, but should be corrected
                         ! later by calling energy subroutine.
                         if (L_elect_table) then
                            velect = velect + qqu(icharge,ii)*qqu(j,jj)*lininter_elect(rij,ntii,ntjj)
                         else if (lewald) then
                            ! compute real space term of velect
                            velect = velect + qqu(icharge,ii)*qqu(j,jj)*erfunc(calp(ibox)*rij)/rij
                         else
                            ! compute all electrostatic interactions
                            velect = velect + qqu(icharge,ii)*qqu(j,jj)/ rij
                         end if
                      end if
                   end do
                end do
             end if
          end do do_nmole

          if (.not.lsami.and..not.lexpsix.and..not.lmmff.and..not.lgenlj.and..not.lninesix.and..not.L_vdW_table) then
             vinter = 4.0d0 * vinter
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
                vext=vext+U_ext(ibox,nchain+2,ii,ntii)
             end do
          end if
       end if

! --------------------------------------------------------------------------
! well potential for thermodynamic integration stages b and c
! --------------------------------------------------------------------------
       vwell = 0.0d0
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
                   if (awell(ii,k,imolty).lt.1.0d-6) cycle
                   rxui = rxp(count,itrial)
                   ryui = ryp(count,itrial)
                   rzui = rzp(count,itrial)
                   rxuij = rxui-rxwell(j,imolty)
                   ryuij = ryui-rywell(j,imolty)
                   rzuij = rzui-rzwell(j,imolty)
                   call mimage(rxuij,ryuij,rzuij,ibox)
                   rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                   vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
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
          my_bfac(my_itrial) = 0.0d0
       else
          if (.not.L_elect_table) then
             velect = velect*qqfact
             vewald = vewald*qqfact
          end if
          v = vinter+vintra+vext+velect+vewald
          if (.not.lnew) then
             my_vipswot(my_itrial) = v
             my_vwellipswot(my_itrial) = vwell
          else
             my_vipswnt(my_itrial) = v
             my_vwellipswnt(my_itrial) = vwell
          end if

          if (lstagea) then
             v = (1.0d0-lambdais*(1.0d0-etais))*v
          else if (lstageb) then
             v = etais*v+lambdais*vwell
          else if (lstagec) then
             v = (etais+(1.0d0-etais)*lambdais)*v+ (1.0d0-lambdais)*vwell
          end if

          my_vtry(my_itrial) = v
          my_vtrintra(my_itrial) = vintra
          my_vtrext(my_itrial)   = vext
          my_vtrinter(my_itrial) = vinter
          my_vtrelect(my_itrial) = velect
          my_vtrewald(my_itrial) = vewald
          ! write(23,*) 'itrial' ,itrial
          ! write(23,*) vtry(itrial), vtrintra(itrial), vtrext(itrial),
          !   &       vtrinter(itrial),vtrelect(itrial), vtrewald(itrial)

          if ((my_vtry(my_itrial)*beta).gt.(2.3d0*softcut))then
             ! write(io_output,*) 'caught by softcut',vtry(itrial)*beta
             my_lovr(my_itrial) = .true.
             my_bfac(my_itrial) = 0.0d0
          else if((my_vtry(my_itrial)*beta).lt.-2.303d0*308)then
             ! write(io_output,*) '### warning: weight too big out of range'
             my_lovr(my_itrial) = .true.
             my_bfac(my_itrial) = 0.0d0
          else
             my_bfac(my_itrial) = dexp ( -(my_vtry(my_itrial)*beta) )
          end if
       end if
    end do

    ! if (ldebug) then
    !    write(100+myid,*) 'boltz for box ',ibox,' molecule ',i,' from myid ',myid
    !    do my_itrial=1,counts(myid+1)
    !       write(100+myid,*) "my_lovr(",my_itrial,") = ",my_lovr(my_itrial), "my_bfac(",my_itrial,") = ",my_bfac(my_itrial)
    !    end do
    ! end if

    call mp_allgather(my_vtry,vtry,rcounts,displs,groupid)
    call mp_allgather(my_vtrintra,vtrintra,rcounts,displs,groupid)
    call mp_allgather(my_vtrext,vtrext,rcounts,displs,groupid)
    call mp_allgather(my_vtrinter,vtrinter,rcounts,displs,groupid)
    call mp_allgather(my_vtrelect,vtrelect,rcounts,displs,groupid)
    call mp_allgather(my_vtrewald,vtrewald,rcounts,displs,groupid)
    call mp_allgather(my_bfac,bfac,rcounts,displs,groupid)
    call mp_allgather(my_vipswot,vipswot,rcounts,displs,groupid)
    call mp_allgather(my_vwellipswot,vwellipswot,rcounts,displs,groupid)
    call mp_allgather(my_vipswnt,vipswnt,rcounts,displs,groupid)
    call mp_allgather(my_vwellipswnt,vwellipswnt,rcounts,displs,groupid)
    call mp_allgather(my_lovr,lovr,rcounts,displs,groupid)
    ovrlap = .true.
    if (ANY(.not.lovr(1:ichoi))) ovrlap=.false.

    ! if (ldebug) then
    !    do my_itrial=1,ichoi
    !       write(100+myid,*)"lovr(",my_itrial,") = ",lovr(my_itrial)
    !       write(100+myid,*)"bfac(",my_itrial,") = ",bfac(my_itrial)
    !    end do
    ! end if
! ----------------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'end BOLTZ in ',myid
#endif
    return
  end subroutine boltz

! **********************************
! tail-corrections in energy ***
! **********************************
  function coru(imolty,jmolty,rho,ibox)

      real::coru,rci3,rho,epsilon2,sigma2
      real::rci1
      integer::imolty,jmolty,ii,jj, ntii, ntjj, ntij ,ibox

      coru = 0.0d0
      do ii = 1, nunit(imolty)
         ntii = ntype(imolty,ii)
         do jj = 1, nunit(jmolty)
            ntjj = ntype(jmolty,jj)
            ntij = type_2body(ntii,ntjj)
            if (lexpsix) then
               coru = coru + rho*consu(ntij)
            else if (lmmff) then
               coru = coru + rho * epsimmff(ntij) * coru_cons(ntij) * sigimmff(ntij)**3.0d0*twopi
            else if (lninesix) then
               coru = coru + 8.0d0*onepi*rho*epsnx(ntij)* rzero(ntij)**3*(rzero(ntij)/rcut(ibox))**3* ((rzero(ntij)/rcut(ibox))**3/3.0d0 - 1.0d0)
            else if (lgenlj) then
               rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
               rci1 = rci3 **(1.0d0/3.0d0)

               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigma_f(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon_f(imolty,ii) *epsilon_f(jmolty,jj))
               else if ( lexpand(imolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigi(ntjj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon_f(imolty,ii)*epsi(ntjj))
               else if ( lexpand(jmolty) ) then
                  sigma2 = (sigma_f(jmolty,jj)+sigi(ntii))**2/4.0d0
                  epsilon2 = dsqrt(epsilon_f(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               end if
               coru = coru + 2.0d0 * onepi * epsilon2 * sigma2 ** (1.50d0) * rho * (  (( (2.0d0**(4.0d0*n1/n0))/(2.0d0*n1-3.0d0)) * rci1 **(2.0d0*n1-3.0d0) ) - ( (2.0d0**((2.0d0*n1/n0)+1.0d0))/(n1-3.0d0)) * rci1 **(n1-3.0d0) )

            else
               rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigma_f(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon_f(imolty,ii) *epsilon_f(jmolty,jj))
               else if ( lexpand(imolty) ) then
                  sigma2 = (sigma_f(imolty,ii)+sigi(ntjj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon_f(imolty,ii)*epsi(ntjj))
               else if ( lexpand(jmolty) ) then
                  sigma2 = (sigma_f(jmolty,jj)+sigi(ntii))**2/4.0d0
                  epsilon2 = dsqrt(epsilon_f(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               end if
               coru = coru +  8.0d0 * onepi * epsilon2 *  sigma2**(1.5d0) *rho *  (rci3 * rci3 * rci3 / 9.0d0 - rci3 / 3.0d0)
            end if
         end do
      end do
      return
  end function coru

  subroutine init_energy_pairwise(file_ff,lmixlb,lmixjo)
    use util_search,only:initiateTable,addToTable
    use util_memory,only:reallocate
    use energy_intramolecular,only:init_energy_bonded
    use energy_external,only:init_energy_external
    use energy_sami,only:susami

    character(LEN=*),INTENT(IN)::file_ff
    logical,intent(in)::lmixlb,lmixjo
    integer,parameter::initial_size=20
    character(LEN=default_string_length)::line_in
    integer::io_ff,jerr,i,j,ij,ji,ibox,nmix
    real::rzeronx(nxatom),epsilonnx(nxatom),rcheck,sr2,sr6,rs1,rs7,sr7,sigmaTmp,epsilonTmp

    if (lshift.or.lmmff.or.lninesix.or.lexpsix.or.lsami) then
       ! Keep the rcut same for each box
       do ibox = 2,nbox
          if (dabs(rcut(1)-rcut(ibox)).gt.1.0d-10) then
             call err_exit(__FILE__,__LINE__,'Keep rcut for each box same',myid+1)
          end if
       end do
    end if

    call initiateTable(atoms,nmolty)

    call init_tabulated_potential_pair()

    io_ff=get_iounit()
    open(unit=io_ff,access='sequential',action='read',file=file_ff,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open forcefield input file',myid+1)
    end if

! Looking for section ATOMS
    CYCLE_READ_ATOMS:DO
       call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) then
          exit cycle_read_atoms
       end if

       if (UPPERCASE(line_in(1:5)).eq.'ATOMS') then
          allocate(atom_type(1:initial_size),sigi(1:initial_size),epsi(1:initial_size),qelect(1:initial_size),mass(1:initial_size),lij(1:initial_size),lqchg(1:initial_size),chemid(1:initial_size),stat=jerr)
          if (jerr.ne.0) then
             call err_exit(__FILE__,__LINE__,'init_pairwise: atoms allocation failed',jerr)
          end if
          sigi = 0.0d0
          epsi = 0.0d0
          qelect = 0.0d0
          mass = 0.0d0
          lij = .true.
          lqchg = .false.

          nntype=0
          do
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) then
                call err_exit(__FILE__,__LINE__,'Reading section ATOMS',jerr)
             end if
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

    ! do j=1,nntype
    ! if (sigi(j).ne.0.or.epsi(j).ne.0.or.mass(j).ne.0.or.qelect(j).ne.0) then
    ! write(104,'(I3,1X,I1,1X,F8.5,1X,F9.4,1X,F7.4,1X,F11.7,1X,A,1X,A)') j,1,sigi(j),epsi(j),qelect(j),mass(j),trim(chemid(j)),'#'//trim(chname(j))
    ! end if
    ! end do

    nmix=nntype*nntype

    allocate(lpl(1:nntype),q1(1:nntype),xiq(1:nntype),jayself(1:nntype),nonbond_type(1:nmix),sig2ij(1:nmix),epsij(1:nmix),rminee(1:nmix),ecut(1:nmix),jayq(1:nmix),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'init_pairwise: nonbond allocation failed',jerr)
    end if

    xiq = 0.0d0
    lpl = .false.
    sig2ij=0.0d0
    epsij=0.0d0
    jayq=0.0d0

    if (lgaro) then
! KEA adding garofalini silica/water potential
! see J. Phys. Chem. 94 5351 (1990)
! form:
! U(2) = A exp(-rij/rho)+[zi zj erfc (rij/beta)]/rij + a/[1+exp(b/(rij-c))]
! U(3) = h3(rij rik thetajik) + h3(rjk rji thetakji) + h3(rki rkj thetaijk)
! h3(rij rik thetajik) = lambda exp[(gamma/(rij-rij*))+(gamma/(rik-rik*))]*
!             [cos(theta)-cos(theta*)]**2     for rij<rij* and rik<rik*
! 0 otherwise
       do i=1,6
          galpha(i) = 0.0d0
          grho(i) = 0.0d0
          gbeta(i) = 0.0d0
          ecut(i) = 0.0d0
          do j=1,3
             ga(i,j) = 0.0d0
             gb(i,j) = 0.0d0
             gc(i,j) = 0.0d0
          end do
       end do
       do i=1,4
          glambda(i) = 0.0d0
          do j=1,2
             grij(i,j) = 0.0d0
             grijsq(i,j) = 0.0d0
             ggamma(i,j) = 0.0d0
          end do
          gtheta(i) = 0.0d0
       end do

! Parameters (galpha,grho,gbeta,ga,gb,gc; lambda,grij,ggamma,gtheta)
! Si-Si
       galpha(1) = 13597175.7d0
       grho(1) = 0.29d0
       gbeta(1) = 2.29d0
       lqchg(1) = .true.
       qelect(1) = 4.0d0
       mass(1) = 28.09d0
       ecut(1) = garofalini(rcut(1)*rcut(1),1,qelect(1),qelect(1) ,1,1)
       chemid(1) = 'Si '

! O-O
       galpha(2) = 5251972.5d0
       grho(2) = 0.29d0
       gbeta(2) = 2.34d0
       lqchg(2) = .true.
       qelect(2) = -2.0d0
       mass(2) = 16.00d0
       ecut(2) = garofalini(rcut(1)*rcut(1),2,qelect(2),qelect(2) ,2,2)
       chemid(2) = 'O  '

! H-H
       galpha(3) = 246299.4d0
       grho(3) = 0.35d0
       gbeta(3) = 2.1d0
       ga(3,1) = -38243.8d0
       gb(3,1) = 6.0d0
       gc(3,1) = 1.51d0
       ga(3,2) = 2515.9d0
       gb(3,2) = 2.0d0
       gc(3,2) = 2.42d0
       lqchg(3) = .true.
       qelect(3) = 1.0d0
       mass(3) = 1.0078d0
       ecut(3) = garofalini(rcut(1)*rcut(1),3,qelect(3),qelect(3) ,3,3)
       chemid(3) = 'H  '

! Si-O
       galpha(4) =  21457024.2d0
       grho(4) = 0.29d0
       gbeta(4) = 2.34d0
       ecut(4) = garofalini(rcut(1)*rcut(1),4,qelect(1),qelect(2) ,1,2)

! Si-H
       galpha(5) = 499842.9d0
       grho(5) = 0.29d0
       gbeta(5) = 2.31d0
       ga(5,1) = -33715.5d0
       gb(5,1) = 6.0d0
       gc(5,1) = 2.2d0
       ecut(5) = garofalini(rcut(1)*rcut(1),5,qelect(1),qelect(3) ,1,3)

! 0-H
       galpha(6) = 2886049.4d0
       grho(6) = 0.29d0
       gbeta(6) = 2.26d0
       ga(6,1) = -15096.7d0
       gb(6,1) = 15.0d0
       gc(6,1) = 1.05d0
       ga(6,2) = 55353.6d0
       gb(6,2) = 3.2d0
       gc(6,2) = 1.50d0
       ga(6,3) = -6038.7d0
       gb(6,3) = 5.0d0
       gc(6,3) = 2.0d0
       ecut(6) = garofalini(rcut(1)*rcut(1),6,qelect(2),qelect(3) ,2,3)

! Si-O-Si
       glambda(1) = 21732.3d0
       ggamma(1,1) = 2.0d0
       grij(1,1) = 2.6d0
       grijsq(1,1) = grij(1,1)*grij(1,1)
       gtheta(1) = dcos(109.5d0*degrad)

! O-Si-O
       glambda(2) = 1376379.0d0
       ggamma(2,1) = 2.8d0
       grij(2,1) = 3.0d0
       grijsq(2,1) = grij(2,1)*grij(2,1)
       gtheta(2) = dcos(109.5d0*degrad)

! H-O-H
       glambda(3) = 2535435.0d0
       ggamma(3,1) = 1.3d0
       grij(3,1) = 1.6d0
       grijsq(3,1) = grij(3,1)*grij(3,1)
       gtheta(3) = dcos(104.5d0*degrad)

! Si-O-H
       glambda(4) = 362205.0d0
       ggamma(4,1) = 2.0d0
       ggamma(4,2) = 1.2d0
       grij(4,1) = grij(1,1)
       grij(4,2) = 1.5d0
       grijsq(4,1) = grij(4,1)*grij(4,1)
       grijsq(4,2) = grij(4,2)*grij(4,2)
       gtheta(4) = dcos(109.5d0*degrad)

       do i=1,6
          write(io_output,*) 'garo ecut',i,ecut(i)
       end do
       return
    else if(lexpsix) then
! Explicit atom carbon Williams Exp-6 potential
! J.Chem.Phys 47 11 4680 (1967) paramter set IV
! note that the combination of C--H has to be the average of C
! and H the way it is set up in the energy subroutines
! natom is set in expsix.inc
! U(r) = A*r^(-6) + B*exp[C*r]
       do i = 1,natom
          aexsix(i) = 0.0d0
          bexsix(i) = 0.0d0
          cexsix(i) = 0.0d0
          qelect(i) = 0.0d0
          xiq(i) = 0.0d0
          lqchg(i) = .false.
          lij(i) = .true.
       end do

! O--- nonbonded interaction
! BKS [O]
       qelect(1) = -1.2
       lqchg(1) = .true.
       aexsix(1) = -175.0000d0*1.602177d-19/1.3806503d-23
       bexsix(1) = 1388.7730d0*1.602177d-19/1.3806503d-23
       cexsix(1) = -2.76000d0
       mass(1) = 15.9994

! O---Si nonbonded interaction
       aexsix(2) = -133.5381d0*1.602177d-19/1.3806503d-23
       bexsix(2) = 18003.7572d0*1.602177d-19/1.3806503d-23
       cexsix(2) = -4.87318d0
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
             consp(i) = (2.0d0/3.0d0)*onepi*(2.0d0*aexsix(i)/(rcut(1) *rcut(1)*rcut(1))+bexsix(i)*dexp(cexsix(i)*rcut(1)) *(-6.0d0/(cexsix(i)*cexsix(i)*cexsix(i))+6.0d0 *rcut(1)/(cexsix(i)*cexsix(i))-3.0d0*rcut(1)* rcut(1)/ cexsix(i)+rcut(1)*rcut(1)*rcut(1)))
             consu(i) = 2.0d0*onepi*(aexsix(i)/(3.0d0*rcut(1)*rcut(1)* rcut(1)) +(-rcut(1)*rcut(1)+2.0d0*rcut(1)/cexsix(i)-2.0d0/ (cexsix(i)* cexsix(i)))*bexsix(i)*dexp(cexsix(i)*rcut(1))/ cexsix(i))
          end do
! write(11,*) 'consp(i)',consp
! write(11,*) 'consu(i)',consu
       end if

       write(io_output,*)  ' i   aexsix       bexsix      cexsix     sexsix'
       do i = 1,natom
          write(io_output,'(i3,2x,4e12.4)')i,aexsix(i),bexsix(i) ,cexsix(i),sexsix(i)
       end do
       return
    else if (lmmff) then
! Merk Molecular Force Field (MMFF)
! J. Am. Chem. Soc. 1992 Vol.114 P7827-7843 (Thomas A. Halgren)
! natomtyp is set in mmff.inc
! U(r) = epsi*(1.07/(rs+0.07))^7 * (1.12/(rs^7+0.12)-2)
! rs = r / sigimmff

! C---C nonbonded interaction ***
       alphammff(1) = 1.050d0
       nmmff(1) = 2.490d0
       ammff(1) = 3.890d0
       gmmff(1) = 1.282d0
       sigimmff(1) = ammff(1)*dsqrt(dsqrt(alphammff(1)))
       epsimmff(1) = 45582.6d0*gmmff(1)*gmmff(1)*dsqrt(nmmff(1)) /(ammff(1)**6.0d0)
! sigimmff(1) = 1.126763255691509d0
! epsimmff(1) = 0.9994354715851470d0
       sigisq(1) = sigimmff(1)*sigimmff(1)
       mass(1) = 12.011d0
! mass(1) = 1.0d0

! H---H nonbonded interaction ***
       alphammff(3) = 0.250d0
       nmmff(3) = 0.800d0
       ammff(3) = 4.200d0
       gmmff(3) = 1.209d0
       sigimmff(3) = ammff(3)*dsqrt(dsqrt(alphammff(3)))
       epsimmff(3) = 45582.6d0*gmmff(3)*gmmff(3)*dsqrt(nmmff(3)) /(ammff(3)**6.0d0)
       sigisq(3) = sigimmff(3)*sigimmff(3)
       mass(3) = 1.0078d0

! C---H nonbonded interaction by using cominbination rule ***
       sigimmff(2) = 0.5d0*(sigimmff(1)+sigimmff(3))*(1.0d0 +  0.2d0*(1.0d0-dexp(-12.d0*(((sigimmff(1)-sigimmff(3))/ (sigimmff(3)+sigimmff(1)))**2.0d0))))
       epsimmff(2) = 91165.1d0*gmmff(1)*gmmff(3)*alphammff(1)* alphammff(3)/((dsqrt(alphammff(1)/nmmff(1))+ dsqrt(alphammff(3)/nmmff(3)))*(sigimmff(2)**6.0d0))
       sigisq(2) = sigimmff(2)*sigimmff(2)

       if (lshift) then
          do i=1,natomtyp
             rs1 = rcut(1)/sigimmff(i)
             rs7 = rs1**7.0d0
             sr7 = (1.07d0/(rs1+0.07d0))**7.0d0
             smmff(i) = epsimmff(i)*sr7*(1.12d0/(rs7+0.12)-2)
          end do
       else
          do i = 1,natomtyp
             smmff(i) = 0.0d0
          end do
          coru_cons(1) = -2.4837937263569310d-02
          ! coru_cons(1) = -5.8244592746534724E-03
          coru_cons(2) = -1.7583010189381791d-02
          coru_cons(3) = -8.3770412792126582d-03
          corp_cons(1) = 0.1696349613545569d0
          ! corp_cons(1) = 4.0098456560842058E-02
          corp_cons(2) = 0.1203650950025348d0
          corp_cons(3) = 5.7576802340310304d-02
       end if

       write(io_output,*) ' i   epsimmff     sigimmff   smmff'
       do i = 1,natom
          write(io_output,'(i3,2x,4e12.4)')i,epsimmff(i),sigimmff(i) ,smmff(i)
       end do
       return
    else if (lninesix) then
! special potential for all-atom formic acid from llnl 4/6/04 jms
! sp2 carbon site H-[C](=O)-OH
       rzeronx(1) = 4.1834161d0
       epsilonnx(1) = 45.224d0
       qelect(1) = 0.44469d0
       lqchg(1) = .true.
       mass(1) = 12.011d0
       chemid(1) = 'C'

! hydroxyl oxygen site H-C(=O)-[O]H
       rzeronx(2) = 3.5694293136d0
       epsilonnx(2) = 47.151d0
       qelect(2) = -0.55296d0
       lqchg(2) = .true.
       mass(2) = 15.9996d0
       chemid(2) = 'OH'

! carbonyl oxygen site H-C[=O]-OH
       rzeronx(3) = 3.0014635172d0
       epsilonnx(3) = 146.008d0
       qelect(3) = -0.43236d0
       lqchg(3) = .true.
       mass(3) = 15.9996d0
       chemid(3) = 'O='

! hydrogen site [H]-C(=O)-OH
       rzeronx(4) = 0.8979696387d0
       epsilonnx(4) = 2.4054d0
       qelect(4) = 0.10732d0
       lqchg(4) = .true.
       mass(4) = 1.00794d0
       chemid(4) = 'HC'

! acidic hydrogen site H-C(=O)-O[H]
       rzeronx(5) = 1.115727276d0
       epsilonnx(5) = 12.027d0
       qelect(5) = 0.43331d0
       lqchg(5) = .true.
       mass(5) = 1.00794d0
       chemid(5) = 'HO'

! calculate all site-site parameters via Lorentz-Berthelot rules
       do i = 1,nxatom
          do j = 1,nxatom
             ij = (i-1)*nxatom + j
             rzero(ij) = 0.5d0*(rzeronx(i) + rzeronx(j))
             epsnx(ij) = dsqrt(epsilonnx(i)*epsilonnx(j))
             if (lshift) then
                shiftnsix(ij) = 4.0d0*epsnx(ij)*(rzero(ij) /rcut(1))**6 * (2.0d0*(rzero(ij)/rcut(1))**3 - 3.0d0)
             end if
          end do
       end do
       return
    end if

! ! --- TraPPE-UA? Methane [CH4] sp3 charged with polarizability
! sigi(28) = 3.73d0
! epsi(28) = 148.0d0
! ! is this correct?
! mass(28) = 16.043d0
! qelect(28) = -0.572d0
! lqchg(28) = .true.
! jayself(28) = 0.5d0*117403d0
! xiq(28) = 9449.3d0
! chname(28) = 'Tr C CH4 chg pol '
! chemid(28)  = 'C  '

! ! --- Methane hydrogen charged with polarizibility
! sigi(29) = 0.0d0
! epsi(29) = 0.0d0
! mass(29) = 1.0078d0
! qelect(29) = 0.143d0
! lqchg(29) = .true.
! jayself(29) = 0.5d0*177700d0
! xiq(29) = 0.0d0
! lij(29) = .false.
! chname(29) = 'Tr H CH4 chg pol '
! chemid(29)  = 'H  '

! ! --- SPC-FQ oxygen [O]   S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(109) = 3.176
! epsi(109) = 148.0d0
! mass(109) = 15.999d0
! qelect(109) = -0.672123708
! lqchg(109) = .true.
! xiq(109) = 36899.0d0
! jayself(109) = (0.5d0)*(503.2d0)*(367.0d0)
! chname(109) = 'SPC-FQ O water   '
! chemid(109)  = '0  '

! ! --- SPC-FQ hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(110) = 0.0d0
! epsi(110) = 0.0d0
! mass(110) = 1.0079d0
! qelect(110) = 0.336061854
! lij(110) = .false.
! lqchg(110) = .true.
! xiq(110) = 0.0d0
! jayself(110) = (0.5d0)*(503.2d0)*(392.2d0)
! chname(110) = 'SPC-FQ H water   '
! chemid(110)  = 'H  '

! ! --- TIP4P-FQ Oxygen [O] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(111) = 3.159d0
! epsi(111) = 144.1d0
! !      epsi(111) = 105.0d0
! mass(111) = 15.999d0
! chname(111) = 'TIP4P-FQ O water '
! chemid(111)  = 'O  '

! ! --- TIP4P-FQ Hydrogen [H] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(112) = 0.0d0
! epsi(112) = 0.0d0
! mass(112) = 1.0079d0
! qelect(112) = 0.35d0
! lij(112) = .false.
! lqchg(112) = .true.
! xiq(112) = 0.0d0
! jayself(112) = (0.5d0)*(503.2d0)*(353.0d0)
! chname(112) = 'TIP4P-FQ H water '
! chemid(112)  = 'H  '

! ! --- TIP4P-FQ Charge [Q] S.W. Rick et al JCP 101 (7), 1 1994 6141
! sigi(113) = 0.0d0
! epsi(113) = 0.0d0
! mass(113) = 0.0d0
! qelect(113) = -0.70d0
! lij(113) = .false.
! lqchg(113) = .true.
! xiq(113) = 34464.0d0
! jayself(113) = (0.5d0)*(503.2d0)*(371.6d0)
! chname(113) = 'TIP4P-FQ M water '
! chemid(113)  = 'M  '

! ! --- TraPPE carbon dioxide carbon in [C]O2-fq (jpotoff 2/15/00)
! sigi(131) = 2.80d0
! epsi(131) = 28.5d0
! mass(131) = 12.011d0
! qelect(131) = 0.6512d0
! lqchg(131) = .true.
! !      xiq(131) = (503.2d0)*123.2d0
! xiq(131) = 0.0d0
! jayself(131) = (0.5d0)*(503.2d0)*(233.5d0)
! chname(131) = 'Tr-FQ C in CO2   '
! chemid(131)  = 'C  '

! ! --- TraPPE carbon dioxide oxygen in C[O]2-fq (jpotoff 2/15/00)
! sigi(132) = 3.06d0
! epsi(132) = 80.5d0
! mass(132) = 15.999d0
! qelect(132) = -0.3256d0
! lqchg(132) = .true.
! !      xiq(132) = (503.2d0)*201.56d0
! xiq(132) = 39430.75d0
! jayself(132) = (0.5d0)*(503.2d0)*(308.17d0)
! chname(132) = 'Tr-FQ O in CO2   '
! chemid(132)  = 'O  '

! ! - CO2-FQ Carbon-Oxygen cross term (JCO)
! i = 131
! j = 132
! djay = (503.2d0)*(133.905d0)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ij) = djay
! jayq(ji) = djay

! ! - CO2-FQ Oxygen-Oxygen cross term (JOO)
! i = 132
! j = 132
! djay = (503.2d0)*(1.09d0)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- SPC-FQ water Oxygen-Hydrogen cross term
! i = 109
! j = 110
! djay = (503.2d0)*(276.0d0)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ij) = djay
! jayq(ji) = djay

! ! --- SPC-FQ water Hydrogen-Hydrogen cross term
! i = 110
! j = 110
! djay = (503.2d0)*(196.0d0)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- TIP4P water Charge-Hydrogen cross term
! i = 112
! j = 113
! djay = (503.2d0)*(286.4d0)
! ij = (i-1)*nntype + j
! ji = (j-1)*nntype + i
! jayq(ji) = djay
! jayq(ij) = djay

! ! --- TIP4P water Hydrogen-Hydrogen cross term
! i = 112
! j = 112
! djay = (503.2d0)*(203.6d0)
! ij = (i-1)*nntype + j
! jayq(ij) = djay

! ! --- Methane C-H cross term
! i = 28
! j = 29
! ij = (i-1)*nntype + j
! ji = (j-i)*nntype + i
! jayq(ji) = 114855.0d0
! jayq(ij) = 114855.0d0

! ! --- Methane H-H cross term
! i = 29
! j = 29
! ij = (i-1)*nntype + j
! jayq(ij) = 112537.0d0

! Computation of un-like interactions

! convert input data to program units ***
    if ( lsami ) then
       call susami
       rcheck = 2.5d0 * 3.527d0
       if ( rcut(1) .ne. rcheck ) then
          write(io_output,*) 'WARNING ### rcut set to 2.5sigma for SAMI'
          rcut(1) = rcheck
       end if
    else
! calculate square sigmas and epsilons for lj-energy subroutines ***
       do i = 1, nntype
          do j = 1, nntype
             ij = (i-1)*nntype + j
             if (lmixlb) then
! Lorentz-Berthelot rules --- sig_ij = 0.5 [ sig_i + sig_j ]
                sig2ij(ij) =(0.5d0*(sigi(i)+sigi(j)))**2
                if (sigi(i).eq.0.0d0.or.sigi(j).eq.0.0d0) sig2ij(ij) = 0.0d0
             else if (lmixjo) then
! Jorgensen mixing rules --- sig_ij = [ sig_i * sig_j ]^(1/2)
                sig2ij(ij) = sigi(i) * sigi(j)
             end if
             epsij(ij) = dsqrt(epsi(i)*epsi(j))

             if (lshift) then
                sr2 = sig2ij(ij)/(rcut(1)*rcut(1))
                sr6 = sr2 * sr2 * sr2
                ecut(ij)= sr6*(sr6-1.0d0)*epsij(ij)
             end if
          end do
       end do
    end if

! Looking for section NONBOND
    REWIND(io_ff)
    CYCLE_READ_NONBOND:DO
       call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) then
          exit cycle_read_nonbond
       end if

       if (UPPERCASE(line_in(1:7)).eq.'NONBOND') then
          nmix=0
          do
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) then
                call err_exit(__FILE__,__LINE__,'Reading section NONBOND',jerr)
             end if
             if (UPPERCASE(line_in(1:11)).eq.'END NONBOND') exit
             nmix=nmix+1
             read(line_in,*) i,j,ji,sigmaTmp,epsilonTmp
             ij=(i-1)*nntype+j
             nonbond_type(ij)=ji
             nonbond_type((j-1)*nntype+i)=ji
             sig2ij(ij)=sigmaTmp*sigmaTmp
             sig2ij((j-1)*nntype+i)=sig2ij(ij)
             epsij(ij)=epsilonTmp
             epsij((j-1)*nntype+i)=epsilonTmp
             if (lshift) then
                sr2 = sig2ij(ij) / (rcut(1)*rcut(1))
                sr6 = sr2 * sr2 * sr2
                ecut(ij)= sr6*(sr6-1.0d0)*epsij(ij)
                ecut((j-1)*nntype+i)=ecut(ij)
             end if
          end do
          exit cycle_read_nonbond
       end if
    END DO CYCLE_READ_NONBOND

! set up the strectching and bending constants
    call init_energy_bonded(io_ff)

    close(io_ff)

    call init_energy_external(nntype)

    return
  end subroutine init_energy_pairwise

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

!     **************************************************************
! calculates the energy using the exp-6 potential       ***
! parameters defined in suijtab.f  M.G. Martin          ***
!     **************************************************************
  function exsix(rijsq,ntij)
      real::rijsq,rij,exsix
      integer::ntij

      rij=dsqrt(rijsq)
      exsix = aexsix(ntij)/(rijsq*rijsq*rijsq) + bexsix(ntij)*dexp(cexsix(ntij)*rij)
      if (lshift) exsix = exsix-sexsix(ntij)
      return
  end function exsix

!     ***************************************************
! calculates the energy of the 9-6 potential ***
! parameters defined in suijtab.f  JMS       ***
!     ***************************************************
  function ninesix(rijsq,ntij)
      real::rijsq,rij,ror,ninesix
      integer::ntij

      rij=dsqrt(rijsq)
      ror = rzero(ntij)/rij
      ninesix = 4.0d0*epsnx(ntij)*ror**6*(2.0d0*ror**3 - 3.0d0)
      if (lshift) then
         ninesix = ninesix - shiftnsix(ntij)
      end if

      return
  end function ninesix

!  ********************************************************************
! calculates the energy of the Generalized Lennard Jones
! potential                                               ***
! parameters defined  poten.inc           ***
!  *******************************************************************
! Also you should make sure that the rcut you choose
! is .gt. sigma*(2**(2/n0)) (because currently tail corrections
! are evaluated with this assumption which is rarely incorrect.)
!  ********************************************************************
  function genlj (rijsq,sr2,epsilon2)
      real::rijsq,rij,srij,sr2,epsilon2,genlj

      rij=(dsqrt(rijsq))
      srij=dsqrt(sr2)

      if ( (rij) .le.(rij*srij)*2.0d0**(2.0d0/n0) ) then
         genlj = 4.0d0*epsilon2*(((srij)**n0)-((srij)**(n0/2.0d0)))
      else
         genlj =epsilon2*(((2.0d0**((4.0d0*n1/n0)))*((srij)** (2.0d0*n1)))-((2.0d0**((2.0d0*(n1/n0))+1.0d0))*((srij)**(n1))))
      end if

! In reduced units
! xij=(1.0d0/dsqrt(sr2))
!
! if ( (xij) .le. 2.0d0**(2.0d0/n0) ) then
! genlj = 4.0d0 * (((1/xij)**n0)-((1/xij)**(n0/2.0d0)))
! else
! genlj =(( (2.0d0**((4.0d0*n1/n0)))*(1/xij)**(2.0d0*n1))-
!     & (2.0d0**((2.0d0*(n1/n0))+1.0d0)*((1/xij)**(n1))))
! end if
! vinter=vinter+e

! not usually used in MC
! if (lshift) then
! genlj = genlj - shiftgenlj()
! end if

      return
  end function genlj

!    *********************************************************
! calculates energy for a polymeric surfactant bead.  **
!    *********************************************************
  function ljpsur ( rijsq, ntij )

      real::ljpsur, rijsq, sr, sr6
      integer::ntij

! --------------------------------------------------------------------
! AT PRESENT: all sigma = 1.0
! epsilon   = 1.0
! --------------------------------------------------------------------

      sr = 1.0d0 / rijsq
      sr6 = sr**3

      if ( ntij .eq. 1 ) then
! nonpolar-nonpolar interaction ( LJ interaction )
         ljpsur = sr6*sr6 - sr6
      else
! polar-polar or polar-nonpolar interaction ( repulsive LJ interaction )
         if ( rijsq .le. 1.259921d0 ) then
            ljpsur = sr6*sr6 - sr6 + 0.25d0
         else
            ljpsur = 0.0d0
         end if
      end if

      return
  end function ljpsur

!   ******************************************************************
! calculates the energy using the Buffered 14-7 potential   ***
! parameters defined in suijtab.f          Bin Chen         ***
!   ******************************************************************
  function mmff(rijsq,ntij)
      real::rijsq,mmff,rs2,rs1,sr1,sr7,rs7
      integer::ntij

        rs2 = rijsq / (sigisq(ntij))
        rs1 = dsqrt(rs2)
        rs7 = rs1*rs2*rs2*rs2
        sr1 = 1.07d0/(rs1+0.07d0)
        sr7 = sr1**7.0d0
        mmff = sr7*(1.12d0/(rs7+0.12d0) - 2.0d0) *epsimmff(ntij)

      if (lshift) mmff = mmff-smmff(ntij)
      return
  end function mmff

!DEC$ ATTRIBUTES FORCEINLINE :: type_2body
  function type_2body(ntii,ntjj)
    integer::type_2body
    integer,intent(in)::ntii,ntjj

    if (lexpsix .or. lmmff) then
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

    U2=0.0d0
    if (L_vdW_table) then
       U2=lininter_vdW(rij,ntii,ntjj)
    else if ( lsami ) then
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
       if ( lexpand(imolty).and.lexpand(jmolty)) then
          sigma2=(sigma_f(imolty,ii)+sigma_f(jmolty,jj))/2.0d0
          sr2 = sigma2*sigma2/rijsq
          epsilon2=dsqrt(epsilon_f(imolty,ii)*epsilon_f(jmolty,jj))
       else if (lexpand(imolty)) then
          sigma2=(sigma_f(imolty,ii)+ sigi(ntjj))/2.0d0
          sr2 = sigma2*sigma2/rijsq
          epsilon2=dsqrt(epsilon_f(imolty,ii)*epsi(ntjj))
       else if (lexpand(jmolty)) then
          sigma2=(sigma_f(jmolty,jj)+ sigi(ntii))/2.0d0
          sr2 = sigma2*sigma2/rijsq
          epsilon2=dsqrt(epsi(ntii)*epsilon_f(jmolty,jj))
       else
          sr2 = sig2ij(ntij) / rijsq
          epsilon2=epsij(ntij)
       end if

       if (lfepsi.and.i.ne.j) then ! NOT for intramolecular interactions
          sr6 = rijsq*rijsq*rijsq
          if ((.not.lqchg(ntii)).and.(.not.lqchg(ntjj))) then
             if (nunit(imolty).eq.4) then
                ! TIP-4P structure (temperary use !??)
                qave=(qqu(i,4)+qqu(j,4))/2.0d0
             else
                qave=(qqu(i,4)+qqu(i,5)+qqu(j,4)+qqu(j,5))*0.85d0
             end if
          else
             qave=(qqu(i,ii)+qqu(j,jj))/2.0d0
          end if
          U2=((aslope*(qave-a0)*(qave-a0)+ashift)/sr6-(bslope*(qave- b0)*(qave-b0)+bshift))/sr6*epsilon2
       else
          sr6 = sr2 * sr2 * sr2
          U2=sr6*(sr6-1.0d0)*epsilon2
          if (lshift) then
             U2=U2-ecut(ntij)
          end if
          if (i.eq.j) then
             ! intramolecular interactions
             U2=U2*ljscale(imolty,ii,jj)
             if (lainclu(imolty,ii,jj)) then
                ! OH 1-5 interaction
               U2=U2+0.25d0*a15(a15type(imolty,ii,jj))/((rijsq**2)*(rijsq**2)*(rijsq**2))
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

    Q2=0.0d0
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
       elseif (lchgall.or.rijsq.lt.rcutsq) then
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
!ccc  Calculates nonbonding van der Waals and electrostatic potential
!ccc  using linear interpolation between two points
!ccc  requires file: fort.43 -- vdW, fort.44 -- Q
!ccc    fort.43: number of tabulated potentials, potential number,
!ccc  number of points per angstrom, tabulated potential
!ccc  (repeat last three parts for each additional potential)
!ccc  for unlike interactions, list 1-2 and 2-1
!ccc  separate potentials with 1000
!ccc  make sure potential does not go up to infinity!
!ccc  bead type numbers should be defined in suijtab, but it doesn't
!ccc  matter what they are (the parameters aren't used)
!ccc  KM 12/03/08
!ccc    fort.44: number of tabulated potentials, potential number,
!ccc  number of points per angstrom, tabulated potential
!ccc  (repeat last three parts for each additional potential)
!ccc  for unlike interactions, list 1-2 and 2-1
!ccc  separate potentials with 1000
!ccc  KM 04/23/09
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine init_tabulated_potential_pair()
    use sim_system,only:L_vdW_table,L_elect_table
    integer,parameter::initial_size=10,grid_size=1500
    integer::jerr

    if (L_vdW_table) then
       allocate(vdWsplits(1:initial_size,1:initial_size),rvdW(1:grid_size,1:initial_size,1:initial_size),tabvdW(1:grid_size,1:initial_size,1:initial_size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_pair: allocation failed for vdW_table',myid+1)
       call read_tabulated_potential_pair('fort.43',ntabvdW,rvdW,tabvdW,vdWsplits,atoms)
       if (myid.eq.0) write(io_output,*) 'using linear interpolation for nonbonded van der Waals interactions'
    end if

    if (L_elect_table) then
       allocate(electsplits(1:initial_size,1:initial_size),relect(1:grid_size,1:initial_size,1:initial_size),tabelect(1:grid_size,1:initial_size,1:initial_size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_pair: allocation failed for elect_table',myid+1)
       call read_tabulated_potential_pair('fort.44',ntabelect,relect,tabelect,electsplits,atoms)
       if (myid.eq.0) write(io_output,*) 'using linear interpolation for electrostatic interactions'
    end if
  end subroutine init_tabulated_potential_pair
end MODULE energy_pairwise
