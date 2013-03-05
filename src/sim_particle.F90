MODULE sim_particle
  use var_type,only:dp,default_string_length
  implicit none
  private
  public::BeadType,AtomType,MoleculeType,init_neighbor_list,rebuild_neighbor_list,update_neighbor_list,check_neighbor_list,lnn,allocate_neighbor_list

  type BeadType
     integer::type
     real::coord(3)
  end type BeadType

  type AtomType
     logical::lqchg,llj
     integer::atomID
     real::mass,charge,eps,sig
     character(LEN=default_string_length)::name,desc !name is the atomic symobl, desc is the full description of the atom
  end type AtomType

  type MoleculeType
     integer::nbead
     type(BeadType),allocatable::bead(:)
  end type MoleculeType

! NEIGH.INC
  logical,allocatable::lnn(:,:)
  real,allocatable::disvec(:,:,:),upnn(:),upnnsq(:),upnndg(:,:)

contains
  subroutine check_neighbor_list(ibox,imol)
    use sim_system,only:rmtrax,rmtray,rmtraz,rmrotx,rmroty,rmrotz,Armtrax,Armtray,Armtraz,io_output
    integer,intent(in)::ibox,imol

    if ( rmtrax(imol,ibox) .gt. upnn(ibox) ) then
       write(io_output,*) ' rmtrax greater than upnn',ibox,imol
       rmtrax(imol,ibox) = upnn(ibox)
    end if
    if ( rmtray(imol,ibox) .gt. upnn(ibox) ) then
       write(io_output,*) ' rmtray greater than upnn',ibox,imol
       rmtray(imol,ibox) = upnn(ibox)
    end if
    if ( rmtraz(imol,ibox) .gt. upnn(ibox) ) then
       write(io_output,*) ' rmtraz greater than upnn',ibox,imol
       rmtraz(imol,ibox) = upnn(ibox)
    end if

    if ( rmrotx(imol,ibox) .gt. upnndg(imol,ibox) ) then
       write(io_output,*) ' rmrotx greater than upnndg',ibox,imol
       rmrotx(imol,ibox) = upnndg(imol,ibox)
    end if
    if ( rmroty(imol,ibox) .gt. upnndg(imol,ibox) ) then
       write(io_output,*) ' rmroty greater than upnndg',ibox,imol
       rmroty(imol,ibox) = upnndg(imol,ibox)
    end if
    if ( rmrotz(imol,ibox) .gt. upnndg(imol,ibox) ) then
       write(io_output,*) ' rmrotz greater than upnndg',ibox,imol
       rmrotz(imol,ibox) = upnndg(imol,ibox)
    end if

    if (Armtrax .gt. upnn(ibox)) then
       Armtrax = upnn(ibox)
       write(io_output,*) '### problem : for target accept ', 'ratio Armtrax should be smaller than upnn'
    end if
    if (Armtray .gt. upnn(ibox)) then
       Armtray = upnn(ibox)
       write(io_output,*) '### problem : for target accept ', 'ratio Armtray should be smaller than upnn'
    end if
    if (Armtraz .gt. upnn(ibox)) then
       Armtraz = upnn(ibox)
       write(io_output,*) '### problem : for target accept ', 'ratio Armtraz should be smaller than upnn'
    end if

  end subroutine check_neighbor_list

  subroutine init_neighbor_list()
    use util_string,only:integer_to_string
    use util_runtime,only:err_exit
    use sim_system,only:nbox,rcut,rcutnn,nmolty,nunit,brvib,io_output
    integer::ibox,imol,j
    real::umatch

! If using neighbour list make sure the rcut & rcutnn is the same
! for all the boxes
    do ibox = 1,nbox
       if ((abs(rcut(1)-rcut(ibox)).gt.1.0E-10_dp).and.(abs(rcutnn(1)-rcutnn(ibox)).gt.1.0E-10_dp)) then
          call err_exit(__FILE__,__LINE__,'Keep rcut and rcutnn for all the boxes same',-1)
       end if

       if (rcut(ibox).ge.rcutnn(ibox)) then
          call err_exit(__FILE__,__LINE__,' rcut greater equal rcutnn for box'//integer_to_string(ibox),-1)
       end if
    end do

! set logical map to .false.
    lnn = .false.
! set displacementvectors to zero
    disvec = 0.0E0_dp

! calculate max. angular displacement that doesn't violate upnn
! calculate max. all-trans chain length ( umatch )
    do ibox=1,nbox
       upnn(ibox) = ( rcutnn(ibox) - rcut(ibox) ) / 3.0E0_dp
       upnnsq(ibox)=upnn(ibox)*upnn(ibox)
       do imol=1,nmolty
          umatch = 0.0E0_dp
          do j = 1, nunit(imol) - 1
             umatch = umatch + brvib(1)
          end do
          upnndg(imol,ibox) = asin( upnn(ibox) / umatch )

          call check_neighbor_list(ibox,imol)
       end do

       call rebuild_neighbor_list(ibox)
    end do

    write(io_output,*)
    write(io_output,*) 'Neighbour list created.'
  end subroutine init_neighbor_list

!> \brief Rebuild neighbour list for ibox
  subroutine rebuild_neighbor_list(ibox)
    use sim_system,only:lpbc,nchain,nboxi,nunit,rxu,ryu,rzu,rcutnn
    use sim_cell,only:setpbc,mimage
    integer,intent(in)::ibox
    integer::i,j,ii,jj
    real::rxui,ryui,rzui,rxuij,ryuij,rzuij,rijsq,rcnnsq

    do i=1,nchain
       if (nboxi(i).eq.ibox) then
! set i-part of logical::map to .false.
          lnn(i,:)=.false.
          lnn(:,i)=.false.
! set displacementvectors to zero
          disvec(:,:,i)=0.0E0_dp
       end if
    end do

    rcnnsq = rcutnn(ibox)**2
    if (lpbc) call setpbc(ibox)
 
! loop over all chains i 
    do i = 1, nchain-1
       if (nboxi(i).ne.ibox) cycle
! loop over all chains j
       molecule2:do j = i+1, nchain
          if (nboxi(j).ne.ibox) cycle
! loop over all beads ii of chain i 
          do ii = 1, nunit(i)
             rxui = rxu(i,ii)
             ryui = ryu(i,ii)
             rzui = rzu(i,ii)
! loop over all beads jj of chain j 
             do jj = 1, nunit(j)
                rxuij = rxui - rxu(j,jj)
                ryuij = ryui - ryu(j,jj)
                rzuij = rzui - rzu(j,jj)
! minimum image the pair separations
                if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
                rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                if ( rijsq .le. rcnnsq ) then
                   lnn(i,j) = .true.
                   lnn(j,i) = .true.
                   cycle molecule2
                end if
             end do
          end do
       end do molecule2
    end do
 
    return
  end subroutine rebuild_neighbor_list

  subroutine update_neighbor_list(i,rx,ry,rz,force_update)
    use sim_system,only:rcutnn,lpbc,nchain,nboxi,nunit,rxu,ryu,rzu
    use sim_cell,only:setpbc,mimage
    implicit none
    integer,intent(in)::i
    real,intent(in)::rx,ry,rz
    logical,intent(in)::force_update
    integer::j,ii,jj,ibox
    real::rxui,ryui,rzui,rxuij,ryuij,rzuij,rijsq,rcnnsq,disvsq1,disvsq2
 
    ibox=nboxi(i)

! check for update of near neighbour bitmap
! check for headgroup
    disvec(1,1,i) = disvec(1,1,i) + rx
    disvec(2,1,i) = disvec(2,1,i) + ry
    disvec(3,1,i) = disvec(3,1,i) + rz
    disvsq1 = disvec(1,1,i) * disvec(1,1,i) + disvec(2,1,i) * disvec(2,1,i) + disvec(3,1,i) * disvec(3,1,i)
! check for last unit
    disvec(1,2,i) = disvec(1,2,i) + rx
    disvec(2,2,i) = disvec(2,2,i) + ry
    disvec(3,2,i) = disvec(3,2,i) + rz
    disvsq2 = disvec(1,2,i) * disvec(1,2,i) + disvec(2,2,i) * disvec(2,2,i) + disvec(3,2,i) * disvec(3,2,i)

    if ((disvsq1.le.upnnsq(ibox)).and.(disvsq2.le.upnnsq(ibox)).and.(.not.force_update)) return

! set i-part of logical::map to .false.
    lnn(i,:)=.false.
    lnn(:,i)=.false.
! set displacementvectors to zero
    disvec(:,:,i)=0.0E0_dp

    rcnnsq = rcutnn(ibox)**2
    if (lpbc) call setpbc(ibox)

! loop over all chains j
    molecule2:do j = 1, nchain
       if ((i.eq.j).or.(nboxi(j).ne.ibox)) cycle
! loop over all beads ii of chain i 
       do ii = 1, nunit(i)
          rxui = rxu(i,ii)
          ryui = ryu(i,ii)
          rzui = rzu(i,ii)
! loop over all beads jj of chain j 
          do jj = 1, nunit(j)
             rxuij = rxui - rxu(j,jj)
             ryuij = ryui - ryu(j,jj)
             rzuij = rzui - rzu(j,jj)
! minimum image the pair separations
             if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox )
             rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
             if ( rijsq .le. rcnnsq ) then
                lnn(i,j) = .true.
                lnn(j,i) = .true.
                cycle molecule2
             end if
          end do
       end do
    end do molecule2
 
    return
  end subroutine update_neighbor_list

  subroutine allocate_neighbor_list
    use util_runtime,only:err_exit
    use sim_system,only:nbxmax,ntmax,nmax,io_output
    integer::jerr
    allocate(lnn(nmax,nmax),disvec(3,2,nmax),upnn(nbxmax),upnnsq(nbxmax),upnndg(ntmax,nbxmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_neighbor_list: allocation failed',jerr)
    end if
  end subroutine allocate_neighbor_list
end MODULE sim_particle
