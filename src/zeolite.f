module zeolite
  use global_data
  use var_type
  use const_math
  use const_phys
  use util_math
  use util_timings, only: time_now
  use util_files
  use util_string
  use sim_cell
  use sim_particle
  implicit none
  include 'common.inc'
  private
  save
  public::zeocoord,suzeo,exzeo

  type ZeoliteUnitCellGridType
     integer::dup(3),ngrid(3) ! dup(:) how many times the unit cell is replicated to form the simulation cell, ngrid(:) number of grid points in each direction
     real(KIND=double_precision)::boxl(3),hmat(3,3),hmati(3,3) ! boxl(:) unit cell length parameters, the angles are the same as for the simulation cell
  end type ZeoliteUnitCellGridType

  type ZeoliteBeadType
     integer::ntype ! number of framework atom types
     integer,allocatable::type(:),num(:) ! index (as in suijtab) and number of atoms of each type
     real(KIND=double_precision),allocatable::radiisq(:) ! square of protective radii of each type
     character(LEN=default_string_length),allocatable::name(:) ! chemical name of each type
  end type ZeoliteBeadType

  type ZeolitePotentialType
     integer::ntype ! number of guest bead types present
     integer,allocatable::table(:) ! bead type of each bead
     real(KIND=double_precision),allocatable::param(:,:,:) ! param(1,i,j)=A=4 eps sig^12, param(2,i,j)=B=4 eps sig^6
  end type ZeolitePotentialType

  real(KIND=double_precision),parameter::dgr=0.2,required_precision=1.0E-2_double_precision,overflow=1.0E+8_double_precision,upperlimit=1.0E+7_double_precision
  integer,parameter::boxZeo=1

  logical,allocatable::lunitcell(:)
  integer(KIND=normal_int)::nlayermax
  real(KIND=double_precision),allocatable::zgrid(:,:,:,:,:),egrid(:,:,:,:)

  type(CellMaskType)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
  type(MoleculeType)::zeo
  type(ZeoliteUnitCellGridType)::zunit
  type(ZeoliteBeadType)::ztype
  type(ZeolitePotentialType)::zpot
  
contains
! Interface to old code using mask data types
  subroutine initZeo()
    integer::i,j
    zcell%solid=>lsolid(boxZeo)
    zcell%ortho=>lrect(boxZeo)
    zcell%pbc(1)=lpbcx
    zcell%pbc(2)=lpbcy
    zcell%pbc(3)=lpbcz
    zcell%cut=>rcut(boxZeo)
    zcell%vol=>cell_vol(boxZeo)
    zcell%calp=>calp(boxZeo)
    zcell%boxl(1)%val=>boxlx(boxZeo)
    zcell%boxl(2)%val=>boxly(boxZeo)
    zcell%boxl(3)%val=>boxlz(boxZeo)
    do i=1,3
       zcell%ang(i)%val=>cell_ang(boxZeo,i)
       zcell%height(i)%val=>min_width(boxZeo,i)
       do j=1,3
          zcell%hmat(j,i)%val=>hmat(boxZeo,3*(i-1)+j)
          zcell%hmati(j,i)%val=>hmati(boxZeo,3*(i-1)+j)
       end do
    end do
  end subroutine initZeo

  subroutine zeocoord(file_zeocoord,lhere)
    character(LEN=*),intent(in)::file_zeocoord
    logical,intent(out)::lhere(:)
    integer(KIND=normal_int)::io_zeocoord,jerr,pos,i,j
    real(KIND=double_precision)::wzeo,scoord(3)
    character(LEN=4)::atom

    call initZeo()

    io_zeocoord=get_iounit()
    open(unit=io_zeocoord,access='sequential',action='read',file=file_zeocoord,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call cleanup('cannot open zeolite coordinate file')
    end if

    read(io_zeocoord,*) zcell%boxl(1)%val,zcell%boxl(2)%val,zcell%boxl(3)%val,zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val,zunit%dup(1),zunit%dup(2),zunit%dup(3)

    if (myid.eq.0) write(iou,"(/,' READING ZEOLITE LATTICE FROM FILE zeolite.cssr:',/&
                                ,' --------------------------------------------------',/&
                                ,' box dimensions                    = ',3f10.3,' Angstrom',/&
                                ,' box angles                        = ',3f10.3,' degrees',/&
                                ,' number of zeolite cells           = ',3i5)") zcell%boxl(1)%val,zcell%boxl(2)%val,zcell%boxl(3)%val,zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val,zunit%dup(1),zunit%dup(2),zunit%dup(3)

    if (abs(zcell%ang(1)%val-90)>eps.or.abs(zcell%ang(2)%val-90)>eps.or.abs(zcell%ang(3)%val-90)>eps) then
       zcell%solid=.true.
       zcell%ortho=.false.
    else
       zcell%solid=.false.
    end if

! Fortran intrinsic trigonometric functions takes radian as input
    forall(i=1:3)
       zcell%ang(i)%val=zcell%ang(i)%val*onepi/180
       zunit%boxl(i)=zcell%boxl(i)%val/zunit%dup(i)
    end forall

     !     --- find closest values for x and y
     zunit%ngrid=int(zunit%boxl/dgr)

! align crystallography axis a with cartesian axis x, b in x-y plane
! hmat is the transformation matrix from cartesian cordinate scoord(2)stem (x,y,z) to (a,b,c)
! hmat=(1,4,7;2,5,8;3,6,9)
    zcell%hmat(1,1)%val=zcell%boxl(1)%val
    zcell%hmat(2,1)%val=0.
    zcell%hmat(3,1)%val=0.
    zcell%hmat(1,2)%val=zcell%boxl(2)%val*cos(zcell%ang(3)%val)
    zcell%hmat(2,2)%val=zcell%boxl(2)%val*sin(zcell%ang(3)%val)
    zcell%hmat(3,2)%val=0.
    zcell%hmat(1,3)%val=zcell%boxl(3)%val*cos(zcell%ang(2)%val)
    zcell%hmat(2,3)%val=zcell%boxl(3)%val*sqrt(cos(zcell%ang(2)%val)**2*cos(zcell%ang(3)%val)**2+cos(zcell%ang(1)%val)**2-2*cos(zcell%ang(1)%val)*cos(zcell%ang(2)%val)*cos(zcell%ang(3)%val))/sin(zcell%ang(3)%val)
    zcell%hmat(3,3)%val=zcell%boxl(3)%val*sqrt(1-cos(zcell%ang(1)%val)**2-cos(zcell%ang(2)%val)**2-cos(zcell%ang(3)%val)**2+2*cos(zcell%ang(1)%val)*cos(zcell%ang(2)%val)*cos(zcell%ang(3)%val))/sin(zcell%ang(3)%val)

! there might be numeric errors in the above trigonometric calculations
    forall(i=1:3,j=1:3,abs(zcell%hmat(j,i)%val).lt.eps) hmat(j,i)=0

    call matops(boxZeo)

    forall(i=1:3,j=1:3)
       zunit%hmat(j,i)=zcell%hmat(j,i)%val/zunit%dup(i)
       zunit%hmati(i,j)=zcell%hmati(i,j)%val*zunit%dup(i)
    end forall

    if (myid.eq.0) write(iou,"(' simulation volume                 = ', f10.1, 20x,' Angst**3')") zcell%vol

    read(io_zeocoord,*) zeo%nbead,ztype%ntype
    allocate(zeo%bead(zeo%nbead),lunitcell(zeo%nbead),ztype%name(ztype%ntype),ztype%radiisq(ztype%ntype),ztype%type(ztype%ntype),ztype%num(ztype%ntype),stat=jerr)
    if (jerr.ne.0) call cleanup('zeocoord: allocation failed')

    if (myid.eq.0) write(iou,"(' number of zeolite atoms           = ',i10,/&
                              ,' number of atomtypes in the lattice= ',i10,/)") zeo%nbead,ztype%ntype
    do i=1,ztype%ntype
       read(io_zeocoord,*) ztype%name(i),ztype%type(i),ztype%radiisq(i)
       lhere(ztype%type(i))=.true.
       ztype%radiisq(i)=ztype%radiisq(i)*ztype%radiisq(i)
    end do

    !     === Converting to absolute coordinates within [0,ri>
    ztype%num=0
    do i = 1,zeo%nbead
       read(io_zeocoord,'(i5,1x,a4,2x,3(f9.5,1x))') pos,atom,zeo%bead(i)%coord(1),zeo%bead(i)%coord(2),zeo%bead(i)%coord(3)

       scoord(1) = zeo%bead(i)%coord(1)*zcell%hmati(1,1)%val+zeo%bead(i)%coord(2)*zcell%hmati(1,2)%val+zeo%bead(i)%coord(3)*zcell%hmati(1,3)%val
       scoord(2) = zeo%bead(i)%coord(1)*zcell%hmati(2,1)%val+zeo%bead(i)%coord(2)*zcell%hmati(2,2)%val+zeo%bead(i)%coord(3)*zcell%hmati(2,3)%val
       scoord(3) = zeo%bead(i)%coord(1)*zcell%hmati(3,1)%val+zeo%bead(i)%coord(2)*zcell%hmati(3,2)%val+zeo%bead(i)%coord(3)*zcell%hmati(3,3)%val
       scoord = scoord - floor(scoord)
       zeo%bead(i)%coord(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
       zeo%bead(i)%coord(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
       zeo%bead(i)%coord(3)=scoord(3)*zcell%hmat(3,3)%val

       if (ALL(scoord*zunit%dup.lt.1)) then
          lunitcell(i)=.true.
       else
          lunitcell(i)=.false.
       end if

       pos=str_search(ztype%name,ztype%ntype,atom)
       if (pos.eq.0) call cleanup('** atomtype: unknown atomtype **')

       zeo%bead(i)%type=pos
       
       ztype%num(pos)=ztype%num(pos)+1

       ! if (myid.eq.0) write(iou,*) zeo%bead(i)%coord(1),zeo%bead(i)%coord(2),zeo%bead(i)%coord(3),zeo%bead(i)%type
    end do

!     do i=1,ztype%ntype
!        if (ztype%num(i).eq.0) then
!           ztype%type(i)=ztype%type(ztype%ntype)
!           ztype%num(i)=ztype%num(ztype%ntype)
!           ztype%radiisq(i)=ztype%radiisq(ztype%ntype)
!           ztype%name(i)=ztype%name(ztype%ntype)
!           ztype%ntype=ztype%ntype-1
!           if (i.ge.ztype%ntype) exit
!        end if
!     end do
    !     === Calculate zeolite density
    wzeo=dot_product(ztype%num(1:ztype%ntype),mass(ztype%type(1:ztype%ntype)))
    if (myid.eq.0) write(iou,"(' Tabulated zeolite potential: ',/&
                              ,' --------------------------------------------------',/&
                              ,<ztype%ntype>(' number of ',A,':',I5,3X),/&
                              ,' mass zeolite                      = ',e12.5,' grams',/ &
                              ,' one adsorbed molecule in sim box  = ',e12.5 ,' mmol/gram',/&
                              ,' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/&
                              ,'         x-dir           : ',i12,'  size: ',f7.4,/&
                              ,'         y-dir           : ',i12,'  size: ',f7.4,/&
                              ,'         z-dir           : ',i12,'  size: ',f7.4,/)") (trim(ztype%name(i)),ztype%num(i),i=1,ztype%ntype),wzeo/6.023e23,1000.0/wzeo,zunit%boxl(1),zunit%boxl(2),zunit%boxl(3),zunit%ngrid(1),zunit%boxl(1)/zunit%ngrid(1),zunit%ngrid(2),zunit%boxl(2)/zunit%ngrid(2),zunit%ngrid(3),zunit%boxl(3)/zunit%ngrid(3)

    close(io_zeocoord)

  end subroutine zeocoord

  subroutine addAllBeadTypes()
    integer::imol,iunit,igtype,i,idi,idj,ntij,jerr,list(nntype)
    real(KIND=double_precision)::sig6

    ! find all bead types and store them in an array
    zpot%ntype=0
    do imol=1,nmolty
       do iunit=1,nunit(imol)
          do igtype=1,zpot%ntype
             if (list(igtype).eq.ntype(imol,iunit)) exit
          end do
          if (igtype.gt.zpot%ntype) then
             zpot%ntype=igtype
             list(igtype)=ntype(imol,iunit)
          end if
       end do
    end do    
    
    allocate(zgrid(3,ztype%ntype,0:zunit%ngrid(1)-1,0:zunit%ngrid(2)-1,0:zunit%ngrid(3)-1),egrid(0:zunit%ngrid(1)-1,0:zunit%ngrid(2)-1,0:zunit%ngrid(3)-1,zpot%ntype),zpot%param(3,ztype%ntype,zpot%ntype),zpot%table(zpot%ntype),stat=jerr)
    if (jerr.ne.0) call cleanup('addAllBeadTypes: allocation failed')

    do igtype=1,zpot%ntype
       idi=list(igtype)
       zpot%table(igtype)=idi
       if (lij(idi).or.lqchg(idi)) then
          do i=1,ztype%ntype
             idj=ztype%type(i)
             if (lij(idi).and.lij(idj)) then
                ntij = (idi-1)*nntype + idj
                ! LJ parameters are scaled to reduce round-off errors, see exzeof()
                sig6=sig2ij(ntij)**3/2e4
                zpot%param(2,i,igtype)=epsij(ntij)*sig6
                zpot%param(1,i,igtype)=zpot%param(2,i,igtype)*sig6
             else
                zpot%param(1:2,i,igtype)=0.
             end if
             if (lqchg(idi).and.lqchg(idj)) then
                zpot%param(3,i,igtype)=qelect(idi)*qelect(idj)
             else
                zpot%param(3,i,igtype)=0.
             end if
          end do
       end if
    end do
    zpot%param(1:2,:,:)=4.*zpot%param(1:2,:,:)
  end subroutine addAllBeadTypes

  subroutine suzeo(file_ztb)
    character(LEN=default_path_length),intent(in)::file_ztb
    character(LEN=default_string_length)::atom
    integer(KIND=normal_int)::io_ztb,igtype,idi,idj,jerr,sw,i,j,k,oldi,oldj,oldk,ngridtmp(3),zntypetmp
    real(KIND=double_precision)::zunittmp(3),zangtmp(3)
    logical::lewaldtmp,ltailctmp,lshifttmp

    ! === tabulation of the zeolite potential
    if (lzgrid) then
       call addAllBeadTypes()
       call setpbc(boxZeo)

       nlayermax=0
       io_ztb=get_iounit()

       open(unit=io_ztb,access='sequential',action='read',file=file_ztb,form='binary',iostat=jerr, status='old')
       if (jerr.eq.0) then
          ! --- read zeolite table from disk
          if (myid.eq.0) write(iou,*) 'read in tabulated potential'
          read(io_ztb) zunittmp,zangtmp,ngridtmp,zntypetmp,lewaldtmp,ltailctmp,lshifttmp
          if (ANY(abs(zunittmp-zunit%boxl).gt.eps).or.ANY(ngridtmp.ne.zunit%ngrid).or.(zntypetmp.ne.ztype%ntype)) call cleanup('problem 1 in zeolite potential table')
          do igtype=1,ztype%ntype
             read(io_ztb) atom
             if (trim(ztype%name(igtype))/=trim(atom)) then
                write(iou,*) igtype,' atom should be ',trim(atom)
                call cleanup('problem 2 in zeolite potential table')
             end if
          end do
          do k=0,zunit%ngrid(3)-1
             do j=0,zunit%ngrid(2)-1
                do i=0,zunit%ngrid(1)-1
                   do igtype=1,ztype%ntype
                      do sw=1,3
                         read(io_ztb,end=100) zgrid(sw,igtype,i,j,k)
                      end do
                   end do
                end do
             end do
          end do
       else
          oldi=0
          oldj=0
          oldk=0
          lewaldtmp=lewald
          ltailctmp=ltailc
          lshifttmp=lshift
       end if

100    if (jerr.eq.0.and.(sw.le.3..or.igtype.le.ztype%ntype.or.i.lt.zunit%ngrid(1).or.j.lt.zunit%ngrid(2).or.k.lt.zunit%ngrid(3))) then
          oldi=i
          oldj=j
          oldk=k
          jerr=1
          close(io_ztb)
          write(iou,*) sw,igtype,oldi,oldj,oldk
       end if

       if (jerr.ne.0) then
          ! make a tabulated potential of the zeolite
          if (myid.eq.0) then
             write(iou,*) 'make tabulated potential'
             open(unit=io_ztb,access='sequential',action='write',file=file_ztb,form='binary',iostat=jerr,status='unknown')
             if (jerr.ne.0) then
                call cleanup('cannot create file for tabulated potential')
             end if
             write(io_ztb) zunit%boxl,zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val,zunit%ngrid,ztype%ntype,lewald,ltailc,lshift
             do igtype=1,ztype%ntype
                write(io_ztb) ztype%name(igtype)
             end do
          end if
          write(iou,*) 'time 1:',time_now()
          do k=0,zunit%ngrid(3)-1
             do j=0,zunit%ngrid(2)-1
!$omp parallel &
!$omp private(i) shared(k,j,zgrid) default(shared)
!$omp do
                do i=0,zunit%ngrid(1)-1
                   ! pass to exzeof arguments in fractional coordinates with respect to the unit cell
                   if (k.gt.oldk.or.(k.eq.oldk.and.j.ge.oldj)) call exzeof(zgrid(:,:,i,j,k),dble(i)/zunit%ngrid(1),dble(j)/zunit%ngrid(2),dble(k)/zunit%ngrid(3),ztype%ntype)
                end do
!$omp end parallel
                if (myid.eq.0) write(io_ztb) zgrid(:,:,:,j,k)
             end do
          end do
          write(iou,*) 'time 2:',time_now()
          if (myid.eq.0.and.ltailc) write(iou,*) 'maxlayer = ',nlayermax
          ! call ztest(idi)
       end if

       do i=1,zpot%ntype
          where(zgrid(1,1,:,:,:).lt.overflow)
             egrid(:,:,:,i)=0.
          elsewhere
             egrid(:,:,:,i)=overflow
          end where
          idi=zpot%table(i)
          if (lij(idi).or.lqchg(idi)) then
             do j=1,ztype%ntype
                idj=ztype%type(j)
                if (lij(idj).and.lij(idi)) then
                   where(egrid(:,:,:,i).lt.overflow) egrid(:,:,:,i)=egrid(:,:,:,i)+zgrid(1,j,:,:,:)*zpot%param(1,j,i)-zgrid(2,j,:,:,:)*zpot%param(2,j,i)
                end if
                if (lqchg(idj).and.lqchg(idi)) then
                   where(egrid(:,:,:,i).lt.overflow) egrid(:,:,:,i)=egrid(:,:,:,i)+qelect(idi)*qelect(idj)*zgrid(3,j,:,:,:)
                end if
             end do
          end if
       end do

       deallocate(lunitcell,zgrid,zeo%bead,ztype%type,ztype%num,ztype%radiisq,ztype%name,zpot%param)
       if (myid.eq.0) then
          close(io_ztb)
          write(iou,*) 'tabulated potential: lewald[',lewaldtmp,'] ltailc[',ltailctmp,'] lshift[',lshifttmp,']'
       end if
    end if

  end subroutine suzeo

! arguments i,j,k are in fractional coordinates with respect to the unit cell
! U_LJ = A/r^12 + B/r^6, scale A by 4*10^8, B by 2*10^4 to reduce round-off error
  subroutine exzeof(tab,i,j,k,zntype)
    integer,intent(in)::zntype
    real(KIND=double_precision),intent(out)::tab(3,zntype)
    real(KIND=double_precision),intent(in)::i,j,k

    integer(KIND=normal_int)::izeo,layer,ii,jj,kk,iztype
    real(KIND=double_precision)::rcutsq,vac,vbc,vb,r,scoord(3),ri(3),dr(3),r2&
     ,vnew(2,ztype%ntype)

    tab=0.
    rcutsq = zcell%cut*zcell%cut
    vbc=2e4/rcutsq**3
    vac=vbc*vbc

    if (lewald.or.(.not.ltailc)) then
       ! further scale i,j,k with respect to the whole cell
       scoord(1)=dble(i)/zunit%dup(1)
       scoord(2)=dble(j)/zunit%dup(2)
       scoord(3)=dble(k)/zunit%dup(3)
       ri(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
       ri(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
       ri(3)=scoord(3)*zcell%hmat(3,3)%val
       do izeo=1,zeo%nbead
          iztype=zeo%bead(izeo)%type
          dr=ri-zeo%bead(izeo)%coord
          call mimage(dr(1),dr(2),dr(3),boxZeo)
          r2=dot_product(dr,dr)
          if (r2.le.ztype%radiisq(iztype)) then
             tab=overflow
             return
          elseif (r2 .lt. rcutsq) then
             if (.not.ltailc) then
                vb=2e4/r2**3
                if (lshift) then
                   tab(2,iztype)=tab(2,iztype)+vb-vbc
                   tab(1,iztype)=tab(1,iztype)+vb*vb-vac
                else
                   tab(2,iztype)=tab(2,iztype)+vb
                   tab(1,iztype)=tab(1,iztype)+vb*vb
                end if
             end if
             r=sqrt(r2)
             if (lewald) then
                tab(3,iztype)=tab(3,iztype)+erfunc(calp(boxZeo)*r)/r
             else
                tab(3,iztype)=tab(3,iztype)+1./r
             end if
          end if
       end do

       if (lewald) call recipzeo(tab(3,:),ri)

    end if

    ! Calculate the Lennard-Jones interactions, include as many layers
    ! of neighboring unit cells as needed for the specified precision
    if (ltailc) then
       layer=0
       do
          vnew=0.
          do izeo=1,zeo%nbead
             if (lunitcell(izeo)) then
                iztype=zeo%bead(izeo)%type
                do ii=-layer,layer
                   do jj=-layer,layer
                      do kk=-layer,layer
                         if (abs(ii).eq.layer .or. abs(jj).eq.layer .or.abs(kk).eq.layer) then
                            scoord(1) = zeo%bead(izeo)%coord(1)*zcell%hmati(1,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(1,2)%val+zeo%bead(izeo)%coord(3)*zcell%hmati(1,3)%val+dble(ii-i)/zunit%dup(1)
                            scoord(2) = zeo%bead(izeo)%coord(1)*zcell%hmati(2,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(2,2)%val+zeo%bead(izeo)%coord(3)*zcell%hmati(2,3)%val+dble(jj-j)/zunit%dup(2)
                            scoord(3) = zeo%bead(izeo)%coord(1)*zcell%hmati(3,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(3,2)%val+zeo%bead(izeo)%coord(3)*zcell%hmati(3,3)%val+dble(kk-k)/zunit%dup(3)
                            dr(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
                            dr(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
                            dr(3)=scoord(3)*zcell%hmat(3,3)%val
                            r2=dot_product(dr,dr)
                            if (r2.le.ztype%radiisq(iztype)) then
                               tab=overflow
                               return
                            end if
                            vb=2e4/r2**3
                            if (lshift) then
                               vnew(2,iztype)=vnew(2,iztype)+vb-vbc
                               vnew(1,iztype)=vnew(1,iztype)+vb*vb-vac
                            else
                               vnew(2,iztype)=vnew(2,iztype)+vb
                               vnew(1,iztype)=vnew(1,iztype)+vb*vb
                            end if
                         end if
                      end do
                   end do
                end do
             end if
          end do
          tab(1:2,:)=tab(1:2,:)+vnew
          layer=layer+1
          if (layer.gt.nlayermax) nlayermax=layer
          if (abs(sum(vnew(1,:)-vnew(2,:))).lt.required_precision) exit
       end do
    end if

  end subroutine exzeof
 
  subroutine recipzeo(tab,ri)
    real(KIND=double_precision),intent(out)::tab(:)
    real(KIND=double_precision),intent(in)::ri(3)

    integer(KIND=normal_int)::kmax(3),i,l,m,n,m_min,n_min
    real(KIND=double_precision)::hmatik(3,3),ki(3),alpsqr4,hmaxsq,ksqr,arg&
     ,sums(ztype%ntype),vrecipz(ztype%ntype)

    ! *** Set up the reciprocal space vectors ***
    vrecipz = 0.0E+0_double_precision

    forall(m=1:3,n=1:3) hmatik(n,m) = twopi*zcell%hmati(n,m)%val
    kmax(1) = int(zcell%hmat(1,1)%val*zcell%calp)+1
    kmax(2) = int(zcell%hmat(2,2)%val*zcell%calp)+1
    kmax(3) = int(zcell%hmat(3,3)%val*zcell%calp)+1
    if (zcell%solid.and..not.zcell%ortho) kmax = kmax+1

    alpsqr4 = 4.0d0*zcell%calp*zcell%calp
    hmaxsq = alpsqr4*onepi*onepi

    ! *** generate the reciprocal-space
    ! here -kmax(1),-kmax(1)+1,...,-1 are skipped, so no need to divide by 2 for the prefactor
    do l = 0,kmax(1) 
       if ( l .eq. 0 ) then
          m_min = 0
       else
          m_min = -kmax(2)
       end if
       do m = m_min, kmax(2)
          if (l .eq. 0 .and. m .eq. 0) then
             n_min = 1
          else
             n_min = -kmax(3)
          end if
          do n = n_min, kmax(3)
             ki(1) = dble(l)*hmatik(1,1)+dble(m)*hmatik(2,1)+ dble(n)*hmatik(3,1)
             ki(2) = dble(l)*hmatik(1,2)+dble(m)*hmatik(2,2)+ dble(n)*hmatik(3,2)
             ki(3) = dble(l)*hmatik(1,3)+dble(m)*hmatik(2,3)+ dble(n)*hmatik(3,3)
             ksqr = dot_product(ki,ki)
             ! if ( ksqr .lt. hmaxsq ) then
             ! --- sometimes these are about equal, which can cause different
             ! --- behavior on 32 and 64 bit machines without this .and. statement
             if ( hmaxsq-ksqr.gt.eps ) then
                sums = 0.0E0_double_precision
                do i = 1,zeo%nbead
                   arg=dot_product(ki,ri-zeo%bead(i)%coord)
                   sums(zeo%bead(i)%type)=sums(zeo%bead(i)%type)+cos(arg)
                end do
                vrecipz=vrecipz+sums*exp(-ksqr/alpsqr4)/ksqr
             end if
          end do
       end do
    end do

    tab=tab+vrecipz*8.*onepi/zcell%vol

  end subroutine recipzeo

  function exzeo(xi,yi,zi,idi)
    real(KIND=double_precision)::exzeo
    real(KIND=double_precision),intent(in)::xi,yi,zi
    integer(KIND=normal_int),intent(in)::idi

    integer(KIND=normal_int),parameter::m=2,mt=2*m+1,mst=-m
    integer(KIND=normal_int)::j,j0,jp,k,k0,kp,l,l0,lp,igtype
    real(KIND=double_precision)::yjtmp(mst:m),yktmp(mst:m),yltmp(mst:m),xt(mst:m),yt(mst:m),zt(mst:m),scoord(3),r(3)

    ! --- fold coordinates into the unit cell, result in fractional coordinates
    scoord(1)=(xi*zunit%hmati(1,1)+yi*zunit%hmati(1,2)+zi*zunit%hmati(1,3))
    scoord(2)=(xi*zunit%hmati(2,1)+yi*zunit%hmati(2,2)+zi*zunit%hmati(2,3))
    scoord(3)=(xi*zunit%hmati(3,1)+yi*zunit%hmati(3,2)+zi*zunit%hmati(3,3))
    scoord(1)=scoord(1)-floor(scoord(1))
    scoord(2)=scoord(2)-floor(scoord(2))
    scoord(3)=scoord(3)-floor(scoord(3))

    if (lzgrid) then
       ! calculation using a grid
       do igtype=1,zpot%ntype
          if (zpot%table(igtype).eq.idi) exit
       end do
       if (igtype.gt.zpot%ntype) then
          call cleanup('exzeo: no such bead type')
       end if

       ! get the Cartesian coordinates of the point in the unit cell
       r(1)=scoord(1)*zunit%hmat(1,1)+scoord(2)*zunit%hmat(1,2)+scoord(3)*zunit%hmat(1,3)
       r(2)=scoord(2)*zunit%hmat(2,2)+scoord(3)*zunit%hmat(2,3)
       r(3)=scoord(3)*zunit%hmat(3,3)
       ! get the index of the grid that the point resides in
       j = scoord(1)*zunit%ngrid(1)
       k = scoord(2)*zunit%ngrid(2)
       l = scoord(3)*zunit%ngrid(3)

       ! --- test if in the reasonable regime
       exzeo=upperlimit
       if (egrid(j,k,l,igtype).ge.upperlimit) return
       ! ---  block m*m*m centered around: j,k,l
       ! ---  set up hulp array: (allow for going beyond unit cell
       !      for polynom fitting)
       do l0=mst,m
          lp=l+l0
          scoord(3)=dble(lp)/zunit%ngrid(3)/zunit%dup(3)
          ! ---    store x,y,z values around xi,yi,zi in arrays
          zt(l0)=scoord(3)*zcell%hmat(3,3)%val
          if (lp.lt.0) lp=lp+zunit%ngrid(3)
          if (lp.ge.zunit%ngrid(3)) lp=lp-zunit%ngrid(3)
          do k0=mst,m
             kp=k+k0
             scoord(2)=dble(kp)/zunit%ngrid(2)/zunit%dup(2)
             yt(k0)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
             if (kp.lt.0) kp=kp+zunit%ngrid(2)
             if (kp.ge.zunit%ngrid(2)) kp=kp-zunit%ngrid(2)
             do j0=mst,m
                jp=j+j0
                scoord(1)=dble(jp)/zunit%ngrid(1)/zunit%dup(1)
                xt(j0)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*hmat(boxZeo ,7)
                if (jp.lt.0) jp=jp+zunit%ngrid(1)
                if (jp.ge.zunit%ngrid(1)) jp=jp-zunit%ngrid(1)
                yjtmp(j0)=egrid(jp,kp,lp,igtype)
                if (yjtmp(j0).ge.upperlimit) return
             end do
             call polint(xt,yjtmp,mt,r(1),yktmp(k0))
          end do
          call polint(yt,yktmp,mt,r(2),yltmp(l0))
       end do
       call polint(zt,yltmp,mt,r(3),exzeo)
    end if

  end function exzeo

end module zeolite
