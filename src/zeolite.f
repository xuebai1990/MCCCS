module zeolite
  use global_data
  use var_type
  use const_math
  use const_phys
  use util_math 
  use util_files
  use util_string
  use sim_cell
  use sim_particle
  implicit none
  private
  save
  public::zeocoord,suzeo,exzeo

  type ZeoliteUnitCellGridType
     integer::dup(3),ngrid(3)
     real(KIND=double_precision)::boxl(3),hmat(3,3),hmati(3,3)
  end type ZeoliteUnitCellGridType

  type ZeoliteBeadType
     integer::ntype
     integer,allocatable::type(:),num(:)
     real(KIND=double_precision),allocatable::radiisq(:)
     character(LEN=default_string_length),allocatable::name(:)
  end type ZeoliteBeadType

  integer,parameter::boxZeo=1,maxPotType=32
  type ZeolitePotentialType
     integer::ntype
     integer::table(0:maxPotType)
  end type ZeolitePotentialType

  logical,allocatable::lunitcell(:)
  integer(KIND=normal_int)::nlayermax
  real(KIND=double_precision),allocatable::egrid(:,:,:,:)
  real(KIND=double_precision),parameter::dgr=0.2,required_precision=1.0E-2_double_precision,overflow=1.0E+8_double_precision,upperlimit=1.0E+7_double_precision

  type(CellMaskType)::zcell
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
    integer(KIND=normal_int)::io_zeocoord,jerr
    integer(KIND=normal_int)::pos,i,j
    real(KIND=double_precision)::scoord(3)
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
    do i=1,3
       zcell%ang(i)%val=zcell%ang(i)%val*onepi/180
    end do

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

    do i=1,3
       do j=1,3
          zunit%hmat(j,i)=zcell%hmat(j,i)%val/zunit%dup(i)
          zunit%hmati(i,j)=zcell%hmati(i,j)%val*zunit%dup(i)
       end do
    end do

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

    close(io_zeocoord)

    return
  end subroutine zeocoord

  subroutine suzeo()
    integer(KIND=normal_int)::io_ztb,imol,iunit,igtype,idi,jerr,i,j,k,oldi,oldj,oldk,ngridtmp(3)
    character(LEN=default_path_length)::file_ztb
    real(KIND=double_precision)::zunittmp(3),wzeo

    ! === tabulation of the zeolite potential
    if (lzgrid) then
       call setpbc(boxZeo)
       ! find all bead types and store them in an array
       zpot%ntype=0
       zpot%table(0)=0
       do imol=1,nmolty
          do iunit=1,nunit(imol)
             do igtype=1,zpot%ntype
                if (zpot%table(igtype).eq.ntype(imol,iunit)) exit
             end do
             if (igtype.gt.zpot%ntype) then
                zpot%ntype=igtype
                zpot%table(zpot%ntype)=ntype(imol,iunit)
             end if
          end do
       end do

       do i=1,3
          zunit%boxl(i) = zcell%boxl(i)%val/zunit%dup(i)
       end do
       !     --- find closest values for x and y
       zunit%ngrid=int(zunit%boxl/dgr)
       allocate(egrid(0:zunit%ngrid(1)-1,0:zunit%ngrid(2)-1,0:zunit%ngrid(3)-1,0:zpot%ntype),stat=jerr)
       if (jerr.ne.0) call cleanup('suzeo: allocation failed')

       !     === Calculate zeolite density
       wzeo=dot_product(ztype%num,mass(ztype%type))
       if (myid.eq.0) write(iou,"(' Tabulated zeolite potential: ',/&
                                 ,' --------------------------------------------------',/&
                                 ,<ztype%ntype>(' number of ',A,':',I5,3X),/&
                                 ,' mass zeolite                      = ',e12.5,' grams',/ &
                                 ,' one adsorbed molecule in sim box  = ',e12.5 ,' mmol/gram',/&
                                 ,' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/&
                                 ,'         x-dir           : ',i12,'  size: ',f7.4,/&
                                 ,'         y-dir           : ',i12,'  size: ',f7.4,/&
                                 ,'         z-dir           : ',i12,'  size: ',f7.4,/)") (trim(ztype%name(i)),ztype%num(i),i=1,ztype%ntype),wzeo/6.023e23,1000.0/wzeo,zunit%boxl(1),zunit%boxl(2),zunit%boxl(3),zunit%ngrid(1),zunit%boxl(1)/zunit%ngrid(1),zunit%ngrid(2),zunit%boxl(2)/zunit%ngrid(2),zunit%ngrid(3),zunit%boxl(3)/zunit%ngrid(3)

       nlayermax=0
       io_ztb=get_iounit()

       do igtype=0,zpot%ntype
          if (igtype.eq.0.and..not.lewald) cycle
          idi=zpot%table(igtype)
          write(file_ztb,'(I3.3,A)'),idi,'.ztb'
          open(unit=io_ztb,access='sequential',action='read',file=file_ztb,form='binary',iostat=jerr, status='old')
          if (jerr.eq.0) then
             ! ---    read zeolite table from disk
             if (myid.eq.0) write(iou,*) 'read in ',trim(file_ztb)
             read(io_ztb) zunittmp(1),zunittmp(2),zunittmp(3),ngridtmp(1),ngridtmp(2),ngridtmp(3)
             if (ANY(abs(zunittmp-zunit%boxl).gt.eps).or. ANY(ngridtmp.ne.zunit%ngrid)) call cleanup('problem in zeolite potential table')
             do i=0,zunit%ngrid(1)-1
                do j=0,zunit%ngrid(2)-1
                   do k=0,zunit%ngrid(3)-1
                      read(io_ztb,end=100) egrid(i,j,k,igtype)
                   end do
                end do
             end do
          end if

100       if (jerr.eq.0.and.i.lt.zunit%ngrid(1)) then
             oldi=i
             oldj=j
             oldk=k
             jerr=1
             close(io_ztb)
             write(iou,*) oldi,oldj,oldk
          else
             oldi=0
             oldj=0
             oldk=0
          end if

          if (jerr.ne.0) then
             ! make a tabulated potential of the zeolite
             if (myid.eq.0) write(iou,*) 'make new ',trim(file_ztb)
             if (myid.eq.0) open(unit=io_ztb,access='sequential',action='write',file=file_ztb,form='binary',iostat=jerr,status='unknown')
             if (jerr.ne.0) then
                call cleanup('cannot create file for tabulated potential')
             end if
             if (myid.eq.0) write(io_ztb) zunit%boxl(1),zunit%boxl(2),zunit%boxl(3),zunit%ngrid(1),zunit%ngrid(2),zunit%ngrid(3)

             do i=0,zunit%ngrid(1)-1
                do j=0,zunit%ngrid(2)-1
                   do k=0,zunit%ngrid(3)-1
                      ! pass to exzeof arguments in fractional coordinates with respect to the unit cell
                      if (i.gt.oldi.or.(i.eq.oldi.and.j.ge.oldj)) egrid(i,j,k,igtype)=exzeof(dble(i)/zunit%ngrid(1),dble(j)/zunit%ngrid(2),dble(k)/zunit%ngrid(3),idi)
                      if (myid.eq.0) write(io_ztb) egrid(i,j,k,igtype)
                   end do
                end do
             end do
             if (myid.eq.0.and.ltailc) write(iou,*) 'maxlayer = ',nlayermax
             !     call ztest(idi)
          end if

          if (idi.ne.0.and.lewald.and.lqchg(idi)) then
             where (egrid(:,:,:,igtype).lt.overflow) egrid(:,:,:,igtype)=egrid(:,:,:,igtype)+qelect(idi)*egrid(:,:,:,0)
          end if
          if (myid.eq.0) then
             close(io_ztb)
          end if
       end do
    end if

    return
  end subroutine suzeo

! arguments i,j,k are in fractional coordinates with respect to the unit cell
  function exzeof(i,j,k,idi)
    real(KIND=double_precision)::exzeof
    real(KIND=double_precision),intent(in)::i,j,k
    integer(KIND=normal_int),intent(in)::idi

    integer(KIND=normal_int)::izeo,idj,ntij,layer,ii,jj,kk
    real(KIND=double_precision)::scoord(3),ri(3),dr(3),r2,rcutsq,r,r2i,r6,vljnew,vqnew

    exzeof=0.

    if (idi.eq.0.or.(.not.ltailc).or.(.not.lewald.and.lqchg(idi))) then
       vljnew=0.
       vqnew=0.
       rcutsq = zcell%cut*zcell%cut
       ! further scale i,j,k with respect to the whole cell
       scoord(1)=dble(i)/zunit%dup(1)
       scoord(2)=dble(j)/zunit%dup(2)
       scoord(3)=dble(k)/zunit%dup(3)
       ri(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
       ri(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
       ri(3)=scoord(3)*zcell%hmat(3,3)%val
       do izeo=1,zeo%nbead
          idj=ztype%type(zeo%bead(izeo)%type)
          ntij = (idi - 1) * nntype + idj
          dr=ri-zeo%bead(izeo)%coord
          call mimage(dr(1),dr(2),dr(3),boxZeo)
          r2=dot_product(dr,dr)
          if (r2.le.ztype%radiisq(zeo%bead(izeo)%type)) then
             exzeof=overflow
             return
          end if
          if (r2 .lt. rcutsq) then
             if (idi.eq.0) then
                if (lqchg(idj)) then
                   r=dsqrt(r2)
                   vqnew=vqnew+qelect(idj)*erfunc(calp(boxZeo)*r)/r
                end if
             else
                if (.not.lewald.and.lqchg(idi).and.lqchg(idj)) then
                   r=dsqrt(r2)
                   vqnew=vqnew+qelect(idi)*qelect(idj)/r
                end if
                if (.not.ltailc.and.lij(idi).and.lij(idj)) then
                   r2i=sig2ij(ntij)/r2
                   r6=r2i*r2i*r2i
                   if (lshift) then     
                      vljnew=vljnew+4.*(epsij(ntij)*(r6-1.0)*r6-ecut(ntij))
                   else
                      vljnew=vljnew+4.*epsij(ntij)*(r6-1.0)*r6
                   end if
                end if
             end if
          end if
       end do

       if (idi.eq.0) vqnew=vqnew+recipzeo(ri,1._double_precision)

       exzeof=vljnew+vqnew*qqfact
    end if

    ! Calculate the Lennard-Jones interactions, include as many layers
    ! of neighboring unit cells as needed for the specified precision
    if (idi.ne.0.and.ltailc) then
       vljnew=required_precision+1.
       layer=0
       do while (abs(vljnew).gt.required_precision)
          if (layer.gt.nlayermax) nlayermax=layer
          vljnew=0.
          do izeo=1,zeo%nbead
             if (lunitcell(izeo)) then
                idj=ztype%type(zeo%bead(izeo)%type)
                ntij = (idi - 1) * nntype + idj
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

                            r2=dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
                            if (r2.le.ztype%radiisq(zeo%bead(izeo)%type)) then
                               exzeof=overflow
                               return
                            end if
                            if (lij(idi).and.lij(idj)) then
                               r2i=sig2ij(ntij)/r2
                               r6=r2i*r2i*r2i
                               if (lshift) then
                                  vljnew=vljnew+4.*(epsij(ntij)*(r6-1.0)*r6-ecut(ntij))
                               else
                                  vljnew=vljnew+4.*epsij(ntij)*(r6-1.0)*r6
                               end if
                            end if
                         end if
                      end do
                   end do
                end do
             end if
          end do
          exzeof=exzeof+vljnew
          layer=layer+1
       end do
    end if

    return
  end function exzeof
 
  function recipzeo(ri,qi)
    real(KIND=double_precision)::recipzeo
    real(KIND=double_precision),intent(in)::ri(3),qi

    integer(KIND=normal_int)::i,l,m,n,m_min,n_min,kmax(3)
    real(KIND=double_precision)::alpsqr4,ksqr,sum,arg,hmatik(3,3),ki(3),hmaxsq

    ! *** Set up the reciprocal space vectors ***
    recipzeo = 0.0E+0_double_precision

    forall(m=1:3,n=1:3) hmatik(n,m) = twopi*zcell%hmati(n,m)%val
    kmax(1) = dint(zcell%hmat(1,1)%val*zcell%calp)+1
    kmax(2) = dint(zcell%hmat(2,2)%val*zcell%calp)+1
    kmax(3) = dint(zcell%hmat(3,3)%val*zcell%calp)+1
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
                sum = 0.0d0
                do i = 1,zeo%nbead
                   arg=ki(1)*(ri(1)-zeo%bead(i)%coord(1))+ki(2)*(ri(2)-zeo%bead(i)%coord(2))+ki(3)*(ri(3)-zeo%bead(i)%coord(3))
                   sum=sum+qelect(ztype%type(zeo%bead(i)%type))*cos(arg)
                end do
                recipzeo=recipzeo+sum*exp(-ksqr/alpsqr4)/ksqr
             end if
          end do
       end do
    end do

    recipzeo=recipzeo*qi*(8.0d0*onepi)/zcell%vol

    return

  end function recipzeo

  function exzeo(xi,yi,zi,idi)
    real(KIND=double_precision)::exzeo
    real(KIND=double_precision),intent(in)::xi,yi,zi
    integer(KIND=normal_int),intent(in)::idi

    integer(KIND=normal_int),parameter::m=2,mt=2*m+1,mst=-m
    integer(KIND=normal_int)::j,j0,jp,k,k0,kp,l,l0,lp,igtype
    real(KIND=double_precision)::yjtmp(mst:m),yktmp(mst:m),yltmp(mst:m)
    real(KIND=double_precision)::xt(mst:m),yt(mst:m),zt(mst:m),scoord(3),r(3)

    ! --- determine cell parameters
    scoord(1)=(xi*zunit%hmati(1,1)+yi*zunit%hmati(1,2)+zi*zunit%hmati(1,3))
    scoord(2)=(xi*zunit%hmati(2,1)+yi*zunit%hmati(2,2)+zi*zunit%hmati(2,3))
    scoord(3)=(xi*zunit%hmati(3,1)+yi*zunit%hmati(3,2)+zi*zunit%hmati(3,3))
    scoord(1)=scoord(1)-floor(scoord(1))
    scoord(2)=scoord(2)-floor(scoord(2))
    scoord(3)=scoord(3)-floor(scoord(3))

    if (.not.lzgrid) then
       exzeo=exzeof(scoord(1),scoord(2),scoord(3),idi)
       if (lewald.and.lqchg(idi)) exzeo=exzeo+exzeof(scoord(1),scoord(2),scoord(3),0)*qelect(idi)
    else
       !     calculation using a grid
       do igtype=1,zpot%ntype
          if (zpot%table(igtype).eq.idi) exit
       end do
       if (igtype.gt.zpot%ntype) then
          call cleanup('exzeo: no such bead type')
       end if

       r(1)=scoord(1)*zunit%hmat(1,1)+scoord(2)*zunit%hmat(1,2)+scoord(3)*zunit%hmat(1,3)
       r(2)=scoord(2)*zunit%hmat(2,2)+scoord(3)*zunit%hmat(2,3)
       r(3)=scoord(3)*zunit%hmat(3,3)
       j = scoord(1)*zunit%ngrid(1)
       k = scoord(2)*zunit%ngrid(2)
       l = scoord(3)*zunit%ngrid(3)

       ! ---    test if in the reasonable regime
       exzeo=upperlimit
       if (egrid(j,k,l,igtype).ge.exzeo) return
       ! --     block m*m*m centered around: j,k,l
       ! ---  set up hulp array: (allow for going beyond unit cell
       !      for polynom fitting)
       do l0=mst,m
          lp=l+l0
          scoord(3)=dble(lp)/zunit%ngrid(3)/zunit%dup(3)
          ! ---    store x,y,z values around xi,yi,zi in arrays
          zt(l0)=scoord(3)*zcell%hmat(3,3)%val
          if (lp.lt.0)    lp=lp+zunit%ngrid(3)
          if (lp.ge.zunit%ngrid(3)) lp=lp-zunit%ngrid(3)
          do k0=mst,m
             kp=k+k0
             scoord(2)=dble(kp)/zunit%ngrid(2)/zunit%dup(2)
             yt(k0)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
             if (kp.lt.0)    kp=kp+zunit%ngrid(2)
             if (kp.ge.zunit%ngrid(2)) kp=kp-zunit%ngrid(2)
             do j0=mst,m
                jp=j+j0
                scoord(1)=dble(jp)/zunit%ngrid(1)/zunit%dup(1)
                xt(j0)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*hmat(boxZeo ,7)
                if (jp.lt.0)    jp=jp+zunit%ngrid(1)
                if (jp.ge.zunit%ngrid(1)) jp=jp-zunit%ngrid(1)
                yjtmp(j0)=egrid(jp,kp,lp,igtype)
                if (yjtmp(j0).ge.exzeo) return
             end do
             call polint(xt,yjtmp,mt,r(1),yktmp(k0))
          end do
          call polint(yt,yktmp,mt,r(2),yltmp(l0))
       end do
       call polint(zt,yltmp,mt,r(3),exzeo)
    end if
    return
  end function exzeo

end module zeolite
