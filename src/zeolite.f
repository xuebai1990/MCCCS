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
  use sim_zeolite
  use parser_pdb
  use parser_cssr
  use parser_cif
  implicit none
  include 'common.inc'
  private
  save
  public::zeocoord,suzeo,exzeo

  real(KIND=double_precision),parameter::requiredPrecision=1.0E-2_double_precision,overlapValue=1.0E+20_double_precision,upperLimit=1.0E+7_double_precision
  integer,parameter::boxZeo=1

  logical,allocatable::lunitcell(:)
  integer(KIND=normal_int)::nlayermax
  real(KIND=double_precision),allocatable::zgrid(:,:,:,:,:),egrid(:,:,:,:)

  type(MoleculeType)::zeo
  type(ZeoliteBeadType)::ztype
  type(CellMaskType)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
  type(ZeoliteUnitCellGridType)::zunit
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

    call initZeo()

    if (index(file_zeocoord,'.cif').gt.0) then
       call readCIF(file_zeocoord,zeo,lunitcell,lhere,ztype,zcell,zunit)
    else if (index(file_zeocoord,'.pdb').gt.0) then
       call readPDB(file_zeocoord,zeo,lunitcell,lhere,ztype,zcell,zunit)
    else if (index(file_zeocoord,'.cssr').gt.0) then
       call readCSSR(file_zeocoord,zeo,lunitcell,lhere,ztype,zcell,zunit)
    end if

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
    real(KIND=double_precision)::wzeo,zunittmp(3),zangtmp(3),rcuttmp,rci3,rci9,rho
    logical::lewaldtmp,ltailczeotmp,lshifttmp

    !     === Calculate zeolite density
    wzeo=dot_product(ztype%num(1:ztype%ntype),mass(ztype%type(1:ztype%ntype)))
    if (myid.eq.0) write(iou,"(' Tabulated zeolite potential: ',/&
                              ,' --------------------------------------------------',/&
                              ,<ztype%ntype>(' number of ',A,':',I5,3X),/&
                              ,' simulation volume                 = ', f10.1, 20x,' Angst**3',/&
                              ,' number of atomtypes in the lattice= ',i10,/&
                              ,' number of zeolite atoms           = ',i10,/&
                              ,' mass zeolite                      = ',e12.5,' grams',/ &
                              ,' one adsorbed molecule in sim box  = ',e12.5 ,' mmol/gram',/&
                              ,' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/&
                              ,'         x-dir           : ',i12,'  size: ',f7.4,/&
                              ,'         y-dir           : ',i12,'  size: ',f7.4,/&
                              ,'         z-dir           : ',i12,'  size: ',f7.4,/)") (trim(ztype%name(i)),ztype%num(i),i=1,ztype%ntype),zcell%vol,ztype%ntype,zeo%nbead,wzeo/6.023e23,1000.0/wzeo,zunit%boxl(1),zunit%boxl(2),zunit%boxl(3),zunit%ngrid(1),zunit%boxl(1)/zunit%ngrid(1),zunit%ngrid(2),zunit%boxl(2)/zunit%ngrid(2),zunit%ngrid(3),zunit%boxl(3)/zunit%ngrid(3)

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
          read(io_ztb) zunittmp,zangtmp,ngridtmp,zntypetmp,lewaldtmp,ltailczeotmp,lshifttmp,rcuttmp
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
!           if (myid.eq.0) then
!              rcuttmp=10.0d0
!              close(io_ztb)
!              open(unit=io_ztb,access='sequential',action='write',file=file_ztb,form='binary',iostat=jerr,status='unknown')
!              if (jerr.ne.0) then
!                 call cleanup('cannot create file for tabulated potential')
!              end if
!              write(io_ztb) zunittmp,zangtmp,ngridtmp,zntypetmp,lewaldtmp,ltailczeotmp,lshifttmp,rcuttmp
!              do igtype=1,ztype%ntype
!                 write(io_ztb) ztype%name(igtype)
!              end do
!              write(io_ztb) zgrid
!           end if
       else
          oldi=0
          oldj=0
          oldk=0
          lewaldtmp=lewald
          ltailczeotmp=ltailcZeo
          lshifttmp=lshift
          rcuttmp=zcell%cut
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
             write(io_ztb) zunit%boxl,zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val,zunit%ngrid,ztype%ntype,lewald,ltailcZeo,lshift,zcell%cut
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
          if (myid.eq.0.and.ltailcZeo) write(iou,*) 'maxlayer = ',nlayermax
          ! call ztest(idi)
       end if

       if (ltailczeo.and..not.ltailcZeotmp) then
          rci3=1./rcuttmp**3
          rci9=4e8*rci3**3
          rci3=2e4*rci3
       end if

       egrid=0.
       forall(k=0:zunit%ngrid(3)-1,j=0:zunit%ngrid(2)-1,i=0:zunit%ngrid(1)-1,zgrid(1,1,i,j,k).ge.overlapValue) egrid(i,j,k,:)=overlapValue

       do j=1,ztype%ntype
          idj=ztype%type(j)
          if (lij(idj).or.lqchg(idj)) then
             rho=ztype%num(j)/zcell%vol
             do i=1,zpot%ntype
                idi=zpot%table(i)
                if (lij(idj).and.lij(idi)) then
                   where(egrid(:,:,:,i).lt.overlapValue) egrid(:,:,:,i)=egrid(:,:,:,i)+zgrid(1,j,:,:,:)*zpot%param(1,j,i)-zgrid(2,j,:,:,:)*zpot%param(2,j,i)
                   if (ltailczeo.and..not.ltailcZeotmp) then
                      where(egrid(:,:,:,i).lt.overlapValue) egrid(:,:,:,i)=egrid(:,:,:,i)+twopi*rho*(rci9*zpot%param(1,j,i)/9-rci3*zpot%param(2,j,i)/3)
                   end if
                end if
                if (lqchg(idj).and.lqchg(idi)) then
                   where(egrid(:,:,:,i).lt.overlapValue) egrid(:,:,:,i)=egrid(:,:,:,i)+qqfact*qelect(idi)*qelect(idj)*zgrid(3,j,:,:,:)
                end if
             end do
          end if
       end do

       deallocate(lunitcell,zgrid,zeo%bead,ztype%type,ztype%num,ztype%radiisq,ztype%name,zpot%param)
       if (myid.eq.0) then
          close(io_ztb)
          write(iou,*) 'tabulated potential: lewald[',lewaldtmp,'] ltailc[',ltailcZeotmp,'] lshift[',lshifttmp,'] rcut[',rcuttmp,'] ltailcZeo[',ltailcZeo,']'
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
    real(KIND=double_precision)::rcutsq,vac,vbc,vb,r,scoord(3),ri(3),dr(3),r2,vnew(2,zntype)

    tab=0.
    rcutsq = zcell%cut*zcell%cut
    vbc=2e4/rcutsq**3
    vac=vbc*vbc

    if (lewald.or.(.not.ltailcZeo)) then
       ! further scale i,j,k with respect to the whole cell
       !!!
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
             tab=overlapValue
             return
          elseif (r2 .lt. rcutsq) then
             if (.not.ltailcZeo) then
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

       if (lewald) call recipzeo(tab(3,:),ri,zntype)

    end if

    ! Calculate the Lennard-Jones interactions, include as many layers
    ! of neighboring unit cells as needed for the specified precision
    if (ltailcZeo) then
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
                            !!!
                            scoord(1) = zeo%bead(izeo)%coord(1)*zcell%hmati(1,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(1,2)%val+zeo%bead(izeo)%coord(3)*zcell%hmati(1,3)%val+dble(ii-i)/zunit%dup(1)
                            scoord(2) = zeo%bead(izeo)%coord(1)*zcell%hmati(2,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(2,2)%val+zeo%bead(izeo)%coord(3)*zcell%hmati(2,3)%val+dble(jj-j)/zunit%dup(2)
                            scoord(3) = zeo%bead(izeo)%coord(1)*zcell%hmati(3,1)%val+zeo%bead(izeo)%coord(2)*zcell%hmati(3,2)%val+zeo%bead(izeo)%coord(3)*zcell%hmati(3,3)%val+dble(kk-k)/zunit%dup(3)
                            dr(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
                            dr(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
                            dr(3)=scoord(3)*zcell%hmat(3,3)%val
                            r2=dot_product(dr,dr)
                            if (r2.le.ztype%radiisq(iztype)) then
                               tab=overlapValue
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
          if (abs(sum(vnew(1,:)-vnew(2,:))).lt.requiredPrecision) exit
       end do
    end if

  end subroutine exzeof
 
  subroutine recipzeo(tab,ri,zntype)
    integer,intent(in)::zntype
    real(KIND=double_precision),intent(out)::tab(zntype)
    real(KIND=double_precision),intent(in)::ri(3)

    integer(KIND=normal_int)::kmax(3),i,j,l,m,n,kmin(2:3)
    real(KIND=double_precision)::hmatik(3,3),ki(3),alpsqr4,hmaxsq,ksqr,arg,sums(zntype),vrecipz(zntype)

    ! *** Set up the reciprocal space vectors ***
    vrecipz = 0.0E+0_double_precision

    forall(i=1:3,j=1:3) hmatik(j,i) = twopi*zcell%hmati(j,i)%val
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
          kmin(2) = 0
       else
          kmin(2) = -kmax(2)
       end if
       do m = kmin(2), kmax(2)
          if (l .eq. 0 .and. m .eq. 0) then
             kmin(3) = 1
          else
             kmin(3) = -kmax(3)
          end if
          do n = kmin(3), kmax(3)
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
    !!!
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
       exzeo=upperLimit
       if (egrid(j,k,l,igtype).ge.upperLimit) return
       ! ---  block m*m*m centered around: j,k,l
       ! ---  set up hulp array: (allow for going beyond unit cell
       !      for polynom fitting)
       do l0=mst,m
          lp=l+l0
          scoord(3)=dble(lp)/zunit%ngrid(3)
          ! ---    store x,y,z values around xi,yi,zi in arrays
          zt(l0)=scoord(3)*zunit%hmat(3,3)
          if (lp.lt.0) lp=lp+zunit%ngrid(3)
          if (lp.ge.zunit%ngrid(3)) lp=lp-zunit%ngrid(3)
          do k0=mst,m
             kp=k+k0
             scoord(2)=dble(kp)/zunit%ngrid(2)
             yt(k0)=scoord(2)*zunit%hmat(2,2)+scoord(3)*zunit%hmat(2,3)
             if (kp.lt.0) kp=kp+zunit%ngrid(2)
             if (kp.ge.zunit%ngrid(2)) kp=kp-zunit%ngrid(2)
             do j0=mst,m
                jp=j+j0
                scoord(1)=dble(jp)/zunit%ngrid(1)
                xt(j0)=scoord(1)*zunit%hmat(1,1)+scoord(2)*zunit%hmat(1,2)+scoord(3)*zunit%hmat(1,3)
                if (jp.lt.0) jp=jp+zunit%ngrid(1)
                if (jp.ge.zunit%ngrid(1)) jp=jp-zunit%ngrid(1)
                yjtmp(j0)=egrid(jp,kp,lp,igtype)
                if (yjtmp(j0).ge.upperLimit) return
             end do
             call polint(xt,yjtmp,mt,r(1),yktmp(k0))
          end do
          call polint(yt,yktmp,mt,r(2),yltmp(l0))
       end do
       call polint(zt,yltmp,mt,r(3),exzeo)
    end if

  end function exzeo

end module zeolite
