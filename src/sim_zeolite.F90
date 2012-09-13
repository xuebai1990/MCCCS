MODULE sim_zeolite
  use global_data,only:myid,iou
  use var_type,only:double_precision,default_string_length
  use const_math,only:onepi
  use util_string,only:str_search
  use sim_cell
  use sim_particle
  implicit none
  private
  public::ZeoliteUnitCellGridType,ZeoliteBeadType,ZeolitePotentialType,setUpAtom,setUpCellStruct,foldToCenterCell,foldToUnitCell,fractionalToAbsolute,absoluteToFractional

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

  integer,parameter::boxZeo=1
  real(KIND=double_precision),parameter::eps=1.0E-4_double_precision,dgr=0.2_double_precision

CONTAINS

  subroutine setUpCellStruct(zcell,zunit)
    type(CellMaskType),intent(inout)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    type(ZeoliteUnitCellGridType),intent(inout)::zunit

    integer::i,j

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
    !bi
    zcell%hmat(1,2)%val=zcell%boxl(2)%val*cos(zcell%ang(3)%val)
    !bj
    zcell%hmat(2,2)%val=zcell%boxl(2)%val*sin(zcell%ang(3)%val)
    zcell%hmat(3,2)%val=0.
    !ci
    zcell%hmat(1,3)%val=zcell%boxl(3)%val*cos(zcell%ang(2)%val)
    !(b*c*cos(alpha)-bi*ci)/bj
    zcell%hmat(2,3)%val=zcell%boxl(3)%val*sqrt(cos(zcell%ang(2)%val)**2*cos(zcell%ang(3)%val)**2+cos(zcell%ang(1)%val)**2-2*cos(zcell%ang(1)%val)*cos(zcell%ang(2)%val)*cos(zcell%ang(3)%val))/sin(zcell%ang(3)%val)
    !sqrt(c**2-ci**2-cj**2)
    zcell%hmat(3,3)%val=zcell%boxl(3)%val*sqrt(1-cos(zcell%ang(1)%val)**2-cos(zcell%ang(2)%val)**2-cos(zcell%ang(3)%val)**2+2*cos(zcell%ang(1)%val)*cos(zcell%ang(2)%val)*cos(zcell%ang(3)%val))/sin(zcell%ang(3)%val)

! there might be numeric errors in the above trigonometric calculations
    forall(i=1:3,j=1:3,abs(zcell%hmat(j,i)%val).lt.eps) zcell%hmat(j,i)%val=0

    call matops(boxZeo)

    forall(i=1:3,j=1:3)
       zunit%hmat(j,i)=zcell%hmat(j,i)%val/zunit%dup(i)
       zunit%hmati(i,j)=zcell%hmati(i,j)%val*zunit%dup(i)
    end forall

  end subroutine setUpCellStruct

  subroutine setUpAtom(atom,i,zeo,lunitcell,ztype,zcell,zunit)
    character(LEN=*),intent(in)::atom
    integer,intent(in)::i
    type(MoleculeType),intent(inout)::zeo
    logical,intent(out)::lunitcell(zeo%nbead)
    type(ZeoliteBeadType),intent(inout)::ztype
    type(CellMaskType),intent(inout)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    type(ZeoliteUnitCellGridType),intent(in)::zunit
    integer::pos
    real(KIND=double_precision)::scoord(3)

    CALL foldToCenterCell(zeo%bead(i)%coord,zcell,scoord)

    if (ALL(scoord*zunit%dup.lt.1)) then
       lunitcell(i)=.true.
    else
       lunitcell(i)=.false.
    end if

    pos=str_search(ztype%name,ztype%ntype,atom)
    if (pos.eq.0) call cleanup('** setUpAtom: unknown atomtype **')

    zeo%bead(i)%type=pos

    ztype%num(pos)=ztype%num(pos)+1

  end subroutine setUpAtom

  subroutine foldToCenterCell(coord,zcell,scoord)
    real(KIND=double_precision),intent(inout)::coord(3)
    type(CellMaskType),intent(in)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    real(KIND=double_precision),intent(out)::scoord(3)
    
    scoord(1) = coord(1)*zcell%hmati(1,1)%val+coord(2)*zcell%hmati(1,2)%val+coord(3)*zcell%hmati(1,3)%val
    scoord(2) = coord(1)*zcell%hmati(2,1)%val+coord(2)*zcell%hmati(2,2)%val+coord(3)*zcell%hmati(2,3)%val
    scoord(3) = coord(1)*zcell%hmati(3,1)%val+coord(2)*zcell%hmati(3,2)%val+coord(3)*zcell%hmati(3,3)%val
    scoord = scoord - floor(scoord)
    coord(1)=scoord(1)*zcell%hmat(1,1)%val+scoord(2)*zcell%hmat(1,2)%val+scoord(3)*zcell%hmat(1,3)%val
    coord(2)=scoord(2)*zcell%hmat(2,2)%val+scoord(3)*zcell%hmat(2,3)%val
    coord(3)=scoord(3)*zcell%hmat(3,3)%val
  end subroutine foldToCenterCell

  subroutine foldToUnitCell(coord,zunit,scoord)
    real(KIND=double_precision),intent(inout)::coord(3)
    type(ZeoliteUnitCellGridType),intent(in)::zunit
    real(KIND=double_precision),intent(out)::scoord(3)
    
    scoord(1) = coord(1)*zunit%hmati(1,1)+coord(2)*zunit%hmati(1,2)+coord(3)*zunit%hmati(1,3)
    scoord(2) = coord(1)*zunit%hmati(2,1)+coord(2)*zunit%hmati(2,2)+coord(3)*zunit%hmati(2,3)
    scoord(3) = coord(1)*zunit%hmati(3,1)+coord(2)*zunit%hmati(3,2)+coord(3)*zunit%hmati(3,3)
    scoord = scoord - floor(scoord)
    coord(1)=scoord(1)*zunit%hmat(1,1)+scoord(2)*zunit%hmat(1,2)+scoord(3)*zunit%hmat(1,3)
    coord(2)=scoord(2)*zunit%hmat(2,2)+scoord(3)*zunit%hmat(2,3)
    coord(3)=scoord(3)*zunit%hmat(3,3)
  end subroutine foldToUnitCell

  subroutine fractionalToAbsolute(newcoord,oldcoord,zcell)
    real(KIND=double_precision),intent(in)::oldcoord(3)
    type(CellMaskType),intent(in)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    real(KIND=double_precision),intent(out)::newcoord(3)
    
    newcoord(1)=oldcoord(1)*zcell%hmat(1,1)%val+oldcoord(2)*zcell%hmat(1,2)%val+oldcoord(3)*zcell%hmat(1,3)%val
    newcoord(2)=oldcoord(2)*zcell%hmat(2,2)%val+oldcoord(3)*zcell%hmat(2,3)%val
    newcoord(3)=oldcoord(3)*zcell%hmat(3,3)%val
  end subroutine fractionalToAbsolute

  subroutine absoluteToFractional(newcoord,oldcoord,zcell)
    real(KIND=double_precision),intent(in)::oldcoord(3)
    type(CellMaskType),intent(in)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    real(KIND=double_precision),intent(out)::newcoord(3)
    
    newcoord(1) = oldcoord(1)*zcell%hmati(1,1)%val+oldcoord(2)*zcell%hmati(1,2)%val+oldcoord(3)*zcell%hmati(1,3)%val
    newcoord(2) = oldcoord(1)*zcell%hmati(2,1)%val+oldcoord(2)*zcell%hmati(2,2)%val+oldcoord(3)*zcell%hmati(2,3)%val
    newcoord(3) = oldcoord(1)*zcell%hmati(3,1)%val+oldcoord(2)*zcell%hmati(3,2)%val+oldcoord(3)*zcell%hmati(3,3)%val
  end subroutine absoluteToFractional
end MODULE sim_zeolite
