MODULE parser_cssr
  use var_type,only:default_string_length
  use util_runtime,only:err_exit
  use util_files,only:get_iounit
  use sim_cell
  use sim_particle
  use sim_zeolite,only:ZeoliteUnitCellGridType,ZeoliteBeadType,setUpAtom,setUpCellStruct
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: readCSSR

CONTAINS

  SUBROUTINE readCSSR(fileCSSR,zeo,lunitcell,ztype,zcell,zunit)
    character(LEN=*),intent(in)::fileCSSR
    type(MoleculeType),intent(out)::zeo
    logical,allocatable,intent(out)::lunitcell(:)
    type(ZeoliteBeadType),intent(out)::ztype
    type(CellMaskType),intent(out)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    type(ZeoliteUnitCellGridType),intent(out)::zunit

    integer::IOCSSR,jerr,i,pos
    character(LEN=default_string_length)::atom

    IOCSSR=get_iounit()
    open(unit=IOCSSR,access='sequential',action='read',file=fileCSSR,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit('cannot open zeolite CSSR file')
    end if

    read(IOCSSR,*) zcell%boxl(1)%val,zcell%boxl(2)%val,zcell%boxl(3)%val,zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val,zunit%dup(1),zunit%dup(2),zunit%dup(3)
    read(IOCSSR,*) zeo%nbead,ztype%ntype
    allocate(zeo%bead(zeo%nbead),lunitcell(zeo%nbead),ztype%name(ztype%ntype),ztype%radiisq(ztype%ntype),ztype%type(ztype%ntype),ztype%num(ztype%ntype),stat=jerr)
    if (jerr.ne.0) call err_exit('readCSSR: allocation failed')

    do i=1,ztype%ntype
       read(IOCSSR,*) ztype%name(i),ztype%type(i),ztype%radiisq(i)
       ztype%radiisq(i)=ztype%radiisq(i)*ztype%radiisq(i)
       ztype%num(i)=0
    end do

    call setUpCellStruct(zcell,zunit)

    do i = 1,zeo%nbead
       read(IOCSSR,'(i5,1x,a4,2x,3(f9.5,1x))') pos,atom,zeo%bead(i)%coord(1),zeo%bead(i)%coord(2),zeo%bead(i)%coord(3)
       call setUpAtom(atom,i,zeo,lunitcell,ztype,zcell,zunit)
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

    close(IOCSSR)

  END SUBROUTINE readCSSR

END MODULE parser_cssr
