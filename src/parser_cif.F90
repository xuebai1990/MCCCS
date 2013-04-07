MODULE parser_cif
! *****************************************************************************
!> \brief  Handles CIF (Crystallographic Information File) files
!>         Format Information implemented:
!>            _cell_length_a
!>            _cell_length_b
!>            _cell_length_c
!>            _cell_angle_alpha
!>            _cell_angle_beta
!>            _cell_angle_gamma
!>            _symmetry_equiv_pos_as_xyz
!>            _atom_site_label
!>            _atom_site_type_symbol
!>            _atom_site_fract_x
!>            _atom_site_fract_y
!>            _atom_site_fract_z
! *****************************************************************************
  use var_type,only:dp,default_string_length
  use util_runtime,only:err_exit
  use util_files,only:get_iounit,readLine
  use util_string,only:splitAndGetNext
  use sim_cell
  use sim_particle
  use sim_zeolite,only:ZeoliteUnitCellGridType,ZeoliteBeadType,setUpAtom,setUpCellStruct,foldToCenterCell,foldToUnitCell,fractionalToAbsolute,absoluteToFractional
  use fparser,only:initf,parsef,evalf,finalizef
  IMPLICIT NONE
  PRIVATE
  PUBLIC::readCIF
  integer,parameter::boxZeo=1
  real,parameter::eps=1.0E-4_dp

CONTAINS

! *****************************************************************************
!> \brief  Performs the real task of reading the proper information from the CIF
!>         file
! *****************************************************************************
  SUBROUTINE readCIF(fileCIF,zeo,lunitcell,ztype,zcell,zunit,lprint)
    character(LEN=*),intent(in)::fileCIF
    type(MoleculeType),intent(out)::zeo
    logical,allocatable,intent(out)::lunitcell(:)
    type(ZeoliteBeadType),intent(out)::ztype
    type(CellMaskType),intent(out)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    type(ZeoliteUnitCellGridType),intent(out)::zunit
    LOGICAL,INTENT(IN)::lprint

    integer,parameter::maxNumField=20,maxNumSymmOp=50
    INTEGER::IOCIF,jerr,i,uninitialized,ia,ib,ic,id,nAtom,nField,nSymm,Field(maxNumField)
    CHARACTER(LEN=default_string_length)::line,atom,SymmOp(maxNumSymmOp,3)
    real::scoord(3),tmpcoord(3),coord(3),dr(3)

    IOCIF=get_iounit()
    open(unit=IOCIF,access='sequential',action='read',file=fileCIF,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open zeolite CIF file',-1)
    end if

    CALL readLine(IOCIF,line,.false.,jerr)
    IF(jerr.ne.0) call err_exit(__FILE__,__LINE__,'wrong CIF file format',-1)

    read(line(2:),*) zeo%nbead,zunit%dup(1),zunit%dup(2),zunit%dup(3),ztype%ntype
    allocate(zeo%bead(zeo%nbead),lunitcell(zeo%nbead),ztype%name(ztype%ntype),ztype%radiisq(ztype%ntype),ztype%type(ztype%ntype),ztype%num(ztype%ntype),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'readCIF: allocation failed',-1)

    do i=1,ztype%ntype
       CALL readLine(IOCIF,line,.false.,jerr)
       IF(jerr.ne.0) call err_exit(__FILE__,__LINE__,'wrong CIF file format',-1)
       read(line(2:),*) ztype%name(i),ztype%type(i),ztype%radiisq(i)
       ztype%radiisq(i)=ztype%radiisq(i)*ztype%radiisq(i)
       ztype%num(i)=0
    end do

    uninitialized=7
    DO
       CALL readLine(IOCIF,line,.true.,jerr)
       if (jerr.ne.0) EXIT

       if (index(line,'_cell_length_a').gt.0) then
          CALL getReal(line(index(line,'_cell_length_a')+14:),zunit%boxl(1),jerr)
          uninitialized=uninitialized-1
       else if (index(line,'_cell_length_b').gt.0) then
          CALL getReal(line(index(line,'_cell_length_b')+14:),zunit%boxl(2),jerr)
          uninitialized=uninitialized-1
       else if (index(line,'_cell_length_c').gt.0) then
          CALL getReal(line(index(line,'_cell_length_c')+14:),zunit%boxl(3),jerr)
          uninitialized=uninitialized-1
       else if (index(line,'_cell_angle_alpha').gt.0) then
          CALL getReal(line(index(line,'_cell_angle_alpha')+17:),zcell%ang(1)%val,jerr)
          uninitialized=uninitialized-1
       else if (index(line,'_cell_angle_beta').gt.0) then
          CALL getReal(line(index(line,'_cell_angle_beta')+16:),zcell%ang(2)%val,jerr)
          uninitialized=uninitialized-1
       else if (index(line,'_cell_angle_gamma').gt.0) then
          CALL getReal(line(index(line,'_cell_angle_gamma')+17:),zcell%ang(3)%val,jerr)
          uninitialized=uninitialized-1
       else if (index(line,'_symmetry_equiv_pos_as_xyz').gt.0) then ! Read in symmetry operations
          i  = 0
          DO
             CALL readLine(IOCIF,line,.true.,jerr)
             if (jerr.ne.0) then
                EXIT
             else if ((INDEX(line,"loop_").ne.0).or.(line(1:1).eq."_")) then
                exit
             end if
             i = i + 1
             ia = INDEX(line,"'")
             ib = INDEX(line(ia+1:),",")+ia
             ic = INDEX(line(ib+1:),",")+ib
             IF (ia.eq.0) THEN
                id = LEN_TRIM(line)+1
             ELSE
                id = INDEX(line(ic+1:),"'")+ic
             END IF
             SymmOp(i,1)=TRIM(line(ia+1:ib-1))
             SymmOp(i,2)=TRIM(line(ib+1:ic-1))
             SymmOp(i,3)=TRIM(line(ic+1:id-1))
          END DO
          nSymm=i
          uninitialized=uninitialized-1
       else if (uninitialized.lt.0.and.index(line,'_atom_site_label').gt.0) then
          i=0
          do
             i = i + 1
             IF (INDEX(line,"_atom_site_label").ne.0) THEN
                Field(i) = 1
             ELSE IF (INDEX(line,"_atom_site_type_symbol").ne.0) THEN
                Field(i) = 2
             ELSE IF (INDEX(line,"_atom_site_fract_x").ne.0) THEN
                Field(i) = 3
             ELSE IF (INDEX(line,"_atom_site_fract_y").ne.0) THEN
                Field(i) = 4
             ELSE IF (INDEX(line,"_atom_site_fract_z").ne.0) THEN
                Field(i) = 5
             ELSE
                Field(i) = 0
             END IF
             CALL readLine(IOCIF,line,.true.,jerr)
             if (jerr.ne.0) then
                EXIT
             else if (index(line,'_atom_site_').eq.0) then
                exit
             end if
          end do
          nField=i
          ! Parse each line entry
          i = 0
          DO WHILE ((INDEX(line,"loop_").eq.0).AND.(line(1:1).ne."_"))
             i = i + 1
             ic=0
             DO id=1,nField
                call splitAndGetNext(line(ic+1:),ia,ib)
                ia=ic+ia
                ib=ic+ib
                SELECT CASE (Field(id))
                CASE (2)
                   read(line(ia:ib),*) atom
                CASE (3)
                   CALL getReal(line(ia:ib),scoord(1),jerr)
                CASE (4)
                   CALL getReal(line(ia:ib),scoord(2),jerr)
                CASE (5)
                   CALL getReal(line(ia:ib),scoord(3),jerr)
                CASE (1,0)
                   ! Skip this field
                CASE DEFAULT
                END SELECT
                ic=ib
             END DO
             scoord=scoord-floor(scoord)
             call fractionalToAbsolute(zeo%bead(i)%coord,scoord/zunit%dup,zcell)
             call setUpAtom(atom,i,zeo,lunitcell,ztype,zcell,zunit)
             CALL readLine(IOCIF,line,.true.,jerr)
             if (jerr.ne.0) EXIT
          END DO
          nAtom=i

          ! Apply symmetry elements and generate the whole set of atoms in the unit cell
          CALL initf(3)
          DO id = 1, nSymm
             CALL parsef(1,SymmOp(id,1),(/'x','y','z'/))
             CALL parsef(2,SymmOp(id,2),(/'x','y','z'/))
             CALL parsef(3,SymmOp(id,3),(/'x','y','z'/))
             DO ic = 1, natom
                tmpcoord=zeo%bead(ic)%coord
                CALL foldToUnitCell(tmpcoord,zunit,scoord)
                tmpcoord(1) = evalf(1,scoord)
                tmpcoord(2) = evalf(2,scoord)
                tmpcoord(3) = evalf(3,scoord)
                tmpcoord = tmpcoord - floor(tmpcoord)
                CALL fractionalToAbsolute(coord,tmpcoord/zunit%dup,zcell)
                DO ib = 1, i
                   dr  = coord-zeo%bead(ib)%coord
                   call mimage(dr(1),dr(2),dr(3),boxZeo)
                   IF (DOT_PRODUCT(dr,dr)<=eps) THEN
                      EXIT
                   END IF
                END DO
                ! If the atom generated is unique let's add to the atom set..
                IF (ib.gt.i) THEN
                   i = i + 1
                   lunitcell(i)=.true.
                   ia=zeo%bead(ic)%type
                   zeo%bead(i)%type=ia
                   ztype%num(ia)=ztype%num(ia)+1
                   zeo%bead(i)%coord=coord
                END IF
             END DO
          END DO
          nAtom = i
          CALL finalizef()
          DO id=1,nAtom
             DO ic=0,zunit%dup(3)-1
                DO ib=0,zunit%dup(2)-1
                   DO ia=0,zunit%dup(1)-1
                      IF (ia.ne.0.or.ib.ne.0.or.ic.ne.0) THEN
                         i=i+1
                         lunitcell(i)=.false.
                         zeo%bead(i)%type=zeo%bead(id)%type
                         ztype%num(zeo%bead(id)%type)=ztype%num(zeo%bead(id)%type)+1
                         call absoluteToFractional(scoord,zeo%bead(id)%coord,zcell)
                         call fractionalToAbsolute(zeo%bead(i)%coord,scoord+real((/ia,ib,ic/),dp)/zunit%dup,zcell)
                      END IF
                   END DO
                END DO
             END DO
          END DO
          nAtom=i
       end if

       if (uninitialized.eq.0) then
          forall(i=1:3) zcell%boxl(i)%val=zunit%boxl(i)*zunit%dup(i)
          call setUpCellStruct(zcell,zunit,lprint)
          uninitialized=-1
       end if
    END DO

    if (nAtom.ne.zeo%nbead) call err_exit(__FILE__,__LINE__,'CIF: Number of atoms incorrect',-1)

  END SUBROUTINE readCIF

! *****************************************************************************
!> \brief  Reads REAL from the CIF file.. This wrapper is needed in order to
!>         treat properly the accuracy specified in the CIF file, i.e. 3.45(6)
! *****************************************************************************
  SUBROUTINE getReal(line,val,err)
    CHARACTER(LEN=*),intent(in)::line
    REAL,intent(out)::val
    INTEGER,INTENT(out)::err
    INTEGER::ib

    ib=INDEX(line,'(')
    IF (ib.gt.0) THEN
       ib = ib - 1
    ELSE
       ib = LEN_TRIM(line)
    END IF
    READ(line(1:ib),*,IOSTAT=err) val

  END SUBROUTINE getReal

END MODULE parser_cif
