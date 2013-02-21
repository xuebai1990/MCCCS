MODULE parser_pdb
! *****************************************************************************
!> \brief Handles PDB (Protein Data Bank) files
!>
!> Protein Data Bank Contents Guide:
!> Atomic Coordinate Entry Format Version 3.3 (July, 2011)
!>   http://www.wwpdb.org/documentation/format33/v3.3.html
!>
!> COLUMNS        DATA  TYPE    FIELD        DEFINITION
!> -------------------------------------------------------------------------------------
!>  1 -  6        Record name   "ATOM  "
!>  7 - 11        Integer       serial       Atom serial number.
!> 13 - 16        Atom          name         Atom name.
!> 17             Character     altLoc       Alternate location indicator.
!> 18 - 20        Residue name  resName      Residue name.
!> 22             Character     chainID      Chain identifier.
!> 23 - 26        Integer       resSeq       Residue sequence number.
!> 27             AChar         iCode        Code for insertion of residues.
!> 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
!> 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
!> 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
!> 55 - 60        Real(6.2)     occupancy    Occupancy.
!> 61 - 66        Real(6.2)     tempFactor   Temperature factor.
!> 77 - 78        LString(2)    element      Element symbol, right-justified.
!> 79 - 80        LString(2)    charge       Charge on the atom.
!> -------------------------------------------------------------------------------------
!>  1 -  6        Record name    "CONECT"
!>  7 - 11        Integer        serial       Atom serial number
!> 12 - 16        Integer        serial       Serial number of bonded atom
!> 17 - 21        Integer        serial       Serial number of bonded atom
!> 22 - 26        Integer        serial       Serial number of bonded atom
!> 27 - 31        Integer        serial       Serial number of bonded atom
! *****************************************************************************
  use var_type,only:default_string_length
  use util_runtime,only:err_exit
  use util_files,only:get_iounit,readLine
  use sim_cell
  use sim_particle
  use sim_zeolite,only:ZeoliteUnitCellGridType,ZeoliteBeadType,setUpAtom,setUpCellStruct
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: readPDB,writePDB

CONTAINS

  SUBROUTINE readPDB(filePDB,zeo,lunitcell,ztype,zcell,zunit)
    character(LEN=*),intent(in)::filePDB
    type(MoleculeType),intent(out)::zeo
    logical,allocatable,intent(out)::lunitcell(:)
    type(ZeoliteBeadType),intent(out)::ztype
    type(CellMaskType),intent(out)::zcell ! the "type" component of (Cell)zcell is the index in (ZeoliteBeadType)ztype
    type(ZeoliteUnitCellGridType),intent(out)::zunit

    INTEGER::IOPDB,jerr,i,resID
    CHARACTER(LEN=default_string_length)::line,atomName,resName,molName,elem
    logical::uninitialized
    real::occup,beta

    IOPDB=get_iounit()
    open(unit=IOPDB,access='sequential',action='read',file=filePDB,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit('cannot open zeolite PDB file')
    end if

    CALL readLine(IOPDB,line,.false.,jerr)
    IF(jerr.ne.0) call err_exit('wrong PDB file format')

    read(line(2:),*) zeo%nbead,zunit%dup(1),zunit%dup(2),zunit%dup(3),ztype%ntype
    allocate(zeo%bead(zeo%nbead),lunitcell(zeo%nbead),ztype%name(ztype%ntype),ztype%radiisq(ztype%ntype),ztype%type(ztype%ntype),ztype%num(ztype%ntype),stat=jerr)
    if (jerr.ne.0) call err_exit('readPDB: allocation failed')

    do i=1,ztype%ntype
       CALL readLine(IOPDB,line,.false.,jerr)
       IF(jerr.ne.0) call err_exit('wrong PDB file format')
       read(line(2:),*) ztype%name(i),ztype%type(i),ztype%radiisq(i)
       ztype%radiisq(i)=ztype%radiisq(i)*ztype%radiisq(i)
       ztype%num(i)=0
    end do

    uninitialized=.true.
    i=0
    DO
       CALL readLine(IOPDB,line,.true.,jerr)
       IF(jerr.ne.0) EXIT

       SELECT CASE (line(1:6))
       CASE ("CRYST1")
          read(line(7:),*) zcell%boxl(1)%val,zcell%boxl(2)%val,zcell%boxl(3)%val,zcell%ang(1)%val,zcell%ang(2)%val,zcell%ang(3)%val
          call setUpCellStruct(zcell,zunit)
          uninitialized=.false.
       CASE ("ATOM","HETATM")
          if (uninitialized) call err_exit('PDB: CRYST1 needs to be defined before ATOM')
          i = i + 1
          READ(line(13:16),*) atomname
          READ(line(18:20),*,IOSTAT=jerr) resName
          READ(line(23:26),*,IOSTAT=jerr) resID
          READ(line(31:38),*,IOSTAT=jerr) zeo%bead(i)%coord(1)
          READ(line(39:46),*,IOSTAT=jerr) zeo%bead(i)%coord(2)
          READ(line(47:54),*,IOSTAT=jerr) zeo%bead(i)%coord(3)
          READ(line(55:60),*,IOSTAT=jerr) occup
          READ(line(61:66),*,IOSTAT=jerr) beta
          READ(line(73:76),*,IOSTAT=jerr) molName
          READ(line(77:78),*,IOSTAT=jerr) elem

          IF (LEN_TRIM(elem).eq.0) THEN
             ! Element is assigned on the basis of the atomName
             elem = atomname
          END IF
          IF (LEN_TRIM(molName).eq.0) THEN
             ! If molname is missing (as in the PDB generated by VMD) let's
             ! use the resname for the molname
             molName =  resname
          END IF

          call setUpAtom(elem,i,zeo,lunitcell,ztype,zcell,zunit)
       CASE ("END","TER")
          EXIT
       CASE ("REMARK")
       CASE DEFAULT
       END SELECT
    END DO

    if (i.ne.zeo%nbead) call err_exit('PDB: Number of atoms incorrect')

  END SUBROUTINE readPDB

! *****************************************************************************
  SUBROUTINE writePDB (filePDB,iBox)
    use sim_system
    character(LEN=*),intent(in)::filePDB
    integer,intent(in)::iBox

    INTEGER::IOPDB,jerr,iChain,iUnit,iAtom,iType

    if (myid.ne.0) return

    IOPDB=get_iounit()
    open(unit=IOPDB,access='sequential',action='write',file=filePDB,form='formatted',iostat=jerr,status='unknown')
    if (jerr.ne.0) then
       call err_exit('cannot open file for writing (PDB)')
    end if

    WRITE(IOPDB,'(A6,3(F9.3),3(F7.2),1X,A10)') "CRYST1",boxlx(ibox),boxly(ibox),boxlz(ibox),cell_ang(ibox,1),cell_ang(ibox,2),cell_ang(ibox,3),"P1        "

    iAtom=0
    DO iChain=1,nChain
       IF (nboxi(iChain).ne.iBox) CYCLE
       iType=molTyp(iChain)
       DO iUnit=1,nUnit(iType)
          iAtom=iAtom+1
          WRITE(IOPDB,'(A6,I5,1X,A4,1X,I3,2X,I4,4X,3(F8.3),2(F6.2),10X,A2)') "ATOM  ",iAtom,ADJUSTL(chemId(nType(iType,iUnit))),iType,iChain,rxu(iChain,iUnit),ryu(iChain,iUnit),rzu(iChain,iUnit),1.0,0.0,ADJUSTL(chemId(nType(iType,iUnit)))
       END DO
    END DO

    WRITE(IOPDB,'(A3,/,A3)') "TER","END"

  END SUBROUTINE writePDB

END MODULE parser_pdb
