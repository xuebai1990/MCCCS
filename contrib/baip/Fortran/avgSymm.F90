!> \file Average values at equivalent points due to space group symmetry

program avgSymm
  use var_type,only:dp,default_string_length,default_path_length
  use util_runtime,only:err_exit
  use util_files,only:get_iounit
  implicit none

  integer,parameter::maxNumSymmOp=200

  type::Coordinate
     real::coord(3)
  end type Coordinate

  type::AtomList
     integer::natom
     type(Coordinate)::atoms(maxNumSymmOp)
  end type AtomList

  CHARACTER(LEN=default_string_length)::SymmOp(maxNumSymmOp,3)
  character(LEN=default_path_length)::fname
  real,allocatable::data(:,:,:),ncount(:,:,:)
  integer::ncell(3),ngr(3),mlen(3),nSymm

  call read_cfg('average-symm.cfg')
  call avg_data()
  call output_data(fname)

contains
  subroutine read_cfg(file_cfg)
    character(LEN=*),intent(in)::file_cfg
    integer::io_cfg,jerr

    io_cfg=get_iounit()
    open(unit=io_cfg,access='sequential',action='read',file=file_cfg,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'read_cfg: cannot open average-symm config file',jerr)
    end if

    read(io_cfg,*) ! ncell ngr mlen path_to_cif
    read(io_cfg,*) ncell,ngr,mlen,fname
    write(*,'(A,3(1X,I0))') 'ncell:',ncell
    write(*,'(A,3(1X,I0))') 'ngr:',ngr
    write(*,'(A,3(1X,I0))') 'mlen:',mlen
    write(*,'(A,A)') 'Symmetry operations from: ',trim(fname)

    if (any(mod(ngr,ncell).ne.0)) call err_exit(__FILE__,__LINE__,'read_cfg: ngr not divisible by ncell',jerr)
    if (any(mod(ngr/ncell,mlen).ne.0)) call err_exit(__FILE__,__LINE__,'read_cfg: ngr/ncell not divisible by mlen',jerr)

    allocate(data(0:ngr(1)-1,0:ngr(2)-1,0:ngr(3)-1),ncount(0:ngr(1)-1,0:ngr(2)-1,0:ngr(3)-1),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'read_cfg: memory allocation',jerr)

    call read_cif(trim(fname))

    read(io_cfg,*) ! path_to_data
    read(io_cfg,*) fname
    write(*,'(A,A)') 'Raw data from: ',trim(fname)

    call read_data(trim(fname))
    write(fname,'(A,A)') 'avg-',trim(fname)
    write(*,'(A,A)') 'Averaged data written to: ',trim(fname)

    close(io_cfg)
  end subroutine read_cfg

  subroutine read_cif(file_cif)
    use util_files,only:readLine
    character(LEN=*),intent(in)::file_cif
    CHARACTER(LEN=default_string_length)::line
    INTEGER::io_cif,jerr,i,ia,ib,ic,id,igr(3)

    io_cif=get_iounit()
    open(unit=io_cif,access='sequential',action='read',file='zeolite.cif',form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'read_cif: cannot open zeolite CIF file',jerr)
    end if

    DO
       CALL readLine(io_cif,line,.true.,jerr)
       if (jerr.ne.0) EXIT

       if (index(line,'_symmetry_equiv_pos_as_xyz').gt.0) then ! Read in symmetry operations
          i  = 0
          DO
             CALL readLine(io_cif,line,.true.,jerr)
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
       end if
    END DO

    close(io_cif)
  end subroutine read_cif

  subroutine read_data(file_data)
    character(LEN=*),intent(in)::file_data
    integer::io_data,jerr,i,j,k

    io_data=get_iounit()
    open(unit=io_data,access='sequential',action='read',file=file_data,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'read_data: cannot open data file',jerr)
    end if

    do k = 0,ngr(3)-1
       do j = 0,ngr(2)-1
          do i = 0,ngr(1)-1
             read(io_data,*) data(i,j,k),ncount(i,j,k)
          end do
       end do
    end do

    close(io_data)
  end subroutine read_data

  subroutine merge_data(halfn,equivalent_length,jump,equivalent_vector)
    real::dataSum,ncountSum
    integer,intent(in)::halfn(3),equivalent_length(3),jump(3),equivalent_vector(3)
    integer::ia,ib,ic,ja,jb,jc,igr(3)

    do ia=0,halfn(1)-1
       do ib=0,halfn(2)-1
          do ic=0,halfn(3)-1
             ncountSum=0._dp
             dataSum=0._dp
             do ja=0,equivalent_length(1)-1
                do jb=0,equivalent_length(2)-1
                   do jc=0,equivalent_length(3)-1
                      igr=(/ia,ib,ic/)*jump+(/ja,jb,jc/)*equivalent_vector
                      if (ncount(igr(1),igr(2),igr(3)).ne.0) then
                         ncountSum=ncountSum+ncount(igr(1),igr(2),igr(3))
                         dataSum=dataSum+data(igr(1),igr(2),igr(3))*ncount(igr(1),igr(2),igr(3))
                      else if (data(igr(1),igr(2),igr(3)).ne.0) then
                         write(*,'(A,9(1X,I0),1X,G16.9)') 'impossible 1:',ia,ib,ic,ja,jb,jc,igr,data(igr(1),igr(2),igr(3))
                      end if
                   end do
                end do
             end do

             if (ncountSum.ne.0) then
                data(ia,ib,ic)=dataSum/ncountSum
                ncount(ia,ib,ic)=ncountSum
             end if
          end do
       end do
    end do
    ngr = halfn
  end subroutine merge_data

  subroutine avg_data()
    use const_math,only:eps
    use fparser,only:initf,parsef,evalf,finalizef
    TYPE(AtomList)::stored
    REAL::scoord(3),tmpcoord(3),dr(3),dataAvg,ncountAvg
    INTEGER::igr(3),ia,ib,ic,i,id
    LOGICAL::processed(0:ngr(1)-1,0:ngr(2)-1,0:ngr(3)-1)

    write(*,'(A,G16.9)') 'Raw data: sum(ncount)=',sum(ncount(0:ngr(1)-1,0:ngr(2)-1,0:ngr(3)-1))
    write(*,'(A,3(1X,I0))') 'ngr:',ngr

    if (any(ncell.gt.1)) then
       igr = ngr / ncell
       call merge_data(igr,ncell,(/1,1,1/),igr)
    end if

    write(*,'(A,G16.9)') 'After super cell folding: sum(ncount)=',sum(ncount(0:ngr(1)-1,0:ngr(2)-1,0:ngr(3)-1))
    write(*,'(A,3(1X,I0))') 'ngr:',ngr

    CALL initf(3)

    ! average 3D data
    processed=.false.
    DO ia=0,ngr(1)-1
       DO ib=0,ngr(2)-1
          DO ic=0,ngr(3)-1
             IF (.NOT. processed(ia,ib,ic)) THEN
                scoord=real((/ia,ib,ic/)+0.5,dp)/ngr

                ncountAvg=0._dp
                dataAvg=0._dp
                stored%natom=0
                ! Apply symmetry elements and generate all equivalent points in the unit cell
                DO id = 1, nSymm
                   CALL parsef(1,SymmOp(id,1),(/'x','y','z'/))
                   CALL parsef(2,SymmOp(id,2),(/'x','y','z'/))
                   CALL parsef(3,SymmOp(id,3),(/'x','y','z'/))
                   tmpcoord(1) = evalf(1,scoord)
                   tmpcoord(2) = evalf(2,scoord)
                   tmpcoord(3) = evalf(3,scoord)
                   tmpcoord = tmpcoord - floor(tmpcoord)
                   DO i = 1, stored%natom
                      dr  = tmpcoord-stored%atoms(i)%coord
                      dr=dr-anint(dr)
                      IF (ALL(abs(dr)<=eps)) THEN
                         EXIT
                      END IF
                   END DO
                   ! If the atom generated is unique let's add to the atom set..
                   IF (i.gt.stored%natom) THEN
                      stored%natom = stored%natom + 1
                      stored%atoms(stored%natom)%coord=tmpcoord
                      igr=tmpcoord*ngr
                      processed(igr(1),igr(2),igr(3))=.true.

                      if (ncount(igr(1),igr(2),igr(3)).ne.0) then
                         ncountAvg=ncountAvg+ncount(igr(1),igr(2),igr(3))
                         dataAvg=dataAvg+data(igr(1),igr(2),igr(3))*ncount(igr(1),igr(2),igr(3))
                      else if (data(igr(1),igr(2),igr(3)).ne.0) then
                         write(*,'(A,6(1X,I0),1X,G16.9)') 'impossible 2:',ia,ib,ic,igr,data(igr(1),igr(2),igr(3))
                      end if
                   END IF
                END DO

                if (ncountAvg.ne.0) then
                   dataAvg=dataAvg/ncountAvg
                   ncountAvg=ncountAvg/real(stored%natom,dp)
                end if

                DO i = 1, stored%natom
                   igr=stored%atoms(i)%coord*ngr
                   ncount(igr(1),igr(2),igr(3))=ncountAvg
                   data(igr(1),igr(2),igr(3))=dataAvg
                END DO
             END IF
          END DO
       END DO
    END DO

    CALL finalizef()

    write(*,'(A,G16.9)') 'After symmetry averaging: sum(ncount)=',sum(ncount(0:ngr(1)-1,0:ngr(2)-1,0:ngr(3)-1))
    write(*,'(A,3(1X,I0))') 'ngr:',ngr

    if (any(mlen.gt.1)) then
       call merge_data(ngr/mlen,mlen,mlen,(/1,1,1/))
    end if

    write(*,'(A,G16.9)') 'After neighbor merging: sum(ncount)=',sum(ncount(0:ngr(1)-1,0:ngr(2)-1,0:ngr(3)-1))
    write(*,'(A,3(1X,I0))') 'ngr:',ngr

    !ncountAvg=sum(ncount)
    !if (ncountAvg.gt.0) ncount=ncount/ncountAvg

  end subroutine avg_data

  subroutine output_data(fname)
    character(LEN=*),intent(in)::fname
    integer::io_output,jerr,i,j,k

    io_output=get_iounit()
    open(unit=io_output,access='sequential',action='write',file=fname,form='formatted',iostat=jerr,status='unknown')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'output_data: cannot open write to disk',jerr)
    end if

    do k=0,ngr(3)-1
       do j=0,ngr(2)-1
          do i=0,ngr(1)-1
             write(io_output,'(G16.9)') data(i,j,k)
          end do
       end do
    end do

    close(io_output)
  end subroutine output_data
END program
