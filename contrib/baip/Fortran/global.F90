module global
  implicit none
  save

  real,parameter::onepi=4.0*atan(1.0),eps=1.0e-6
    
  integer::nframe& !< number frames or snapshots from the simulation
  ,nchain& !< number chains (molecules) in the simulation
  ,nmolty& !< number of molecule types
  ,nbox& !< number of boxes in the simulation
  ,nhere& !< number of bead types
  ,ncycles !< number of cycles at frame was taken

  integer,dimension(:),allocatable::temphere& !< type of bead i
  ,nunit& !< number of beads in the molecule type i
  ,moltyp& !< molecule type of chain i
  ,molbox& !< box of chain i
  ,actualBox !< actualBox(i) always store the box number of a certain identity even the boxes change identities
  integer,dimension(:,:),allocatable::ncmt& !< number of chains j in box i
  ,beadtype& !< type of bead for molecule type i, bead j
  ,nvib& !< number of vibration for chain i bead j
  ,ntor !< number of torsions for chain i bead j
  integer,dimension(:,:,:),allocatable::ijvib& !< bead involved in vibration with chain i, bead j, vib k
  ,itor2& !< 1st bead involved in torsion with chain i, bead j, tors k
  ,itor3& !< 2nd bead involved in torsion with chain i, bead j, tors k
  ,itor4 !< 3rd bead involved in torsion with chain i, bead j, tors k
  real,dimension(:),allocatable::rcut& !< potential cutoff in box i
  ,boxlx& !< x egdelength of box i
  ,boxly& !< y egdelength of box i
  ,boxlz& !< z egdelength of box i
  ,xcom& !< x center of mass position for chain i
  ,ycom& !< y center of mass position for chain i
  ,zcom !< z center of mass position for chain i
  real,dimension(:,:),allocatable::rxu& !< x bead position for bead j of chain i
  ,ryu& !< y bead position for bead j of chain i
  ,rzu& !< z bead position for bead j of chain i
  ,qbead !< charge for bead j of chain i

contains
  function square_distance(dr,boxLength) result(rsq)
    real,intent(in)::dr(3),boxLength(3)
    real::rsq

    real::r(3)

    r=dr-anint(dr/boxLength)*boxLength
    rsq=dot_product(r,r)
  end function square_distance

!> \brief Read movie file header
  subroutine readPreamble(io_movie)
    integer,intent(in)::io_movie

    logical,save::initialized=.false.
    integer::i,imolty,iunit,invib,intor&
    ,most_units& !< most beads in any molecule
    ,most_vibs& !< most torsions in any molecule
    ,most_tors !< most vibrations in any molecule

#ifdef OLD_CODE
    if (.not.initialized) then
       allocate(temphere(1:nhere),rcut(1:nbox))
    end if

    read(io_movie,*) nframe, nchain, nmolty, (rcut(i), i=1,nbox)
    read(io_movie,*) nhere, (temphere(i), i=1,nhere)

    if (.not.initialized) then
       allocate(nunit(1:nmolty),moltyp(1:nchain),molbox(1:nchain),ncmt(1:nbox,1:nmolty),boxlx(1:nbox),boxly(1:nbox),boxlz(1:nbox),xcom(1:nchain),ycom(1:nchain),zcom(1:nchain),actualBox(nbox))

#else
    read(io_movie,*) nframe, nchain, nmolty, nbox, nhere

    if (.not.initialized) then
       allocate(temphere(1:nhere),nunit(1:nmolty),moltyp(1:nchain),molbox(1:nchain),ncmt(1:nbox,1:nmolty),rcut(1:nbox),boxlx(1:nbox),boxly(1:nbox),boxlz(1:nbox),xcom(1:nchain),ycom(1:nchain),zcom(1:nchain),actualBox(nbox))
#endif
       most_units = -1
       most_vibs  = -1
       most_tors  = -1
    end if

#ifndef OLD_CODE
    read(io_movie,*) (rcut(i), i=1,nbox)
    read(io_movie,*) (temphere(i), i=1,nhere)
#endif

    do imolty = 1,nmolty
       read(io_movie,*) nunit(imolty)
       if (.not.initialized) most_units = max(most_units,nunit(imolty))

       do iunit = 1,nunit(imolty)
          if (.not.initialized) then
             read(io_movie,*) invib
             most_vibs = max(most_vibs,invib)
          else
             read(io_movie,*) nvib(imolty,iunit), (ijvib(imolty,iunit,i), i=1,nvib(imolty,iunit))
          end if
       end do

       do iunit=1,nunit(imolty)
          if (.not.initialized) then
             read(io_movie,*) intor
             most_tors = max(most_tors,intor)
#ifdef OLD_CODE
             do i=1,intor
                read(io_movie,*)
             end do
#endif
          else
#ifdef OLD_CODE
             read(io_movie,*) ntor(imolty,iunit)
             do i=1,ntor(imolty,iunit)
                read(io_movie,*) itor2(imolty,iunit,i),itor3(imolty,iunit,i),itor4(imolty,iunit,i)
             end do
#else
             read(io_movie,*) ntor(imolty,iunit),(itor2(imolty,iunit,i),itor3(imolty,iunit,i),itor4(imolty,iunit,i),i=1,ntor(imolty,iunit))
#endif
          end if
       end do
    end do

    if (.not.initialized) then
         allocate(beadtype(1:nmolty,1:most_units),nvib(1:nmolty,1:most_units),ntor(1:nmolty,1:most_units),ijvib(1:nmolty,1:most_units,1:most_vibs),itor2(1:nmolty,1:most_units,1:most_tors),itor3(1:nmolty,1:most_units,1:most_tors),itor4(1:nmolty,1:most_units,1:most_tors),rxu(1:nchain,1:most_units),ryu(1:nchain,1:most_units),rzu(1:nchain,1:most_units),qbead(1:nchain,1:most_units))
         initialized=.true.
    end if
  end subroutine readPreamble

!> \brief Read the iframe-th frame
  subroutine readFrame(io_movie,iframe,incomplete)
    logical,intent(out)::incomplete
    integer,intent(in)::io_movie,iframe

    integer::ibox,imolty,ichain,jchain,iunit

    incomplete=.false.
    
    read(io_movie,*,END=100) ncycles

    do ibox = 1,nbox
       read(io_movie,*,END=100) (ncmt(ibox,imolty), imolty = 1,nmolty)
       read(io_movie,*,END=100) boxlx(ibox), boxly(ibox), boxlz(ibox)
    end do

    do ichain = 1,nchain
       read(io_movie,*,END=100) jchain,moltyp(jchain),nunit(moltyp(jchain)),molbox(jchain),xcom(jchain),ycom(jchain),zcom(jchain)
       if (ichain.ne.jchain) then
          write(6,*) 'ichain: ',ichain,' .ne. jchain: ',jchain,' in frame ',iframe
       end if
       do iunit = 1,nunit(moltyp(jchain))
          read(io_movie,*,END=100) rxu(jchain,iunit),ryu(jchain,iunit)&
           ,rzu(jchain,iunit),qbead(jchain,iunit),beadtype(moltyp(jchain),iunit)
       end do
    end do

    return

100 incomplete=.true.
  end subroutine readFrame
end module global
