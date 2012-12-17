program analysis
  use global
  use rdf
  use hbond
  implicit none

  ! Local variables
  logical::incomplete
  integer::ierr,iframe

  ! Read simulation parameters from the movie file and allocate arrays
  open(unit=10,access='sequential',action='read',file='fort.10',form='formatted',iostat=ierr,status='old')
  call readPreamble(10)

  ! Read analysis input file
  call readConfig()

  ! Read and analyze frame by frame
  rewind(10)
  call readPreamble(10)
  actualBox=(/1,2,3/)
  do iframe=1,nframe
     call readFrame(10,iframe,incomplete)
     !>>> If simulation boxes can change identities, uncomment the following lines so that the mole fraction of molecule type 1 in box 1* is always larger than its mole fraction in box 2*; i.e., say if molecule type 1 is water, then box 1* always contain the results for the aqueous phase <<<
     ! if (real(ncmt(1,1),8)/sum(ncmt(1,1:nmolty)).ge.real(ncmt(2,1),8)/sum(ncmt(2,1:nmolty))) then
     !    actualBox(1)=1
     !    actualBox(2)=2
     !    actualBox(3)=3
     ! else
     !    actualBox(1)=2
     !    actualBox(2)=1
     !    actualBox(3)=3
     ! end if

     if (incomplete) exit
     if (iframe.eq.1) then
        if (lrdf) call initializeRDF()
        if (lhbond) call initializeHBond()
     end if
     if (lrdf) then
        call calcRdf()
        call calcRdfCom()
     end if
     if (lhbond) call calcHBond()
  end do
  close(10)

  if (incomplete) then
     nframe=iframe-1
     write(*,*) 'simulation not complete, nframe = ',nframe
     if (nframe.eq.0) stop
  end if

  if (lrdf) then
     call outputRdf()
     call outputRdfCom()
  end if
  if (lhbond) call outputHBond('hbond.dat')

contains
!> \brief Read analysis input file
  subroutine readConfig()
    integer::ierr,imolty,ielem

    open(unit=100,access='sequential',action='read',file='analysis.cfg',form='formatted',iostat=ierr,status='old')
    !variables for radial distribution function
    allocate(rdf_width(1:nbox),hbatomdef(1:nmolty))
    read(100,*)
    read(100,*) lrdf,rdf_block,rdf_nblock,rdf_width
    read(100,*) 
    read(100,*) mtype(1),nbtype(1),(btype(1,imolty),imolty=1,nbtype(1))
    read(100,*) 
    read(100,*) mtype(2),nbtype(2),(btype(2,imolty),imolty=1,nbtype(2))
    !variables for hbond profile
    read(100,*) 
    read(100,*) lhbond,rOOsq,rOHsq,cosOHO
    read(100,*) 

    rOOsq=rOOsq*rOOsq
    rOHsq=rOHsq*rOHsq

    do imolty=1,nmolty
       read(100,*) hbatomdef(imolty)%nO
       allocate(hbatomdef(imolty)%OId(1:hbatomdef(imolty)%nO),hbatomdef(imolty)%hlist(1:hbatomdef(imolty)%nO))
       read(100,*) hbatomdef(imolty)%OId,(hbatomdef(imolty)%hlist(ielem)%nelem,ielem=1,hbatomdef(imolty)%nO)
       do ielem=1,hbatomdef(imolty)%nO
          allocate(hbatomdef(imolty)%hlist(ielem)%elem(1:hbatomdef(imolty)%hlist(ielem)%nelem))
       end do
       read(100,*) (hbatomdef(imolty)%hlist(ielem)%elem,ielem=1,hbatomdef(imolty)%nO)
    end do
  end subroutine readConfig
end program analysis
