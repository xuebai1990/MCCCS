subroutine read10() !
  !*****************************************************************************
  ! Reads all variables from fort.10 movie file
  !*****************************************************************************

  use fort10_vars !stores all variables from fort.10

  implicit none

  !Counters and array sizes
  integer::ibox,ihere,imolty,iunit,ivib,itor,iframe,ichain,invib,intor,jchain
  logical::ltest
  ltest = .false.

  !-----READ IN VARIABLES FROM FORT.10 NEEDED FOR ARRAY ALLOCATION 
  write(6,*)
  write(6,*) '   Reading fort.10 ...'

  allocate(rcut(1:nbox))
  allocate(maxl(nbox))

  read(10,*) nframe, nchain, nmolty, (rcut(ibox), ibox=1,nbox)
  read(10,*) nhere
  close(10)

  allocate(temphere(1:nhere))
  allocate(nunit(1:nmolty))

  allocate(ncycles(1:nframe))

  !number count stuff
  allocate(ncmt(1:nframe,1:nbox,1:nmolty))

  !box_lenghts
  allocate(boxlx(1:nframe,1:nbox))
  allocate(boxly(1:nframe,1:nbox))
  allocate(boxlz(1:nframe,1:nbox))

  !-----Chain data
  !box location
  allocate(molbox(1:nframe, 1:nchain))

  !COM coordinate arrays
  allocate(xcom(1:nframe,1:nchain))
  allocate(ycom(1:nframe,1:nchain))
  allocate(zcom(1:nframe,1:nchain))
  allocate(moltyp(1:nchain))

  most_units = -1
  most_vibs  = -1
  most_tors  = -1

  read(10,*) nframe, nchain, nmolty, (rcut(ibox), ibox=1,nbox)
  read(10,*) nhere, (temphere(ihere), ihere=1,nhere)

  loop_molty1: do imolty = 1,nmolty
     read(10,*) nunit(imolty)
     most_units = max(most_units,nunit(imolty))

     loop_unit1A: do iunit = 1,nunit(imolty)
        read(10,*) invib
        most_vibs = max(most_vibs,invib)
     enddo loop_unit1A

     loop_unit1B: do iunit=1,nunit(imolty)
        read(10,*) intor
        most_tors = max(most_tors,intor)
        loop_tor1: do itor = 1,intor
           read(10,*) !itor2....
        enddo loop_tor1
     enddo loop_unit1B
  enddo loop_molty1

  close(10)

  !----FINISH ALLOCATING THE REST OF THE ARRAYS SO REST OF FORT.10
  !----MAY BE READ IN
  allocate(beadtype(1:nmolty, 1:most_units))

  !vibration stuff
  allocate(nvib(1:nmolty,1:most_units))
  allocate(ijvib(1:nmolty,1:most_units,1:most_vibs))

  !torsion stuff
  allocate(ntor(1:nmolty,1:most_units))
  allocate(itor2(1:nmolty,1:most_units,1:most_tors))
  allocate(itor3(1:nmolty,1:most_units,1:most_tors))
  allocate(itor4(1:nmolty,1:most_units,1:most_tors))

  !-----BEAD DATA
  !bead coordinate arrays etc
  allocate(xbead(1:nframe, 1:nchain, 1:most_units))
  allocate(ybead(1:nframe, 1:nchain, 1:most_units))
  allocate(zbead(1:nframe, 1:nchain, 1:most_units))
  allocate(qbead(1:nframe, 1:nchain, 1:most_units))

  !-----FINISH READING FORT.10 (from the top again!)
  read(10,*) nframe, nchain, nmolty, (rcut(ibox), ibox=1,nbox)
  read(10,*) nhere, (temphere(ihere), ihere=1,nhere)

  loop_molty2: do imolty = 1,nmolty
     read(10,*) nunit(imolty)

     loop_units2A: do iunit = 1,nunit(imolty)
        read(10,*) nvib(imolty,iunit)  &
             , (ijvib(imolty,iunit,ivib), ivib=1,nvib(imolty,iunit))
     enddo loop_units2A

     loop_units2B: do iunit = 1,nunit(imolty)
        read(10,*) ntor(imolty,iunit)
        loop_tors2: do itor = 1,ntor(imolty,iunit)
           read(10,*) itor2(imolty,iunit,itor)  &
                ,itor3(imolty,iunit,itor)  &
                ,itor4(imolty,iunit,itor)  
        enddo loop_tors2
     enddo loop_units2B
  enddo loop_molty2

  !-----READ IN INFORMATION FROM EACH OF THE FRAMES
  maxl=0.
  loop_frame4: do iframe = 1,nframe
     read(10,*,END=100) ncycles(iframe)
     loop_box4: do ibox = 1,nbox
        read(10,*,END=100) (ncmt(iframe,ibox,imolty), imolty = 1,nmolty)
        read(10,*,END=100) boxlx(iframe,ibox), boxly(iframe,ibox), boxlz(iframe,ibox)
        maxl(ibox)=max(maxl(ibox),boxlx(iframe,ibox))
        maxl(ibox)=max(maxl(ibox),boxly(iframe,ibox))
        maxl(ibox)=max(maxl(ibox),boxlz(iframe,ibox))
     enddo loop_box4

     loop_chains4: do ichain = 1,nchain
        read(10,*,END=100) jchain, moltyp(jchain),nunit(moltyp(jchain))  &
             ,molbox(iframe,jchain),xcom(iframe,jchain)  &
             ,ycom(iframe,jchain), zcom(iframe,jchain)  
        if (ichain.ne.jchain) then
           write(6,*) 'ichain: ',ichain,' .ne. jchain: ',jchain,' in frame ',iframe
           stop
        end if
        loop_units4: do iunit = 1,nunit(moltyp(jchain))
           read(10,*,END=100) xbead(iframe,jchain,iunit)  &
                ,ybead(iframe,jchain,iunit)  &
                ,zbead(iframe,jchain,iunit)  &
                ,qbead(iframe,jchain,iunit)  &
                ,beadtype(moltyp(jchain),iunit)
        enddo loop_units4
     enddo loop_chains4
  enddo loop_frame4
  close(10)

100 if (iframe.le.nframe) then
     nframe=iframe-1
     write(6,*) 'simulation not complete, nframe = ',nframe
  end if

end subroutine read10
