module HBond
  use global
  implicit none
  save

  type::List
     integer,allocatable::elem(:)
     integer::nelem=0,max
  end type List

  type::HydrogenBond
     integer,allocatable::OId(:)
     integer::nO=0
     type(List),allocatable::hlist(:)
  end type HydrogenBond

  logical::lhbond
  integer,allocatable::nframe_eff(:,:)
  real::rOOsq,rOHsq,cosOHO
  real,allocatable::hbondavg(:,:,:)
  type(HydrogenBond),allocatable::hbatomdef(:)

contains
  subroutine initializeHBond()
    allocate(nframe_eff(nmolty,nbox),hbondavg(4,nmolty,nbox)) !< 1: total, 2: donor, 3: acceptor, 4: intra
    nframe_eff=0
    hbondavg=0.
  end subroutine initializeHBond

!> \brief Analyze hydrogen bonding
  subroutine calcHBond()
    integer::ichain,imolty,ibox,iOId,iO,iH,iHId,jchain,jmolty,jbox,jOId,jO,jH,jHId
    real::drOOsq,drOHsq,drHOsq,hbondtmp(4,nmolty,nbox)

    hbondtmp=0.
    do ichain = 1,nchain
       imolty = moltyp(ichain)
       ibox = actualBox(molbox(ichain))
       loop_jchain: do jchain = ichain,nchain
          jmolty = moltyp(jchain)
          jbox = actualBox(molbox(jchain))
          if (jbox.ne.ibox) cycle loop_jchain
          jbox=molbox(jchain) ! to reference boxl{x,y,z}

          do iOId=1,hbatomdef(imolty)%nO
             iO=hbatomdef(imolty)%OId(iOId)
             do jOId=1,hbatomdef(jmolty)%nO
                if (ichain.eq.jchain.and.iOId.eq.jOId) cycle
                jO=hbatomdef(jmolty)%OId(jOId)

                drOOsq = square_distance((/rxu(ichain,iO)-rxu(jchain,jO),ryu(ichain,iO)-ryu(jchain,jO),rzu(ichain,iO)-rzu(jchain,jO)/),(/boxlx(jbox),boxly(jbox),boxlz(jbox)/))

                if (drOOsq.lt.rOOsq) then
                   ! ichain-O to jchain-H
                   do jHId=1,hbatomdef(jmolty)%hlist(jOId)%nelem
                      jH=hbatomdef(jmolty)%hlist(jOId)%elem(jHId)
                      drOHsq=square_distance((/rxu(ichain,iO)-rxu(jchain,jH),ryu(ichain,iO)-ryu(jchain,jH),rzu(ichain,iO)-rzu(jchain,jH)/),(/boxlx(jbox),boxly(jbox),boxlz(jbox)/))
                      drHOsq=square_distance((/rxu(jchain,jH)-rxu(jchain,jO),ryu(jchain,jH)-ryu(jchain,jO),rzu(jchain,jH)-rzu(jchain,jO)/),(/boxlx(jbox),boxly(jbox),boxlz(jbox)/))
                      if (drOHsq.lt.rOHsq.and.(drOHsq+drHOsq-drOOsq)/sqrt(drOHsq*drHOsq).lt.2*cosOHO) then
                         if (ichain.ne.jchain) then
                            hbondtmp(3,imolty,ibox)=hbondtmp(3,imolty,ibox)+1.
                            hbondtmp(2,jmolty,ibox)=hbondtmp(2,jmolty,ibox)+1.
                         else
                            hbondtmp(4,imolty,ibox)=hbondtmp(4,imolty,ibox)+1.
                         end if
                      end if
                   end do

                   if (ichain.eq.jchain) cycle

                   ! ichain-H to jchain-O
                   do iHId=1,hbatomdef(imolty)%hlist(iOId)%nelem
                      iH=hbatomdef(imolty)%hlist(iOId)%elem(iHId)
                      drHOsq=square_distance((/rxu(ichain,iH)-rxu(jchain,jO),ryu(ichain,iH)-ryu(jchain,jO),rzu(ichain,iH)-rzu(jchain,jO)/),(/boxlx(jbox),boxly(jbox),boxlz(jbox)/))
                      drOHsq=square_distance((/rxu(ichain,iO)-rxu(ichain,iH),ryu(ichain,iO)-ryu(ichain,iH),rzu(ichain,iO)-rzu(ichain,iH)/),(/boxlx(jbox),boxly(jbox),boxlz(jbox)/))
                      if (drHOsq.lt.rOHsq.and.(drHOsq+drOHsq-drOOsq)/sqrt(drHOsq*drOHsq).lt.2*cosOHO) then
                         hbondtmp(2,imolty,ibox)=hbondtmp(2,imolty,ibox)+1.
                         hbondtmp(3,jmolty,ibox)=hbondtmp(3,jmolty,ibox)+1.
                      end if
                   end do
                end if
             end do
          end do
       end do loop_jchain
    end do

    forall (ibox=1:nbox,imolty=1:nmolty,ncmt(actualBox(ibox),imolty).ne.0)
       nframe_eff(imolty,ibox)=nframe_eff(imolty,ibox)+1
       hbondtmp(2:4,imolty,ibox)=hbondtmp(2:4,imolty,ibox)/ncmt(actualBox(ibox),imolty)
       hbondavg(2:4,imolty,ibox)=real(nframe_eff(imolty,ibox)-1)/real(nframe_eff(imolty,ibox))*hbondavg(2:4,imolty,ibox)+hbondtmp(2:4,imolty,ibox)/nframe_eff(imolty,ibox)
    end forall
  end subroutine calcHBond

  subroutine outputHBond(fname)
    character(len=*),intent(in)::fname

    integer::ierr,ibox,imolty,i
    hbondavg(1,:,:)=hbondavg(2,:,:)+hbondavg(3,:,:)+hbondavg(4,:,:)

    open(unit=13,access='sequential',action='write',file=fname,form='formatted',iostat=ierr,status='unknown')
!    write(13,'("# ibox-imolty-(1,2,3) numHBond stdev nframe (rOO = ",F5.2,", rOH = ",F5.2,", cosOHO = ",F5.2,")")') sqrt(rOOsq),sqrt(rOHsq),cosOHO
    do ibox=1,nbox
       do imolty=1,nmolty
          do i=1,4
             write(13,'(I1,"-",I2.2,"-",I1,1X,G16.9,1X,G16.9,I5)') ibox,imolty,i,hbondavg(i,imolty,ibox),0.0,nframe_eff(imolty,ibox)
          end do
       end do
    end do
    close(13)
  end subroutine outputHBond
end module HBond
