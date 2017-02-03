module Rdf
  use global
  implicit none
  save

  logical::lrdf
  integer::rdf_block,rdf_nblock,mtype(2),nbtype(2),btype(2,10)
  integer,allocatable::nframe_eff(:),hbn(:),n(:,:,:)
  real::beadcount(2)
  real,allocatable::rdf_width(:)
  real,dimension(:,:,:),allocatable::numint,g,navg,gavg,navgCom,gavgCom

contains
  subroutine initializeRDF()
    integer::ibox,i,imolty,iunit

    allocate(nframe_eff(nbox),hbn(nbox),n(rdf_block,rdf_nblock,nbox),numint(0:rdf_block,rdf_nblock,nbox),g(rdf_block,rdf_nblock,nbox),navg(rdf_block,rdf_nblock,nbox),gavg(rdf_block,rdf_nblock,nbox),navgCom(rdf_block,rdf_nblock,nbox),gavgCom(rdf_block,rdf_nblock,nbox))

    ! determine number of beads of each type
    ! for each molecule type (SPC/E has 2 beads of type 108, etc)
    beadcount=0
    do i=1,2
       imolty=mtype(i)
       do iunit=1,nunit(imolty)
          if (ANY(beadtype(imolty,iunit).eq.btype(i,1:nbtype(i)))) then
             beadcount(i) = beadcount(i)+1
          end if
       end do
    end do
    if (ANY(beadcount.lt.0.5)) then
       write(6,*) 'Specified bead not found!'
       stop
    end if

    nframe_eff=0
    hbn=0

    navg=0.
    gavg=0.

    navgCom=0.
    gavgCom=0.
  end subroutine initializeRDF

!> \brief Calculate site-site radial distribution function
  subroutine calcRdf()
    integer::i,imolty,iunit,ichain,ibox,jmolty,junit,jchain,jbox,ibin,iblock,nct(2,nbox)
    real::dr,const,nideal

    nct=0
    n=0
    ! Calculate specific pair distances
    ! Loop over chains - btype1 becomes the origin bead
    ! rdf is directional btype1 to btype2, 
    ! although btype2 to btype1 should be same
    ! number integrals will be different depending on direction
    do ichain=1,nchain
       imolty=moltyp(ichain)
       ibox=actualBox(molbox(ichain))

       ! Count number of molecules that could interact (mtype2)
       do i=1,2
          if (imolty.eq.mtype(i)) then
             nct(i,ibox)=nct(i,ibox)+1
          end if
       end do

       ! First find molecule type 1, then set btype1
       if (imolty.ne.mtype(1) ) cycle

       do iunit=1,nunit(imolty)
          if (ALL(beadtype(imolty,iunit).ne.btype(1,1:nbtype(1)))) cycle

          ! Loop over interacting chains - only btype2 can interact
          do jchain=1,nchain
             jmolty=moltyp(jchain)
             jbox=actualBox(molbox(jchain))

             ! Check that the interacting molecule is in the same box and is mtype2
             if (ichain.eq.jchain.or.ibox.ne.jbox.or.jmolty.ne.mtype(2)) cycle
             jbox=molbox(jchain) ! to reference boxl{x,y,z}

             do junit=1,nunit(jmolty)
                if(ALL(beadtype(jmolty,junit).ne.btype(2,1:nbtype(2)))) cycle

                ! Bin distances
                dr=sqrt(square_distance((/rxu(ichain,iunit)-rxu(jchain,junit),ryu(ichain,iunit)-ryu(jchain,junit),rzu(ichain,iunit)-rzu(jchain,junit)/),(/boxlx(jbox),boxly(jbox),boxlz(jbox)/)))
                ibin=int(dr/rdf_width(ibox))+1
                hbn(ibox)=max(hbn(ibox),ibin)
                iblock=ceiling(log(dble(ibin))/log(dble(rdf_block)))
                if (iblock.eq.0) iblock=1
                if (iblock>1) ibin=int(dr/(rdf_width(ibox)*rdf_block**(iblock-1)))+1

                if (iblock.le.rdf_nblock.and.ibin.le.rdf_block) then
                   n(ibin,iblock,ibox)=n(ibin,iblock,ibox)+1
                else
                   write(6,*) 'ichain: ',ichain,', iunit: ',iunit,'; jchain: ',jchain,', junit: ',junit
                   write(6,*) 'ibox: ',jbox,', ibin: ',ibin,', iblock: ',iblock
                end if
             end do ! end loop over units junit
          end do ! end loop over molecules jchain
       end do ! end loop over units iunit
    end do ! end loop over molecules ichain

    numint=0.
    g=0.

    do ibox=1,nbox
       if (nct(1,ibox).eq.0) cycle
       jbox=actualBox(ibox)
       ! Normalize to ideal gas
       const=4.0*onepi*(nct(2,ibox)*dble(beadcount(2))/(boxlx(jbox)*boxly(jbox)*boxlz(jbox)))/3.0

       ! calculate number integral            
       do iblock=1,rdf_nblock
          if (iblock.gt.1) n(1,iblock,ibox)=numint(rdf_block,iblock-1,ibox)
          do ibin=1,rdf_block
             numint(ibin,iblock,ibox)=numint(ibin-1,iblock,ibox)+dble(n(ibin,iblock,ibox))/nct(1,ibox)
          end do
       end do

       if (nct(2,ibox).eq.0) cycle
       nframe_eff(ibox)=nframe_eff(ibox)+1
       ! calculate rdf
       do iblock=1,rdf_nblock
          do ibin=1,rdf_block
             nideal=const*((rdf_width(ibox)*rdf_block**(iblock-1)*ibin)**3-(rdf_width(ibox)*rdf_block**(iblock-1)*(ibin-1))**3)
             g(ibin,iblock,ibox)=n(ibin,iblock,ibox)/dble(beadcount(1))/nct(1,ibox)/nideal
          end do
       end do

       ! accumulate averages
       navg(:,:,ibox)=real(nframe_eff(ibox)-1)/real(nframe_eff(ibox))*navg(:,:,ibox)+numint(:,:,ibox)/nframe_eff(ibox)
       gavg(:,:,ibox)=real(nframe_eff(ibox)-1)/real(nframe_eff(ibox))*gavg(:,:,ibox)+g(:,:,ibox)/nframe_eff(ibox)
    end do ! end loop over boxes

  end subroutine calcRdf

!> \brief Average and output site-site RDF
  subroutine outputRdf()
    integer::ibox,ierr,nblock,iblock,nbin,ibin
    real::rx
    character(len=64)::fname1,fname2

    navg=navg/beadcount(1)
    do ibox=1,nbox
       write(fname1,'("nint",I1,"-",I1,"-",I3.3,"-",I1,"-",I3.3,".dat")') ibox,mtype(1),btype(1,1),mtype(2),btype(2,1)
       write(fname2,'("rdf",I1,"-",I1,"-",I3.3,"-",I1,"-",I3.3,".dat")') ibox,mtype(1),btype(1,1),mtype(2),btype(2,1)
       open(unit=11,access='sequential',action='write',file=fname1,form='formatted',iostat=ierr,status='unknown')
       open(unit=12,access='sequential',action='write',file=fname2,form='formatted',iostat=ierr,status='unknown')

       nblock=ceiling(log(dble(hbn(ibox)))/log(dble(rdf_block)))
       if (nblock.le.0) nblock=1
       nbin=int(hbn(ibox)/(rdf_block**(nblock-1)))+1
       do iblock=1,nblock
          do ibin=1,rdf_block
             if (iblock.eq.nblock.and.ibin.gt.nbin) exit
             if (iblock.gt.1.and.ibin.eq.1) cycle
             rx=rdf_width(ibox)*rdf_block**(iblock-1)*(dble(ibin)-0.5d0)
             write(11,*) rx,navg(ibin,iblock,ibox),' 0.0 ',nframe_eff(ibox)
             write(12,*) rx,gavg(ibin,iblock,ibox),' 0.0 ',nframe_eff(ibox)
          end do
       end do
       close(11)
       close(12)
    end do ! end loop over boxes
  end subroutine outputRdf

!> \brief Calculate COM-COM radial distribution function
  subroutine calcRdfCom()
    integer::i,imolty,iunit,ichain,ibox,jmolty,junit,jchain,jbox,ibin,iblock,nct(2,nbox)
    real::dr,const,nideal

    nct=0
    n=0
    ! Calculate specific pair distances
    ! Loop over chains - btype1 becomes the origin bead
    ! rdf is directional btype1 to btype2, 
    ! although btype2 to btype1 should be same
    ! number integrals will be different depending on direction
    do ichain=1,nchain
       imolty=moltyp(ichain)
       ibox=actualBox(molbox(ichain))

       ! Count number of molecules that could interact (mtype2)
       do i=1,2
          if (imolty.eq.mtype(i)) then
             nct(i,ibox)=nct(i,ibox)+1
          end if
       end do

       ! First find molecule type 1, then set btype1
       if (imolty.ne.mtype(1)) cycle

       ! Loop over interacting chains - only btype2 can interact
       do jchain=1,nchain
          jmolty=moltyp(jchain)
          jbox=actualBox(molbox(jchain))

          ! Check that the interacting molecule is in the same box and is mtype2 
          if (ichain.eq.jchain.or.ibox.ne.jbox.or.jmolty.ne.mtype(2)) cycle
          jbox=molbox(jchain) ! to reference boxl{x,y,z}

          ! Bin distances
          dr=sqrt(square_distance((/xcom(ichain)-xcom(jchain),ycom(ichain)-ycom(jchain),zcom(ichain)-zcom(jchain)/),(/boxlx(jbox),boxly(jbox),boxlz(jbox)/)))
          ibin=int(dr/rdf_width(ibox))+1
          hbn(ibox)=max(hbn(ibox),ibin)
          iblock=ceiling(log(dble(ibin))/log(dble(rdf_block)))
          if (iblock.eq.0) iblock=1
          if (iblock>1) ibin=int(dr/(rdf_width(ibox)*rdf_block**(iblock-1)))+1

          if (iblock.le.rdf_nblock.and.ibin.le.rdf_block) then
             n(ibin,iblock,ibox)=n(ibin,iblock,ibox)+1
          else
             write(6,*) 'ichain: ',ichain,'; jchain: ',jchain
             write(6,*) 'ibox: ',jbox,', ibin: ',ibin,', iblock: ',iblock
          end if
       end do ! end loop over molecules jchain
    end do ! end loop over molecules ichain

    numint=0.
    g=0.

    do ibox=1,nbox
       if (nct(1,ibox).eq.0) cycle
       jbox=actualBox(ibox)
       ! Normalize to ideal gas
       const=4.0*onepi*(nct(2,ibox)/(boxlx(jbox)*boxly(jbox)*boxlz(jbox)))/3.0

       ! calculate number integral            
       do iblock=1,rdf_nblock
          if (iblock.gt.1) n(1,iblock,ibox)=numint(rdf_block,iblock-1,ibox)
          do ibin=1,rdf_block
             numint(ibin,iblock,ibox)=numint(ibin-1,iblock,ibox)+dble(n(ibin,iblock,ibox))/nct(1,ibox)
          end do
       end do

       if (nct(2,ibox).eq.0) cycle
       ! calculate rdf
       do iblock=1,rdf_nblock
          do ibin=1,rdf_block
             nideal=const*((rdf_width(ibox)*rdf_block**(iblock-1)*ibin)**3-(rdf_width(ibox)*rdf_block**(iblock-1)*(ibin-1))**3)
             g(ibin,iblock,ibox)=n(ibin,iblock,ibox)/dble(nct(1,ibox))/nideal
          end do
       end do
    end do ! end loop over boxes

    ! accumulate averages
    forall(ibox=1:nbox,nframe_eff(ibox).gt.0)
       navgCom(:,:,ibox)=real(nframe_eff(ibox)-1)/real(nframe_eff(ibox))*navgCom(:,:,ibox)+numint(:,:,ibox)/nframe_eff(ibox)
       gavgCom(:,:,ibox)=real(nframe_eff(ibox)-1)/real(nframe_eff(ibox))*gavgCom(:,:,ibox)+g(:,:,ibox)/nframe_eff(ibox)
    end forall
  end subroutine calcRdfCom

!> \brief Average and output COM-COM RDF
  subroutine outputRdfCom()
    integer::ibox,ierr,nblock,iblock,nbin,ibin
    real::rx
    character(len=64)::fname1,fname2

    do ibox=1,nbox
       write(fname1,'("nint",I1,"-",I1,"-",I1,".dat")') ibox,mtype(1),mtype(2)
       write(fname2,'("rdf",I1,"-",I1,"-",I1,".dat")') ibox,mtype(1),mtype(2)
       open(unit=11,access='sequential',action='write',file=fname1,form='formatted',iostat=ierr,status='unknown')
       open(unit=12,access='sequential',action='write',file=fname2,form='formatted',iostat=ierr,status='unknown')

       nblock=ceiling(log(dble(hbn(ibox)))/log(dble(rdf_block)))
       if (nblock.le.0) nblock=1
       nbin=int(hbn(ibox)/(rdf_block**(nblock-1)))+1
       do iblock=1,nblock
          do ibin=1,rdf_block
             if (iblock.eq.nblock.and.ibin.gt.nbin) exit
             if (iblock.gt.1.and.ibin.eq.1) cycle
             rx=rdf_width(ibox)*rdf_block**(iblock-1)*(dble(ibin)-0.5d0)
             write(11,*) rx,navgCom(ibin,iblock,ibox),' 0.0 ',nframe_eff(ibox)
             write(12,*) rx,gavgCom(ibin,iblock,ibox),' 0.0 ',nframe_eff(ibox)
          end do
       end do
       close(11)
       close(12)
    end do ! end loop over boxes
  end subroutine outputRdfCom
end module Rdf
