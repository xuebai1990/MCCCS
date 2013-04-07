MODULE moves_volume
  use util_random,only:random
  use util_runtime,only:err_exit
  use sim_system
  use sim_particle,only:rebuild_neighbor_list
  use sim_cell
  use energy_kspace,only:recip,calp,save_kvector,restore_kvector
  use energy_pairwise,only:sumup
  implicit none
  private
  save
  public::volume_1box,volume_2box,init_moves_volume,update_volume_max_displacement,output_volume_stats,read_checkpoint_volume,write_checkpoint_volume

  real,allocatable,public::acsvol(:),acnvol(:),acshmat(:,:),acnhmat(:,:),bsvol(:),bnvol(:),bshmat(:,:),bnhmat(:,:)
  real,allocatable::vboxn(:,:),vboxo(:,:),bxo(:),byo(:),bzo(:),xcmo(:),ycmo(:),zcmo(:),rxuo(:,:),ryuo(:,:),rzuo(:,:),qquo(:,:)
  integer,allocatable::neigho_cnt(:),neigho(:,:)
  real::hmato(9),hmatio(9)

contains
!***********************************************************
! makes an isotropic volume change for NVT-Gibbs ensemble **
! the maximum change is controlled by rmtrax and the      **
! number of successful trial moves is stored in bsvol.    **
!***********************************************************
!
! perform change of the volume: random walk in ln(V1/V2) with V1+V2=const
!
  subroutine volume_2box()
    real::rpair,rm,rbox,volo(nbxmax),volt,voln(nbxmax),rbcut(nbxmax),dfac(nbxmax),df,dx,dy,dz,expdv,min_boxl,v(nEnergy),dele
    integer::ipair,ipairb,boxa,boxb,ibox,i,hbox,jbox,jhmat,imolty,j,ichoiq
    logical::lncubic,lx(nbxmax),ly(nbxmax),lz(nbxmax),ovrlap
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start VOLUME_2BOX in ',myid
#endif

    ! select pair of boxes to do the volume move
    if ( nvolb .gt. 1 ) then
       rpair = random(-1)
       do ipair = 1, nvolb
          if ( rpair .lt. pmvolb(ipair) ) then
             ipairb = ipair
             exit
          end if
       end do
    else
       ipairb = 1
    end if
    boxa = box5(ipairb)
    boxb = box6(ipairb)

    bnvol(ipairb) = bnvol(ipairb) + 1.0E0_dp

    call save_box(boxa)
    call save_box(boxb)
    call save_configuration((/boxa,boxb/))

    lncubic = .false.
    lx = .false.
    ly = .false.
    lz = .false.
    do ibox = 1, 2
       if (ibox .eq. 1) i = boxa
       if (ibox .eq. 2) i = boxb

       if (lsolid(i)) then
          ! volume move independently in x, y, z directions
          rm = random(-1)
          if ( rm .le. pmvolx ) then
             lx(i) = .true.
          else if ( rm .le. pmvoly ) then
             ly(i) = .true.
          else
             lz(i) = .true.
          end if

          if (.not. lrect(i)) then
             lncubic = .true.
             hbox = i
             volo(i) = cell_vol(i)

             ! select one of the cell edge
             if ( lx(i) ) then
                rbox = 3.0_dp*random(-1)
                if ( rbox .lt. 1.0E0_dp ) then
                   jhmat = 1
                else if (rbox .lt. 2.0E0_dp ) then
                   jhmat = 4
                else
                   jhmat = 7
                end if
             else if ( ly(i) ) then
                rbox = 2.0_dp*random(-1)
                if ( rbox .lt. 1.0E0_dp ) then
                   jhmat = 5
                else
                   jhmat = 8
                end if
             else
                jhmat = 9
             end if
          end if
       end if

       if (.not.lsolid(i).or.lrect(i)) then
          jbox = i
          if ( lpbcz ) then
             volo(i) = bxo(i)*byo(i)*bzo(i)
          else
             volo(i) = bxo(i)*byo(i)
          end if
       end if
    end do

    ! calculate total volume
    volt = volo(boxa) + volo(boxb)

    if ( lncubic ) then
       ! hbox is the non-orthorhombic box, jbox is the other box
       bnhmat(hbox,jhmat) = bnhmat(hbox,jhmat) + 1.0E0_dp
       hmat(hbox,jhmat) = hmat(hbox,jhmat) + rmhmat(hbox,jhmat)*( 2.0E0_dp*random(-1) - 1.0E0_dp )
       call matops(hbox)

       voln(hbox) = cell_vol(hbox)
       voln(jbox) = volt-voln(hbox)
       rbcut(hbox) = 2.0_dp*rcut(hbox)
       rbcut(jbox) = 2.0_dp*rcut(jbox)

       if (lsolid(jbox)) then
          ! volume move independently in x, y, z directions
          dfac(jbox)=voln(jbox)/volo(jbox)
          if (lx(jbox)) boxlx(jbox) = boxlx(jbox) * dfac(jbox)
          if (ly(jbox)) boxly(jbox) = boxly(jbox) * dfac(jbox)
          if (lz(jbox)) boxlz(jbox) = boxlz(jbox) * dfac(jbox)
       else
          if ( lpbcz ) then
             dfac(jbox)= (voln(jbox)/volo(jbox))**(1.0E0_dp/3.0E0_dp)
             boxlz(jbox) = boxlz(jbox)*dfac(jbox)
          else
             dfac(jbox)= sqrt(voln(jbox)/volo(jbox))
          end if
          boxlx(jbox) = boxlx(jbox)*dfac(jbox)
          boxly(jbox) = boxly(jbox)*dfac(jbox)
       end if

       if (ANY(rbcut(hbox) .gt. min_width(hbox,:)) .or. boxlx(jbox) .lt. rbcut(jbox) .or. boxly(jbox) .lt. rbcut(jbox) .or. (lpbcz .and. boxlz(jbox) .lt. rbcut(jbox))) then
          hmat(hbox,jhmat) = hmato(jhmat)
          boxlx(jbox) = bxo(jbox)
          boxly(jbox) = byo(jbox)
          if ( lpbcz ) then
             boxlz(jbox) = bzo(jbox)
          end if
          call dump('final-config')
          write(io_output,*) 'w1:',min_width(hbox,1),'w2:',min_width(hbox,2),'w3:',min_width(hbox,3)
          call err_exit(__FILE__,__LINE__,'non-orthorhombic volume_2box move rejected. box width below cutoff size',myid+1)
       end if

       ! determine the displacement of the COM
       df = dfac(jbox) - 1.0E0_dp
       do i = 1,nchain
          ibox = nboxi(i)
          imolty = moltyp(i)
          if (ibox .eq. hbox) then
             if ( lx(ibox) ) then
                dx = sxcm(i)*(hmat(hbox,1)-hmato(1))+sycm(i)*(hmat(hbox,4)-hmato(4))+szcm(i)*(hmat(hbox,7)-hmato(7))
                xcm(i) = xcm(i) + dx
                do j = 1, nunit(imolty)
                   rxu(i,j) = rxu(i,j) + dx
                end do
             else if ( ly(ibox) ) then
                dy = sycm(i)*(hmat(hbox,5)-hmato(5))+szcm(i)*(hmat(hbox,8)-hmato(8))
                ycm(i) = ycm(i) + dy
                do j = 1, nunit(imolty)
                   ryu(i,j) = ryu(i,j) + dy
                end do
             else
                dz = szcm(i)*(hmat(hbox,9)-hmato(9))
                zcm(i) = zcm(i) + dz
                do j = 1, nunit(imolty)
                   rzu(i,j) = rzu(i,j) + dz
                end do
             end if
          else if (ibox .eq. jbox) then
             if (lsolid(jbox)) then
                if ( lx(ibox) ) then
                   dx = xcm(i) * df
                   xcm(i) = xcm(i) + dx
                   do j = 1, nunit(imolty)
                      rxu(i,j) = rxu(i,j) + dx
                   end do
                else if ( ly(ibox) ) then
                   dy = ycm(i) * df
                   ycm(i) = ycm(i) + dy
                   do j = 1, nunit(imolty)
                      ryu(i,j) = ryu(i,j) + dy
                   end do
                else
                   dz = zcm(i) * df
                   zcm(i) = zcm(i) + dz
                   do j = 1, nunit(imolty)
                      rzu(i,j) = rzu(i,j) + dz
                   end do
                end if
             else
                dx = xcm(i) * df
                dy = ycm(i) * df
                if ( lpbcz ) dz = zcm(i) * df
                xcm(i) = xcm(i) + dx
                ycm(i) = ycm(i) + dy
                if ( lpbcz ) zcm(i) = zcm(i) + dz
                do j = 1, nunit(imolty)
                   rxu(i,j) = rxu(i,j) + dx
                   ryu(i,j) = ryu(i,j) + dy
                   if ( lpbcz ) rzu(i,j) = rzu(i,j) + dz
                end do
             end if
          end if
       end do
    else
       ! calculate new volume
       expdv=rmvol(ipairb)*(2.0E0_dp*random(-1)-1.0E0_dp)
       expdv = volo(boxa)/volo(boxb)*exp(expdv)
       voln(boxa)= expdv*volt/(1+expdv)
       voln(boxb)= volt-voln(boxa)
       rbcut(boxa) = 2.0_dp*rcut(boxa)
       rbcut(boxb) = 2.0_dp*rcut(boxb)

       do i=1,2
          if (i.eq.1) ibox=boxa
          if (i.eq.2) ibox=boxb

          if (lsolid(ibox)) then
             ! volume move independently in x, y, z directions
             dfac(ibox)=voln(ibox)/volo(ibox)
             if (lx(ibox)) boxlx(ibox) = boxlx(ibox) * dfac(ibox)
             if (ly(ibox)) boxly(ibox) = boxly(ibox) * dfac(ibox)
             if (lz(ibox)) boxlz(ibox) = boxlz(ibox) * dfac(ibox)
          else
             if ( lpbcz ) then
                dfac(ibox)= (voln(ibox)/volo(ibox))**(1.0E0_dp/3.0E0_dp)
                boxlz(ibox) = boxlz(ibox) * dfac(ibox)
             else
                dfac(ibox)= sqrt(voln(ibox)/volo(ibox))
             end if
             boxlx(ibox) = boxlx(ibox) * dfac(ibox)
             boxly(ibox) = boxly(ibox) * dfac(ibox)
          end if
       end do

       if ( boxlx(boxa) .lt. rbcut(boxa) .or.  boxly(boxa) .lt. rbcut(boxa) .or.  (lpbcz .and. boxlz(boxa) .lt. rbcut(boxa)) .or. boxlx(boxb) .lt. rbcut(boxb) .or.  boxly(boxb) .lt. rbcut(boxb) .or.  (lpbcz .and. boxlz(boxb) .lt. rbcut(boxb)) ) then
          boxlx(boxa) = bxo(boxa)
          boxlx(boxb) = bxo(boxb)
          boxly(boxa) = byo(boxa)
          boxly(boxb) = byo(boxb)
          if ( lpbcz ) then
             boxlz(boxa) = bzo(boxa)
             boxlz(boxb) = bzo(boxb)
          end if
          call dump('final-config')
          call err_exit(__FILE__,__LINE__,'A move was attempted that would lead to a boxlength less than twice rcut',myid+1)
       end if

       ! determine new positions of the molecules
       ! calculate centre of mass and its displacement
       do i = 1, nchain
          ibox = nboxi(i)
          if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
             imolty = moltyp(i)
             df = dfac(ibox) - 1.0E0_dp
             if (lsolid(ibox)) then
                if ( lx(ibox) ) then
                   dx = xcm(i) * df
                   xcm(i) = xcm(i) + dx
                   do j = 1, nunit(imolty)
                      rxu(i,j) = rxu(i,j) + dx
                   end do
                else if ( ly(ibox) ) then
                   dy = ycm(i) * df
                   ycm(i) = ycm(i) + dy
                   do j = 1, nunit(imolty)
                      ryu(i,j) = ryu(i,j) + dy
                   end do
                else
                   dz = zcm(i) * df
                   zcm(i) = zcm(i) + dz
                   do j = 1, nunit(imolty)
                      rzu(i,j) = rzu(i,j) + dz
                   end do
                end if
             else
                dx = xcm(i) * df
                dy = ycm(i) * df
                if ( lpbcz ) dz = zcm(i) * df
                xcm(i) = xcm(i) + dx
                ycm(i) = ycm(i) + dy
                if ( lpbcz ) zcm(i) = zcm(i) + dz
                do j = 1, nunit(imolty)
                   rxu(i,j) = rxu(i,j) + dx
                   ryu(i,j) = ryu(i,j) + dy
                   if ( lpbcz ) rzu(i,j) = rzu(i,j) + dz
                end do
             end if
          end if
       end do
    end if

    if ( lchgall ) then
       if (lsolid(boxa).and.(.not.lrect(boxa))) then
          min_boxl = min(min_width(boxa,1),min_width(boxa,2),min_width(boxa,3))
       else
          min_boxl = min(boxlx(boxa),boxly(boxa),boxlz(boxa))
       end if
       calp(boxa) = kalp(boxa)/min_boxl
       if (lsolid(boxb).and.(.not.lrect(boxb))) then
          min_boxl = min(min_width(boxb,1),min_width(boxb,2), min_width(boxb,3))
       else
          min_boxl = min(boxlx(boxb),boxly(boxb),boxlz(boxb))
       end if
       calp(boxb) = kalp(boxb)/min_boxl
    end if

    do i = 1,2
       if ( i .eq. 1 ) ibox = boxa
       if ( i .eq. 2 ) ibox = boxb
       call sumup(ovrlap,v,ibox,.true.)
       if ( ovrlap ) goto 500
       vboxn(:,ibox) = v
       vboxn(10,ibox)= v3garo
       vboxn(1,ibox) = vboxo(1,ibox) + (vboxn(2,ibox)-vboxo(2,ibox)) + (vboxn(9,ibox)-vboxo(9,ibox)) + (vboxn(8,ibox)-vboxo(8,ibox)) + (vboxn(10,ibox)-vboxo(10,ibox)) ! inter, ext, elect, garo
    end do

    if ( lanes ) then
       ! for ANES algorithm, optimize the charge configuration
       ! on the new coordinates, continue to use the fluctuating charge
       ! algorithm to optimize the charge configurations, update the
       ! energy, coordinates and the ewald sum
       do i = 1,2
          if ( i .eq. 1 ) ibox = boxa
          if ( i .eq. 2 ) ibox = boxb
          vbox(1,ibox) = vbox(1,ibox) + (vboxn(1,ibox) - vboxo(1,ibox))
          vbox(2,ibox)  = vbox(2,ibox) +  (vboxn(2,ibox) - vboxo(2,ibox))
          vbox(3,ibox) = vbox(3,ibox) + (vboxn(3,ibox) - vboxo(3,ibox))
          vbox(9,ibox) = vbox(9,ibox) + (vboxn(9,ibox) - vboxo(9,ibox))
          vbox(8,ibox) = vbox(8,ibox) +  (vboxn(8,ibox) - vboxo(8,ibox))
          do ichoiq = 1,nchoiq(ibox)
             call flucq(0,ibox)
          end do
       end do
       dele = (vbox(1,boxa) - vboxo(1,boxa))+( vbox(1,boxb)- vboxo(1,boxb)) - ((nchbox(boxa)+1+ghost_particles(boxa)) *log(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+1+ghost_particles(boxb)) *log(voln(boxb)/volo(boxb))/beta)
    else if (lncubic) then
       dele = (vboxn(1,boxa)-vboxo(1,boxa)) + (vboxn(1,boxb)-vboxo(1,boxb)) - ((nchbox(boxa)+ghost_particles(boxa))*log(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+ghost_particles(boxb))*log(voln(boxb)/volo(boxb))/beta)
    else
       dele = (vboxn(1,boxa)-vboxo(1,boxa)) + (vboxn(1,boxb)-vboxo(1,boxb)) - ((nchbox(boxa)+1+ghost_particles(boxa))*log(voln(boxa)/volo(boxa))/beta) - ((nchbox(boxb)+1+ghost_particles(boxb))*log(voln(boxb)/volo(boxb))/beta)
    end if

    ! acceptance test
    if (random(-1) .lt. exp(-beta*dele) ) then
       ! accepted
       bsvol(ipairb) = bsvol(ipairb) + 1.0E0_dp
       if ( lncubic ) then
          bshmat(hbox,jhmat) = bshmat(hbox,jhmat) + 1.0E0_dp
       end if
       call update_box(boxa)
       call update_box(boxb)
       return
    end if

    ! rejected
500 call restore_box(boxa)
    call restore_box(boxb)
    call restore_configuration((/boxa,boxb/))

#ifdef __DEBUG__
    write(io_output,*) 'end VOLUME_2BOX in ',myid,boxa,boxb
#endif
    return
  end subroutine volume_2box

!***********************************************************
! makes an isotropic volume change under const. pressure  **
! the maximum change is controlled by rmtrax and the      **
! number of successful trial moves is stored in bsvol.    **
!***********************************************************
!
! perform change of the volume: random walk in V
!
  subroutine volume_1box()
    real::rbox,volo,voln,rbcut,dx,dy,dz,dfac,df,v(nEnergy),dele,min_boxl
    integer::ibox,boxvch,jhmat,i,imolty,j,ichoiq
    logical::lx,ly,lz,ovrlap
! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'start VOLUME_1BOX in ',myid
#endif

    ! Select a box at  random to change the volume of box
    rbox = random(-1)
    do ibox = 1,nbox
       if (rbox .lt. pmvlmt(ibox) ) then
          boxvch=ibox
          exit
       end if
    end do

    bnvol(boxvch) = bnvol(boxvch) + 1.0E0_dp

    call save_box(boxvch)
    call save_configuration((/boxvch/))

    lx = .false.
    ly = .false.
    lz = .false.
    if ( lsolid(boxvch) ) then
       ! volume move independently in x, y, z directions
       rbox = random(-1)
       if ( rbox .le. pmvolx ) then
          lx = .true.
       else if ( rbox .le. pmvoly ) then
          ly = .true.
       else
          lz = .true.
       end if

       if (.not.lrect(boxvch)) then
          volo = cell_vol(boxvch)
          ! select one of the cell edge
          if ( lx ) then
             rbox = 3.0_dp*random(-1)
             if ( rbox .le. 1.0E0_dp ) then
                jhmat = 1
             else if (rbox .le. 2.0E0_dp ) then
                jhmat = 4
             else
                jhmat = 7
             end if
          else if ( ly ) then
             rbox = 2.0_dp*random(-1)
             if ( rbox .le. 1.0E0_dp ) then
                jhmat = 5
             else
                jhmat = 8
             end if
          else
             jhmat = 9
          end if
       end if
    end if

    if (.not.lsolid(boxvch).or.lrect(boxvch)) then
       if ( lpbcz ) then
          volo = bxo(boxvch)*byo(boxvch)*bzo(boxvch)
       else
          volo = bxo(boxvch)*byo(boxvch)
       end if
    end if

    if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then
       bnhmat(boxvch,jhmat) = bnhmat(boxvch,jhmat) + 1.0E0_dp
       hmat(boxvch,jhmat) = hmat(boxvch,jhmat) + rmhmat(boxvch,jhmat)* ( 2.0E0_dp*random(-1) - 1.0E0_dp )
       call matops(boxvch)

       voln = cell_vol(boxvch)
       rbcut = 2.0_dp*rcut(boxvch)

       if (ANY(rbcut .gt. min_width(boxvch,:))) then
          hmat(boxvch,jhmat) = hmato(jhmat)
          call dump('final-config')
          write(io_output,*) 'w1:',min_width(boxvch,1),'w2:',min_width(boxvch,2),'w3:',min_width(boxvch,3)
          call err_exit(__FILE__,__LINE__,'non-rectangular volume move rejected. box width below cutoff size',myid+1)
       end if

       ! determine the displacement of the COM
       do i = 1,nchain
          if (nboxi(i) .eq. boxvch) then
             imolty = moltyp(i)
             if ( lx ) then
                dx = sxcm(i)*(hmat(boxvch,1)-hmato(1))+sycm(i)*(hmat(boxvch,4)-hmato(4))+szcm(i)*(hmat(boxvch,7)-hmato(7))
                xcm(i) = xcm(i) + dx
                do j = 1, nunit(imolty)
                   rxu(i,j) = rxu(i,j) + dx
                end do
             else if ( ly ) then
                dy = sycm(i)*(hmat(boxvch,5)-hmato(5))+szcm(i)*(hmat(boxvch,8)-hmato(8))
                ycm(i) = ycm(i) + dy
                do j = 1, nunit(imolty)
                   ryu(i,j) = ryu(i,j) + dy
                end do
             else
                dz = szcm(i)*(hmat(boxvch,9)-hmato(9))
                zcm(i) = zcm(i) + dz
                do j = 1, nunit(imolty)
                   rzu(i,j) = rzu(i,j) + dz
                end do
             end if
          end if
       end do
    else
       ! calculate new volume
       voln = volo + rmvol(boxvch) * ( 2.0E0_dp*random(-1) - 1.0E0_dp )
       rbcut = 2.0_dp*rcut(boxvch)

       if (lsolid(boxvch)) then
          ! volume move independently in x, y, z directions
          dfac=voln/volo
          if (lx) boxlx(boxvch) = boxlx(boxvch) * dfac
          if (ly) boxly(boxvch) = boxly(boxvch) * dfac
          if (lz) boxlz(boxvch) = boxlz(boxvch) * dfac
       else
          if ( lpbcz ) then
             dfac = (voln/volo)**(1.0E0_dp/3.0E0_dp)
             boxlz(boxvch) = boxlz(boxvch) * dfac
          else
             dfac= sqrt(voln/volo)
          end if
          boxlx(boxvch) = boxlx(boxvch) * dfac
          boxly(boxvch) = boxly(boxvch) * dfac
       end if

       if (boxlx(boxvch) .lt. rbcut .or. boxly(boxvch) .lt. rbcut .or. (lpbcz .and. boxlz(boxvch) .lt. rbcut) ) then
          boxlx(boxvch) = bxo(boxvch)
          boxly(boxvch) = byo(boxvch)
          if ( lpbcz ) then
             boxlz(boxvch) = bzo(boxvch)
          end if
          call dump('final-config')
          write(io_output,*) 'boxvch',boxvch
          call err_exit(__FILE__,__LINE__,'A move was attempted that would lead to a boxlength less than twice rcut',myid+1)
       end if

       ! determine new positions of the molecules
       ! calculate centre of mass and its displacement
       df = dfac - 1.0E0_dp
       do i = 1, nchain
          ! Check if the chain i is in the correct box
          if (nboxi(i) .eq. boxvch) then
             imolty = moltyp(i)
             if (lsolid(boxvch)) then
                if ( lx ) then
                   dx = xcm(i) * df
                   xcm(i) = xcm(i) + dx
                   do j = 1, nunit(imolty)
                      rxu(i,j) = rxu(i,j) + dx
                   end do
                else if ( ly ) then
                   dy = ycm(i) * df
                   ycm(i) = ycm(i) + dy
                   do j = 1, nunit(imolty)
                      ryu(i,j) = ryu(i,j) + dy
                   end do
                else
                   dz = zcm(i) * df
                   zcm(i) = zcm(i) + dz
                   do j = 1, nunit(imolty)
                      rzu(i,j) = rzu(i,j) + dz
                   end do
                end if
             else
                dx = xcm(i) * df
                dy = ycm(i) * df
                if ( lpbcz ) dz = zcm(i) * df
                xcm(i) = xcm(i) + dx
                ycm(i) = ycm(i) + dy
                if ( lpbcz ) zcm(i) = zcm(i) + dz
                do j = 1, nunit(imolty)
                   rxu(i,j) = rxu(i,j) + dx
                   ryu(i,j) = ryu(i,j) + dy
                   if ( lpbcz ) rzu(i,j) = rzu(i,j) + dz
                end do
             end if
          end if
       end do
    end if

    if ( lchgall ) then
       if (lsolid(boxvch).and.(.not.lrect(boxvch))) then
          min_boxl = min(min_width(boxvch,1),min_width(boxvch,2),min_width(boxvch,3))
       else
          min_boxl = min(boxlx(boxvch),boxly(boxvch),boxlz(boxvch))
       end if
       calp(boxvch) = kalp(boxvch)/boxlx(boxvch)
    end if

    call sumup(ovrlap,v,boxvch,.true.)
    if ( ovrlap ) goto 500
    vboxn(:,boxvch) = v
    vboxn(10,boxvch)= v3garo
    vboxn(1,boxvch) = vboxo(1,boxvch) + (vboxn(2,boxvch)-vboxo(2,boxvch)) + (vboxn(9,boxvch)-vboxo(9,boxvch)) + (vboxn(8,boxvch)-vboxo(8,boxvch)) + (vboxn(10,boxvch)-vboxo(10,boxvch)) !inter, ext, elect, garo

    if ( lanes ) then
       ! for ANES algorithm, optimize the charge configuration
       ! on the new coordinates, continue to use the fluctuating charge
       ! algorithm to optimize the charge configurations, update the
       ! energy, coordinates and the ewald sum
       vbox(1,boxvch)=vbox(1,boxvch)+(vboxn(1,boxvch)-vboxo(1,boxvch))
       vbox(2,boxvch) = vbox(2,boxvch) + (vboxn(2,boxvch)-vboxo(2,boxvch))
       vbox(3,boxvch)  = vbox(3,boxvch) + (vboxn(3,boxvch)-vboxo(3,boxvch))
       vbox(9,boxvch)   = vbox(9,boxvch) + (vboxn(9,boxvch)-vboxo(9,boxvch))
       vbox(8,boxvch) = vbox(8,boxvch) + (vboxn(8,boxvch)-vboxo(8,boxvch))
       do ichoiq = 1,nchoiq(boxvch)
          call flucq(0,boxvch)
       end do
       dele = (vbox(1,boxvch) - vboxo(1,boxvch)) + express(boxvch)*(voln-volo) - ((nchbox(boxvch)+ghost_particles(boxvch)) * log(voln/volo) / beta )
    else
       dele = ( vboxn(1,boxvch) - vboxo(1,boxvch) ) + express(boxvch)*(voln-volo) - ((nchbox(boxvch)+ghost_particles(boxvch))*log(voln/volo)/beta)
    end if

    ! acceptance test
    if (random(-1) .lt. exp(-(beta*dele)) ) then
       ! accepted
       bsvol(boxvch) = bsvol(boxvch) + 1.0E0_dp
       if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then
          bshmat(boxvch,jhmat) = bshmat(boxvch,jhmat) + 1.0E0_dp
       end if
       call update_box(boxvch)
       return
    end if

    ! rejected
500 call restore_box(boxvch)
    call restore_configuration((/boxvch/))

#ifdef __DEBUG__
    write(io_output,*) 'end VOLUME_1BOX in ',myid,boxvch
#endif
    return
  end subroutine volume_1box

  subroutine init_moves_volume
    integer::jerr
    allocate(acsvol(nbxmax),acnvol(nbxmax),acshmat(nbxmax,9),acnhmat(nbxmax,9),bsvol(nbxmax),bnvol(nbxmax),bshmat(nbxmax,9),bnhmat(nbxmax,9),vboxn(nEnergy,nbxmax),vboxo(nEnergy,nbxmax),bxo(nbxmax),byo(nbxmax),bzo(nbxmax),xcmo(nmax),ycmo(nmax),zcmo(nmax),rxuo(nmax,numax),ryuo(nmax,numax),rzuo(nmax,numax),qquo(nmax,numax),neigho_cnt(nmax),neigho(100,nmax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'init_moves_volume: allocation failed',jerr)
    end if

    acsvol = 0.E0_dp
    acnvol = 0.E0_dp
    acshmat = 0.0E0_dp
    acnhmat = 0.0E0_dp
    bsvol = 0.0E0_dp
    bnvol = 0.0E0_dp
    bshmat = 0.0E0_dp
    bnhmat = 0.0E0_dp
  end subroutine init_moves_volume

!> \brief Adjust maximum volume displacement
  subroutine update_volume_max_displacement(io_output)
    integer,intent(in)::io_output
    integer::ibox,j
    real::ratvol

    do ibox = 1, nbox
       if (lsolid(ibox) .and. .not. lrect(ibox)) then
          do j = 1,9
             if ( bnhmat(ibox,j) .gt. 0.5E0_dp ) then
                ratvol = bshmat(ibox,j) / bnhmat(ibox,j)
                if (ratvol .eq. 0.0E0_dp) then
                   rmhmat(ibox,j) = rmhmat(ibox,j) * 0.1E0_dp
                else
                   rmhmat(ibox,j) = rmhmat(ibox,j)*ratvol/tavol
                end if
             end if
          end do
       else
          if ( bnvol(ibox) .gt. 0.5E0_dp ) then
             ratvol = bsvol(ibox) / bnvol(ibox)
             if ( ratvol .eq. 0.0E0_dp ) then
                rmvol(ibox) = rmvol(ibox) * 0.1E0_dp
             else
                rmvol(ibox) = rmvol(ibox) * ratvol / tavol
                if (rmvol(ibox).gt.(0.10E0_dp*boxlx(ibox)*boxly(ibox)*boxlz(ibox))) then
                   rmvol(ibox)=0.1E0_dp*(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
                end if
             end if
          end if
       end if
    end do

    if (myid.eq.rootid) then
       do ibox = 1, nbox
          if (lsolid(ibox) .and. .not. lrect(ibox)) then
             do j = 1,9
                write(io_output,"(' h-matrix change:  bn =',f8.1, '   bs =',f8.1,'   max.displ. =',e12.5)") bnhmat(ibox,j),bshmat(ibox,j), rmhmat(ibox,j)
             end do
          else
             write(io_output,"(' volume change:  bn =',f8.1, '   bs =',f8.1,'   max.displ. =',e12.5)") bnvol(ibox),bsvol(ibox),rmvol(ibox)
          end if
       end do
    end if

    do ibox = 1, nbox
       if (lsolid(ibox) .and. .not. lrect(ibox)) then
          do j = 1,9
             acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
             acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
             bshmat(ibox,j) = 0.0E0_dp
             bnhmat(ibox,j) = 0.0E0_dp
          end do
       else
          acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
          acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
          bnvol(ibox) = 0.0E0_dp
          bsvol(ibox) = 0.0E0_dp
       end if
    end do
  end subroutine update_volume_max_displacement

  subroutine output_volume_stats(io_output)
    integer,intent(in)::io_output
    integer::ibox,j
    real::ratvol
    character(LEN=default_path_length)::fmt="(' h-matrix attempts =',f8.1,'   ratio =',f6.3, '   max.displ. =',e11.4)"

    write(io_output,*)
    write(io_output,*) '### Volume change       ###'
    do ibox = 1,nbox
       if (lsolid(ibox) .and. .not. lrect(ibox)) then
          do j = 1,9
             acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
             acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
             if ( acshmat(ibox,j) .gt. 0.5E0_dp) then
                write(io_output,fmt) acnhmat(ibox,j), acshmat(ibox,j)/acnhmat(ibox,j),rmhmat(ibox,j)
             else
                write(io_output,fmt) acnhmat(ibox,j),0.0E0_dp,rmhmat(ibox,j)
             end if
          end do
       else
          acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
          acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
          if ( acnvol(ibox) .ne. 0.0E0_dp ) then
             ratvol = acsvol(ibox) / acnvol(ibox)
          else
             ratvol = 0.0E0_dp
          end if
          write(io_output,"(' attempts =',f8.1,'   ratio =',f6.3, '   max.displ. =',e11.4)") acnvol(ibox),ratvol,rmvol(ibox)
       end if
    end do
  end subroutine output_volume_stats

! store old box lengths and energy
  subroutine save_box(box)
    integer,intent(in)::box
    integer::j
    real::vdum

    vboxo(1,box)   = vbox(1,box)
    vboxo(2,box) = vbox(2,box)
    vboxo(3,box)  = vbox(3,box)
    vboxo(9,box)   = vbox(9,box)
    vboxo(8,box) = vbox(8,box)
    vboxo(11,box) = vbox(11,box)
    vboxo(10,box)     = vbox(10,box)

    bxo(box) = boxlx(box)
    byo(box) = boxly(box)
    if ( lpbcz ) bzo(box) = boxlz(box)

    if (lsolid(box) .and. .not. lrect(box)) then
       do j = 1,9
          hmato(j) = hmat(box,j)
          hmatio(j) = hmati(box,j)
       end do
    end if

    if ( lewald ) then
       call save_kvector(box)
       call recip(box,vdum,vdum,3)
    end if
  end subroutine save_box

! store old chain configuration
  subroutine save_configuration(boxes)
    integer,intent(in)::boxes(:)
    integer::i,j

    do i = 1, nchain
       if (ANY(nboxi(i).eq.boxes)) then
          xcmo(i) = xcm(i)
          ycmo(i) = ycm(i)
          if (lpbcz) zcmo(i) = zcm(i)
          do j = 1, nunit(moltyp(i))
             rxuo(i,j) = rxu(i,j)
             ryuo(i,j) = ryu(i,j)
             if ( lpbcz ) rzuo(i,j) = rzu(i,j)
             qquo(i,j) = qqu(i,j)
          end do
          if (lneighbor) then
             neigho_cnt(i) = neigh_cnt(i)
             do j = 1,neigho_cnt(i)
                neigho(j,i)=neighbor(j,i)
             end do
          end if

          ! store neighbor list for garofalini --- KEA
          if (lgaro) then
             neigh_o(i) = neigh_cnt(i)
             do j=1,neigh_o(i)
                neighboro(j,i) = neighbor(j,i)
                ndijo(j,i) = ndij(j,i)
                nxijo(j,i) = nxij(j,i)
                nyijo(j,i) = nyij(j,i)
                nzijo(j,i) = nzij(j,i)
             end do
          end if
       end if
    end do
  end subroutine save_configuration

! restore old energy, box lengths
  subroutine restore_box(box)
    integer,intent(in)::box
    integer::j
    real::vdum

    vbox(1,box)    = vboxo(1,box)
    vbox(2,box) = vboxo(2,box)
    vbox(3,box)  = vboxo(3,box)
    vbox(9,box)   = vboxo(9,box)
    vbox(8,box) = vboxo(8,box)
    vbox(11,box) = vboxo(11,box)
    vbox(10,box) = vboxo(10,box)

    boxlx(box)   = bxo(box)
    boxly(box)   = byo(box)
    if ( lpbcz ) boxlz(box)   = bzo(box)

    if (lsolid(box) .and. .not. lrect(box)) then
       do j = 1,9
          hmat(box,j) = hmato(j)
          hmati(box,j) = hmatio(j)
       end do
       call matops(box)
    end if

    if ( lewald ) then
       call restore_kvector(box)
       call recip(box,vdum,vdum,4)
    end if
  end subroutine restore_box

! restore old energy, box lengths
  subroutine restore_configuration(boxes)
    integer,intent(in)::boxes(:)
    integer::i,j

    do i = 1, nchain
       if (ANY(nboxi(i).eq.boxes)) then
          xcm(i) = xcmo(i)
          ycm(i) = ycmo(i)
          if ( lpbcz ) zcm(i) = zcmo(i)
          do j = 1, nunit(moltyp(i))
             rxu(i,j) = rxuo(i,j)
             ryu(i,j) = ryuo(i,j)
             if ( lpbcz ) rzu(i,j) = rzuo(i,j)
             qqu(i,j) = qquo(i,j)
          end do
          if (lneighbor) then
             neigh_cnt(i) = neigho_cnt(i)
             do j = 1,neigh_cnt(i)
                neighbor(j,i)=neigho(j,i)
             end do
          end if

          ! restore old neighbor list for garofalini --- KEA
          if (lgaro) then
             neigh_cnt(i) = neigh_o(i)
             do j=1,neigh_cnt(i)
                neighbor(j,i) = neighboro(j,i)
                ndij(j,i) = ndijo(j,i)
                nxij(j,i) = nxijo(j,i)
                nyij(j,i) = nyijo(j,i)
                nzij(j,i) = nzijo(j,i)
             end do
          end if
       end if
    end do
  end subroutine restore_configuration

  subroutine update_box(box)
    integer,intent(in)::box

    if ( .not. lanes ) then
       vbox(1,box)    = vbox(1,box) + (vboxn(1,box) - vboxo(1,box))
       vbox(2,box) = vbox(2,box) + (vboxn(2,box) - vboxo(2,box))
       vbox(3,box)  = vbox(3,box) + (vboxn(3,box) - vboxo(3,box))
       vbox(9,box)   = vbox(9,box) + (vboxn(9,box) - vboxo(9,box))
       vbox(8,box) = vbox(8,box) + (vboxn(8,box) - vboxo(8,box))
       vbox(10,box) = vbox(10,box) + (vboxn(10,box)-vboxo(10,box))
    end if

    ! update centers of mass
    call ctrmas(.true.,box,0,5)
    ! update linkcell, if applicable
    if (licell .and. (box .eq. boxlink)) then
       call build_linked_cell()
    end if
    if (lneigh) then
       call rebuild_neighbor_list(box)
    end if
  end subroutine update_box

  subroutine read_checkpoint_volume(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) bnvol,bsvol,acnvol,acsvol,bnhmat,bshmat,acnhmat,acshmat
    call mp_bcast(bnvol,nbxmax,rootid,groupid)
    call mp_bcast(bsvol,nbxmax,rootid,groupid)
    call mp_bcast(acnvol,nbxmax,rootid,groupid)
    call mp_bcast(acsvol,nbxmax,rootid,groupid)
    call mp_bcast(bnhmat,nbxmax*9,rootid,groupid)
    call mp_bcast(bshmat,nbxmax*9,rootid,groupid)
    call mp_bcast(acnhmat,nbxmax*9,rootid,groupid)
    call mp_bcast(acshmat,nbxmax*9,rootid,groupid)
  end subroutine read_checkpoint_volume

  subroutine write_checkpoint_volume(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) bnvol,bsvol,acnvol,acsvol,bnhmat,bshmat,acnhmat,acshmat
  end subroutine write_checkpoint_volume
end MODULE moves_volume
