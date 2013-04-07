module transfer_swap
  use util_random,only:random
  use util_runtime,only:err_exit
  use sim_system
  use sim_cell
  use energy_kspace,only:recip
  use energy_pairwise,only:energy,boltz,coru
  use energy_intramolecular,only:U_bonded
  use moves_cbmc,only:rosenbluth,schedule,explct,safeschedule
  use transfer_shared,only:lopt_bias,update_bias
  implicit none
  private
  save
  public::swap,init_swap,cnt,output_swap_stats,read_checkpoint_swap,write_checkpoint_swap

  real,allocatable,public::acchem(:,:)
  integer,allocatable,public::bnchem(:,:)
  integer,allocatable::bsswap(:,:,:),bnswap(:,:,:),bnswap_in(:,:),bnswap_out(:,:)
  integer::cnt_wf1(0:6,0:6,4),cnt_wf2(0:6,0:6,4),cnt_wra1(1000,4),cnt_wra2(1000,4)

contains
!******************************************************************
! removes a molecule from one box and inserts it into the other  **
! using CBMC insertion techniques.  Works for linear or branched **
! molecules and for DC-CBMC and Explicit atom                    **
! Rewritten from old swap and swapbr subroutines by M.G. Martin  **
! 9-18-97                                                        **
!******************************************************************
  subroutine swap()
    use sim_particle,only:update_neighbor_list

    logical::ovrlap,lterm,lnew,lempty,ldone,ltors,lovrh(nchmax),lfavor,laccept,lswapinter,lrem_out,lins_in,lneighij,linsk_in,lremk_in,lfixnow

    integer::boxins,boxrem,imol,ichoi,ip,iwalk,idum,iins1,imolty1
    integer::istt,iett,itype,ipair,ipairb,beg

    integer::iutry,icbu,ifrom,irem,iins,glist(numax),findex,iii,j,ibox,iunit,ic,pointp,imolty,imt,jmt,igrow,pointp2,jins,jmolty,neighj_num,neighk_num,joffset,koffset,kmolty,kins,target,neigh_old

    real::sx,sy,sz

    real::v(nEnergy),delen,deleo,rpair
    real::vnewflucq,voldflucq,qion(numax),ctorfo,ctorfn
    real::rxuold(numax),ryuold(numax),rzuold(numax)
    real::rmol,rbf,bsum
    real::waddnew,waddold

    real::total_NBE,vtgn,vbendn,vvibn

    real::v1insext,v1remext,v1ins,w1ins,v1rem,w1rem,v1insint,v1remint,v1insewd,v1remewd,wnlog,wolog,wdlog,wratio,vinsta,vremta,volins,volrem,rho,arg,v1inselc,v1remelc
    real::rvol,x,y,z,rijsq,wbias_ins,wbias_rem,r,xi1,xi2,xisq
    real::vrecipn,vrecipo,vdum,whins,whrem
    real::rxuh(numax,nchmax),ryuh(numax,nchmax),rzuh(numax,nchmax),delenh(nchmax),vtrhext(nchmax),vtrhintra(nchmax),vtrhinter(nchmax),vtrhelect(nchmax),vtrhewald(nchmax),vtrhtg(nchmax),bfach(nchmax)

! --------------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'START SWAP in ',myid
#endif
! write(11,*) '1:',neigh_cnt(18)

    lempty = .false.
    lfixnow = .false.
    lins_in = .false.
    linsk_in = .false.

    if (lexpee.and.leemove) then
       imolty = ee_moltyp(nstate)
       imolty1 = imolty
       irem = eeirem
       pointp = eepointp
       boxrem = boxrem1
       boxins = boxins1
       if (boxins.eq.boxrem) then
          lswapinter = .false.
       else
          lswapinter = .true.
       end if
! write(io_output,*) 'ee val', imolty, irem, pointp, boxrem, boxins
    else
       wee_ratio = 1.0E0_dp

! select a molecule typ with probabilities given in pmswmt
       rmol = random(-1)
       do imol = 1, nmolty
          if ( rmol .lt. pmswmt(imol) ) then
             imolty = imol
             exit
          end if
       end do

       if (temtyp(imolty).eq.0) return
       imolty1 = imolty

! select a box given in pmswatyp
       if ( nswapb(imolty) .gt. 1 ) then
          rpair = random(-1)
          do ipair = 1, nswapb(imolty)
             if ( rpair .lt. pmswapb(imolty,ipair) ) then
                ipairb = ipair
                exit
             end if
          end do
       else
          ipairb = 1
       end if

       if (random(-1).lt.0.5E0_dp) then
          boxins=box1(imolty,ipairb)
          boxrem=box2(imolty,ipairb)
       else
          boxins=box2(imolty,ipairb)
          boxrem=box1(imolty,ipairb)
       end if

       if ( boxins .eq. boxrem ) then
          lswapinter = .false.
       else
          lswapinter = .true.
       end if
! write(io_output,*) 'boxins:',boxins,'boxrem:',boxrem
       if ( .not. (lgibbs .or. lgrand) .and. lswapinter ) then
          call err_exit(__FILE__,__LINE__,'no interbox swap if not gibbs/grand ensemble!',myid+1)
       end if

! select a chain in BOXREM at random ***

       if ( ncmt(boxrem,imolty) .eq. 0 ) then
          lempty = .true.
          if ( .not. lswapinter ) return
       else if ( lswapinter .or. lavbmc1(imolty) .and. .not. (lavbmc2(imolty) .or. lavbmc3(imolty)) ) then
! for the advanced AVBMC algorithm, this particle will be selected in
! sub-regions defined by Vin
197       pointp = int( dble(ncmt(boxrem,imolty))*random(-1) ) + 1
          if (lexpee) then
             if ((pointp.eq.eepointp).and. (boxrem.eq.box_state(mstate)).and. (ncmt(boxrem,imolty).gt.1)) then
                goto 197
             else
                return
             end if
          end if

          irem = parbox(pointp,boxrem,imolty)

          if ( moltyp(irem) .ne. imolty ) write(io_output,*) 'screwup swap, irem:',irem,moltyp(irem),imolty
          ibox = nboxi(irem)
          if ( ibox .ne. boxrem ) then
             call err_exit(__FILE__,__LINE__,'problem in swap',myid+1)
          end if
       end if

!$$$      write(io_output,*) 'particle ',irem,' is being removed, imolty is:',
!$$$     &     imolty,' and the box is:',boxrem


! ===>  for both gibbs and grand-canonical we have:
! insert a chain in box: boxins
! remove one in box: boxrem
! write(io_output,*) 'boxrem',boxrem,' imolty',imolty,' lempty',lempty
! bnswap(imolty,ipair,X) decoder X = 1 is # non-empty attempts into box 1
! X = 2 is # non-empty attempts into box 2, ..., X = nbox
! X = nbox + 1 is # success swaps into box 1, ... X = nbox + nbox
! bsswap is the same thing but keeps track of all attempts (X = 1 to nbox)
! and successful growths (X = nbox + 1 to nbox + nbox)

       bsswap(imolty,ipairb,boxins)=bsswap(imolty,ipairb,boxins)+1
       if (.not. lempty) then
          bnswap(imolty,ipairb,boxins) = bnswap(imolty,ipairb,boxins) + 1
       else if (lrigid(imolty)) then
          ! old configuration for irem must exist to grow the new configuration
          return
       end if
    end if

! store number of units in iunit ***
    iunit = nunit(imolty)
    igrow = nugrow(imolty)
! give i a phony number ***
    if ( lswapinter ) then
       iins = nchain + 1
       iins1 = iins
       moltyp(iins) = imolty
! give charges to phony number
       if ((.not.lexpee).and.(.not.leemove)) then
          if ( lempty ) then
             do icbu = 1, iunit
                qqu(iins,icbu) = qelect(ntype(imolty,icbu))
             end do
          else
             do icbu = 1, iunit
                qqu(iins,icbu) = qqu(irem,icbu)
             end do
          end if
       else
          do icbu = 1, iunit
             qqu(iins,icbu) = ee_qqu(icbu,nstate)
          end do
       end if
    else if ( lavbmc2(imolty) .or. lavbmc3(imolty) ) then
       iins = 0
    else
       iins = irem !??? irem not defined
       if (lexpee.and.leemove) then
          iins1 = iins
          do icbu = 1, iunit
             qqu(iins,icbu) = ee_qqu(icbu,nstate)
          end do
       end if
    end if

! select a position of the first/starting unit at RANDOM ***
! and calculate the boltzmann weight                     ***
! for the chain to be INSERTED                           ***

    if (lrigid(imolty)) then
       beg = riutry(imolty,1)
    else
       beg = 1
    end if

    wbias_ins = 1.0E0_dp
    wbias_rem = 1.0E0_dp

    ichoi = nchoi1(imolty)

    if ( .not. lswapinter .and. lbias(imolty)) then

       if (boxins .ne. boxrem) then
          write(io_output,*) 'avbmc, boxins, boxrem:',boxins,boxrem
          call err_exit(__FILE__,__LINE__,'',myid+1)
       end if

! *********************
! First AVBMC Steps *
! *********************
       rmol = random(-1)
       ldone = .false.
       do imol = 1, nmolty
          if ( rmol .lt. pmbsmt(imol) ) then
             if ( .not. ldone ) then
                jmolty = imol
                ldone = .true.
             end if
          end if
       end do

       if ( ncmt(boxins,jmolty) .eq. 0 .or. (jmolty .eq. imolty .and. ncmt(boxins,jmolty) .eq. 1) ) then

! ??? what shall we do now?
          lempty = .true.
          return

       else
111       pointp2=int( dble(ncmt(boxins,jmolty))*random(-1))+1
          jins = parbox(pointp2,boxins,jmolty)
          if ( jins .eq. iins ) goto 111
          if ( moltyp(jins) .ne. jmolty )  write(io_output,*) 'screwup swap, jins:' ,jins,moltyp(jins),jmolty
          if ( nboxi(jins) .ne. boxins ) then
             write(io_output,*) 'problem in swap with jins'
          end if
       end if

       if ( lavbmc3(imolty) ) then
! define a second bonding region bounded by kins with kmolty type
          rmol = random(-1)
          ldone = .false.
          do imol = 1, nmolty
             if ( rmol .lt. pmbsmt(imol) ) then
                if ( .not. ldone ) then
                   kmolty = imol
                   ldone = .true.
                end if
             end if
          end do
          if ( ncmt(boxins,kmolty) .eq. 0 .or. (kmolty .eq. imolty .and. jmolty .eq. imolty .and.  ncmt(boxins,kmolty) .eq. 2) .or. (kmolty .eq. imolty .and. ncmt(boxins,kmolty) .eq. 1) ) then

! ??? what shall we do now?
             lempty = .true.
             return
          else
112          pointp2=int( dble(ncmt(boxins,kmolty))*random(-1))+1
             kins = parbox(pointp2,boxins,kmolty)

! make sure the two regions bounded by jins and kins do not
! overlap (sampling from cluster to cluster)
! *** ??? Problems: potentially infinite loop if there is only
! one region in the whole system

             if ( jins .eq. kins ) goto 112
! do ip = 1,neighj_num
! do ic = 1,neighk_num
! if ( neighbor(ip,jins) .eq. kins .or.
!     &                    neighbor(ic,kins) .eq. jins .or.
!     &                    neighbor(ip,jins) .eq. neighbor(ic,kins))
!     &                    goto 112
! end do
! end do
             x = rxu(jins,1) - rxu(kins,1)
             y = ryu(jins,1) - ryu(kins,1)
             z = rzu(jins,1) - rzu(kins,1)
             if ( lpbc ) call mimage(x,y,z,boxins)
             rijsq = x*x + y*y + z*z
             if ( rijsq .lt. (2.0E0_dp*rbsmax)**2 ) goto 112

             if ( moltyp(kins) .ne. kmolty )  write(io_output,*) 'screwup swap, kins:' ,kins,moltyp(kins),kmolty
             if ( nboxi(kins) .ne. boxins ) then
                call err_exit(__FILE__,__LINE__,'problem in swap with kins',myid+1)
             end if
          end if
          neighk_num = neigh_cnt(kins)
       end if

       if ( lavbmc2(imolty) .or. lavbmc3(imolty) ) then
! neighj_num = neigh_cnt(jins,imolty)
          neighj_num = neigh_cnt(jins)
          if (jmolty .eq. imolty) then
             joffset = 1
          else
             joffset = 0
          end if
          if (lavbmc3(imolty) .and. kmolty .eq. imolty) then
             koffset = 1
          else
             koffset = 0
          end if
       end if

       if ( lpbc ) call setpbc(boxins)

       if ( random(-1) .lt. pmbias(imolty) ) then
          wbias_rem = pmbias(imolty) * (1.0E0_dp/vol_eff)
          lins_in = .true.
          if (lavbmc2(imolty) .or. (lavbmc3(imolty)  .and. random(-1) .gt. pmbias2(imolty)) ) then
! select a particle in the out region and move this particle into
! the in region defined by the particle jins's bonding region
! out -> j case
             wbias_rem=wbias_rem /dble(ncmt(boxrem,imolty)-neighj_num-joffset)

             lrem_out = .true.
119          pointp=int(dble(ncmt(boxrem,imolty))*random(-1))+1
             irem = parbox(pointp,boxrem,imolty)
             if ( irem .eq. jins ) goto 119
             x = rxu(irem,1) - rxu(jins,1)
             y = ryu(irem,1) - ryu(jins,1)
             z = rzu(irem,1) - rzu(jins,1)
             if ( lpbc ) call mimage(x,y,z,boxins)
             rijsq = x*x + y*y + z*z
             if ( rijsq .lt. rbsmax**2 .and. rijsq .gt. rbsmin**2)  goto 119

             if ( moltyp(irem) .ne. imolty )  write(io_output,*) 'screwup swap1, irem:',irem, moltyp(irem),imolty
             ibox = nboxi(irem)
             if ( ibox .ne. boxrem ) then
                write(io_output,*) 'problem in swap'
                call err_exit(__FILE__,__LINE__,'',myid+1)
             end if
             iins = irem
             lremk_in = .false.
          else if ( lavbmc3(imolty) ) then
! select a particle in the region bounded by kins and move this particle
! into the in region defined by the particle jins's bonding region
! k -> j case
             lrem_out = .false.
             lremk_in = .true.

             if ( neighk_num .eq. 0 .or. neighk_num .eq. 1 .and.  neighbor(1,kins) .eq. jins) then
! ??? what shall we do now?
                lempty = .true.
                return
             else
113             pointp=int(dble(neighk_num)*random(-1))+1
! irem = neighbor(pointp,kins,imolty)
! write(io_output,*) 'kins,irem:',kins,irem,neighk_num
                irem = neighbor(pointp,kins)
                if ( irem .eq. jins ) goto 113
                iins = irem
             end if
             wbias_rem=wbias_rem/dble(neighk_num)
          end if

! write(io_output,*) 'move in'
! write(io_output,*) '3:',wbias_rem,boxlx(boxins),vol_eff

          do icbu = 1,ichoi
! choose a random association distance
             rvol = random(-1)
             r = (rbsmax*rbsmax*rbsmax*rvol +  (1.0E0_dp-rvol)*rbsmin*rbsmin*rbsmin)**(1.0E0_dp/3.0E0_dp)

! calculate random vector on the unit sphere ---
109          xi1 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
             xi2 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
             xisq = xi1**2 + xi2**2
             if ( xisq .lt. 1.0E0_dp ) then
                x = r * 2.0E0_dp * xi1 * sqrt( 1.0E0_dp - xisq )
                y = r * 2.0E0_dp * xi2 * sqrt( 1.0E0_dp - xisq )
                z = r * ( 1.0E0_dp - 2.0E0_dp * xisq )
             else
                goto 109
             end if
             rxp(1,icbu) = rxu(jins,1) + x
             ryp(1,icbu) = ryu(jins,1) + y
             rzp(1,icbu) = rzu(jins,1) + z
          end do
       else
          lins_in = .false.
          if ( lavbmc2(imolty) .or. lavbmc3(imolty) ) then
! select a particle in the in region defined by the particle
! jins's bonding region through neighbor list
! and move to the out region or the region bounded by kins in AVBMC3
             if ( neighj_num .eq. 0 .or. (neighj_num .eq. 1 .and. neighbor(1,jins) .eq. kins) ) then
! ??? what shall we do now?
                lempty = .true.
                return
             else
114             pointp=int(dble(neighj_num)*random(-1))+1
! irem = neighbor(pointp,jins,imolty)
! write(io_output,*) 'jins,irem:',jins,irem,neighj_num
                irem = neighbor(pointp,jins)
                if ( irem .eq. kins ) goto 114
                if ( moltyp(irem) .ne. imolty )  write(io_output,*) 'screwup swap2, irem:',irem, moltyp(irem),imolty,neighj_num,pointp,jins
                ibox = nboxi(irem)
                if ( ibox .ne. boxrem ) then
                   write(io_output,*) 'problem in swap'
                   call err_exit(__FILE__,__LINE__,'',myid+1)
                end if
                iins = irem
             end if
             lrem_out = .false.
             lremk_in = .false.
          end if
! write(io_output,*) 'move out'
! write(io_output,*) '4:',wbias_rem,boxlx(boxins),vol_eff

          if ( lavbmc3(imolty) .and.  random(-1) .lt. pmbias2(imolty) ) then
! move to the region bounded by kins in AVBMC3
! j -> k case
             linsk_in = .true.
             wbias_rem = (1-pmbias(imolty))/vol_eff/dble(neighj_num)

             do icbu = 1,ichoi
! choose a random association distance
                rvol = random(-1)
                r = (rbsmax*rbsmax*rbsmax*rvol +  (1-rvol)*rbsmin*rbsmin*rbsmin)**(1.0/3.0E0_dp)

! calculate random vector on the unit sphere ---
139             xi1 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
                xi2 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
                xisq = xi1**2 + xi2**2
                if ( xisq .lt. 1.0E0_dp ) then
                   x = r * 2.0E0_dp * xi1 * sqrt( 1.0E0_dp - xisq )
                   y = r * 2.0E0_dp * xi2 * sqrt( 1.0E0_dp - xisq )
                   z = r * ( 1.0E0_dp - 2.0E0_dp * xisq )
                else
                   goto 139
                end if
                rxp(1,icbu) = rxu(kins,1) + x
                ryp(1,icbu) = ryu(kins,1) + y
                rzp(1,icbu) = rzu(kins,1) + z
             end do
          else
             linsk_in = .false.
! move to the out region
! j -> out case
             wbias_rem = (1-pmbias(imolty))/ (boxlx(boxins)*boxly(boxins)* boxlz(boxins)-vol_eff)
             if ( lavbmc2(imolty) .or. lavbmc3(imolty) )  wbias_rem = wbias_rem/dble(neighj_num)

             do icbu = 1,ichoi
222             rxp(1,icbu) = boxlx(boxins) * random(-1)
                ryp(1,icbu) = boxly(boxins) * random(-1)
                if (lpbcz .or. lslit) then
                   rzp(1,icbu) = boxlz(boxins) * random(-1)
                else if ( lsami .or. lmuir .or. ljoe ) then
                   if ( lempty ) then
                      rzp(1,icbu) = 20*random(-1)-10
                   else
                      rzp(1,icbu) = rzu(irem,1)
                   end if
                else
                   rzp(1,icbu) = 0.0E0_dp
                end if

! determine whether it is inside the sphere with particle jins
                x = rxp(1,icbu) - rxu(jins,1)
                y = ryp(1,icbu) - ryu(jins,1)
                z = rzp(1,icbu) - rzu(jins,1)
                if ( lpbc ) call mimage(x,y,z,boxins)
                rijsq = x*x + y*y + z*z
                if (rijsq .lt. rbsmax**2 .and. rijsq .gt. rbsmin**2)  goto 222
             end do
          end if
       end if
! *************************
! end first Avbmc calcs *
! *************************
    else
       if (lsolid(boxins) .and. .not. lrect(boxins)) then
          ibox = boxins
          do icbu = 1,ichoi
!$$$               if (lbnow.and.boxins.eq.boxsite) then
!$$$                  x = random(-1)
!$$$                  y = random(-1)
!$$$                  z = random(-1)
!$$$
!$$$c     --- gaussian picking probabilities
!$$$c                  x = adev * sqrt(- 2.0E0_dp * log(x*adev*twopiroot))
!$$$c                  y = adev * sqrt(- 2.0E0_dp * log(y*adev*twopiroot))
!$$$c                  z = adev * sqrt(- 2.0E0_dp * log(z*adev*twopiroot))
!$$$
!$$$c     --- old picking probabilities
!$$$                  rvol = random(-1)
!$$$                  r = (rbsmax2**3*rvol +
!$$$     &                 (1-rvol)*rbsmin*rbsmin*rbsmin)**(1.0/3.0E0_dp)
!$$$
!$$$c --- calculate random vector on the unit sphere ---
!$$$ 129              xi1 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
!$$$                  xi2 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
!$$$                  xisq = xi1**2 + xi2**2
!$$$                  if ( xisq .lt. 1.0E0_dp ) then
!$$$                     x = r * 2.0E0_dp * xi1 * sqrt( 1.0E0_dp - xisq )
!$$$                     y = r * 2.0E0_dp * xi2 * sqrt( 1.0E0_dp - xisq )
!$$$                     z = r * ( 1.0E0_dp - 2.0E0_dp * xisq )
!$$$                  else
!$$$                     goto 129
!$$$                  end if
!$$$
!$$$                  rxp(1,icbu) = rxu(jins,1) + x
!$$$                  ryp(1,icbu) = ryu(jins,1) + y
!$$$                  rzp(1,icbu) = rzu(jins,1) + z
!$$$
!$$$               else
             sx = random(-1)
             sy = random(-1)
             sz = random(-1)

             rxp(1,icbu) = sx*hmat(ibox,1)+sy*hmat(ibox,4) +sz*hmat(ibox,7)
             ryp(1,icbu) = sx*hmat(ibox,2)+sy*hmat(ibox,5) +sz*hmat(ibox,8)
             rzp(1,icbu) = sx*hmat(ibox,3)+sy*hmat(ibox,6) +sz*hmat(ibox,9)
!$$$               end if
          end do
       else
          do icbu = 1,ichoi
             rxp(1,icbu) = boxlx(boxins) * random(-1)
             ryp(1,icbu) = boxly(boxins) * random(-1)
             if (lpbcz .or. lslit) then
                rzp(1,icbu) = boxlz(boxins) * random(-1)
             else if ( lsami .or. lmuir .or. ljoe ) then
                if ( lempty ) then
                   rzp(1,icbu) = 20*random(-1)-10
                else
                   rzp(1,icbu) = rzu(irem,1)
                end if
             else
                rzp(1,icbu) = 0.0E0_dp
             end if
          end do
       end if
    end if

! select starting unit ***
! always using bead 1 as the starting unit
    iutry = beg
    lnew = .true.
    glist(1) = beg

! insert the first atom
    call boltz(lnew,.true.,ovrlap,iins,iins,imolty,boxins,ichoi,idum,1,glist,0.0E0_dp)

    if ( lanes .and. lflucq(imolty) ) then
! lfavor is used to set up the favor and favor2 to preferentially sample
! the electronic degrees of freedom
       lfavor = .true.
    else
       lfavor = .false.
       if ( lswapinter ) bnchem(boxins,imolty) = bnchem(boxins,imolty) + 1
    end if

    if ( ovrlap ) return

! perform the walk according to the availibility of the choices ***
! and calculate the correct weight for the trial walk           ***
    w1ins = 0.0E0_dp
    do ip = 1, ichoi
       w1ins = w1ins + bfac(ip)
    end do

! check for termination of walk ---
    if ( w1ins .lt. softlog ) then
       write(io_output,*) 'caught in swap'
       return
    end if

! select one position at random ---
    if ( ichoi .gt. 1 ) then
       rbf = w1ins * random(-1)
       bsum = 0.0E0_dp
       do ip = 1, ichoi
          if ( .not. lovr(ip) ) then
             bsum = bsum + bfac(ip)
             if ( rbf .lt. bsum ) then
! select ip position ---
                iwalk = ip
                exit
             end if
          end if
       end do
       if (ip.gt.ichoi) then
          write(io_output,*) 'w1ins:',w1ins,'rbf:',rbf
          call err_exit(__FILE__,__LINE__,'big time screwup -- w1ins',myid+1)
       end if
    else
       iwalk = 1
    end if

    v1ins =  vtr(1,iwalk)
    v1insext = vtr(9,iwalk)
    v1insint = vtr(2,iwalk)
    v1inselc = vtr(8,iwalk)
    v1insewd = vtr(14,iwalk)
! KM
! if (myid.eq.rootid) write(io_output,*) 'vtry swap ', iwalk, vtr(1,iwalk)

! neigh_icnt = ntr_icnt(iwalk)
! do ip = 1,neigh_icnt
! neighi(ip) = ntr(ip,iwalk)
! end do
    rxnew(beg) = rxp(1,iwalk)
    rynew(beg) = ryp(1,iwalk)
    rznew(beg) = rzp(1,iwalk)

    if (lrigid(imolty)) then
! calculate new vector from initial bead
       do j = beg,iunit
          rxnew(j) = rxnew(beg)  - (rxu(irem,beg) - rxu(irem,j))
          rynew(j) = rynew(beg)  - (ryu(irem,beg) - ryu(irem,j))
          rznew(j) = rznew(beg)  - (rzu(irem,beg) - rzu(irem,j))
          if ((.not.lexpee).and.(.not.leemove)) then
             qqu(iins,j) = qqu(irem,j)
          else
             qqu(iins,j) = ee_qqu(j,nstate)
          end if
       end do
       call schedule(igrow,imolty,ifrom,iutry,0,4)
    else if (lring(imolty)) then
       lfixnow = .true.
       call safeschedule(igrow,imolty,ifrom,iutry,findex,2)
    else
       call schedule(igrow,imolty,ifrom,iutry,0,2)
    end if

!     ------------------------------------------------------------------
    waddnew = 1.0E0_dp
    lterm = .false.

    call rosenbluth(.true.,lterm,iins,iins,imolty,ifrom,boxins,igrow,waddnew,lfixnow,ctorfn,2)

! termination of cbmc attempt due to walk termination ---
    if ( lterm ) return

    if ( ldual .or. lewald .or. iunit .ne. igrow  .or. ((.not. lchgall) .and. lelect(imolty)) ) then
! Put on hydrogens for explicit AA model for calculation of COM
! and assign all of the grown new and old beads to rxuion
! with rxuion: new = 2
       iii = 2
       do j=1,igrow
          rxuion(j,iii) = rxnew(j)
          ryuion(j,iii) = rynew(j)
          rzuion(j,iii) = rznew(j)
          qquion(j,iii) = qqu(iins,j)
       end do
       if ( igrow .ne. iunit .and. .not. lrigid(imolty) ) then
          idum = nchain+1
          do j=1,igrow
             rxu(idum,j) = rxnew(j)
             ryu(idum,j) = rynew(j)
             rzu(idum,j) = rznew(j)
             ! write(io_output,*) rxu(idum,j),ryu(idum,j),rzu(idum,j)
          end do
          moltyp(idum) = imolty
          call explct(idum,v(7),.false.,.false.)
          do j=igrow+1,iunit
             rxuion(j,iii) = rxu(idum,j)
             ryuion(j,iii) = ryu(idum,j)
             rzuion(j,iii) = rzu(idum,j)
             qquion(j,iii) = qqu(iins,j)
             ! write(io_output,*) rxu(idum,j),ryu(idum,j),rzu(idum,j)
          end do
       end if
       moltion(iii) = imolty

       ibox=boxins
       nboxi(iins) = ibox
    end if

! Begin DC-CBMC Corrections for NEW configuration
    if (ldual .or. ((.not. lchgall) .and. lelect(imolty)) .or. (lchgall .and. lewald .and. (.not. ldual))) then
! calculate the true site-site energy
       istt = 1
       iett = igrow
       ! write(io_output,*) igrow

       call energy(iins,imolty,v,iii,ibox,istt,iett,.true.,ovrlap,.false.,.false.,lfavor,.false.)

       if (ovrlap) then
          write(io_output,*) 'iins',iins,'irem',irem
          call err_exit(__FILE__,__LINE__,'strange screwup in DC-CBMC swap',myid+1)
       end if
! v1insewd, vnewewald and vnewintra now accounted for in v from energy
!$$$         delen = v(1) - ( vnew(2) + vnew(9) +vnew(8)) - (v1ins - v1insewd)
       delen = v(1) - ( vnew(2) + vnew(9) +vnew(8) + vnew(4) + vnew(14) + v1ins)
       waddnew = waddnew*exp(-beta*delen)
       vnew(1)     = vnew(1) + delen
       vnew(2) = v(2) - v1insint
       vnew(9)   = v(9) - v1insext
       vnew(8) = v(8) - v1inselc
       vnew(14) = v(14) - v1insewd
       vnew(4) = v(4)
    end if
! End DC-CBMC Corrections for NEW configuration

! Begin Explicit Atom Corrections for NEW configuration
    if ( iunit .ne. igrow ) then
! calculate the true Lennard-Jones energy for the hydrogens
! iii=2 new conformation
       istt = igrow+1
       iett = iunit
       ltors = .false.
       whins = 0.0E0_dp
       ichoi = nchoih(imolty)

       do ip = 1, ichoi
! place the explicit hydrogens nchoih times and compute the
! full intermolecular and intramolecular energy.  Then select
! one with the usual rosenbluth scheme
          if ( ip .eq. 1) then
! ip .eq. 1, uses the hydrogen positions generated above
             do j = igrow+1, iunit
                rxuh(j,ip) = rxuion(j,iii)
                ryuh(j,ip) = ryuion(j,iii)
                rzuh(j,ip) = rzuion(j,iii)
             end do
          else
             idum = nchain+1
! generate positions for the hydrogens
             call explct(idum,v(7),.false.,.false.)
             do j = igrow+1, iunit
                rxuion(j,iii) = rxu(idum,j)
                ryuion(j,iii) = ryu(idum,j)
                rzuion(j,iii) = rzu(idum,j)
                rxuh(j,ip) = rxu(idum,j)
                ryuh(j,ip) = ryu(idum,j)
                rzuh(j,ip) = rzu(idum,j)
             end do
          end if

! ??? problem here on calculation of favor and favor2 when ichoi > 1

          call energy(iins,imolty,v,iii,ibox,istt,iett,.true.,ovrlap,ltors,.true.,lfavor,.false.)
          ! write(io_output,*) 'ovrlap:',ovrlap
          ! if ( iins .eq. 118) write(io_output,*) 'vinter:',v(2)

          if (ovrlap) then
             lovrh(ip) = .true.
          else
             delenh(ip) = v(1) + v(7)

             if ( delenh(ip)*beta .gt. (3.3E0_dp*softcut) ) then
                lovrh(ip) = .true.
             else
                ! calculate the boltzman and rosenbluth weight
                bfach(ip) = exp(-beta*delenh(ip))
                whins = whins + bfach(ip)
                lovrh(ip) = .false.
                vtrhintra(ip) = v(4)
                vtrhinter(ip) = v(2)
                vtrhext(ip) = v(9)
                vtrhtg(ip) = v(7)
                vtrhelect(ip) = v(8)
                vtrhewald(ip) = v(14)
             end if
          end if
       end do

! check for termination of walk
       if ( whins .lt. softlog ) return

       if ( ichoi .gt. 1 ) then
          rbf = whins * random(-1)
          bsum = 0.0E0_dp
          do ip = 1, ichoi
             if ( .not. lovrh(ip) ) then
                bsum = bsum + bfach(ip)
                if ( rbf .lt. bsum ) then
                   iwalk = ip
                   goto 190
                end if
             end if
          end do
          call err_exit(__FILE__,__LINE__,'screw up in explicit hydrogen',myid+1)
       else
          iwalk = 1
       end if

190    waddnew = waddnew * whins
       vnew(1)     = vnew(1) + delenh(iwalk)
       vnew(4) = vnew(4) + vtrhintra(iwalk)
       vnew(2) = vnew(2) + vtrhinter(iwalk)
       vnew(9)   = vnew(9) + vtrhext(iwalk)
       vnew(7) = vnew(7) + vtrhtg(iwalk)
       vnew(8) = vnew(8) + vtrhelect(iwalk)
       vnew(14) = vnew(14) + vtrhewald(iwalk)
       do j = igrow+1, iunit
          rxuion(j,iii) = rxuh(j,iwalk)
          ryuion(j,iii) = ryuh(j,iwalk)
          rzuion(j,iii) = rzuh(j,iwalk)
       end do
    end if
! End Explicit Atom Corrections for NEW configuration

! Begin Ewald-sum Corrections
    if (lewald.and..not.lideal(boxins)) then
! reciprocal space sum
! prepare qquion(jj,1) etc
       moltion(1) = imolty
       if ( lswapinter ) then
          do j = 1,iunit
             qquion(j,1) = 0.0E0_dp
          end do
          call recip(boxins,vrecipn,vrecipo,1)
          delen = vrecipn - vrecipo
          waddnew = waddnew*exp(-beta*delen)
          vnew(8) = vnew(8) + delen
          vnew(1) = vnew(1) + delen
       end if
    end if
! End Ewald-sum Corrections

    if (lpbcz) then
       if (lsolid(boxins) .and. .not. lrect(boxins)) then
          ibox = boxins
          volins = cell_vol(ibox)
       else
          volins=boxlx(boxins)*boxly(boxins)*boxlz(boxins)
       end if
    else
       volins=boxlx(boxins)*boxly(boxins)
    end if

! Begin Fluctuating Charge corrections for NEW configuration
    vnewflucq = 0.0E0_dp
    if ( lflucq(imolty) ) then
       do j = 1,iunit
          qion(j) = qqu(iins,j)
       end do
       call charge(iins, qion, v(11),vdum)
! once we go to fully flexible will need to compute this in boltz
       vnewflucq = v(11)
       vnew(1) = vnew(1) + v(11)
       waddnew = waddnew*exp(-beta*v(11))
    end if
! End Fluctuating Charge corrections for NEW configuration

! Begin Tail corrections for BOXINS with inserted particle
! JLR 11-24-09 don't compute tail correction if lideal(boxins)
    if (ltailc.and.lswapinter.and..not.lideal(boxins)) then
! END JLR 11-24-09
       vinsta = 0.0E0_dp
       do jmt = 1, nmolty
          if ( jmt .eq. imolty ) then
             rho = dble( ncmt(boxins,jmt) + 1 ) / volins
          else
             rho = dble( ncmt(boxins,jmt) ) / volins
          end if
          do imt = 1, nmolty
             if ( imt .eq. imolty ) then
                vinsta = vinsta + dble(ncmt(boxins,imt) + 1)* coru(imt,jmt,rho,boxins)
             else
                vinsta = vinsta + dble(ncmt(boxins,imt))* coru(imt,jmt,rho,boxins)
             end if
          end do
       end do

!$$$         if (boxins .eq. 1 .and. lexzeo) then
!$$$            do jmt = 1,zntype
!$$$               rho = znum(jmt)/volins
!$$$               do imt = 1, nmolty
!$$$                  if (imt.eq.imolty) then
!$$$                     vinsta = vinsta + dble(ncmt(boxins,imt)+1)*
!$$$     &                    coruz(imt,jmt,rho,boxins)
!$$$                  else
!$$$                     vinsta = vinsta + dble(ncmt(boxins,imt))*
!$$$     &                    coruz(imt,jmt,rho,boxins)
!$$$                  end if
!$$$               end do
!$$$            end do
!$$$         end if

! if(LSOLPAR.and.(boxins.eq.2))then
! vinsta = 0.0E0_dp
! end if

       vinsta = vinsta - vbox(3, boxins )
       waddnew = waddnew*exp(-beta*vinsta)
       vnew(1) = vnew(1) + vinsta
       vnew(2) = vnew(2) + vinsta

! this approximate method of tail correction was used for chem. pot.
! until 1-26-98 MGM
! dtest = 0.0E0_dp
! do jmt = 1, nmolty
! if ( jmt .eq. imolty ) then
! rho = dble(ncmt(boxins,jmt)+1)/volins
! else
! rho = dble(ncmt(boxins,jmt))/volins
! end if
! dtest = dtest + 2.0E0_dp*coru(imolty,jmt,rho)
! arg=arg*exp(-beta*2.0E0_dp*coru(imolty,jmt,rho))
! end do
! write(io_output,*) 'dtest',dtest
    else
       vinsta = 0.0E0_dp
    end if
! End Tail corrections for BOXINS with inserted particle

    if ( .not. lanes ) then
       if ( lswapinter ) then
          arg = w1ins * waddnew * weight * volins  / real(ncmt(boxins,imolty)+1,dp)
          acchem(boxins,imolty) = acchem(boxins,imolty)+arg
       end if
    end if

    if (lexpee.and.leemove) then
       imolty = ee_moltyp(mstate)
       moltion(1) = imolty
       do icbu = 1, iunit
          qqu(irem,icbu) = ee_qqu(icbu,mstate)
       end do
       goto 500
    end if

    bsswap(imolty,ipairb,boxins+nbox) =  bsswap(imolty,ipairb,boxins+nbox) + 1

! Compute weights for the molecule to be removed from boxrem

! check that there is at least one molecule in BOXREM ***
    if ( lempty ) then
! write(io_output,*) 'no molecule in BOXREM'
       if (lgrand) then
          if (boxrem.eq.2) write(io_output,*) ' ERROR ***** array too low !'
       end if
       return
    end if

500 continue

! select a position of the first/starting unit at RANDOM ***
! and calculate the boltzmann weight                     ***
! for the chain to be REMOVED                           ***
    rxp(1,1) = rxu(irem,beg)
    ryp(1,1) = ryu(irem,beg)
    rzp(1,1) = rzu(irem,beg)

    ichoi = nchoi1(imolty)

    if ( .not. lswapinter ) then
! *********************
! second AVMBC part *
! *********************
       if (boxins .ne. boxrem) then
          write(io_output,*) 'avbmc, boxins, boxrem:',boxins,boxrem
          call err_exit(__FILE__,__LINE__,'',myid+1)
       end if

       if ( lbias(imolty) ) then
          if (lavbmc1(imolty)) then
! determine whether it is in or out
             x = rxu(irem,1) - rxu(jins,1)
             y = ryu(irem,1) - ryu(jins,1)
             z = rzu(irem,1) - rzu(jins,1)
             if ( lpbc ) call mimage(x,y,z,boxins)
             rijsq = x*x + y*y + z*z
             if ( rijsq .gt. rbsmax*rbsmax .or.  rijsq .lt. rbsmin*rbsmin ) then
                lrem_out = .true.
             else
                lrem_out = .false.
             end if
          end if
          if ( lrem_out ) then
! out-> j case
             wbias_ins = (1-pmbias(imolty))/ (boxlx(boxins)*boxly(boxins)*boxlz(boxins)- vol_eff)

             if (lavbmc2(imolty) .or. lavbmc3(imolty) )  wbias_ins=wbias_ins/dble(neighj_num+1)

! write(io_output,*) 'originally out'
! write(io_output,*) '1:',wbias_ins,boxlx(boxins),vol_eff

             do icbu = 2,ichoi
232             rxp(1,icbu) = boxlx(boxins) * random(-1)
                ryp(1,icbu) = boxly(boxins) * random(-1)
                if (lpbcz .or. lslit) then
                   rzp(1,icbu) = boxlz(boxins) * random(-1)
                else if ( lsami .or. lmuir .or. ljoe ) then
                   if ( lempty ) then
                      rzp(1,icbu) = 20*random(-1)-10
                   else
                      rzp(1,icbu) = rzu(irem,1)
                   end if
                else
                   rzp(1,icbu) = 0.0E0_dp
                end if

! determine whether it is inside the sphere with particle jins
                x = rxp(1,icbu) - rxu(jins,1)
                y = ryp(1,icbu) - ryu(jins,1)
                z = rzp(1,icbu) - rzu(jins,1)
                if ( lpbc ) call mimage(x,y,z,boxins)
                rijsq = x*x + y*y + z*z
                if ( rijsq .lt. rbsmax**2 .and. rijsq .gt. rbsmin**2 )  goto 232
             end do
          else
             if (lavbmc3(imolty) .and. lremk_in) then
! selected from the region bounded by kins
! k -> j case
                target = kins
                wbias_ins=(1-pmbias(imolty))/vol_eff /dble(neighj_num+1)
             else
! selected from the region bounded by jins
                target = jins

                wbias_ins = pmbias(imolty)/vol_eff

                if ( lavbmc2(imolty) )  wbias_ins=wbias_ins/dble( ncmt(boxrem,imolty)-neighj_num-joffset+1)

                if ( lavbmc3(imolty) ) then
                   if ( linsk_in ) then
! j -> k case
                      wbias_ins = wbias_ins/(neighk_num+1)
                   else
! j -> out case
                      wbias_ins=wbias_ins/dble( ncmt(boxrem,imolty)-neighj_num-joffset+1)
                   end if
                end if
             end if

! write(io_output,*) 'originally in'
! write(io_output,*) '2:',wbias_ins,boxlx(boxins),vol_eff

             do icbu = 2,ichoi
! choose a random association distance
                rvol = random(-1)
                r = (rbsmax*rbsmax*rbsmax*rvol +  (1-rvol)*rbsmin*rbsmin*rbsmin)**(1.0/3.0E0_dp)
! calculate random vector on the unit sphere ---
129             xi1 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
                xi2 = ( 2.0E0_dp * random(-1) ) - 1.0E0_dp
                xisq = xi1**2 + xi2**2
                if ( xisq .lt. 1.0E0_dp ) then
                   x = r * 2.0E0_dp * xi1 * sqrt( 1.0E0_dp - xisq )
                   y = r * 2.0E0_dp * xi2 * sqrt( 1.0E0_dp - xisq )
                   z = r * ( 1.0E0_dp - 2.0E0_dp * xisq )
                else
                   goto 129
                end if
                rxp(1,icbu) = rxu(target,1) + x
                ryp(1,icbu) = ryu(target,1) + y
                rzp(1,icbu) = rzu(target,1) + z
             end do

          end if
       end if
! *************************
! end second AVMBC part *
! *************************
    else
       if (lsolid(boxrem) .and. .not. lrect(boxrem)) then
          ibox = boxrem
          do icbu = 2,ichoi
             sx = random(-1)
             sy = random(-1)
             sz = random(-1)
             rxp(1,icbu) = sx*hmat(ibox,1)+sy*hmat(ibox,4) +sz*hmat(ibox,7)
             ryp(1,icbu) = sx*hmat(ibox,2)+sy*hmat(ibox,5) +sz*hmat(ibox,8)
             rzp(1,icbu) = sx*hmat(ibox,3)+sy*hmat(ibox,6) +sz*hmat(ibox,9)
          end do
       else
          do icbu = 2,ichoi
             rxp(1,icbu) = boxlx(boxrem) * random(-1)
             ryp(1,icbu) = boxly(boxrem) * random(-1)
             if (lpbcz .or. lslit) then
                rzp(1,icbu) = boxlz(boxrem) * random(-1)
             else if ( lsami .or. lmuir .or. ljoe ) then
                if ( lempty ) then
                   rzp(1,icbu) = 20*random(-1)-10
                else
                   rzp(1,icbu) = rzu(irem,1)
                end if
             else
                rzp(1,icbu) = 0.0E0_dp
             end if
          end do
       end if
    end if

    if ( .not. lswapinter .and. lbias(imolty) ) then
       if ( lrem_out .and. lins_in ) then
          bnswap_in(imolty,1) = bnswap_in(imolty,1) + 1
       else if ( (.not. lrem_out) .and. (.not. lins_in ) ) then
          bnswap_out(imolty,1) = bnswap_out(imolty,1) + 1
       end if
    end if

! calculate the boltzmann weight of first bead          ***
    lnew = .false.
    call boltz(lnew,.true.,ovrlap,irem,irem,imolty,boxrem,ichoi,idum ,1,glist,0.0E0_dp)

    if ( ovrlap ) then
       write(io_output,*) 'SWAP:1st bead overlap in rembox',boxrem ,' for moltyp',imolty
    end if

! calculate the correct weight for the  old  walk ***
    w1rem = 0.0E0_dp
    do ip = 1, ichoi
       w1rem = w1rem + bfac(ip)
    end do

! check for termination of walk ---
    if ( w1rem .lt. softlog ) then
       write(io_output,*) 'SWAP:soft overlap in rembox',boxrem,' for moltyp' ,imolty
    end if

    v1rem = vtr(1,1)
    v1remint = vtr(2,1)
    v1remext = vtr(9,1)
    v1remelc = vtr(8,1)
    v1remewd = vtr(14,1)

    waddold = 1.0E0_dp

! call rosenbluth for old conformation
    call rosenbluth(.false.,lterm,irem,irem,imolty,ifrom ,boxrem,igrow,waddold,lfixnow,ctorfo,2)

    if ( lterm ) then
! write(io_output,*) 'SWAP: rosenbluth old rejected'
       return
    end if

    if ( ldual .or. lewald .or. igrow .ne. iunit  .or. ((.not. lchgall) .and. lelect(imolty)) ) then
! store the old grown beads and explict placed beads positions
! 1 = old conformation
       iii = 1
       moltion(1) = imolty
       do j = 1,iunit
          rxuion(j,1) = rxu(irem,j)
          ryuion(j,1) = ryu(irem,j)
          rzuion(j,1) = rzu(irem,j)
          qquion(j,1) = qqu(irem,j)
       end do
    end if

! Begin Correction for DC-CBMC for OLD configuration
    if (ldual .or. ((.not. lchgall) .and. lelect(imolty)) .or. (lchgall .and. lewald .and. (.not. ldual))) then
! correct the acceptance rules
! calculate the Full rcut site-site energy
       istt=1
       iett = igrow

       call energy(irem,imolty,v,iii,boxrem,istt,iett,.true.,ovrlap,.false.,.false.,lfavor,.false.)

       if (ovrlap) then
          write(io_output,*) 'disaster ovrlap in old conf SWAP'
          call err_exit(__FILE__,__LINE__,'',myid+1)
       end if
! v now includes vnewintra,v1remewd and voldewald, take out
!$$$         deleo = v(1) - ( vold(2) + vold(9) +vold(8))
!$$$     &        - (v1rem - v1remewd)
       deleo = v(1) - ( vold(2) + vold(9) +vold(8) + vold(4)  + vold(14) + v1rem)
       waddold = waddold*exp(-beta*deleo)
       vold(1)     = vold(1) + deleo
       vold(2) = v(2) - v1remint
       vold(9)   = v(9) - v1remext
       vold(8) = v(8) - v1remelc
       vold(14) = v(14) - v1remewd
       vold(4) = v(4)
    end if
! End Correction for DC-CBMC for OLD configuration

! Begin Correction for Explicit Atom for OLD configuration
    if ( iunit .ne. igrow ) then
! calculate the true Lennard-Jones energy for the hydrogens
! Calculate the energy of the non-backbone beads
! iii=1 old conformation
       ibox = boxrem
       ltors = .true.
       istt = igrow + 1
       iett = iunit

       call energy(irem,imolty,v,iii,ibox,istt,iett,.true.,ovrlap,ltors,.true.,lfavor,.false.)

       ! if (irem .eq. 118)  write(io_output,*) 'for old',v(2)
       if (ovrlap) then
          write(io_output,*) 'disaster ovrlap in old conf SWAP'
          call err_exit(__FILE__,__LINE__,'',myid+1)
       end if
       deleo = v(1) + v(7)
       vold(1)     = vold(1) + deleo
       vold(4) = vold(4) + v(4)
       vold(2) = vold(2) + v(2)
       vold(9)   = vold(9) + v(9)
       vold(7)    = vold(7) + v(7)
       vold(8) = vold(8) + v(8)
       vold(14) = vold(14) + v(14)

       whrem = exp(-beta*deleo)
       ichoi = nchoih(imolty)

       if ( ichoi .ne. 1 ) then
! torsion is calculated in explicit for the placed atoms
          ltors = .false.
! store the old coordinates for the hydrogens
          do j = igrow+1, iunit
             rxuold(j) = rxuion(j,iii)
             ryuold(j) = ryuion(j,iii)
             rzuold(j) = rzuion(j,iii)
          end do

          do ip = 2,ichoi
! rosenbluth weight for multiple placement of explicit hydrogens
             call explct(irem,v(7),.false.,.false.)
             do j = igrow + 1, iunit
                rxuion(j,iii) = rxu(irem,j)
                ryuion(j,iii) = ryu(irem,j)
                rzuion(j,iii) = rzu(irem,j)
             end do

! ??? problem here on calculation of favor and favor2
             call energy(irem,imolty,v,iii,ibox,istt,iett,.true.,ovrlap,ltors,.true.,lfavor,.false.)

             deleo = v(1) + v(7)
             if ( .not. ovrlap ) then
                whrem = whrem + exp(-beta*deleo)
             end if
          end do
          if ( ichoi .gt. 1 ) then
! restore the old hydrogen positions
             do j = igrow+1, iunit
                rxuion(j,iii) = rxuold(j)
                ryuion(j,iii) = ryuold(j)
                rzuion(j,iii) = rzuold(j)
                rxu(irem,j) = rxuold(j)
                ryu(irem,j) = ryuold(j)
                rzu(irem,j) = rzuold(j)
             end do
          end if
       end if

       waddold = waddold*whrem
    end if
! End Correction for Explicit Atom for OLD configuration

! Begin Ewald-sum Corrections for OLD configuration
    if (lewald.and..not.lideal(boxrem)) then
! reciprocal space sum on r*uion
! prepare qquion(jj,1) etc
       if ( lswapinter ) then
          do j = 1,iunit
             qquion(j,2) = 0.0E0_dp
          end do
       end if
       call recip(boxrem,vrecipn,vrecipo,1)
       deleo = vrecipo - vrecipn
       vold(1) = vold(1) + deleo
       vold(8) = vold(8) + deleo
       waddold = waddold * exp(-beta*deleo)
    end if
! End Ewald-sum Corrections for OLD configuration

! Begin Fluctuating Charge corrections for OLD configuration
    voldflucq = 0.0E0_dp
    if ( lflucq(imolty) ) then
! this is a temporary fix - eventually should implement fully
! flexible fluctuating charge calculation into boltz
! right now the old flucq energy is the same as the new flucq
       voldflucq = vnewflucq
       vold(1) = vold(1) + voldflucq
       waddold = waddold*exp(-beta*voldflucq)
    end if
! End Fluctuating Charge corrections for OLD configuration

    if (lpbcz) then
       if (lsolid(boxrem) .and. .not. lrect(boxrem)) then
          ibox = boxrem
          volrem = cell_vol(ibox)
       else
          volrem=boxlx(boxrem)*boxly(boxrem)*boxlz(boxrem)
       end if
    else
       volrem=boxlx(boxrem)*boxly(boxrem)
    end if

! Start of intermolecular tail correction for boxrem
! JLR 11-24-09 don't compute tail correction if lideal(boxrem)
    if (ltailc.and.lswapinter.and..not.lideal(boxrem)) then
! END JLR 11-24-09
! BOXREM without removed particle
       vremta = 0.0E0_dp
       do jmt = 1, nmolty
          if ( jmt .eq. imolty ) then
             rho = dble( ncmt(boxrem,jmt) - 1 ) / volrem
          else
             rho = dble( ncmt(boxrem,jmt) ) / volrem
          end if
          do imt = 1, nmolty
             if ( imt .eq. imolty ) then
                vremta = vremta + dble(ncmt(boxrem,imt) - 1)* coru(imt,jmt,rho,boxrem)
             else
                vremta = vremta + dble(ncmt(boxrem,imt))* coru(imt,jmt,rho,boxrem)
             end if
          end do
       end do

! if(LSOLPAR.and.(boxrem.eq.2)) then
! vremta = 0.0E0_dp
! end if

!$$$         if (boxrem .eq. 1 .and. lexzeo) then
!$$$            do jmt = 1,zntype
!$$$               rho = znum(jmt)/volrem
!$$$               do imt = 1, nmolty
!$$$                  if (imt.eq.imolty) then
!$$$                     vremta = vremta + dble(ncmt(boxrem,imt)-1)*
!$$$     &                    coruz(imt,jmt,rho,boxrem)
!$$$                  else
!$$$                     vremta = vremta + dble(ncmt(boxrem,imt))*
!$$$     &                    coruz(imt,jmt,rho,boxrem)
!$$$                  end if
!$$$               end do
!$$$            end do
!$$$         end if

       vremta = - vremta + vbox(3, boxrem )
       waddold=waddold*exp(-beta*vremta)
       vold(1) = vold(1) + vremta
       vold(2) = vold(2) + vremta
    else
       vremta = 0.0E0_dp
    end if
! End of intermolecular tail correction for boxrem

! Add contributions of the first bead and additional beads:
    vnew(1)     = vnew(1)  + v1ins
    vnew(2) = vnew(2) + v1insint
    vnew(9)   = vnew(9) + v1insext
    vnew(8) = vnew(8) + v1inselc
    vnew(14) = vnew(14) + v1insewd

    vold(1)     = vold(1)  + v1rem
    vold(2) = vold(2)+(v1remint)
    vold(9)   = vold(9)+v1remext
    vold(8) = vold(8) + v1remelc
    vold(14) = vold(14) + v1remewd

    weight= w1ins * waddnew * weight
    weiold= w1rem * waddold * weiold

    wnlog = log10 ( weight )
    wolog = log10 ( weiold )
    wdlog = wnlog - wolog

    if ( wdlog .lt. -softcut ) then
! write(io_output,*) '### underflow in wratio calculation ###'
       return
    end if

! For ANES algorithm, do the Fluctuating charge moves.
    if ( lanes ) then
       if ( lswapinter ) then
          nboxi(irem) = boxins
          parbox(ncmt(boxins,imolty)+1,boxins,imolty)= irem
          parbox(pointp,boxrem,imolty)= parbox(ncmt(boxrem,imolty),boxrem,imolty)
          parbox(ncmt(boxrem,imolty),boxrem,imolty)=0
          nchbox(boxins) = nchbox(boxins) + 1
          nchbox(boxrem) = nchbox(boxrem) - 1
          ncmt(boxins,imolty) = ncmt(boxins,imolty) + 1
          ncmt(boxrem,imolty) = ncmt(boxrem,imolty) - 1
       end if

       favor(irem) = 1.0E0_dp
       favor2(irem) = 1.0E0_dp

       call anes(irem,boxins,boxrem,3,laccept,v,vinsta,vremta,vnewflucq,voldflucq,lswapinter)

       if ( lswapinter ) then
          arg = weight * volins  / dble( ncmt(boxins,imolty) )
          acchem(boxins,imolty) = acchem(boxins,imolty)+arg
          bnchem(boxins,imolty) = bnchem(boxins,imolty) + 1
       end if

       if ( laccept ) then
          bnswap(imolty,ipairb,boxins+nbox)= bnswap(imolty,ipairb,boxins+nbox)+1
       else if ( lswapinter ) then
          nboxi(irem) = boxrem
          parbox(ncmt(boxrem,imolty)+1,boxrem,imolty)=irem
          parbox(ncmt(boxins,imolty),boxins,imolty)=0
          nchbox(boxins) = nchbox(boxins) - 1
          nchbox(boxrem) = nchbox(boxrem) + 1
          ncmt(boxins,imolty) = ncmt(boxins,imolty) - 1
          ncmt(boxrem,imolty) = ncmt(boxrem,imolty) + 1
       end if
       return
    end if

    if ( lswapinter ) then
       if (lgibbs.and.((.not.leemove).and.(.not.lexpee))) then
! Note: acceptance based on only molecules of type imolty
          wratio = ( weight / weiold ) * wee_ratio * ( volins * dble( ncmt(boxrem,imolty) ) / ( volrem * dble( ncmt(boxins,imolty1) + 1 ) ) ) * exp(beta*(eta2(boxrem,imolty)- eta2(boxins,imolty1)))
       else if (lgibbs.and.(leemove.and.lexpee)) then
          wratio = ( weight / weiold ) * wee_ratio * volins/volrem * exp(beta*(eta2(boxrem,imolty)- eta2(boxins,imolty1)))
       else if (lgrand) then
          if (boxins.eq.1) then
! molecule added to box 1
             wratio = (weight /  weiold ) * volins * B(imolty) / (ncmt(boxins,imolty)+1)
          else
! molecule removed from box 1
             wratio = (weight /weiold)* dble(ncmt(boxrem,imolty))/(volrem*B(imolty))
          end if
       end if
    else
       wratio = (weight*wbias_ins)/(weiold*wbias_rem)*wee_ratio
    end if

    if (.not. lswapinter) then
       neigh_old = neigh_cnt(irem)
       if ( neigh_old .le. 6 .and. neigh_icnt .le. 6 ) then
          if ( lavbmc1(imolty) ) then
             if ( lrem_out ) then
                if ( lins_in ) then
                   ip = 1
                else
                   ip = 2
                end if
             else
                if ( lins_in ) then
                   ip = 3
                else
                   ip = 4
                end if
             end if
          else if ( lavbmc2(imolty) ) then
             if ( lins_in ) then
                ip = 1
             else
                ip = 2
             end if
          else
             if ( lins_in ) then
                if ( lremk_in ) then
                   ip = 1
                else
                   ip = 2
                end if
             else if ( linsk_in ) then
                ip =3
             else
                ip =4
             end if
          end if
! write(io_output,*) neigh_old,neigh_icnt,ip
! if ( .not. lins_in .and. neigh_old .eq. 0 ) then
! write(io_output,*) '####error:',irem,jins
! end if
          cnt_wf1(neigh_old,neigh_icnt,ip) = cnt_wf1(neigh_old,neigh_icnt,ip)+1
          if (neigh_old .eq. 0 .and. neigh_icnt .eq. 1) then
             wdlog = log10 (wratio)
             ic = aint((wdlog+95.0E0_dp)/0.1E0_dp)+1
             if ( ic .lt. 1 ) ic = 1
             if ( ic .gt. 1000 ) ic = 1000
             cnt_wra1(ic,ip) = cnt_wra1(ic,ip) + 1
          else if (neigh_old .eq. 1 .and. neigh_icnt .eq. 0) then
             wdlog = log10 (wratio)
             ic = aint((wdlog+95.0E0_dp)/0.1E0_dp) + 1
             if ( ic .lt. 1 ) ic = 1
             if ( ic .gt. 1000 ) ic = 1000
             cnt_wra2(ic,ip) = cnt_wra2(ic,ip) + 1
          end if
       end if
    end if

    if (lfixnow) then
       wratio = wratio * ctorfo / ctorfn
    end if

    if (lopt_bias(imolty)) call update_bias(log(wratio*2.0)/beta,boxrem,boxins,imolty)

    if ( random(-1) .le. wratio ) then
! write(io_output,*) 'SWAP MOVE ACCEPTED',irem
! we can now accept !!!!! ***
       if ((.not.leemove).and.(.not.lexpee)) then
          bnswap(imolty,ipairb,boxins+nbox) =  bnswap(imolty,ipairb,boxins+nbox) + 1
       end if
       if ( .not. lswapinter .and. lbias(imolty) ) then
          if ( lrem_out .and. lins_in ) then
             bnswap_in(imolty,2) = bnswap_in(imolty,2) + 1
          else if ( (.not. lrem_out) .and. (.not. lins_in ) ) then
             bnswap_out(imolty,2) = bnswap_out(imolty,2) + 1
          end if
       end if

! update the position, it will be used to get the bonded energy
       do ic = 1,igrow
          rxu(irem,ic) = rxnew(ic)
          ryu(irem,ic) = rynew(ic)
          rzu(irem,ic) = rznew(ic)
          if (leemove.and.lexpee) qqu(irem,ic) = qqu(iins1,ic)
       end do
       do ic = igrow+1,iunit
          rxu(irem,ic) = rxuion(ic,2)
          ryu(irem,ic) = ryuion(ic,2)
          rzu(irem,ic) = rzuion(ic,2)
       end do

       vtgn = 0.0_dp
       vbendn = 0.0_dp
       vvibn = 0.0_dp
       if (lrigid(imolty).and.(pm_atom_tra.gt.0.000001)) then
          call U_bonded(irem,imolty,vvibn,vbendn,vtgn)
       end if
! total_NBE = vintran+velectn+vewaldn+vtgn+vbendn+vvibn
       total_NBE = vtgn+vbendn+vvibn
! write(io_output,*) vintran,velectn,vewaldn

! write(io_output,*) 'irem', irem

! update energies:
       vbox(1,boxrem)     = vbox(1,boxrem)     - vold(1) - total_NBE
       vbox(2,boxrem)  = vbox(2,boxrem)  - vold(2)
       vbox(3,boxrem)   = vbox(3,boxrem)   - vremta
       vbox(4,boxrem)  = vbox(4,boxrem)  - vold(4)
       vbox(5,boxrem)    = vbox(5,boxrem)    - vold(5) - vvibn
       vbox(7,boxrem)     = vbox(7,boxrem)     - vold(7) - vtgn
       vbox(9,boxrem)    = vbox(9,boxrem)    - vold(9)
       vbox(6,boxrem)   = vbox(6,boxrem)   - vold(6) - vbendn
       vbox(8,boxrem)  = vbox(8,boxrem)  - (vold(8)+vold(14))
       vbox(11,boxrem)  = vbox(11,boxrem)  - voldflucq

       vbox(1,boxins)     = vbox(1,boxins)     + vnew(1) + total_NBE
       vbox(2,boxins)  = vbox(2,boxins)  + vnew(2)
       vbox(3,boxins)   = vbox(3,boxins)   + vinsta
       vbox(4,boxins)  = vbox(4,boxins)  + vnew(4)
       vbox(5,boxins)    =  vbox(5,boxins)   + vnew(5) + vvibn
       vbox(7,boxins)     = vbox(7,boxins)     + vnew(7) + vtgn
       vbox(9,boxins)    = vbox(9,boxins)    + vnew(9)
       vbox(6,boxins)   = vbox(6,boxins)   + vnew(6) + vbendn
       vbox(8,boxins)  = vbox(8,boxins)  + (vnew(8)+vnew(14))
       vbox(11,boxins)  = vbox(11,boxins)  + vnewflucq

! update book keeping
       if ( lswapinter ) then
          nboxi(irem) = boxins
          if ((.not.leemove).and.(.not.lexpee)) then
             parbox(ncmt(boxins,imolty)+1,boxins,imolty)= irem
             parbox(pointp,boxrem,imolty)= parbox(ncmt(boxrem,imolty),boxrem,imolty)
             parbox(ncmt(boxrem,imolty),boxrem,imolty)=0

             nchbox(boxins) = nchbox(boxins) + 1
             nchbox(boxrem) = nchbox(boxrem) - 1
             ncmt(boxins,imolty) = ncmt(boxins,imolty) + 1
             ncmt(boxrem,imolty) = ncmt(boxrem,imolty) - 1
          else
             parbox(ncmt(boxins,imolty1)+1,boxins,imolty1)= irem
             parbox(pointp,boxrem,imolty)= parbox(ncmt(boxrem,imolty),boxrem,imolty)
             parbox(ncmt(boxrem,imolty),boxrem,imolty)=0

             nchbox(boxins) = nchbox(boxins) + 1
             nchbox(boxrem) = nchbox(boxrem) - 1
             ncmt(boxins,imolty1) = ncmt(boxins,imolty1) + 1
             ncmt(boxrem,imolty) = ncmt(boxrem,imolty) - 1
             eepointp = ncmt(boxins,imolty1)
          end if

          if ( lexpand(imolty) ) then
             itype = eetype(imolty)
             ncmt2(boxins,imolty,itype) =  ncmt2(boxins,imolty,itype) + 1
             ncmt2(boxrem,imolty,itype) =  ncmt2(boxrem,imolty,itype) + 1
          end if
       end if

! do ic = 1,igrow
! rxu(irem,ic) = rxnew(ic)
! ryu(irem,ic) = rynew(ic)
! rzu(irem,ic) = rznew(ic)
! end do
! do ic = igrow+1,iunit
! rxu(irem,ic) = rxuion(ic,2)
! ryu(irem,ic) = ryuion(ic,2)
! rzu(irem,ic) = rzuion(ic,2)
! end do

       if ( lewald ) then
! update reciprocal-space sum
          if ( lswapinter ) then
             if (.not.lideal(boxins)) call recip(boxins,vdum,vdum,2)
             if (.not.lideal(boxrem)) call recip(boxrem,vdum,vdum,2)
          else if (.not.lideal(boxins)) then
             call recip(boxins,vdum,vdum,2)
          end if
       end if

! update center of mass
       call ctrmas(.false.,boxins,irem,3)
! update linkcell, if applicable
       if ( licell .and. ((boxins .eq. boxlink) .or. (boxrem .eq. boxlink))) then
          call update_linked_cell(irem)
       end if

! KM for HD, 01/2010
! for computing dielectric constant with swap/AVBMC moves
! unclear as yet if this is correct!
       if (ldielect) then
! update dipole term
          call dipole(ibox,1)
       end if

       if ( lneigh ) call update_neighbor_list(irem,0.,0.,0.,.true.)
       if ( lneighbor ) then
          neigh_old = neigh_cnt(irem)
          if ( neigh_old .le. 6 .and. neigh_icnt .le. 6 ) then
             cnt_wf2(neigh_old,neigh_icnt,ip) = cnt_wf2(neigh_old,neigh_icnt,ip)+1
          end if

          do ic = 1, neigh_old
             j = neighbor(ic,irem)
             do ip = 1,neigh_cnt(j)
                if ( neighbor(ip,j) .eq. irem ) then
                   neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                   neigh_cnt(j) = neigh_cnt(j)-1
                   exit
                end if
             end do
          end do
          neigh_cnt(irem) = neigh_icnt
! write(io_output,*) 'irem:',irem,neigh_icnt
          do ic = 1,neigh_icnt
             j = neighi(ic)
             neighbor(ic,irem)=j
             lneighij = .false.
             do ip = 1,neigh_cnt(j)
                if ( neighbor(ip,j) .eq. irem ) then
                   lneighij = .true.
                end if
             end do
             if ( .not. lneighij ) then
                neigh_cnt(j) = neigh_cnt(j)+1
                neighbor(neigh_cnt(j),j) = irem
             end if
          end do
       end if
    end if
! -----------------------------------------------------------------
#ifdef __DEBUG__
    write(io_output,*) 'end SWAP in ',myid,irem
#endif
    return
  end subroutine swap

  subroutine init_swap
    integer::jerr
    allocate(acchem(nbxmax,ntmax),bnchem(nbxmax,ntmax),bsswap(ntmax,npabmax,nbxmax*2),bnswap(ntmax,npabmax,nbxmax*2),bnswap_in(ntmax,2),bnswap_out(ntmax,2),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'init_swap: allocation failed',jerr)
    end if

    bnswap=0
    bsswap=0
    bnswap_in=0
    bnswap_out=0
    cnt_wf1=0
    cnt_wf2=0
    cnt_wra1=0
    cnt_wra2=0
    acchem=0.0_dp
    bnchem=0
  end subroutine init_swap

  subroutine output_swap_stats(io_output)
    integer,intent(in)::io_output
    integer::i,j,ibox,jbox,jbox_max

    write(io_output,*)
    write(io_output,*) '### Molecule swap       ###'
    write(io_output,*)
    do i = 1, nmolty
       write(io_output,*) 'molecule typ =',i
       do j=1,nswapb(i)
          if ( box1(i,j) .eq. box2(i,j) ) then
             jbox_max = 1
          else
             jbox_max = 2
          end if
          do jbox = 1,jbox_max
             if ( jbox .eq. 1 ) ibox = box1(i,j)
             if ( jbox .eq. 2 ) ibox = box2(i,j)
             write(io_output,"('between box ',i0,' and ',i0,' into box ',i0, '   uattempts = ',I0,' attempts = ',I0 ,'   accepted = ',I0)") box1(i,j),box2(i,j),ibox,bsswap(i,j,ibox),bnswap(i,j,ibox),bnswap(i,j,ibox+nbox)
             if (bnswap(i,j,ibox) .gt. 0) then
                write(io_output,"(' suc.growth % =',f7.3,'   accepted % =',f7.3)") bsswap(i,j,ibox+nbox)*100.0_dp/bnswap(i,j,ibox),bnswap(i,j,ibox+nbox)*100.0_dp/bsswap(i,j,ibox)
             end if
          end do
       end do
       write(io_output,"('number of times move in: ', I0,  '  accepted = ',I0)") bnswap_in(i,1), bnswap_in(i,2)
       write(io_output,"('number of times move out: ', I0,  '  accepted = ',I0)") bnswap_out(i,1), bnswap_out(i,2)
    end do
  end subroutine output_swap_stats

  subroutine cnt()
    logical::lpr
    integer::i,j,nnn

    lpr = .false.
    do i = 1,ntmax
       if (lbias(i)) then
          lpr = .true.
       end if
    end do

    if (lpr.and.myid.eq.rootid) then
       do nnn = 1,4
          write(31,*)
          write(31,*) 'nnn:',nnn
          do i = 0,6
             do j = 0,6
                if (cnt_wf1(i,j,nnn) .gt. 0 ) then
                   write(31,*) i,j,cnt_wf1(i,j,nnn),cnt_wf2(i,j,nnn)
                end if
             end do
          end do
          write(32,*)
          write(32,*) 'nnn:', nnn
          write(33,*)
          write(33,*) 'nnn:', nnn

          do i = 1,1000
             if ( cnt_wra1(i,nnn) .gt. 0 )  write(32,*) (dble(i)-0.5E0_dp)* 0.1E0_dp-95.0E0_dp,cnt_wra1(i,nnn)
             if ( cnt_wra2(i,nnn) .gt. 0 )  write(33,*) (dble(i)-0.5E0_dp)* 0.1E0_dp-95.0E0_dp,cnt_wra2(i,nnn)
          end do
       end do
    end if
  end subroutine cnt

  subroutine read_checkpoint_swap(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) bnswap,bsswap,bnchem,acchem,bnswap_in,bnswap_out
    call mp_bcast(bnswap,ntmax*npabmax*nbxmax*2,rootid,groupid)
    call mp_bcast(bsswap,ntmax*npabmax*nbxmax*2,rootid,groupid)
    call mp_bcast(bnchem,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(acchem,nbxmax*ntmax,rootid,groupid)
    call mp_bcast(bnswap_in,ntmax*2,rootid,groupid)
    call mp_bcast(bnswap_out,ntmax*2,rootid,groupid)
  end subroutine read_checkpoint_swap

  subroutine write_checkpoint_swap(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) bnswap,bsswap,bnchem,acchem,bnswap_in,bnswap_out
  end subroutine write_checkpoint_swap
end module transfer_swap
