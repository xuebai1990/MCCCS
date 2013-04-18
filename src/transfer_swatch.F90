module transfer_swatch
  use util_runtime,only:err_exit
  use util_random,only:random
  use sim_system
  use sim_cell
  use energy_kspace,only:recip
  use energy_pairwise,only:energy,coru
  use moves_cbmc,only:rosenbluth,schedule,explct,lexshed
  use transfer_shared,only:lopt_bias,update_bias
  implicit none
  private
  save
  public::swatch,init_swatch,output_swatch_stats,read_checkpoint_swatch,write_checkpoint_swatch

  real,allocatable::bnswat(:,:),bnswat_empty(:,:),bsswat(:,:) !< accumulators for swatch performance
contains
!> Added intrabox move for two particles within one box
!> in combined move that shares the same parameters.
!> Will also accept rigid (lrigid) molecules.  Contains
!> several critical bug fixes as well.
!> \since 9-25-02 JMS
  subroutine swatch()
    use sim_particle,only:update_neighbor_list,ctrmas

      logical::lempty,lterm

      integer::type_a,type_b,from(2*numax),prev(2*numax),self,iboxnew,iboxold,imolty,igrow,new,old,islen,ifirst,iprev,iii,j
      integer::oldchain,newchain,oldunit,newunit,iunit,iins

      integer::ic,ibox,icbu,jj,mm,imt,jmt,imolin,imolrm
      integer::boxa,boxb,ipair,imolta,imoltb,iboxa,iboxb,iboxal,iboxbl,iboxia,iboxib,iunita,iunitb,orgaia,orgaib,orgbia,orgbib,ipairb

      real::tweight,tweiold,rxut(4,numax),ryut(4,numax),rzut(4,numax),dvol,vola,volb,rho,dinsta,rpair

      real::vnbox(nbxmax),vninte(nbxmax),vnintr(nbxmax),vnvibb(nbxmax),vntgb(nbxmax),vnextb(nbxmax),vnbend(nbxmax),vntail(nbxmax),vnelect(nbxmax),vnewald(nbxmax),wnlog,wolog,wdlog,wswat

      real::v(nEnergy),vdum,delen,deleo,dicount,vrecipn,vrecipo

      ! additions from iswatch
      integer::izz,box,iboxi,bdmol_a,bdmol_b,iparty
      integer::imola,imolb,moltaid,moltbid
      integer::s_type,o_type,thisbox,otherbox
      real::rx_1(numax),ry_1(numax),rz_1(numax)
      ! end additions

      real::waddold,waddnew,vdum2

      ! JLR 11-24-09
      integer::icallrose
      ! END JLR 11-24-09

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef __DEBUG__
      write(io_output,*) 'start SWATCH in ',myid
#endif

      lempty = .false.
! randomly select chains to switch between boxa boxb
! select a pair type to switch
      if ( nswaty .gt. 1 ) then
         rpair = random(-1)
         do 96 ipair = 1, nswaty
            if ( rpair .lt. pmsatc(ipair) ) then
               iparty = ipair
               goto 97
            end if
 96      continue
      else
         iparty = 1
      end if

! randomly select the molecules from the pair ***

      if (random(-1) .lt. 0.5E0_dp) then
         type_a = 1
         type_b = 2
      else
         type_a = 2
         type_b = 1
      end if

! select the molecules from the pair
 97   imolta = nswatb(iparty,1)
      imoltb = nswatb(iparty,2)
!> \bug no longer random
      type_a = 1
      type_b = 2

! choose box A and box B at random
      if ( nswtcb(iparty) .gt. 1 ) then
         rpair = random(-1)
         do 98 ipair = 1, nswtcb(iparty)
            if ( rpair .lt. pmswtcb(iparty,ipair) ) then
               ipairb = ipair
               goto 99
            end if
 98      continue
      else
         ipairb = 1
      end if

 99   if (random(-1) .lt. 0.5E0_dp) then
         boxa=box3(iparty,ipairb)
         boxb=box4(iparty,ipairb)
      else
         boxa=box4(iparty,ipairb)
         boxb=box3(iparty,ipairb)
      end if

      if (boxa .eq. boxb) then

! ****************************************************************
! INTRABOX SWATCH ADD *
! ****************************************************************
! ****************************************************************

! liswatch prevents non-grown beads from being included in
! the new growth in boltz.f

! write(io_output,*) 'start iSWATCH'

!$$$c *** randomly select chains to switch ***
!$$$c *** select a pair type to switch ***
!$$$      if ( niswaty .gt. 1 ) then
!$$$         rpair = random(-1)
!$$$         do 96 ipair = 1, niswaty
!$$$            if ( rpair .lt. pmisatc(ipair) ) then
!$$$               iparty = ipair
!$$$               goto 97
!$$$            end if
!$$$ 96      continue
!$$$      else
!$$$         iparty = 1
!$$$      end if

!$$$c *** box determinations ***
!$$$
!$$$      if ( lgibbs ) then
!$$$         if (random(-1) .lt. 0.5E0_dp) then
!$$$            thisbox = 1
!$$$            otherbox = 2
!$$$         else
!$$$            thisbox = 2
!$$$            otherbox = 1
!$$$         end if
!$$$      else
!$$$         thisbox = 1
!$$$         otherbox = 2
!$$$      end if
!$$$
!$$$      box = thisbox

! box determinations ***

         if ( (ncmt(boxa,imolta) .eq. 0) .or. (ncmt(boxb,imoltb) .eq. 0) ) then
            lempty = .true.
         else

           thisbox = boxa
           box = thisbox

           if ( lgibbs ) then
              if (thisbox .lt. nbxmax) then
                 otherbox = thisbox + 1
              else
                 otherbox = thisbox - 1
              end if
           else
              otherbox = 2
           end if

! get particle of type a ***

           moltaid = int(dble (ncmt(box,imolta))*random(-1)) + 1

! imola is the overall molecule id number (runs from 1 to total #) *
           imola = parbox( moltaid, box, imolta)

           if (moltyp(imola) .ne. imolta)  write(io_output,*) 'screwup in iswatch'

           iboxi = nboxi(imola)

           if (iboxi .ne. box) call err_exit(__FILE__,__LINE__,'problem in iswatch',myid+1)

! get particle of type b ***

           moltbid = int(dble (ncmt(box,imoltb))*random(-1)) + 1
           imolb = parbox( moltbid, box, imoltb)

           if (moltyp(imolb) .ne. imoltb)  write(io_output,*) 'screwup in iswatch'
           iboxi = nboxi(imolb)

           if (iboxi .ne. box) call err_exit(__FILE__,__LINE__,'problem in iswatch',myid+1)
         end if
! add one attempt to the count for iparty
!$$$  bniswat(iparty,box) = bniswat(iparty,box) + 1.0E0_dp
         bnswat(iparty,ipairb) = bnswat(iparty,ipairb) + 1.0E0_dp
! JLR 12-1-09 count the empty attempts
         if (lempty) then
            bnswat_empty(iparty,ipairb) = bnswat_empty(iparty,ipairb) + 1.0E0_dp
            return
         end if
! END JLR 12-1-09

! write out the molecule numbers of the pair *
! write(io_output,*) imola,imolb

!     ***************************
! Begin Growth Setups ***
!     ***************************

! assign from and prev for each moltyp

         do izz = 1,ncut(iparty,type_a)

            from(type_a+2*(izz-1)) = gswatc(iparty,type_a,1+2*(izz-1))
            prev(type_a+2*(izz-1)) = gswatc(iparty,type_a,2+2*(izz-1))

         end do

         do izz = 1,ncut(iparty,type_b)

            from(type_b+2*(izz-1)) = gswatc(iparty,type_b,1+2*(izz-1))
            prev(type_b+2*(izz-1)) = gswatc(iparty,type_b,2+2*(izz-1))

         end do

!$$$         write(io_output,*) 'mol a'
!$$$
!$$$         do izz = 1,ncut(type_a)
!$$$            write(io_output,*) izz,'from',from(type_a+2*(izz-1)),' prev',
!$$$     &        prev(type_a+2*(izz-1))
!$$$         end do
!$$$
!$$$         write(io_output,*) 'mol b'
!$$$
!$$$         write(io_output,*) gswatc(iparty,type_b,3),gswatc(iparty,type_b,4)
!$$$         write(io_output,*) from(3),prev(3),type_b
!$$$
!$$$         do izz = 1,ncut(type_b)
!$$$            write(io_output,*) izz,'from',from(type_b+2*(izz-1)),' prev',
!$$$     &        prev(type_b+2*(izz-1))
!$$$         end do
!$$$
!$$$         call err_exit(__FILE__,__LINE__,'',myid+1)

! store number of units in iunita and iunitb
         iunita = nunit(imolta)
         iunitb = nunit(imoltb)

! initialize trial weights ***
         tweight = 1.0E0_dp
         tweiold = 1.0E0_dp

! set the trial energies to zero ***
         vnbox(box)   = 0.0E0_dp
         vninte(box)  = 0.0E0_dp
         vnintr(box)  = 0.0E0_dp
         vnvibb(box)  = 0.0E0_dp
         vntgb(box)   = 0.0E0_dp
         vnextb(box)  = 0.0E0_dp
         vnbend(box)  = 0.0E0_dp
         vnelect(box) = 0.0E0_dp
         vnewald(box) = 0.0E0_dp

! store position 1, a's original site ***

         do izz = 1, nunit(imolta)

            rx_1(izz) = rxu(imola,izz)
            ry_1(izz) = ryu(imola,izz)
            rz_1(izz) = rzu(imola,izz)

         end do

! store same bead coordinates for molecules a and b ***

         do icbu = 1, nsampos(iparty)

            bdmol_a = splist(iparty,icbu,type_a)

            rxut(3,icbu) = rxu(imola,bdmol_a)
            ryut(3,icbu) = ryu(imola,bdmol_a)
            rzut(3,icbu) = rzu(imola,bdmol_a)

            bdmol_b = splist(iparty,icbu,type_b)

            rxut(4,icbu) = rxu(imolb,bdmol_b)
            ryut(4,icbu) = ryu(imolb,bdmol_b)
            rzut(4,icbu) = rzu(imolb,bdmol_b)

         end do

!     *****************************
! Liswinc Determination ***
!     *****************************

! only need this for molecule b, since b is grown after a
! but a is grown when some of b doesn't exist

         self = imolb
         other = imola
         imolty = imoltb
         igrow = nugrow(imoltb)

         ifirst = from(type_b)
         iprev = prev(type_b)

! determine which beads aren't in the same positions *

         call schedule(igrow,imolty,islen,ifirst,iprev,3)

         if (ncut(iparty,type_b) .gt. 1) then
            do izz = 2,ncut(iparty,type_b)
               ifirst = from(type_b+2*(izz-1))
               iprev = prev(type_b+2*(izz-1))

               call schedule(igrow,imolty,islen,ifirst,iprev,5)

            end do
         end if

! assign growth schedule for molecule b *

         do izz = 1,nunit(imolty)
            liswinc(izz,imolty) = lexshed(izz)
            ! write(io_output,*) izz, liswinc(izz,imolty)
         end do

!     ******************
! ROSENBLUTH ***
!     ******************

! start rosenbluth weight calculations ***
         do ic = 1,2
            if ( ic .eq. 1 ) then
               liswatch = .true.

               self = imola
               other = imolb

               s_type = type_a
               o_type = type_b

               imolty = imolta
               igrow = nugrow(imolta)
               iunit = nunit(imolta)

            else
               liswatch = .false.

               self = imolb
               other = imola

               s_type = type_b
               o_type = type_a

               imolty = imoltb
               igrow = nugrow(imoltb)
               iunit = nunit(imoltb)

            end if

! assign same bead coordinates for new growth ***
            do icbu = 1, nsampos(iparty)

               new = splist(iparty,icbu,s_type)
               old = splist(iparty,icbu,o_type)

               rxnew(new) = rxu(other,old)
               rynew(new) = ryu(other,old)
               rznew(new) = rzu(other,old)

            end do

! set up growth schedule
            ifirst = from(s_type)
            iprev = prev(s_type)

! lrigid add on

            if (lrigid(imolty)) then
!cc--- JLR 11-24-09
!cc--- adding some if statements for rigid swaps:
!cc--- if we have swapped the whole rigid part we will
!cc--- not compute the vectors from ifirst in old position
               if ( nsampos(iparty).lt.iunit) then

                  if ( (rindex(imolty).eq.0) .or. (ifirst.lt.riutry(imolty,1))  ) then

                     if (nsampos(iparty).ge.3) then

                        call align_planes(iparty,self,other, s_type,o_type,rxnew,rynew,rznew)

                     else

! calculate new vector from initial bead
                        do j = 1,iunit
                           rxnew(j) = rxnew(ifirst) - (rxu(self,ifirst) - rxu(self,j))
                           rynew(j) = rynew(ifirst) - (ryu(self,ifirst) - ryu(self,j))
                           rznew(j) = rznew(ifirst) - (rzu(self,ifirst) - rzu(self,j))
                        end do

                     end if      !nsampos.eq.3
                  end if         !ifirst.lt.riutry
               end if            !nsampos.lt.iunit
               call schedule(igrow,imolty,islen,ifirst,iprev,4)
!cc---END JLR 11-24-09
            else
               call schedule(igrow,imolty,islen,ifirst,iprev,3)
            end if

!> \bug Need to check: I wonder if this works with lrigid?
            if (ncut(iparty,s_type) .gt. 1) then
               do izz = 2,ncut(iparty,s_type)
                  ifirst = from(s_type+2*(izz-1))
                  iprev = prev(s_type+2*(izz-1))

                  call schedule(igrow,imolty,islen,ifirst,iprev,5)

               end do
            end if



! Commenting out for new combined code JMS 6-20-00 *****
! Call qqcheck to setup the group based qq cutoff ***
!$$$  if ( lelect(imolty) ) then
!$$$  call qqcheck(other,box,rxnew(1),rynew(1),rznew(1))
!$$$  end if

!     ******************
! new growth ***
!     ******************

! moving molecules for rosenbluth ***
            if (ic .eq. 1) then

! putting molecule b in position 1 *

               do izz = 1, nsampos(iparty)

                  bdmol_b = splist(iparty,izz,type_b)

                  rxu(other,bdmol_b) = rxut(3,izz)
                  ryu(other,bdmol_b) = ryut(3,izz)
                  rzu(other,bdmol_b) = rzut(3,izz)

               end do

            else

! putting molecule a into its (fully grown) trial position 2 *

               do izz = 1, nunit(moltyp(other))

                  rxu(other,izz) = rxut(1,izz)
                  ryu(other,izz) = ryut(1,izz)
                  rzu(other,izz) = rzut(1,izz)

               end do

            end if

! grow molecules

! changing for lrigid to include waddnew
            waddnew = 1.0E0_dp

! grow new chain conformation
! JLR 11-24-09
! Different logic/calls to rosenbluth for rigid molecules
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there
                  !don't regrow anything in rosenbluth
                  icallrose = 4
               else if( (nsampos(iparty).ge.3) .or. (ifirst.ge.riutry(imolty,1) )) then
                  ! rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
                  ! rigid part is not grown, do rigrot
                  icallrose = 2
               end if
            else
               ! flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            end if

            call rosenbluth(.true.,lterm,self,self,imolty,islen ,box ,igrow,waddnew,.false.,vdum2,icallrose)
! END JLR 11-24-09

! propagate new rosenbluth weight
            tweight = tweight*weight*waddnew

! end rigid add on

! moving molecules back *
            if (ic .eq. 1) then

               do izz = 1, nsampos(iparty)

                  bdmol_b = splist(iparty,izz,type_b)

                  rxu(other,bdmol_b) = rxut(4,izz)
                  ryu(other,bdmol_b) = ryut(4,izz)
                  rzu(other,bdmol_b) = rzut(4,izz)

               end do

            else

               do izz = 1, nunit(moltyp(other))
                  rxu(other,izz) = rx_1(izz)
                  ryu(other,izz) = ry_1(izz)
                  rzu(other,izz) = rz_1(izz)
               end do

            end if

            if ( lterm ) then
! termination of cbmc attempt due to walk termination ***
! write(io_output,*) 'iSWATCH:new growth rejected',ic
! reset liswatch
               liswatch = .false.
               return
            else
! write(io_output,*) 'iSWATCH:new growth accepted',ic
            end if

! save the new coordinates ***
            do jj = 1,igrow
               rxut(ic,jj) = rxnew(jj)
               ryut(ic,jj) = rynew(jj)
               rzut(ic,jj) = rznew(jj)
            end do

! Corrections for switched beads, and DC-CBMC
! Assign all of the grown new and old beads to rxuion
! with rxuion: new = 2
            iii = 2
            do j=1,igrow
               rxuion(j,iii) = rxnew(j)
               ryuion(j,iii) = rynew(j)
               rzuion(j,iii) = rznew(j)
               qquion(j,iii) = qqu(self,j)
            end do

! added from new-combined code ***
            if ( iunit .ne. igrow ) then

! for explicit-hydrogen model, put on the hydrogens
! use phony number iins and call explct to add constrained hydrogens

               iins = nchain + 1
               moltyp(iins) = imolty

               do j=1,igrow
                  rxu(iins,j) = rxnew(j)
                  ryu(iins,j) = rynew(j)
                  rzu(iins,j) = rznew(j)
               end do

               call explct(iins,vdum,.false.,.false.)

               do j = igrow + 1, iunit

                  rxuion(j,iii) = rxu(iins,j)
                  ryuion(j,iii) = ryu(iins,j)
                  rzuion(j,iii) = rzu(iins,j)
                  qquion(j,iii) = qqu(self,j)

                  rxut(ic,j) = rxu(iins,j)
                  ryut(ic,j) = ryu(iins,j)
                  rzut(ic,j) = rzu(iins,j)

               end do
            end if


! Begin DC-CBMC and switched bead Corrections for NEW configuration
! calculate the true site-site energy

            if (ic .eq. 1) then
! exclude molecule b from energy calculation, put in other box *
               nboxi(imolb) = otherbox

            else
! put molecule a into position 2 (fully grown trial position)
! for the energy of b's new position (second time around)
               do izz = 1,nunit(imolta)
                  rxu(imola,izz) = rxut(1,izz)
                  ryu(imola,izz) = ryut(1,izz)
                  rzu(imola,izz) = rzut(1,izz)
               end do
            end if

! get energy of configuration ***

            call ctrmas(.false.,box,self,8)
            call energy(self,imolty,v,iii,box,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)

! return to normal ***
            if (ic .eq. 1) then

! return b to original box *
               nboxi(imolb) = thisbox
            else

! return a to position 1 *
               do izz = 1, nunit(imolta)
                  rxu(imola,izz) = rx_1(izz)
                  ryu(imola,izz) = ry_1(izz)
                  rzu(imola,izz) = rz_1(izz)

               end do

            end if

            if (lterm) then
               write(io_output,*) 'other ',other,' self ',self,ic
               call err_exit(__FILE__,__LINE__,'interesting screwup in CBMC iswatch',myid+1)
            end if

! add on the changes in energy ***

            delen = v(1) - ( vnew(1) - (vnew(5) + vnew(6) + vnew(7) ))

            tweight = tweight*exp(-beta*delen)
            vnew(1)     = vnew(1) + delen
            vnew(2) = v(2)
            vnew(4) = v(4)
            vnew(9)   = v(9)
            vnew(8) = v(8)
            vnew(14) = v(14)

! End DC-CBMC and switched bead Corrections for NEW configuration

! save the trial energies ***
            vnbox(box)   = vnbox(box) + vnew(1)
            vninte(box)  = vninte(box) + vnew(2)
            vnintr(box)  = vnintr(box) + vnew(4)
            vnvibb(box)  = vnvibb(box) + vnew(5)
            vntgb(box)   = vntgb(box) + vnew(7)
            vnextb(box)  = vnextb(box) + vnew(9)
            vnbend(box)  = vnbend(box) + vnew(6)
            vnelect(box) = vnelect(box) + vnew(8)
            vnewald(box) = vnewald(box) + vnew(14)

!     ******************
! old growth ***
!     ******************

! lrigid add on

            waddold = 1.0E0_dp

! grow old chain conformation
! JLR 11-24-09
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
! molecule is all there
! don't regrow anything in rosenbluth
                  icallrose = 4
               else if( (nsampos(iparty).ge.3)      .or. (ifirst.ge.riutry(imolty,1)) ) then
! rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
! rigid part is not grown, do rigrot
                  icallrose = 2
               end if
            else
! flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            end if

            call rosenbluth(.false.,lterm,self,self,imolty, islen,box,igrow,waddold,.false.,vdum2,icallrose)
! END JLR 11-24-09

            if ( lterm ) then
! termination of old walk due to problems generating orientations ***
! write(io_output,*) 'iSWATCH:old growth rejected',ic
! reset liswatch
               liswatch = .false.
               return
            else
! write(io_output,*) 'iSWATCH:old growth accepted',ic
            end if

! propagate old rosenbluth weight
            tweiold = tweiold*weiold*waddold

! end rigid add on

! store the old grown beads and explict placed beads positions
! 1 = old conformation
            iii = 1
            do j = 1,iunit
               rxuion(j,1) = rxu(self,j)
               ryuion(j,1) = ryu(self,j)
               rzuion(j,1) = rzu(self,j)
               qquion(j,1) = qqu(self,j)
            end do

! Begin Correction for DC-CBMC and switched beads for OLD configuration
! correct the acceptance rules
! calculate the Full rcut site-site energy

! excluding molecule b for first loop ***
            if (ic .eq. 1) then
               nboxi(imolb) = otherbox
            end if

! get total energy *

            call energy(self,imolty,v,iii,box,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)

            if (ic .eq. 1) then
! return b to current box *
               nboxi(imolb) = thisbox
            end if

            if (lterm) call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf iSWATCH',myid+1)
            deleo = v(1) - ( vold(1) - (vold(5) + vold(6) + vold(7)) )

            tweiold = tweiold*exp(-beta*deleo)
            vold(1)     = vold(1) + deleo
            vold(4) = v(4)
            vold(2) = v(2)
            vold(9)   = v(9)
            vold(8) = v(8)
            vold(14) = v(14)

! End Correction for DC-CBMC and switched beads for OLD configuration

! save the trial energies ***
            vnbox(box)   = vnbox(box)  - vold(1)
            vninte(box)  = vninte(box) - vold(2)
            vnintr(box)  = vnintr(box) - vold(4)
            vnvibb(box)  = vnvibb(box) - vold(5)
            vntgb(box)   = vntgb(box)  - vold(7)
            vnextb(box)  = vnextb(box) - vold(9)
            vnbend(box)  = vnbend(box) - vold(6)
            vnelect(box) = vnelect(box) - vold(8)
            vnewald(box) = vnewald(box) - vold(14)

         end do

!     *****************************
! Ewald sum corrections ***
!     *****************************

! Perform the Ewald sum reciprical space corrections
         if (lewald.and..not.lideal(box)) then
! added into tweight even though it really contains new-old

! store the reciprocal space vector
            call recip(box,vdum,vdum,3)

! Position 1 ***
            oldchain = imola
            newchain = imolb
            oldunit = nunit(imolta)
            newunit = nunit(imoltb)

            do j = 1,oldunit
               rxuion(j,1) = rxu(oldchain,j)
               ryuion(j,1) = ryu(oldchain,j)
               rzuion(j,1) = rzu(oldchain,j)
               qquion(j,1) = qqu(oldchain,j)
            end do
            moltion(1) = imolta
            do j = 1,newunit
               rxuion(j,2) = rxut(2,j)
               ryuion(j,2) = ryut(2,j)
               rzuion(j,2) = rzut(2,j)
               qquion(j,2) = qqu(newchain,j)
            end do
            moltion(2) = imoltb

            call recip(box,vrecipn,vrecipo,1)

            delen = vrecipn - vrecipo
            tweight = tweight * exp(-beta*delen)

            vnewald(box) = vnewald(box) + delen
            vnbox(box) = vnbox(box) + delen

! update the reciprocal space terms *
            call recip(box,vdum,vdum,2)

! Position 2 ***
            oldchain = imolb
            newchain = imola
            oldunit = nunit(imoltb)
            newunit = nunit(imolta)

            do j = 1,oldunit
               rxuion(j,1) = rxu(oldchain,j)
               ryuion(j,1) = ryu(oldchain,j)
               rzuion(j,1) = rzu(oldchain,j)
               qquion(j,1) = qqu(oldchain,j)
            end do
            moltion(1) = imoltb
            do j = 1,newunit
               rxuion(j,2) = rxut(1,j)
               ryuion(j,2) = ryut(1,j)
               rzuion(j,2) = rzut(1,j)
               qquion(j,2) = qqu(newchain,j)
            end do
            moltion(2) = imolta

            call recip(box,vrecipn,vrecipo,1)

            delen = vrecipn - vrecipo
            tweight = tweight * exp(-beta*delen)

            vnewald(box) = vnewald(box) + delen
            vnbox(box) = vnbox(box) + delen

         end if
! End Ewald-sum Corrections

!     ----------------------------------------------------------------------

!     *************************
! Ratio Calculation ***
!     *************************

         wnlog = log10( tweight )
         wolog = log10( tweiold )
         wdlog = wnlog - wolog
         if (wdlog.lt.-softcut.and..not.lideal(box)) then
! write(io_output,*) '### underflow in wratio calculation ###'
            call recip(box,vdum,vdum,4)
            return
         end if

         wswat = ( tweight / tweiold )

         if ( random(-1) .le. wswat ) then
! we can now accept !!!!! ***
            bsswat(iparty,ipairb) = bsswat(iparty,ipairb) + 1.0E0_dp
! write(io_output,*) 'SWATCH ACCEPTED',imola,imolb

            vbox(1,box)     = vbox(1,box)    + vnbox(box)
            vbox(2,box)  = vbox(2,box) + vninte(box)
            vbox(4,box)  = vbox(4,box) + vnintr(box)
            vbox(5,box)    = vbox(5,box)   + vnvibb(box)
            vbox(7,box)     = vbox(7,box)    + vntgb(box)
            vbox(9,box)    = vbox(9,box)   + vnextb(box)
            vbox(6,box)   = vbox(6,box)  + vnbend(box)
            vbox(8,box)  = vbox(8,box) + vnelect(box) + vnewald(box)

! update book keeping ***
! assign new geometries ***
            do ic = 1,iunita
               rxu(imola,ic) = rxut(1,ic)
               ryu(imola,ic) = ryut(1,ic)
               rzu(imola,ic) = rzut(1,ic)
! write(io_output,*) 'imola:',imola
! write(io_output,*) rxu(imola,ic),ryu(imola,ic),rzu(imola,ic)
            end do

            do ic = 1,iunitb
               rxu(imolb,ic) = rxut(2,ic)
               ryu(imolb,ic) = ryut(2,ic)
               rzu(imolb,ic) = rzut(2,ic)
! write(io_output,*) 'imolb:',imolb
! write(io_output,*) rxu(imolb,ic),ryu(imolb,ic),rzu(imolb,ic)
            end do
            if (lewald.and..not.lideal(box)) then
! update reciprocal-space sum
               call recip(box,vdum,vdum,2)
            end if

! update center of mass
            call ctrmas(.false.,box,imolb,8)
            call ctrmas(.false.,box,imola,8)

            if (licell .and. (box .eq. boxlink))  call err_exit(__FILE__,__LINE__,'not yet implemented!',myid+1)

! call nearest neighbor list
            if ( lneigh ) then
               call update_neighbor_list(imola,0.,0.,0.,.true.)
               call update_neighbor_list(imolb,0.,0.,0.,.true.)
            end if
         else if (lewald.and..not.lideal(box)) then
! recover the reciprocal space vectors
! if the move is not accepted ***
            call recip(box,vdum,vdum,4)
         end if

! ****************************************************************
! END INTRABOX SWATCH ADD *
! ****************************************************************

      else

! check if the particle types are in their boxes
         if ( (ncmt(boxa,imolta) .eq. 0) .or.  (ncmt(boxb,imoltb) .eq. 0) ) then
            lempty = .true.
! write(io_output,*) 'one box out of swatch particle'
         else
! get particle from box a

            iboxal = int( real(ncmt(boxa,imolta),dp)*random(-1) ) + 1
            iboxa = parbox(iboxal,boxa,imolta)
            if ( moltyp(iboxa) .ne. imolta ) write(io_output,*) 'screwup'
            iboxia = nboxi(iboxa)
            if (iboxia .ne. boxa) call err_exit(__FILE__,__LINE__,'problem in swatch',myid+1)

! get particle from box b

            iboxbl = int( real(ncmt(boxb,imoltb),dp)*random(-1) ) + 1
            iboxb = parbox(iboxbl,boxb,imoltb)
            if ( moltyp(iboxb) .ne. imoltb ) write(io_output,*) 'screwup'
            iboxib = nboxi(iboxb)
            if (iboxib .ne. boxb) call err_exit(__FILE__,__LINE__,'problem in swatch',myid+1)

!cc--!!!JLR - for test write coordinates of a and b
! open(91,file='a_init.xyz',status='unknown')
! write(91,*) nunit(imolta)
! write(91,*)
! do izz = 1,nunit(imolta)
! write(91,*) 'C ', rxu(iboxa,izz),ryu(iboxa,izz),
!     &              rzu(iboxa,izz)
! end do
! close(91)
! open(92,file='b_init.xyz',status='unknown')
! write(92,*) nunit(imoltb)
! write(92,*)
! do izz = 1,nunit(imoltb)
! write(92,*) 'O ', rxu(iboxb,izz),ryu(iboxb,izz),
!     &              rzu(iboxb,izz)
! end do
! close(92)
!cc--!!!JLR - end of coordinate test

         end if

! add one attempt to the count for iparty
! write(io_output,*) 'iparty:',iparty,'boxa:',boxa
         bnswat(iparty,ipairb) = bnswat(iparty,ipairb) + 1.0E0_dp

! JLR 12-1-09, Count the empty attempts
! if (lempty) return
         if (lempty) then
            bnswat_empty(iparty,ipairb) = bnswat_empty(iparty,ipairb) + 1.0E0_dp
            return
         end if
! END JLR 12-1-09 ---

!$$$c     --- assign from and prev for each moltyp
!$$$         from(1) = gswatc(iparty,type_a,1)
!$$$         prev(1) = gswatc(iparty,type_a,2)
!$$$         from(2) = gswatc(iparty,type_b,1)
!$$$         prev(2) = gswatc(iparty,type_b,2)

         do izz = 1,ncut(iparty,type_a)

            from(type_a+2*(izz-1)) = gswatc(iparty,type_a,1+2*(izz-1))
            prev(type_a+2*(izz-1)) = gswatc(iparty,type_a,2+2*(izz-1))

         end do

         do izz = 1,ncut(iparty,type_b)

            from(type_b+2*(izz-1)) = gswatc(iparty,type_b,1+2*(izz-1))
            prev(type_b+2*(izz-1)) = gswatc(iparty,type_b,2+2*(izz-1))

         end do
! store number of units in iunita and iunitb
         iunita = nunit(imolta)
         iunitb = nunit(imoltb)

! store number of each type in the boxes
         orgaia = ncmt(boxa, imolta)
         orgbia = ncmt(boxa, imoltb)
         orgaib = ncmt(boxb, imolta)
         orgbib = ncmt(boxb, imoltb)
         tweight = 1.0E0_dp
         tweiold = 1.0E0_dp

! set the trial energies to zero
         vnbox(boxa)   = 0.0E0_dp
         vninte(boxa)  = 0.0E0_dp
         vnintr(boxa)  = 0.0E0_dp
         vnvibb(boxa)  = 0.0E0_dp
         vntgb(boxa)   = 0.0E0_dp
         vnextb(boxa)  = 0.0E0_dp
         vnbend(boxa)  = 0.0E0_dp
         vnelect(boxa) = 0.0E0_dp
         vnewald(boxa) = 0.0E0_dp

         vnbox(boxb)   = 0.0E0_dp
         vninte(boxb)  = 0.0E0_dp
         vnintr(boxb)  = 0.0E0_dp
         vnvibb(boxb)  = 0.0E0_dp
         vntgb(boxb)   = 0.0E0_dp
         vnextb(boxb)  = 0.0E0_dp
         vnbend(boxb)  = 0.0E0_dp
         vnelect(boxb) = 0.0E0_dp
         vnewald(boxb) = 0.0E0_dp

! compute the rosenbluth weights for each molecule type in each box
         do ic = 1,2
            if ( ic .eq. 1 ) then
               self = iboxa
               other = iboxb
! s_type = imolta
! o_type = imoltb
               s_type = type_a
               o_type = type_b
               iboxnew = iboxib
               iboxold = iboxia
               imolty = imolta
               igrow = nugrow(imolta)
               iunit = nunit(imolta)
            else
               self = iboxb
               other = iboxa
! s_type = imoltb
! o_type = imolta
               s_type = type_b
               o_type = type_a
               iboxnew = iboxia
               iboxold = iboxib
               imolty = imoltb
               igrow = nugrow(imoltb)
               iunit = nunit(imoltb)
            end if

! store the beads that are identical
            do icbu = 1, nsampos(iparty)
               if ( ic .eq. 1 ) then
! new conformation is for type a
                  new = splist(iparty,icbu,type_a)
                  old = splist(iparty,icbu,type_b)
               else
! new conformation is for type b
                  new = splist(iparty,icbu,type_b)
                  old = splist(iparty,icbu,type_a)
               end if

               rxnew(new) = rxu(other,old)
               rynew(new) = ryu(other,old)
               rznew(new) = rzu(other,old)
            end do

! set up growth schedule
            ifirst = from(ic)
            iprev = prev(ic)

! rigid molecule add on

            if (lrigid(imolty)) then
!cc--- JLR 11-24-09
!cc---adding two if statements for rigid swaps
               if (nsampos(iparty).lt.iunit) then

                  if ( (rindex(imolty).eq.0) .or. (ifirst.lt.riutry(imolty,1))  ) then

                     if (nsampos(iparty).ge.3) then

                        call align_planes(iparty,self,other, s_type,o_type,rxnew,rynew,rznew)

                     else

! calculate new vector from initial bead
                        do j = 1,iunit

                           rxnew(j) = rxnew(ifirst) - (rxu(self,ifirst) - rxu(self,j))
                           rynew(j) = rynew(ifirst) - (ryu(self,ifirst) - ryu(self,j))
                           rznew(j) = rznew(ifirst) - (rzu(self,ifirst) - rzu(self,j))
                           end do

                     end if      !nsampos.lt.3
                  end if         !ifirst.lt.riutry...
               end if            !nsampos.lt.iunit
               call schedule(igrow,imolty,islen,ifirst,iprev,4)
!cc---END JLR 11-24-09
            else
               call schedule(igrow,imolty,islen,ifirst,iprev,3)
            end if !lrigid

! Adding in multiple end regrowths

            if (ncut(iparty,s_type) .gt. 1) then
               do izz = 2,ncut(iparty,s_type)
                  ifirst = from(s_type+2*(izz-1))
                  iprev = prev(s_type+2*(izz-1))

                  call schedule(igrow,imolty,islen,ifirst,iprev,5)

               end do
            end if

            waddnew = 1.0E0_dp
! JLR 11-24-09 New stuff for rigid swatch
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there
                  !don't regrow anything in rosenbluth
                  icallrose = 4
               else if( (nsampos(iparty).ge.3)      .or. (ifirst.ge.riutry(imolty,1)) ) then
! rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
! rigid part is not grown, do rigrot
                  icallrose = 2
               end if
            else
! flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            end if


            call rosenbluth(.true.,lterm,other,self,imolty,islen, iboxnew,igrow,waddnew,.false.,vdum2,icallrose)
! END JLR 11-24-09
! termination of cbmc attempt due to walk termination ---
            if ( lterm ) return

! propagate new rosenbluth weight
            tweight = tweight*weight*waddnew

! end rigid add on

! save the new coordinates
! write(io_output,*) 'new:',moltyp(self),iunit
            do jj = 1,igrow
               rxut(ic,jj) = rxnew(jj)
               ryut(ic,jj) = rynew(jj)
               rzut(ic,jj) = rznew(jj)
! write(io_output,*) rxnew(jj),rynew(jj),rznew(jj)
            end do
! write(io_output,*) 'old:',moltyp(other),nunit(moltyp(other))
! do jj = 1,nunit(moltyp(other))
! write(io_output,*) rxu(other,jj),ryu(other,jj),rzu(other,jj)
! end do

! Corrections for switched beads, and DC-CBMC
! Assign all of the grown new and old beads to rxuion
! with rxuion: new = 2
            iii = 2
            do j=1,igrow
               rxuion(j,iii) = rxnew(j)
               ryuion(j,iii) = rynew(j)
               rzuion(j,iii) = rznew(j)
               qquion(j,iii) = qqu(self,j)
            end do

            if ( iunit .ne. igrow ) then

! for explicit-hydrogen model, put on the hydrogens
! use phony number iins and call explct to add constrained hydrogens

               iins = nchain + 1
               moltyp(iins) = imolty
               do j=1,igrow
                  rxu(iins,j) = rxnew(j)
                  ryu(iins,j) = rynew(j)
                  rzu(iins,j) = rznew(j)
               end do
               call explct(iins,vdum,.false.,.false.)
               do j = igrow + 1, iunit
                  rxuion(j,iii) = rxu(iins,j)
                  ryuion(j,iii) = ryu(iins,j)
                  rzuion(j,iii) = rzu(iins,j)
                  qquion(j,iii) = qqu(self,j)
                  rxut(ic,j) = rxu(iins,j)
                  ryut(ic,j) = ryu(iins,j)
                  rzut(ic,j) = rzu(iins,j)
! write(io_output,*) 'new:',rxut(ic,j),ryut(ic,j),rzut(ic,j)
               end do
            end if

! Begin DC-CBMC, explicit-hydrogen and
! switched bead Corrections for NEW configuration
! calculate the true site-site energy

!> \bug energy problem for cases which involve the change of the bending type
!> and torsional type for those units swatched!!!!!!

            call ctrmas(.false.,iboxnew,other,8)
            call energy(other,imolty,v,iii,iboxnew,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)

            if (lterm) then
! write(io_output,*) 'other ',other,' self ',self
! call err_exit(__FILE__,__LINE__,'interesting screwup in CBMC swatch',myid+1)
               return
            end if
            delen = v(1) - ( vnew(1) - (vnew(5) + vnew(6) + vnew(7) ))
            tweight = tweight*exp(-beta*delen)

            vnew(1)     = vnew(1) + delen
            vnew(2) = v(2)
            vnew(4) = v(4)
            vnew(9)   = v(9)
            vnew(8) = v(8)
            vnew(14) = v(14)

! End DC-CBMC and switched bead Corrections for NEW configuration

! save the trial energies
            vnbox(iboxnew)   = vnbox(iboxnew) + vnew(1)
            vninte(iboxnew)  = vninte(iboxnew) + vnew(2)
            vnintr(iboxnew)  = vnintr(iboxnew) + vnew(4)
            vnvibb(iboxnew)  = vnvibb(iboxnew) + vnew(5)
            vntgb(iboxnew)   = vntgb(iboxnew) + vnew(7)
            vnextb(iboxnew)  = vnextb(iboxnew) + vnew(9)
            vnbend(iboxnew)  = vnbend(iboxnew) + vnew(6)
            vnelect(iboxnew) = vnelect(iboxnew) + vnew(8)
            vnewald(iboxnew) = vnewald(iboxnew) + vnew(14)
! rigid add on

            waddold = 1.0E0_dp

! grow old chain conformation
! JLR 11-24-09 New stuff for rigid swatch
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there
                  !don't regrow anything in rosenbluth
                  icallrose = 4
               else if( (nsampos(iparty).ge.3)      .or. (ifirst.ge.riutry(imolty,1) )) then
! rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
! rigid part is not grown, do rigrot
                  icallrose = 2
               end if
            else
! flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            end if

            call rosenbluth(.false.,lterm,self,self,imolty, islen,iboxold,igrow,waddold,.false.,vdum2,icallrose)
! END JLR 11-24-09

! termination of old walk due to problems generating orientations
            if ( lterm ) then
               write(io_output,*) 'SWATCH: old growth rejected'
               return
            end if

! propagate old rosenbluth weight
            tweiold = tweiold*weiold*waddold

! store the old grown beads and explict placed beads positions
! 1 = old conformation
            iii = 1
            do j = 1,iunit
               rxuion(j,1) = rxu(self,j)
               ryuion(j,1) = ryu(self,j)
               rzuion(j,1) = rzu(self,j)
               qquion(j,1) = qqu(self,j)
            end do

! Begin Correction for DC-CBMC and switched beads for OLD configuration
! correct the acceptance rules
! calculate the Full rcut site-site energy

            call energy(self,imolty,v,iii,iboxold,1,iunit,.true.,lterm,.false.,.false.,.false.,.false.)

            if (lterm) call err_exit(__FILE__,__LINE__,'disaster ovrlap in old conf SWATCH',myid+1)
            deleo = v(1) - ( vold(1) - (vold(5) + vold(6) + vold(7)) )

            tweiold = tweiold*exp(-beta*deleo)

            vold(1)     = vold(1) + deleo
            vold(4) = v(4)
            vold(2) = v(2)
            vold(9)   = v(9)
            vold(8) = v(8)
            vold(14) = v(14)

! End Correction for DC-CBMC and switched beads for OLD configuration

! save the trial energies
            vnbox(iboxold)   = vnbox(iboxold)  - vold(1)
            vninte(iboxold)  = vninte(iboxold) - vold(2)
            vnintr(iboxold)  = vnintr(iboxold) - vold(4)
            vnvibb(iboxold)  = vnvibb(iboxold) - vold(5)
            vntgb(iboxold)   = vntgb(iboxold)  - vold(7)
            vnextb(iboxold)  = vnextb(iboxold) - vold(9)
            vnbend(iboxold)  = vnbend(iboxold) - vold(6)
            vnelect(iboxold) = vnelect(iboxold) - vold(8)
            vnewald(iboxold) = vnewald(iboxold) - vold(14)

         end do

! Perform the Ewald sum reciprical space corrections
         if ( lewald ) then
! added into tweight even though it really contains new-old
! Box A
            if (.not.lideal(iboxia)) then
               oldchain = iboxa
               newchain = iboxb
               oldunit = nunit(imolta)
               newunit = nunit(imoltb)

               do j = 1,oldunit
                  rxuion(j,1) = rxu(oldchain,j)
                  ryuion(j,1) = ryu(oldchain,j)
                  rzuion(j,1) = rzu(oldchain,j)
                  qquion(j,1) = qqu(oldchain,j)
               end do
               moltion(1) = imolta
               do j = 1,newunit
                  rxuion(j,2) = rxut(2,j)
                  ryuion(j,2) = ryut(2,j)
                  rzuion(j,2) = rzut(2,j)
                  qquion(j,2) = qqu(newchain,j)
               end do
               moltion(2) = imoltb

               call recip(iboxia,vrecipn,vrecipo,1)

               delen = vrecipn - vrecipo

               tweight = tweight * exp(-beta*delen)

               vnewald(iboxia) = vnewald(iboxia) + delen
               vnbox(iboxia) = vnbox(iboxia) + delen
            end if
            if (.not.lideal(iboxib)) then
! Box B
               oldchain = iboxb
               newchain = iboxa
               oldunit = nunit(imoltb)
               newunit = nunit(imolta)

               do j = 1,oldunit
                  rxuion(j,1) = rxu(oldchain,j)
                  ryuion(j,1) = ryu(oldchain,j)
                  rzuion(j,1) = rzu(oldchain,j)
                  qquion(j,1) = qqu(oldchain,j)
               end do
               moltion(1) = imoltb
               do j = 1,newunit
                  rxuion(j,2) = rxut(1,j)
                  ryuion(j,2) = ryut(1,j)
                  rzuion(j,2) = rzut(1,j)
                  qquion(j,2) = qqu(newchain,j)
               end do
               moltion(2) = imolta

               call recip(iboxib,vrecipn,vrecipo,1)

               delen = vrecipn - vrecipo
               tweight = tweight * exp(-beta*delen)

               vnewald(iboxib) = vnewald(iboxib) + delen
               vnbox(iboxib) = vnbox(iboxib) + delen
            end if
         end if
! End Ewald-sum Corrections

!     ----------------------------------------------------------------------

         if (ltailc) then
! add tail corrections
            if (lpbcz) then
               if (lsolid(boxa) .and. .not. lrect(boxa)) then
                  vola = (hmat(boxa,1) * (hmat(boxa,5) * hmat(boxa,9) - hmat(boxa,8)*hmat(boxa,6))+ hmat(boxa,4)*(hmat(boxa,8) * hmat(boxa,3)-hmat(boxa,2)* hmat(boxa,9))+hmat(boxa,7) * (hmat(boxa,2)*hmat(boxa,6)- hmat(boxa,5)*hmat(boxa,3)))
               else
                  vola=boxlx(boxa)*boxly(boxa)*boxlz(boxa)
               end if

               if (lsolid(boxb) .and. .not. lrect(boxb)) then
                  volb = (hmat(boxb,1) * (hmat(boxb,5) * hmat(boxb,9) - hmat(boxb,8)*hmat(boxb,6))+ hmat(boxb,4)*(hmat(boxb,8) * hmat(boxb,3)-hmat(boxb,2)* hmat(boxb,9))+hmat(boxb,7) * (hmat(boxb,2)*hmat(boxb,6)- hmat(boxb,5)*hmat(boxb,3)))
               else
                  volb=boxlx(boxb)*boxly(boxb)*boxlz(boxb)
               end if
            else
               vola=boxlx(boxa)*boxly(boxa)
               volb=boxlx(boxb)*boxly(boxb)
            end if

! for new BOXINS with inserted particle
            do mm = 1,2
               dinsta = 0.0E0_dp
               if ( mm .eq. 1 ) then
                  ibox   = boxa
                  imolin = imoltb
                  dvol   = vola
                  imolrm = imolta
               else
                  ibox   = boxb
                  imolin = imolta
                  dvol   = volb
                  imolrm = imoltb
               end if
! JLR 11-24-09 don't do tail corrections for ideal box
               if (.not.lideal(ibox)) then
! new logic for tail correction (same answer) MGM 3-25-98
                  do jmt = 1, nmolty
                     rho = real( ncmt(ibox,jmt) ,dp)
                     if ( jmt .eq. imolin ) rho = rho + 1.0E0_dp
                     if ( jmt .eq. imolrm ) rho = rho - 1.0E0_dp
                     rho = rho / dvol
                     do imt = 1, nmolty
                        dicount = ncmt(ibox,imt)
                        if ( imt .eq. imolin ) dicount = dicount + 1
                        if ( imt .eq. imolrm ) dicount = dicount - 1
                        dinsta = dinsta +  dicount * coru(imt,jmt,rho,ibox)
                     end do
                  end do
!$$$                  if (ibox .eq. 1 .and. lexzeo) then
!$$$                     do jmt = 1,zntype
!$$$                        rho = znum(jmt)/dvol
!$$$                        do imt = 1, nmolty
!$$$                           dicount=ncmt(ibox,imt)
!$$$                           if ( imt .eq. imolin ) dicount=dicount+1
!$$$                           if ( imt .eq. imolrm ) dicount=dicount-1
!$$$                           dinsta=dinsta+dicount*coruz(imt,jmt,rho,ibox)
!$$$                        end do
!$$$                     end do
!$$$                  end if
               else
                  dinsta = 0.0E0_dp
               end if
! END JLR 11-24-09
               dinsta = dinsta - vbox(3,ibox)

               tweight=tweight*exp(-beta*dinsta)

               vntail(ibox) = dinsta
            end do
         else
            vntail(boxa) = 0.0E0_dp
            vntail(boxb) = 0.0E0_dp
         end if

         wnlog = log10( tweight )
         wolog = log10( tweiold )
         wdlog = wnlog - wolog
         if ( wdlog .lt. -softcut ) then
! write(io_output,*) '### underflow in wratio calculation ###'
            return
         end if
         wswat = ( tweight / tweiold ) * ( real(orgaia*orgbib,dp) / real((orgbia+1)*(orgaib+1),dp) ) * exp(beta*(eta2(boxa,imolta) +eta2(boxb,imoltb)-eta2(boxa,imoltb) -eta2(boxb,imolta)))

         if (lopt_bias(imolta)) call update_bias(log(wswat*2.0)/beta/2.0,boxa,boxb,imolta)
         if (lopt_bias(imoltb)) call update_bias(log(wswat*2.0)/beta/2.0,boxb,boxa,imoltb)

! write(io_output,*) 'imolta,imoltb',imolta,imoltb
! write(io_output,*) 'wswat,tweight,tweiold',wswat,tweight,tweiold
         if ( random(-1) .le. wswat ) then
! we can now accept !!!!! ***
            bsswat(iparty,ipairb) = bsswat(iparty,ipairb) + 1.0E0_dp
! write(io_output,*) 'SWATCH ACCEPTED',iboxa,iboxb
            do jj = 1,2
               if ( jj .eq. 1 ) ic = boxa
               if ( jj .eq. 2 ) ic = boxb
               vbox(1,ic)     = vbox(1,ic)    + vnbox(ic) + vntail(ic)
               vbox(2,ic)  = vbox(2,ic) + vninte(ic) + vntail(ic)
               vbox(4,ic)  = vbox(4,ic) + vnintr(ic)
               vbox(5,ic)    = vbox(5,ic)   + vnvibb(ic)
               vbox(7,ic)     = vbox(7,ic)    + vntgb(ic)
               vbox(9,ic)    = vbox(9,ic)   + vnextb(ic)
               vbox(6,ic)   = vbox(6,ic)  + vnbend(ic)
               vbox(3,ic)   = vbox(3,ic)  + vntail(ic)
               vbox(8,ic)  = vbox(8,ic) + vnelect(ic) + vnewald(ic)
            end do

! update book keeping
            nboxi(iboxa) = boxb
            nboxi(iboxb) = boxa

            parbox(orgbia+1,boxa,imoltb)= iboxb
            parbox(orgaib+1,boxb,imolta)= iboxa
            parbox(iboxal,boxa,imolta) = parbox(orgaia,boxa,imolta)
            parbox(iboxbl,boxb,imoltb) = parbox(orgbib,boxb,imoltb)
            parbox(orgaia,boxa,imolta) = 0
            parbox(orgbib,boxb,imoltb) = 0

            ncmt(boxa,imolta) = orgaia - 1
            ncmt(boxa,imoltb) = orgbia + 1
            ncmt(boxb,imolta) = orgaib + 1
            ncmt(boxb,imoltb) = orgbib - 1

            do ic = 1,iunita
               rxu(iboxa,ic) = rxut(1,ic)
               ryu(iboxa,ic) = ryut(1,ic)
               rzu(iboxa,ic) = rzut(1,ic)
            end do
            do ic = 1,iunitb
               rxu(iboxb,ic) = rxut(2,ic)
               ryu(iboxb,ic) = ryut(2,ic)
               rzu(iboxb,ic) = rzut(2,ic)
            end do
            if ( lewald ) then
! update reciprocal-space sum
               if (.not.lideal(boxa)) call recip(boxa,vdum,vdum,2)
               if (.not.lideal(boxb)) call recip(boxb,vdum,vdum,2)
            end if

! update center of mass
            call ctrmas(.false.,boxa,iboxb,8)
            call ctrmas(.false.,boxb,iboxa,8)

            if (licell .and. (boxa .eq. boxlink .or. boxb .eq. boxlink)) call err_exit(__FILE__,__LINE__,'not yet implemented!',myid+1)

! call nearest neighbor list
            if ( lneigh ) then
               call update_neighbor_list(iboxa,0.,0.,0.,.true.)
               call update_neighbor_list(iboxb,0.,0.,0.,.true.)
            end if

!cc--!!!JLR - for test 2 write final coordinates of a and b
! open(93,file='a_final.xyz',status='unknown')
! write(93,*) iunita
! write(93,*)
! do izz = 1,iunita
! write(93,*) 'C ', rxu(iboxa,izz),ryu(iboxa,izz),rzu(iboxa,izz)
! end do
!
! open(93,file='b_final.xyz',status='unknown')
! write(93,*) iunitb
! write(93,*)
! do izz = 1,iunitb
! write(93,*) 'C ', rxu(iboxb,izz),ryu(iboxb,izz),rzu(iboxb,izz)
! end do
! call err_exit(__FILE__,__LINE__,'END OF SWATCH TEST',myid+1)
!cc--!!!JLR - end of coordinate test

         end if
! -----------------------------------------------------------------

      end if

#ifdef __DEBUG__
      write(io_output,*) 'end SWATCH in ',myid
#endif

      return
  end subroutine swatch

  subroutine init_swatch(io_input,lprint)
    use var_type,only:default_string_length
    use util_string,only:uppercase
    use util_files,only:readLine
    integer,intent(in)::io_input
    LOGICAL,INTENT(IN)::lprint
    character(LEN=default_string_length)::line_in
    integer::jerr,i,j,k
    namelist /mc_swatch/ pmswat,nswaty,pmsatc

    allocate(bnswat(npamax,npabmax),bnswat_empty(npamax,npabmax),bsswat(npamax,npabmax),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_swatch: allocation failed',jerr)

    bnswat = 0.0E0_dp
    bsswat = 0.0E0_dp
    bnswat_empty = 0.0E0_dp

    !> read namelist mc_swatch
    nswaty=nmolty*(nmolty-1)/2
    do i=1,nswaty
       pmsatc(i)=real(i,dp)/nswaty
    end do

    rewind(io_input)
    read(UNIT=io_input,NML=mc_swatch,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) call err_exit(__FILE__,__LINE__,'reading namelist: mc_swatch',jerr)

    if (nswaty.gt.npamax) call err_exit(__FILE__,__LINE__,'nswaty gt npamax',myid+1)

    if (lprint) then
       write(io_output,'(/,A,/,A)') 'NAMELIST MC_SWATCH','------------------------------------------'
       write(io_output,'(A,G16.9)') 'pmswat: ',pmswat
       write(io_output,'(A,I0)') '   number of swatch pairs (nswaty): ',nswaty
       do i=1,nswaty
          write(io_output,'(A,G16.9)') '   probability of each swatch pair: ',pmsatc(i)
       end do
    end if

    ! Looking for section MC_SWATCH
    REWIND(io_input)
    CYCLE_READ_SWATCH:DO
       call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Section MC_SWATCH not found',jerr)

       if (UPPERCASE(line_in(1:9)).eq.'MC_SWATCH') then
          do i=1,nswaty+1
             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
             if (UPPERCASE(line_in(1:13)).eq.'END MC_SWATCH') then
                if (i.ne.nswaty+1) call err_exit(__FILE__,__LINE__,'Section MC_SWATCH not complete!',myid+1)
                exit
             else if (i.eq.nswaty+1) then
                call err_exit(__FILE__,__LINE__,'Section MC_SWATCH has more than nswaty records!',myid+1)
             end if

             ! moltyp1<->moltyp2 nsampos 2xncut
             read(line_in,*) nswatb(i,1:2),nsampos(i),ncut(i,1:2)
             if (nswatb(i,1).eq.nswatb(i,2)) then
                ! safety checks on swatch
                write(io_output,*) 'nswaty ',i,' has identical moltyp'
                call err_exit(__FILE__,__LINE__,'cannot swatch identical moltyp',myid+1)
             end if

             if (lprint) then
                write(io_output,'(/,A,2(4X,I0))') '   swatch molecule type pairs:',nswatb(i,1:2)
                write(io_output,'(A,I0,A,2(2X,I0))') '   nsampos: ',nsampos(i),', ncut:',ncut(i,1:2)
             end if

             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
             ! gswatc 2x(ifrom, iprev)
             read(line_in,*) (gswatc(i,j,1:2*ncut(i,j)),j=1,2)

             if (lprint) then
                do j=1,2
                   write(io_output,FMT='(A,I0)') '   molecule ',j
                   do k = 1,ncut(i,j)
                      write(io_output,'(3(A,I0))') '   ncut ',k,': grom from ',gswatc(i,j,2*k-1),', prev ',gswatc(i,j,2*k-1)
                   end do
                end do
             end if

             do j = 1,nsampos(i)
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                ! splist
                read(line_in,*) splist(i,j,1:2)
                if (lprint) then
                   write(io_output,'(A,2(4X,I0))') '   splist:',splist(i,j,1:2)
                end if
             end do

             call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
             ! nswtcb pmswtcb
             read(line_in,*) nswtcb(i),pmswtcb(i,1:nswtcb(i))
             if (lprint) then
                write(io_output,'(A,I0)') '   number of swatch box pairs: ',nswtcb(i)
                do j=1,nswtcb(i)
                   write(io_output,'(A,G16.9)') '   probability of the swatch box pair: ',pmswtcb(i,j)
                end do
             end if

             do j = 1,nswtcb(i)
                call readLine(io_input,line_in,skipComment=.true.,iostat=jerr)
                if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'Reading section MC_SWATCH',jerr)
                ! box numbers
                read(line_in,*) box3(i,j),box4(i,j)
                if (lprint) then
                   write(io_output,'(A,2(4X,I0))') '   box pair:',box3(i,j),box4(i,j)
                end if
             end do
          end do
          exit cycle_read_swatch
       end if
    END DO CYCLE_READ_SWATCH
  end subroutine init_swatch

  subroutine output_swatch_stats(io_output)
    integer,intent(in)::io_output
    integer::i,j

    write(io_output,*)
    write(io_output,*) '### Molecule swatch     ###'
    write(io_output,*)
    do i = 1, nswaty
       write(io_output,*) 'pair typ =',i
       write(io_output,*) 'moltyps = ',nswatb(i,1),' and',nswatb(i,2)
       do j = 1, nswtcb(i)
          ! JLR 12-1-09 changing to exclude empty box attempts from swatch rate
          write(io_output,"('between box ',i2,' and ',i2, '   uattempts =',f12.1,   '  attempts =',f9.1, '  accepted =',f8.1)") box3(i,j),box4(i,j), bnswat(i,j),bnswat(i,j)-bnswat_empty(i,j),bsswat(i,j)
          if (bnswat(i,j) .gt. 0.5E0_dp ) then
             write(io_output,"(' accepted % =',f7.3)") 100.0E0_dp * bsswat(i,j)/ (bnswat(i,j)-bnswat_empty(i,j))
          end if
          ! EN JLR 12-1-09
       end do
    end do
  end subroutine output_swatch_stats

  subroutine read_checkpoint_swatch(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) bnswat,bnswat_empty,bsswat
    call mp_bcast(bnswat,npamax*npabmax,rootid,groupid)
    call mp_bcast(bnswat_empty,npamax*npabmax,rootid,groupid)
    call mp_bcast(bsswat,npamax*npabmax,rootid,groupid)
  end subroutine read_checkpoint_swatch

  subroutine write_checkpoint_swatch(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) bnswat,bnswat_empty,bsswat
  end subroutine write_checkpoint_swatch
end module transfer_swatch
