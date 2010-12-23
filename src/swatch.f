      subroutine swatch 

!      **********************************************************
!    *** Added intrabox move for two particles within one box   ***
!    *** in combined move that shares the same parameters.      ***
!    *** Will also accept rigid (lrigid) molecules.  Contains   ***
!    *** several critical bug fixes as well.                    ***
!    ***                                     9-25-02 JMS        ***
!      **********************************************************
 
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
 
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'system.inc'
!$$$      include 'external.inc'
!$$$      include 'zeopoten.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'rosen.inc'
!$$$      include 'swtcmove.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'cell.inc'

      logical::lempty,lterm

      integer(KIND=normal_int)::type_a,type_b,from,prev,self,iboxnew
     & ,iboxold,imolty,igrow,new,old,islen,ifirst,iprev,iii,j
      integer(KIND=normal_int)::oldchain,newchain,oldunit,newunit,ncount
     & ,iunit,iins

      integer(KIND=normal_int)::ic,ibox,icbu,jj,mm,imt,jmt,imolin,imolrm
      integer(KIND=normal_int)::boxa,boxb,ipair,imolta,imoltb,iboxa
     & ,iboxb,iboxal,iboxbl,iboxia,iboxib,iunita,iunitb,orgaia ,orgaib
     & ,orgbia,orgbib,ipairb 

      real(KIND=double_precision)::random,tweight,tweiold,rxut,ryut,rzut
     & ,dvol,vola,volb,rho,coru,coruz,dinsta,rpair

      real(KIND=double_precision)::vnbox,vninte,vnintr,vnvibb,vntgb
     & ,vnextb,vnbend,vntail,vnelect,vnewald,wnlog,wolog,wdlog,wswat
         
      real(KIND=double_precision)::v,vintra,vinter,vext,velect,vewald
     & ,vdum,delen,deleo,dicount,vrecipn,vrecipo
! * additions from iswatch

      integer(KIND=normal_int)::izz,zzz,box,iboxi,bdmol_a,bdmol_b
      integer(KIND=normal_int)::imola,imolb,moltaid,moltbid,i,ii,iu

      integer(KIND=normal_int)::s_type, o_type, thisbox, otherbox

      real(KIND=double_precision)::rx_1(numax),ry_1(numax),rz_1(numax)
     & ,dummy

!      dimension from(2),prev(2)
      dimension from(2*numax),prev(2*numax)

!      dimension rxut(2,numax),ryut(2,numax),rzut(2,numax)
      dimension rxut(4,numax),ryut(4,numax),rzut(4,numax)

! * end additions

      real(KIND=double_precision)::waddold,waddnew,vdum2

      dimension vnbox(nbxmax),vninte(nbxmax),vnintr(nbxmax)
     & ,vnvibb(nbxmax),vntgb(nbxmax),vnextb(nbxmax),vnbend(nbxmax)
     & ,vntail(nbxmax),vnelect(nbxmax),vnewald(nbxmax)

! --- JLR 11-24-09
      integer(KIND=normal_int)::icallrose
! --- END JLR 11-24-09

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!      write(iou,*) 'start SWATCH'

      lempty = .false.
!     ---randomly select chains to switch between boxa boxb
!     --- select a pair type to switch 
      if ( nswaty .gt. 1 ) then
         rpair = random()
         do 96 ipair = 1, nswaty
            if ( rpair .lt. pmsatc(ipair) ) then
               iparty = ipair
               goto 97
            end if
 96      continue
      else
         iparty = 1 
      end if

! *** randomly select the molecules from the pair ***

      if (random() .lt. 0.5d0) then
         type_a = 1
         type_b = 2
      else
         type_a = 2
         type_b = 1
      end if

!     ---select the molecules from the pair
 97   imolta = nswatb(iparty,1)
      imoltb = nswatb(iparty,2)
! ???!!!
      type_a = 1
      type_b = 2

!     ---choose box A and box B at random
      if ( nswtcb(iparty) .gt. 1 ) then
         rpair = random()
         do 98 ipair = 1, nswtcb(iparty)
            if ( rpair .lt. pmswtcb(iparty,ipair) ) then
               ipairb = ipair
               goto 99
            end if
 98      continue
      else
         ipairb = 1
      end if
         
 99   if (random() .lt. 0.5d0) then
         boxa=box3(iparty,ipairb)
         boxb=box4(iparty,ipairb)
      else
         boxa=box4(iparty,ipairb)
         boxb=box3(iparty,ipairb)
      end if

      if (boxa .eq. boxb) then

! ****************************************************************
! * INTRABOX SWATCH ADD *
! ****************************************************************
! ****************************************************************
 
! * liswatch prevents non-grown beads from being included in 
! * the new growth in boltz.f

!      write(iou,*) 'start iSWATCH'

!$$$c *** randomly select chains to switch ***
!$$$c *** select a pair type to switch *** 
!$$$      if ( niswaty .gt. 1 ) then
!$$$         rpair = random()
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
!$$$         if (random() .lt. 0.5d0) then
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

! *** box determinations ***

         if ( (ncmt(boxa,imolta) .eq. 0) .or.
     &        (ncmt(boxb,imoltb) .eq. 0) ) then
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

!     *** get particle of type a ***

           moltaid = idint(dble (ncmt(box,imolta))*random()) + 1

!     * imola is the overall molecule id number (runs from 1 to total #) *
           imola = parbox( moltaid, box, imolta) 

           if (moltyp(imola) .ne. imolta) 
     &                   write(iou,*) 'screwup in iswatch'

           iboxi = nboxi(imola)

           if (iboxi .ne. box) call cleanup('problem in iswatch')

!     *** get particle of type b ***

           moltbid = idint(dble (ncmt(box,imoltb))*random()) + 1
           imolb = parbox( moltbid, box, imoltb) 

           if (moltyp(imolb) .ne. imoltb) 
     &                write(iou,*) 'screwup in iswatch'
           iboxi = nboxi(imolb)

           if (iboxi .ne. box) call cleanup('problem in iswatch')
         end if  
!     *** add one attempt to the count for iparty
!$$$  bniswat(iparty,box) = bniswat(iparty,box) + 1.0d0   
         bnswat(iparty,ipairb) = bnswat(iparty,ipairb) + 1.0d0   
! --- JLR 12-1-09 count the empty attempts
!         if (lempty) return 
         if (lempty) then
            bnswat_empty(iparty,ipairb) =
     &           bnswat_empty(iparty,ipairb) + 1.0d0
            return
         end if
! --- END JLR 12-1-09

!     * write out the molecule numbers of the pair *
!     write(iou,*) imola,imolb

!     ***************************
!     *** Begin Growth Setups ***
!     ***************************

!     *** assign from and prev for each moltyp

         do izz = 1,ncut(iparty,type_a)

            from(type_a+2*(izz-1)) = gswatc(iparty,type_a,1+2*(izz-1)) 
            prev(type_a+2*(izz-1)) = gswatc(iparty,type_a,2+2*(izz-1)) 

         end do

         do izz = 1,ncut(iparty,type_b)

            from(type_b+2*(izz-1)) = gswatc(iparty,type_b,1+2*(izz-1)) 
            prev(type_b+2*(izz-1)) = gswatc(iparty,type_b,2+2*(izz-1)) 

         end do

!$$$         write(iou,*) 'mol a'
!$$$  
!$$$         do izz = 1,ncut(type_a)
!$$$            write(iou,*) izz,'from',from(type_a+2*(izz-1)),' prev',
!$$$     &        prev(type_a+2*(izz-1))
!$$$         end do
!$$$  
!$$$         write(iou,*) 'mol b'
!$$$  
!$$$         write(iou,*) gswatc(iparty,type_b,3),gswatc(iparty,type_b,4)
!$$$         write(iou,*) from(3),prev(3),type_b
!$$$  
!$$$         do izz = 1,ncut(type_b)
!$$$            write(iou,*) izz,'from',from(type_b+2*(izz-1)),' prev',
!$$$     &        prev(type_b+2*(izz-1))
!$$$         end do
!$$$         
!$$$         call cleanup('')
      
!     *** store number of units in iunita and iunitb
         iunita = nunit(imolta)
         iunitb = nunit(imoltb)

!     *** initialize trial weights ***
         tweight = 1.0d0
         tweiold = 1.0d0

!     *** set the trial energies to zero ***
         vnbox(box)   = 0.0d0
         vninte(box)  = 0.0d0
         vnintr(box)  = 0.0d0
         vnvibb(box)  = 0.0d0
         vntgb(box)   = 0.0d0
         vnextb(box)  = 0.0d0
         vnbend(box)  = 0.0d0
         vnelect(box) = 0.0d0
         vnewald(box) = 0.0d0

!     *** store position 1, a's original site ***

         do izz = 1, nunit(imolta)

            rx_1(izz) = rxu(imola,izz)
            ry_1(izz) = ryu(imola,izz)
            rz_1(izz) = rzu(imola,izz)

         end do

!     *** store same bead coordinates for molecules a and b ***

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
!     *** Liswinc Determination ***
!     *****************************

!     *** only need this for molecule b, since b is grown after a
!     *** but a is grown when some of b doesn't exist

         self = imolb
         other = imola
         imolty = imoltb
         igrow = nugrow(imoltb)

         ifirst = from(type_b)
         iprev = prev(type_b)

!     * determine which beads aren't in the same positions *

         call schedule(igrow,imolty,islen,ifirst,iprev,3)

         if (ncut(iparty,type_b) .gt. 1) then
            do izz = 2,ncut(iparty,type_b)
               ifirst = from(type_b+2*(izz-1))
               iprev = prev(type_b+2*(izz-1))

               call schedule(igrow,imolty,islen,ifirst,iprev,5)

            end do
         end if

!     * assign growth schedule for molecule b *

         do izz = 1,nunit(imolty)

            liswinc(izz,imolty) = lexshed(izz)

!     write(iou,*) izz, liswinc(izz,imolty)
            
         end do

!     ******************
!     *** ROSENBLUTH ***
!     ******************

!     *** start rosenbluth weight calculations ***
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

!     *** assign same bead coordinates for new growth ***
            do icbu = 1, nsampos(iparty)

               new = splist(iparty,icbu,s_type)
               old = splist(iparty,icbu,o_type)

               rxnew(new) = rxu(other,old)
               rynew(new) = ryu(other,old)
               rznew(new) = rzu(other,old)

            end do

!        --- set up growth schedule
            ifirst = from(s_type)
            iprev = prev(s_type)

!     --- lrigid add on 

            if (lrigid(imolty)) then
!cc--- JLR 11-24-09
!cc--- adding some if statements for rigid swaps:
!cc--- if we have swapped the whole rigid part we will 
!cc--- not compute the vectors from ifirst in old position                                   
               if ( nsampos(iparty).lt.iunit) then

                  if ( (rindex(imolty).eq.0) .or.
     &                 (ifirst.lt.riutry(imolty,1))  ) then

                     if (nsampos(iparty).ge.3) then

                        call align_planes(iparty,self,other,
     &                       s_type,o_type,rxnew,rynew,rznew)

                     else

!     --- calculate new vector from initial bead                                                    
                        do j = 1,iunit
                           rxnew(j) = rxnew(ifirst)
     &                          - (rxu(self,ifirst) - rxu(self,j))
                           rynew(j) = rynew(ifirst)
     &                          - (ryu(self,ifirst) - ryu(self,j))
                           rznew(j) = rznew(ifirst)
     &                          - (rzu(self,ifirst) - rzu(self,j))
                        end do

                     end if      !nsampos.eq.3                                                       
                  end if         !ifirst.lt.riutry...                                                
               end if            !nsampos.lt.iunit                                                   
               call schedule(igrow,imolty,islen,ifirst,iprev,4)
!cc---END JLR 11-24-09
            else
               call schedule(igrow,imolty,islen,ifirst,iprev,3)
            end if

! * I wonder if this works with lrigid?                                                             
            if (ncut(iparty,s_type) .gt. 1) then
               do izz = 2,ncut(iparty,s_type)
                  ifirst = from(s_type+2*(izz-1))
                  iprev = prev(s_type+2*(izz-1))

                  call schedule(igrow,imolty,islen,ifirst,iprev,5)

               end do
            end if



!     *********** Commenting out for new combined code JMS 6-20-00 *****
!     *** Call qqcheck to setup the group based qq cutoff ***
!$$$  if ( lelect(imolty) ) then
!$$$  call qqcheck(other,box,rxnew(1),rynew(1),rznew(1))
!$$$  end if

!     ******************
!     *** new growth ***
!     ******************

!     *** moving molecules for rosenbluth ***
            if (ic .eq. 1) then

!     * putting molecule b in position 1 *

               do izz = 1, nsampos(iparty)

                  bdmol_b = splist(iparty,izz,type_b)

                  rxu(other,bdmol_b) = rxut(3,izz)
                  ryu(other,bdmol_b) = ryut(3,izz)
                  rzu(other,bdmol_b) = rzut(3,izz)

               end do

            else

!     * putting molecule a into its (fully grown) trial position 2 *

               do izz = 1, nunit(moltyp(other))

                  rxu(other,izz) = rxut(1,izz)
                  ryu(other,izz) = ryut(1,izz)
                  rzu(other,izz) = rzut(1,izz)

               end do

            end if

!     --- grow molecules

!     --- changing for lrigid to include waddnew
            waddnew = 1.0d0

!     --- grow new chain conformation
! --- JLR 11-24-09
! --- Different logic/calls to rosenbluth for rigid molecules 
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there                                                      
                  !don't regrow anything in rosenbluth  
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1) )) then
!                 rigid part is grown, don't do rigrot in rosebluth 
                  icallrose = 3
               else
!                 rigid part is not grown, do rigrot  
                  icallrose = 2
               end if
            else
!              flexible molecule call rosenbluth in normal fashion  
               icallrose = 2
            end if

            call rosenbluth(.true.,lterm,self,self,imolty,islen
     &           ,box ,igrow,waddnew,.false.,vdum2,icallrose)
! ---END JLR 11-24-09

!     --- propagate new rosenbluth weight
            tweight = tweight*weight*waddnew

! --- end rigid add on

!     * moving molecules back *
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
!     *** termination of cbmc attempt due to walk termination ***
!     write(iou,*) 'iSWATCH:new growth rejected',ic
! * reset liswatch
               liswatch = .false.
               return
            else
!     write(iou,*) 'iSWATCH:new growth accepted',ic
            end if

!     *** save the new coordinates ***
            do jj = 1,igrow
               rxut(ic,jj) = rxnew(jj)
               ryut(ic,jj) = rynew(jj)
               rzut(ic,jj) = rznew(jj)
            end do

!     --- Corrections for switched beads, and DC-CBMC
!     --- Assign all of the grown new and old beads to rxuion
!     --- with rxuion: new = 2
            iii = 2
            do j=1,igrow
               rxuion(j,iii) = rxnew(j)
               ryuion(j,iii) = rynew(j)
               rzuion(j,iii) = rznew(j)
               qquion(j,iii) = qqu(self,j)
            end do

!     *** added from new-combined code ***
            if ( iunit .ne. igrow ) then

!     --- for explicit-hydrogen model, put on the hydrogens
!     --- use phony number iins and call explct to add constrained hydrogens

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


!     --- Begin DC-CBMC and switched bead Corrections for NEW configuration
!     --- calculate the true site-site energy

            if (ic .eq. 1) then
!     * exclude molecule b from energy calculation, put in other box *
               nboxi(imolb) = otherbox

            else
!     * put molecule a into position 2 (fully grown trial position)
!     * for the energy of b's new position (second time around)
               do izz = 1,nunit(imolta)
                  rxu(imola,izz) = rxut(1,izz)
                  ryu(imola,izz) = ryut(1,izz)
                  rzu(imola,izz) = rzut(1,izz)
               end do
            end if

!     *** get energy of configuration ***

            call ctrmas(.false.,box,self,8)
            call energy (self,imolty, v, vintra,vinter,vext
     &           ,velect,vewald,iii,box,1,iunit,.true.,lterm,.false.
     &           ,vdum,.false.,.false.)

!     *** return to normal ***
            if (ic .eq. 1) then

!     * return b to original box *
               nboxi(imolb) = thisbox
            else

!     * return a to position 1 *
               do izz = 1, nunit(imolta)
                  rxu(imola,izz) = rx_1(izz)
                  ryu(imola,izz) = ry_1(izz)
                  rzu(imola,izz) = rz_1(izz)

               end do

            end if
            
            if (lterm) then
               write(iou,*) 'other ',other,' self ',self,ic
               call cleanup('interesting screwup in CBMC iswatch')
            end if

!     *** add on the changes in energy ***

            delen = v - ( vnewt - (vnewbvib + vnewbb + vnewtg )) 

            tweight = tweight*dexp(-beta*delen)
            vnewt     = vnewt + delen
            vnewinter = vinter 
            vnewintra = vintra
            vnewext   = vext 
            vnewelect = velect
            vnewewald = vewald

!     End DC-CBMC and switched bead Corrections for NEW configuration

!     *** save the trial energies ***
            vnbox(box)   = vnbox(box) + vnewt
            vninte(box)  = vninte(box) + vnewinter
            vnintr(box)  = vnintr(box) + vnewintra
            vnvibb(box)  = vnvibb(box) + vnewbvib
            vntgb(box)   = vntgb(box) + vnewtg
            vnextb(box)  = vnextb(box) + vnewext
            vnbend(box)  = vnbend(box) + vnewbb
            vnelect(box) = vnelect(box) + vnewelect
            vnewald(box) = vnewald(box) + vnewewald

!     ******************
!     *** old growth ***
!     ******************

!     --- lrigid add on

            waddold = 1.0d0

!        --- grow old chain conformation
! --- JLR 11-24-09
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
!                 molecule is all there
!                 don't regrow anything in rosenbluth 
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1)) ) then
!                 rigid part is grown, don't do rigrot in rosebluth 
                  icallrose = 3
               else
!                 rigid part is not grown, do rigrot
                  icallrose = 2
               end if
            else
!              flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            end if

            call rosenbluth(.false.,lterm,self,self,imolty,
     &           islen,box,igrow,waddold,.false.,vdum2,icallrose)
! --- END JLR 11-24-09

            if ( lterm ) then
!     *** termination of old walk due to problems generating orientations ***
!     write(iou,*) 'iSWATCH:old growth rejected',ic
! * reset liswatch
               liswatch = .false.
               return
            else 
!     write(iou,*) 'iSWATCH:old growth accepted',ic
            end if

!        --- propagate old rosenbluth weight
            tweiold = tweiold*weiold*waddold

! --- end rigid add on

!     --- store the old grown beads and explict placed beads positions
!     --- 1 = old conformation
            iii = 1
            do j = 1,iunit
               rxuion(j,1) = rxu(self,j)
               ryuion(j,1) = ryu(self,j)
               rzuion(j,1) = rzu(self,j)
               qquion(j,1) = qqu(self,j)
            end do

!     Begin Correction for DC-CBMC and switched beads for OLD configuration
!     --- correct the acceptance rules 
!     --- calculate the Full rcut site-site energy

!     *** excluding molecule b for first loop ***
            if (ic .eq. 1) then
               nboxi(imolb) = otherbox
            end if

!     * get total energy *

            call energy (self,imolty, v, vintra,vinter,vext,velect
     &           ,vewald,iii,box, 1,iunit,.true.,lterm,.false.,vdum
     &           ,.false.,.false.)

            if (ic .eq. 1) then
!     * return b to current box *
               nboxi(imolb) = thisbox
            end if

            if (lterm) call cleanup('disaster ovrlap in old
     &                          conf iSWATCH')
            deleo = v - ( voldt - (voldbvib + voldbb + voldtg) ) 

            tweiold = tweiold*dexp(-beta*deleo)
            voldt     = voldt + deleo
            voldintra = vintra
            voldinter = vinter 
            voldext   = vext 
            voldelect = velect
            voldewald = vewald

!     End Correction for DC-CBMC and switched beads for OLD configuration

!     *** save the trial energies ***
            vnbox(box)   = vnbox(box)  - voldt
            vninte(box)  = vninte(box) - voldinter
            vnintr(box)  = vnintr(box) - voldintra
            vnvibb(box)  = vnvibb(box) - voldbvib 
            vntgb(box)   = vntgb(box)  - voldtg 
            vnextb(box)  = vnextb(box) - voldext
            vnbend(box)  = vnbend(box) - voldbb 
            vnelect(box) = vnelect(box) - voldelect
            vnewald(box) = vnewald(box) - voldewald

         end do

!     *****************************
!     *** Ewald sum corrections ***
!     *****************************

!     --- Perform the Ewald sum reciprical space corrections
         if ( lewald ) then
!     --- added into tweight even though it really contains new-old

!     *** store the reciprocal space vector
            call recip(box,vdum,vdum,3)

!     *** Position 1 ***
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
            tweight = tweight * dexp(-beta*delen)

            vnewald(box) = vnewald(box) + delen
            vnbox(box) = vnbox(box) + delen

!     * update the reciprocal space terms *
            call recip(box,vdum,vdum,2)

!     *** Position 2 ***
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
            tweight = tweight * dexp(-beta*delen)

            vnewald(box) = vnewald(box) + delen
            vnbox(box) = vnbox(box) + delen

         end if
!     End Ewald-sum Corrections

!     ----------------------------------------------------------------------

!     *************************
!     *** Ratio Calculation ***
!     *************************

         wnlog = dlog10( tweight )
         wolog = dlog10( tweiold )
         wdlog = wnlog - wolog
         if ( wdlog .lt. -softcut ) then
!     write(iou,*) '### underflow in wratio calculation ###'
            call recip(box,vdum,vdum,4)
            return
         end if
         
         wswat = ( tweight / tweiold ) 

         if ( random() .le. wswat ) then
!     *** we can now accept !!!!! ***
            bsswat(iparty,ipairb) = bsswat(iparty,ipairb) + 1.0d0
!     write(iou,*) 'SWATCH ACCEPTED',imola,imolb

            vbox(box)     = vbox(box)    + vnbox(box)
            vinterb(box)  = vinterb(box) + vninte(box)
            vintrab(box)  = vintrab(box) + vnintr(box)
            vvibb(box)    = vvibb(box)   + vnvibb(box)
            vtgb(box)     = vtgb(box)    + vntgb(box)
            vextb(box)    = vextb(box)   + vnextb(box)
            vbendb(box)   = vbendb(box)  + vnbend(box)
            velectb(box)  = velectb(box) + vnelect(box) + vnewald(box)

!     *** update book keeping ***
!     *** assign new geometries ***
            do ic = 1,iunita
               rxu(imola,ic) = rxut(1,ic)
               ryu(imola,ic) = ryut(1,ic)
               rzu(imola,ic) = rzut(1,ic)
!     write(iou,*) 'imola:',imola
!     write(iou,*) rxu(imola,ic),ryu(imola,ic),rzu(imola,ic)
            end do

            do ic = 1,iunitb
               rxu(imolb,ic) = rxut(2,ic)
               ryu(imolb,ic) = ryut(2,ic)
               rzu(imolb,ic) = rzut(2,ic)
!     write(iou,*) 'imolb:',imolb
!     write(iou,*) rxu(imolb,ic),ryu(imolb,ic),rzu(imolb,ic)
            end do
            if ( lewald ) then
!     -- update reciprocal-space sum
               call recip(iboxi,vdum,vdum,2)
            end if

!     *** update center of mass
            call ctrmas(.false.,box,imolb,8)
            call ctrmas(.false.,box,imola,8)

            if (licell .and. (box .eq. boxlink)) 
     &           call cleanup('not yet implemented!')

!     --- call nearest neighbor list
            if ( lneigh ) then
               call updnn( imola )
               call updnn( imolb )
            end if
         else
!     *** recover the reciprocal space vectors
!     if the move is not accepted ***

            call recip(box,vdum,vdum,4)

         end if

! ****************************************************************
! * END INTRABOX SWATCH ADD *
! ****************************************************************

      else

!     --- check if the particle types are in their boxes
         if ( (ncmt(boxa,imolta) .eq. 0) .or. 
     &        (ncmt(boxb,imoltb) .eq. 0) ) then
            lempty = .true.
!     write(iou,*) 'one box out of swatch particle'
         else
!     ---get particle from box a

            iboxal = idint( dble(ncmt(boxa,imolta))*random() ) + 1
            iboxa = parbox(iboxal,boxa,imolta)
            if ( moltyp(iboxa) .ne. imolta ) write(iou,*) 'screwup'
            iboxia = nboxi(iboxa)
            if (iboxia .ne. boxa) call cleanup('problem in swatch')

!     ---get particle from box b

            iboxbl = idint( dble(ncmt(boxb,imoltb))*random() ) + 1
            iboxb = parbox(iboxbl,boxb,imoltb)
            if ( moltyp(iboxb) .ne. imoltb ) write(iou,*) 'screwup'
            iboxib = nboxi(iboxb)
            if (iboxib .ne. boxb) call cleanup('problem in swatch')

!cc--!!!JLR - for test write coordinates of a and b                                              
!            open(unit=91,file='a_init.xyz',status='unknown')                                 
!            write(91,*) nunit(imolta)                                                        
!            write(91,*)                                                                      
!            do izz = 1,nunit(imolta)                                                          
!               write(91,*) 'C ', rxu(iboxa,izz),ryu(iboxa,izz),                                
!     &              rzu(iboxa,izz)                                                            
!            end do                                                                            
!            close(91)                                                                        
!            open(unit=92,file='b_init.xyz',status='unknown')                                 
!            write(92,*) nunit(imoltb)                                                        
!            write(92,*)                                                                      
!            do izz = 1,nunit(imoltb)                                                          
!               write(92,*) 'O ', rxu(iboxb,izz),ryu(iboxb,izz),                                
!     &              rzu(iboxb,izz)                                                            
!            end do                                                                            
!            close(92)                                                                        
!cc--!!!JLR - end of coordinate test            

         end if

!     ---add one attempt to the count for iparty
!     write(iou,*) 'iparty:',iparty,'boxa:',boxa
         bnswat(iparty,ipairb) = bnswat(iparty,ipairb) + 1.0d0

! --- JLR 12-1-09, Count the empty attempts
!         if (lempty) return
         if (lempty) then
            bnswat_empty(iparty,ipairb) =
     &           bnswat_empty(iparty,ipairb) + 1.0d0
            return
         end if
! --- END JLR 12-1-09 ---

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
!     ---store number of units in iunita and iunitb
         iunita = nunit(imolta)
         iunitb = nunit(imoltb)

!     ---store number of each type in the boxes
         orgaia = ncmt(boxa, imolta)
         orgbia = ncmt(boxa, imoltb)
         orgaib = ncmt(boxb, imolta)
         orgbib = ncmt(boxb, imoltb)
         tweight = 1.0d0
         tweiold = 1.0d0

!     --- set the trial energies to zero
         vnbox(boxa)   = 0.0d0
         vninte(boxa)  = 0.0d0
         vnintr(boxa)  = 0.0d0
         vnvibb(boxa)  = 0.0d0
         vntgb(boxa)   = 0.0d0
         vnextb(boxa)  = 0.0d0
         vnbend(boxa)  = 0.0d0
         vnelect(boxa) = 0.0d0
         vnewald(boxa) = 0.0d0

         vnbox(boxb)   = 0.0d0
         vninte(boxb)  = 0.0d0
         vnintr(boxb)  = 0.0d0
         vnvibb(boxb)  = 0.0d0
         vntgb(boxb)   = 0.0d0
         vnextb(boxb)  = 0.0d0
         vnbend(boxb)  = 0.0d0
         vnelect(boxb) = 0.0d0
         vnewald(boxb) = 0.0d0

!     --- compute the rosenbluth weights for each molecule type in each box
         do ic = 1,2
            if ( ic .eq. 1 ) then
               self = iboxa
               other = iboxb
!               s_type = imolta
!               o_type = imoltb
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
!               s_type = imoltb
!               o_type = imolta
               s_type = type_b
               o_type = type_a
               iboxnew = iboxia
               iboxold = iboxib
               imolty = imoltb
               igrow = nugrow(imoltb)
               iunit = nunit(imoltb)
            end if

!     --- store the beads that are identical
            do icbu = 1, nsampos(iparty)
               if ( ic .eq. 1 ) then
!     --- new conformation is for type a
                  new = splist(iparty,icbu,type_a)
                  old = splist(iparty,icbu,type_b)
               else
!     --- new conformation is for type b
                  new = splist(iparty,icbu,type_b)
                  old = splist(iparty,icbu,type_a)
               end if

               rxnew(new) = rxu(other,old)
               rynew(new) = ryu(other,old)
               rznew(new) = rzu(other,old)
            end do

!        --- set up growth schedule
            ifirst = from(ic)
            iprev = prev(ic)

! --- rigid molecule add on

            if (lrigid(imolty)) then
!cc--- JLR 11-24-09
!cc---adding two if statements for rigid swaps
               if (nsampos(iparty).lt.iunit) then

                  if ( (rindex(imolty).eq.0) .or.
     &                 (ifirst.lt.riutry(imolty,1))  ) then

                     if (nsampos(iparty).ge.3) then

                        call align_planes(iparty,self,other,
     &                       s_type,o_type,rxnew,rynew,rznew)

                     else

!     --- calculate new vector from initial bead
                        do j = 1,iunit

                           rxnew(j) = rxnew(ifirst)
     &                          - (rxu(self,ifirst) - rxu(self,j))
                           rynew(j) = rynew(ifirst)
     &                          - (ryu(self,ifirst) - ryu(self,j))
                           rznew(j) = rznew(ifirst)
     &                          - (rzu(self,ifirst) - rzu(self,j))
                           end do

                     end if      !nsampos.lt.3
                  end if         !ifirst.lt.riutry...
               end if            !nsampos.lt.iunit 
               call schedule(igrow,imolty,islen,ifirst,iprev,4)
!cc---END JLR 11-24-09
            else
               call schedule(igrow,imolty,islen,ifirst,iprev,3)
            end if !lrigid
         
!           --- Adding in multiple end regrowths

            if (ncut(iparty,s_type) .gt. 1) then
               do izz = 2,ncut(iparty,s_type)
                  ifirst = from(s_type+2*(izz-1))
                  iprev = prev(s_type+2*(izz-1))
                  
                  call schedule(igrow,imolty,islen,ifirst,iprev,5)
                  
               end do
            end if

            waddnew = 1.0d0
! --- JLR 11-24-09 New stuff for rigid swatch
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there
                  !don't regrow anything in rosenbluth 
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1)) ) then
!                 rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
!                 rigid part is not grown, do rigrot
                  icallrose = 2
               end if
            else
!              flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            end if


            call rosenbluth(.true.,lterm,other,self,imolty,islen,
     &           iboxnew,igrow,waddnew,.false.,vdum2,icallrose)
! --- END JLR 11-24-09
!        --- termination of cbmc attempt due to walk termination ---
            if ( lterm ) return

!        --- propagate new rosenbluth weight
            tweight = tweight*weight*waddnew

! --- end rigid add on 

!     --- save the new coordinates
!     write(iou,*) 'new:',moltyp(self),iunit
            do jj = 1,igrow
               rxut(ic,jj) = rxnew(jj)
               ryut(ic,jj) = rynew(jj)
               rzut(ic,jj) = rznew(jj)
!     write(iou,*) rxnew(jj),rynew(jj),rznew(jj)
            end do
!     write(iou,*) 'old:',moltyp(other),nunit(moltyp(other))
!     do jj = 1,nunit(moltyp(other))
!     write(iou,*) rxu(other,jj),ryu(other,jj),rzu(other,jj)
!     end do

!     --- Corrections for switched beads, and DC-CBMC
!     --- Assign all of the grown new and old beads to rxuion
!     --- with rxuion: new = 2
            iii = 2
            do j=1,igrow
               rxuion(j,iii) = rxnew(j)
               ryuion(j,iii) = rynew(j)
               rzuion(j,iii) = rznew(j)
               qquion(j,iii) = qqu(self,j)
            end do

            if ( iunit .ne. igrow ) then

!     --- for explicit-hydrogen model, put on the hydrogens
!     --- use phony number iins and call explct to add constrained hydrogens

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
!     write(iou,*) 'new:',rxut(ic,j),ryut(ic,j),rzut(ic,j)
               end do
            end if

!     --- Begin DC-CBMC, explicit-hydrogen and
!     --- switched bead Corrections for NEW configuration
!     --- calculate the true site-site energy

!     ??? energy problem for cases which involve the change of the bending type
!     and torsional type for those units swatched!!!!!!

            call ctrmas(.false.,iboxnew,other,8)
            call energy (other,imolty, v, vintra,vinter,vext
     &           ,velect,vewald,iii,iboxnew,1,iunit,.true.,lterm,.false.
     &           ,vdum,.false.,.false.)
            
            if (lterm) then
!     write(iou,*) 'other ',other,' self ',self
!     call cleanup('interesting screwup in CBMC swatch')
               return
            end if
            delen = v - ( vnewt - (vnewbvib + vnewbb + vnewtg )) 
            tweight = tweight*dexp(-beta*delen)

            vnewt     = vnewt + delen
            vnewinter = vinter 
            vnewintra = vintra
            vnewext   = vext 
            vnewelect = velect
            vnewewald = vewald
            
!     End DC-CBMC and switched bead Corrections for NEW configuration

!     --- save the trial energies
            vnbox(iboxnew)   = vnbox(iboxnew) + vnewt
            vninte(iboxnew)  = vninte(iboxnew) + vnewinter
            vnintr(iboxnew)  = vnintr(iboxnew) + vnewintra
            vnvibb(iboxnew)  = vnvibb(iboxnew) + vnewbvib
            vntgb(iboxnew)   = vntgb(iboxnew) + vnewtg
            vnextb(iboxnew)  = vnextb(iboxnew) + vnewext
            vnbend(iboxnew)  = vnbend(iboxnew) + vnewbb
            vnelect(iboxnew) = vnelect(iboxnew) + vnewelect
            vnewald(iboxnew) = vnewald(iboxnew) + vnewewald 
! --- rigid add on

            waddold = 1.0d0

!        --- grow old chain conformation
! --- JLR 11-24-09 New stuff for rigid swatch
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there
                  !don't regrow anything in rosenbluth
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1) )) then
!                 rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
!                 rigid part is not grown, do rigrot
                  icallrose = 2
               end if
            else
!              flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            end if

            call rosenbluth(.false.,lterm,self,self,imolty,
     &           islen,iboxold,igrow,waddold,.false.,vdum2,icallrose)
! --- END JLR 11-24-09

!        --- termination of old walk due to problems generating orientations
            if ( lterm ) then
               write(iou,*) 'SWATCH: old growth rejected'
               return
            end if

!        --- propagate old rosenbluth weight
            tweiold = tweiold*weiold*waddold

!     --- store the old grown beads and explict placed beads positions
!     --- 1 = old conformation
            iii = 1
            do j = 1,iunit
               rxuion(j,1) = rxu(self,j)
               ryuion(j,1) = ryu(self,j)
               rzuion(j,1) = rzu(self,j)
               qquion(j,1) = qqu(self,j)
            end do

!     Begin Correction for DC-CBMC and switched beads for OLD configuration
!     --- correct the acceptance rules 
!     --- calculate the Full rcut site-site energy

            call energy (self,imolty, v, vintra,vinter,vext,velect
     &           ,vewald,iii,iboxold, 1,iunit,.true.,lterm,.false.,vdum
     &           ,.false.,.false.)

            if (lterm) call cleanup('disaster ovrlap in old
     &                           conf SWATCH')
            deleo = v - ( voldt - (voldbvib + voldbb + voldtg) ) 
            
            tweiold = tweiold*dexp(-beta*deleo)

            voldt     = voldt + deleo
            voldintra = vintra
            voldinter = vinter 
            voldext   = vext 
            voldelect = velect
            voldewald = vewald

!     End Correction for DC-CBMC and switched beads for OLD configuration

!     --- save the trial energies
            vnbox(iboxold)   = vnbox(iboxold)  - voldt
            vninte(iboxold)  = vninte(iboxold) - voldinter
            vnintr(iboxold)  = vnintr(iboxold) - voldintra
            vnvibb(iboxold)  = vnvibb(iboxold) - voldbvib 
            vntgb(iboxold)   = vntgb(iboxold)  - voldtg 
            vnextb(iboxold)  = vnextb(iboxold) - voldext
            vnbend(iboxold)  = vnbend(iboxold) - voldbb 
            vnelect(iboxold) = vnelect(iboxold) - voldelect
            vnewald(iboxold) = vnewald(iboxold) - voldewald

         end do

!     --- Perform the Ewald sum reciprical space corrections
         if ( lewald ) then
!     --- added into tweight even though it really contains new-old
!     --- Box A
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

            tweight = tweight * dexp(-beta*delen)

            vnewald(iboxia) = vnewald(iboxia) + delen
            vnbox(iboxia) = vnbox(iboxia) + delen

!     --- Box B
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
            tweight = tweight * dexp(-beta*delen)

            vnewald(iboxib) = vnewald(iboxib) + delen
            vnbox(iboxib) = vnbox(iboxib) + delen

         end if
!     End Ewald-sum Corrections

!     ----------------------------------------------------------------------

         if (ltailc) then
!     ---    add tail corrections
            if (lpbcz) then
               if (lsolid(boxa) .and. .not. lrect(boxa)) then
                  vola = (hmat(boxa,1) * (hmat(boxa,5) * hmat(boxa,9) -
     &                 hmat(boxa,8)*hmat(boxa,6))+
     &                 hmat(boxa,4)*(hmat(boxa,8)
     &                 * hmat(boxa,3)-hmat(boxa,2)*
     &                 hmat(boxa,9))+hmat(boxa,7)
     &                 * (hmat(boxa,2)*hmat(boxa,6)-
     &                 hmat(boxa,5)*hmat(boxa,3)))
               else
                  vola=boxlx(boxa)*boxly(boxa)*boxlz(boxa)
               end if

               if (lsolid(boxb) .and. .not. lrect(boxb)) then
                  volb = (hmat(boxb,1) * (hmat(boxb,5) * hmat(boxb,9) -
     &                 hmat(boxb,8)*hmat(boxb,6))+
     &                 hmat(boxb,4)*(hmat(boxb,8)
     &                 * hmat(boxb,3)-hmat(boxb,2)*
     &                 hmat(boxb,9))+hmat(boxb,7)
     &                 * (hmat(boxb,2)*hmat(boxb,6)-
     &                 hmat(boxb,5)*hmat(boxb,3)))
               else
                  volb=boxlx(boxb)*boxly(boxb)*boxlz(boxb)
               end if
            else
               vola=boxlx(boxa)*boxly(boxa)
               volb=boxlx(boxb)*boxly(boxb)
            end if

!     - for new BOXINS with inserted particle
            do mm = 1,2
               dinsta = 0.0d0
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
! --- JLR 11-24-09 don't do tail corrections for ideal box
               if (.not.lideal(ibox)) then 
!     --- new logic for tail correction (same answer) MGM 3-25-98
                  do jmt = 1, nmolty
                     rho = dble( ncmt(ibox,jmt) )
                     if ( jmt .eq. imolin ) rho = rho + 1.0d0
                     if ( jmt .eq. imolrm ) rho = rho - 1.0d0
                     rho = rho / dvol
                     do imt = 1, nmolty
                        dicount = ncmt(ibox,imt)
                        if ( imt .eq. imolin ) dicount = dicount + 1
                        if ( imt .eq. imolrm ) dicount = dicount - 1
                        dinsta = dinsta + 
     &                       dicount * coru(imt,jmt,rho,ibox)
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
                  dinsta = 0.0d0
               end if
! --- END JLR 11-24-09              
               dinsta = dinsta - vtailb( ibox )

               tweight=tweight*dexp(-beta*dinsta)

               vntail(ibox) = dinsta 
            end do
         else
            vntail(boxa) = 0.0d0
            vntail(boxb) = 0.0d0
         end if

         wnlog = dlog10( tweight )
         wolog = dlog10( tweiold )
         wdlog = wnlog - wolog
         if ( wdlog .lt. -softcut ) then
!     write(iou,*) '### underflow in wratio calculation ###'
            return
         end if
         wswat = ( tweight / tweiold ) * ( dble(orgaia*orgbib) /
     &        dble((orgbia+1)*(orgaib+1)) ) *
     &        dexp(beta*(eta2(boxa,imolta)
     &        +eta2(boxb,imoltb)-eta2(boxa,imoltb)
     &        -eta2(boxb,imolta)))

!     write(iou,*) 'imolta,imoltb',imolta,imoltb
!     write(iou,*) 'wswat,tweight,tweiold',wswat,tweight,tweiold
         if ( random() .le. wswat ) then
!     *** we can now accept !!!!! ***
            bsswat(iparty,ipairb) = bsswat(iparty,ipairb) + 1.0d0
!     write(iou,*) 'SWATCH ACCEPTED',iboxa,iboxb
            do jj = 1,2
               if ( jj .eq. 1 ) ic = boxa
               if ( jj .eq. 2 ) ic = boxb
               vbox(ic)     = vbox(ic)    + vnbox(ic) + vntail(ic)
               vinterb(ic)  = vinterb(ic) + vninte(ic) + vntail(ic)
               vintrab(ic)  = vintrab(ic) + vnintr(ic)
               vvibb(ic)    = vvibb(ic)   + vnvibb(ic)
               vtgb(ic)     = vtgb(ic)    + vntgb(ic)
               vextb(ic)    = vextb(ic)   + vnextb(ic)
               vbendb(ic)   = vbendb(ic)  + vnbend(ic)
               vtailb(ic)   = vtailb(ic)  + vntail(ic)
               velectb(ic)  = velectb(ic) + vnelect(ic) + vnewald(ic)
            end do

!     ---update book keeping
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
!     -- update reciprocal-space sum
               call recip(boxa,vdum,vdum,2)
               call recip(boxb,vdum,vdum,2)
            end if

!     ---update center of mass
            call ctrmas(.false.,boxa,iboxb,8)
            call ctrmas(.false.,boxb,iboxa,8)

            if (licell .and. (boxa .eq. boxlink .or. boxb .eq. boxlink))
     &           call cleanup('not yet implemented!')

!     --- call nearest neighbor list
            if ( lneigh ) then
               call updnn( iboxa )
               call updnn( iboxb )
            end if

!cc--!!!JLR - for test 2 write final coordinates of a and b                                      
!         open(unit=93,file='a_final.xyz',status='unknown')                                   
!         write(93,*) iunita                                                                  
!         write(93,*)                                                                         
!         do izz = 1,iunita                                                                    
!            write(93,*) 'C ', rxu(iboxa,izz),ryu(iboxa,izz),rzu(iboxa,izz)                      
!         end do                                                                               
!                                                                                             
!         open(unit=93,file='b_final.xyz',status='unknown')                                   
!         write(93,*) iunitb                                                                  
!         write(93,*)                                                                         
!         do izz = 1,iunitb                                                                    
!            write(93,*) 'C ', rxu(iboxb,izz),ryu(iboxb,izz),rzu(iboxb,izz)                      
!         end do                                                                               
!         call cleanup('END OF SWATCH TEST')                                                           
!cc--!!!JLR - end of coordinate test         

         end if
! -----------------------------------------------------------------

      end if

!      write(iou,*) 'end SWATCH'

      return
      end




