      subroutine swatch 

c swatch
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
c John Stubbs, and Collin Wick and Ilja Siepmann  
c                     
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to 
c
c Free Software Foundation, Inc. 
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      **********************************************************
c    *** Added intrabox move for two particles within one box   ***
c    *** in combined move that shares the same parameters.      ***
c    *** Will also accept rigid (lrigid) molecules.  Contains   ***
c    *** several critical bug fixes as well.                    ***
c    ***                                     9-25-02 JMS        ***
c      **********************************************************
 
      implicit none
 
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'cbmc.inc'
      include 'inputdata.inc'
      include 'rosen.inc'
      include 'swtcmove.inc'
      include 'ewaldsum.inc'
      include 'cell.inc'

      logical lempty,lterm

      integer type_a,type_b,from,prev,self,iboxnew,iboxold
     &     ,imolty,igrow,new,old,islen,ifirst,iprev,iii,j
      integer oldchain,newchain,oldunit,newunit,ncount,iunit,iins

      integer ic,ibox,icbu,jj,mm,imt,jmt,imolin,imolrm
      integer boxa,boxb,ipair,imolta,imoltb,iboxa
     &  ,iboxb,iboxal,iboxbl,iboxia,iboxib,iunita,iunitb,orgaia
     &  ,orgaib,orgbia,orgbib,ipairb 

      double precision random,tweight,tweiold,rxut,ryut,rzut
     &  ,dvol,vola,volb,rho,coru,dinsta,rpair

      double precision vnbox,vninte,vnintr,vnvibb,vntgb,vnextb
     &  ,vnbend,vntail,vnelect,vnewald,wnlog,wolog,wdlog,wswat
         
      double precision v,vintra,vinter,vext,velect,vewald,vdum,delen
     &     ,deleo,dicount,vrecipn,vrecipo
c * additions from iswatch

      integer zz,zzz,box,iboxi,bdmol_a,bdmol_b
      integer imola,imolb,moltaid,moltbid,i,ii,iu

      integer s_type, o_type, thisbox, otherbox

      double precision rx_1(numax),ry_1(numax),rz_1(numax),dummy

c      dimension from(2),prev(2)
      dimension from(2*numax),prev(2*numax)

c      dimension rxut(2,numax),ryut(2,numax),rzut(2,numax)
      dimension rxut(4,numax),ryut(4,numax),rzut(4,numax)

c * end additions

      double precision waddold,waddnew,vdum2

      dimension vnbox(nbxmax),vninte(nbxmax),vnintr(nbxmax)
     &  ,vnvibb(nbxmax),vntgb(nbxmax),vnextb(nbxmax),vnbend(nbxmax)
     &  ,vntail(nbxmax),vnelect(nbxmax),vnewald(nbxmax)

c --- JLR 11-24-09
      integer icallrose
c --- END JLR 11-24-09

c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c      write(iou,*) 'start SWATCH'

      lempty = .false.
c     ---randomly select chains to switch between boxa boxb
c     --- select a pair type to switch 
      if ( nswaty .gt. 1 ) then
         rpair = random()
         do 96 ipair = 1, nswaty
            if ( rpair .lt. pmsatc(ipair) ) then
               iparty = ipair
               goto 97
            endif
 96      continue
      else
         iparty = 1 
      endif

c *** randomly select the molecules from the pair ***

      if (random() .lt. 0.5d0) then
         type_a = 1
         type_b = 2
      else
         type_a = 2
         type_b = 1
      endif

c     ---select the molecules from the pair
 97   imolta = nswatb(iparty,1)
      imoltb = nswatb(iparty,2)
c ???!!!
      type_a = 1
      type_b = 2

c     ---choose box A and box B at random
      if ( nswtcb(iparty) .gt. 1 ) then
         rpair = random()
         do 98 ipair = 1, nswtcb(iparty)
            if ( rpair .lt. pmswtcb(iparty,ipair) ) then
               ipairb = ipair
               goto 99
            endif
 98      continue
      else
         ipairb = 1
      endif
         
 99   if (random() .lt. 0.5d0) then
         boxa=box3(iparty,ipairb)
         boxb=box4(iparty,ipairb)
      else
         boxa=box4(iparty,ipairb)
         boxb=box3(iparty,ipairb)
      endif

      if (boxa .eq. boxb) then

c ****************************************************************
c * INTRABOX SWATCH ADD *
c ****************************************************************
c ****************************************************************
 
c * liswatch prevents non-grown beads from being included in 
c * the new growth in boltz.f

c      write(iou,*) 'start iSWATCH'

c$$$c *** randomly select chains to switch ***
c$$$c *** select a pair type to switch *** 
c$$$      if ( niswaty .gt. 1 ) then
c$$$         rpair = random()
c$$$         do 96 ipair = 1, niswaty
c$$$            if ( rpair .lt. pmisatc(ipair) ) then
c$$$               iparty = ipair
c$$$               goto 97
c$$$            endif
c$$$ 96      continue
c$$$      else
c$$$         iparty = 1 
c$$$      endif

c$$$c *** box determinations ***
c$$$
c$$$      if ( lgibbs ) then
c$$$         if (random() .lt. 0.5d0) then
c$$$            thisbox = 1
c$$$            otherbox = 2
c$$$         else
c$$$            thisbox = 2
c$$$            otherbox = 1
c$$$         endif
c$$$      else
c$$$         thisbox = 1
c$$$         otherbox = 2
c$$$      endif
c$$$            
c$$$      box = thisbox

c *** box determinations ***

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
              endif
           else
              otherbox = 2
           endif

c     *** get particle of type a ***

           moltaid = idint(dble (ncmt(box,imolta))*random()) + 1

c     * imola is the overall molecule id number (runs from 1 to total #) *
           imola = parbox( moltaid, box, imolta) 

           if (moltyp(imola) .ne. imolta) 
     &                   write(iou,*) 'screwup in iswatch'

           iboxi = nboxi(imola)

           if (iboxi .ne. box) stop 'problem in iswatch'

c     *** get particle of type b ***

           moltbid = idint(dble (ncmt(box,imoltb))*random()) + 1
           imolb = parbox( moltbid, box, imoltb) 

           if (moltyp(imolb) .ne. imoltb) 
     &                write(iou,*) 'screwup in iswatch'
           iboxi = nboxi(imolb)

           if (iboxi .ne. box) stop 'problem in iswatch'
         endif  
c     *** add one attempt to the count for iparty
c$$$  bniswat(iparty,box) = bniswat(iparty,box) + 1.0d0   
         bnswat(iparty,ipairb) = bnswat(iparty,ipairb) + 1.0d0   
c --- JLR 12-1-09 count the empty attempts
c         if (lempty) return 
         if (lempty) then
            bnswat_empty(iparty,ipairb) =
     &           bnswat_empty(iparty,ipairb) + 1.0d0
            return
         endif
c --- END JLR 12-1-09

c     * write out the molecule numbers of the pair *
c     write(iou,*) imola,imolb

c     ***************************
c     *** Begin Growth Setups ***
c     ***************************

c     *** assign from and prev for each moltyp

         do zz = 1,ncut(iparty,type_a)

            from(type_a+2*(zz-1)) = gswatc(iparty,type_a,1+2*(zz-1)) 
            prev(type_a+2*(zz-1)) = gswatc(iparty,type_a,2+2*(zz-1)) 

         enddo

         do zz = 1,ncut(iparty,type_b)

            from(type_b+2*(zz-1)) = gswatc(iparty,type_b,1+2*(zz-1)) 
            prev(type_b+2*(zz-1)) = gswatc(iparty,type_b,2+2*(zz-1)) 

         enddo

c$$$         write(iou,*) 'mol a'
c$$$  
c$$$         do zz = 1,ncut(type_a)
c$$$            write(iou,*) zz,'from',from(type_a+2*(zz-1)),' prev',
c$$$     &        prev(type_a+2*(zz-1))
c$$$         enddo
c$$$  
c$$$         write(iou,*) 'mol b'
c$$$  
c$$$         write(iou,*) gswatc(iparty,type_b,3),gswatc(iparty,type_b,4)
c$$$         write(iou,*) from(3),prev(3),type_b
c$$$  
c$$$         do zz = 1,ncut(type_b)
c$$$            write(iou,*) zz,'from',from(type_b+2*(zz-1)),' prev',
c$$$     &        prev(type_b+2*(zz-1))
c$$$         enddo
c$$$         
c$$$         stop
      
c     *** store number of units in iunita and iunitb
         iunita = nunit(imolta)
         iunitb = nunit(imoltb)

c     *** initialize trial weights ***
         tweight = 1.0d0
         tweiold = 1.0d0

c     *** set the trial energies to zero ***
         vnbox(box)   = 0.0d0
         vninte(box)  = 0.0d0
         vnintr(box)  = 0.0d0
         vnvibb(box)  = 0.0d0
         vntgb(box)   = 0.0d0
         vnextb(box)  = 0.0d0
         vnbend(box)  = 0.0d0
         vnelect(box) = 0.0d0
         vnewald(box) = 0.0d0

c     *** store position 1, a's original site ***

         do zz = 1, nunit(imolta)

            rx_1(zz) = rxu(imola,zz)
            ry_1(zz) = ryu(imola,zz)
            rz_1(zz) = rzu(imola,zz)

         enddo

c     *** store same bead coordinates for molecules a and b ***

         do icbu = 1, nsampos(iparty)

            bdmol_a = splist(iparty,icbu,type_a)

            rxut(3,icbu) = rxu(imola,bdmol_a)
            ryut(3,icbu) = ryu(imola,bdmol_a)
            rzut(3,icbu) = rzu(imola,bdmol_a)

            bdmol_b = splist(iparty,icbu,type_b)

            rxut(4,icbu) = rxu(imolb,bdmol_b)
            ryut(4,icbu) = ryu(imolb,bdmol_b)
            rzut(4,icbu) = rzu(imolb,bdmol_b)

         enddo

c     *****************************
c     *** Liswinc Determination ***
c     *****************************

c     *** only need this for molecule b, since b is grown after a
c     *** but a is grown when some of b doesn't exist

         self = imolb
         other = imola
         imolty = imoltb
         igrow = nugrow(imoltb)

         ifirst = from(type_b)
         iprev = prev(type_b)

c     * determine which beads aren't in the same positions *

         call schedule(igrow,imolty,islen,ifirst,iprev,3)

         if (ncut(iparty,type_b) .gt. 1) then
            do zz = 2,ncut(iparty,type_b)
               ifirst = from(type_b+2*(zz-1))
               iprev = prev(type_b+2*(zz-1))

               call schedule(igrow,imolty,islen,ifirst,iprev,5)

            enddo
         endif

c     * assign growth schedule for molecule b *

         do zz = 1,nunit(imolty)

            liswinc(zz,imolty) = lexshed(zz)

c     write(iou,*) zz, liswinc(zz,imolty)
            
         enddo

c     ******************
c     *** ROSENBLUTH ***
c     ******************

c     *** start rosenbluth weight calculations ***
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

            endif

c     *** assign same bead coordinates for new growth ***
            do icbu = 1, nsampos(iparty)

               new = splist(iparty,icbu,s_type)
               old = splist(iparty,icbu,o_type)

               rxnew(new) = rxu(other,old)
               rynew(new) = ryu(other,old)
               rznew(new) = rzu(other,old)

            enddo

c        --- set up growth schedule
            ifirst = from(s_type)
            iprev = prev(s_type)

c     --- lrigid add on 

            if (lrigid(imolty)) then
ccc--- JLR 11-24-09
ccc--- adding some if statements for rigid swaps:
ccc--- if we have swapped the whole rigid part we will 
ccc--- not compute the vectors from ifirst in old position                                   
               if ( nsampos(iparty).lt.iunit) then

                  if ( (rindex(imolty).eq.0) .or.
     &                 (ifirst.lt.riutry(imolty,1))  ) then

                     if (nsampos(iparty).ge.3) then

                        call align_planes(iparty,self,other,
     &                       s_type,o_type,rxnew,rynew,rznew)

                     else

c     --- calculate new vector from initial bead                                                    
                        do j = 1,iunit
                           rxnew(j) = rxnew(ifirst)
     &                          - (rxu(self,ifirst) - rxu(self,j))
                           rynew(j) = rynew(ifirst)
     &                          - (ryu(self,ifirst) - ryu(self,j))
                           rznew(j) = rznew(ifirst)
     &                          - (rzu(self,ifirst) - rzu(self,j))
                        enddo

                     endif      !nsampos.eq.3                                                       
                  endif         !ifirst.lt.riutry...                                                
               endif            !nsampos.lt.iunit                                                   
               call schedule(igrow,imolty,islen,ifirst,iprev,4)
ccc---END JLR 11-24-09
            else
               call schedule(igrow,imolty,islen,ifirst,iprev,3)
            endif

c * I wonder if this works with lrigid?                                                             
            if (ncut(iparty,s_type) .gt. 1) then
               do zz = 2,ncut(iparty,s_type)
                  ifirst = from(s_type+2*(zz-1))
                  iprev = prev(s_type+2*(zz-1))

                  call schedule(igrow,imolty,islen,ifirst,iprev,5)

               enddo
            endif



c     *********** Commenting out for new combined code JMS 6-20-00 *****
c     *** Call qqcheck to setup the group based qq cutoff ***
c$$$  if ( lelect(imolty) ) then
c$$$  call qqcheck(other,box,rxnew(1),rynew(1),rznew(1))
c$$$  endif

c     ******************
c     *** new growth ***
c     ******************

c     *** moving molecules for rosenbluth ***
            if (ic .eq. 1) then

c     * putting molecule b in position 1 *

               do zz = 1, nsampos(iparty)

                  bdmol_b = splist(iparty,zz,type_b)

                  rxu(other,bdmol_b) = rxut(3,zz)
                  ryu(other,bdmol_b) = ryut(3,zz)
                  rzu(other,bdmol_b) = rzut(3,zz)

               enddo

            else

c     * putting molecule a into its (fully grown) trial position 2 *

               do zz = 1, nunit(moltyp(other))

                  rxu(other,zz) = rxut(1,zz)
                  ryu(other,zz) = ryut(1,zz)
                  rzu(other,zz) = rzut(1,zz)

               enddo

            endif

c     --- grow molecules

c     --- changing for lrigid to include waddnew
            waddnew = 1.0d0

c     --- grow new chain conformation
c --- JLR 11-24-09
c --- Different logic/calls to rosenbluth for rigid molecules 
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there                                                      
                  !don't regrow anything in rosenbluth  
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1) )) then
c                 rigid part is grown, don't do rigrot in rosebluth 
                  icallrose = 3
               else
c                 rigid part is not grown, do rigrot  
                  icallrose = 2
               endif
            else
c              flexible molecule call rosenbluth in normal fashion  
               icallrose = 2
            endif

            call rosenbluth(.true.,lterm,self,self,imolty,islen
     &           ,box ,igrow,waddnew,.false.,vdum2,icallrose)
c ---END JLR 11-24-09

c     --- propagate new rosenbluth weight
            tweight = tweight*weight*waddnew

c --- end rigid add on

c     * moving molecules back *
            if (ic .eq. 1) then

               do zz = 1, nsampos(iparty)

                  bdmol_b = splist(iparty,zz,type_b)

                  rxu(other,bdmol_b) = rxut(4,zz)
                  ryu(other,bdmol_b) = ryut(4,zz)
                  rzu(other,bdmol_b) = rzut(4,zz)

               enddo

            else

               do zz = 1, nunit(moltyp(other))
                  rxu(other,zz) = rx_1(zz)
                  ryu(other,zz) = ry_1(zz)
                  rzu(other,zz) = rz_1(zz)
               enddo

            endif

            if ( lterm ) then
c     *** termination of cbmc attempt due to walk termination ***
c     write(iou,*) 'iSWATCH:new growth rejected',ic
c * reset liswatch
               liswatch = .false.
               return
            else
c     write(iou,*) 'iSWATCH:new growth accepted',ic
            endif

c     *** save the new coordinates ***
            do jj = 1,igrow
               rxut(ic,jj) = rxnew(jj)
               ryut(ic,jj) = rynew(jj)
               rzut(ic,jj) = rznew(jj)
            enddo

c     --- Corrections for switched beads, and DC-CBMC
c     --- Assign all of the grown new and old beads to rxuion
c     --- with rxuion: new = 2
            iii = 2
            do j=1,igrow
               rxuion(j,iii) = rxnew(j)
               ryuion(j,iii) = rynew(j)
               rzuion(j,iii) = rznew(j)
               qquion(j,iii) = qqu(self,j)
            enddo

c     *** added from new-combined code ***
            if ( iunit .ne. igrow ) then

c     --- for explicit-hydrogen model, put on the hydrogens
c     --- use phony number iins and call explct to add constrained hydrogens

               iins = nchain + 1
               moltyp(iins) = imolty

               do j=1,igrow
                  rxu(iins,j) = rxnew(j)
                  ryu(iins,j) = rynew(j)
                  rzu(iins,j) = rznew(j)
               enddo

               call explct(iins,vdum,.false.,.false.)

               do j = igrow + 1, iunit

                  rxuion(j,iii) = rxu(iins,j)
                  ryuion(j,iii) = ryu(iins,j)
                  rzuion(j,iii) = rzu(iins,j)
                  qquion(j,iii) = qqu(self,j)

                  rxut(ic,j) = rxu(iins,j)
                  ryut(ic,j) = ryu(iins,j)
                  rzut(ic,j) = rzu(iins,j)

               enddo
            endif


c     --- Begin DC-CBMC and switched bead Corrections for NEW configuration
c     --- calculate the true site-site energy

            if (ic .eq. 1) then
c     * exclude molecule b from energy calculation, put in other box *
               nboxi(imolb) = otherbox

            else
c     * put molecule a into position 2 (fully grown trial position)
c     * for the energy of b's new position (second time around)
               do zz = 1,nunit(imolta)
                  rxu(imola,zz) = rxut(1,zz)
                  ryu(imola,zz) = ryut(1,zz)
                  rzu(imola,zz) = rzut(1,zz)
               enddo
            endif

c     *** get energy of configuration ***

            call ctrmas(.false.,box,self,8)
            call energy (self,imolty, v, vintra,vinter,vext
     &           ,velect,vewald,iii,box,1,iunit,.true.,lterm,.false.
     &           ,vdum,.false.,.false.)

c     *** return to normal ***
            if (ic .eq. 1) then

c     * return b to original box *
               nboxi(imolb) = thisbox
            else

c     * return a to position 1 *
               do zz = 1, nunit(imolta)
                  rxu(imola,zz) = rx_1(zz)
                  ryu(imola,zz) = ry_1(zz)
                  rzu(imola,zz) = rz_1(zz)

               enddo

            endif
            
            if (lterm) then
               write(iou,*) 'other ',other,' self ',self,ic
               stop 'interesting screwup in CBMC iswatch'
            endif

c     *** add on the changes in energy ***

            delen = v - ( vnewt - (vnewbvib + vnewbb + vnewtg )) 

            tweight = tweight*dexp(-beta*delen)
            vnewt     = vnewt + delen
            vnewinter = vinter 
            vnewintra = vintra
            vnewext   = vext 
            vnewelect = velect
            vnewewald = vewald

c     End DC-CBMC and switched bead Corrections for NEW configuration

c     *** save the trial energies ***
            vnbox(box)   = vnbox(box) + vnewt
            vninte(box)  = vninte(box) + vnewinter
            vnintr(box)  = vnintr(box) + vnewintra
            vnvibb(box)  = vnvibb(box) + vnewbvib
            vntgb(box)   = vntgb(box) + vnewtg
            vnextb(box)  = vnextb(box) + vnewext
            vnbend(box)  = vnbend(box) + vnewbb
            vnelect(box) = vnelect(box) + vnewelect
            vnewald(box) = vnewald(box) + vnewewald

c     ******************
c     *** old growth ***
c     ******************

c     --- lrigid add on

            waddold = 1.0d0

c        --- grow old chain conformation
c --- JLR 11-24-09
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
c                 molecule is all there
c                 don't regrow anything in rosenbluth 
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1)) ) then
c                 rigid part is grown, don't do rigrot in rosebluth 
                  icallrose = 3
               else
c                 rigid part is not grown, do rigrot
                  icallrose = 2
               endif
            else
c              flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            endif

            call rosenbluth(.false.,lterm,self,self,imolty,
     &           islen,box,igrow,waddold,.false.,vdum2,icallrose)
c --- END JLR 11-24-09

            if ( lterm ) then
c     *** termination of old walk due to problems generating orientations ***
c     write(iou,*) 'iSWATCH:old growth rejected',ic
c * reset liswatch
               liswatch = .false.
               return
            else 
c     write(iou,*) 'iSWATCH:old growth accepted',ic
            endif

c        --- propagate old rosenbluth weight
            tweiold = tweiold*weiold*waddold

c --- end rigid add on

c     --- store the old grown beads and explict placed beads positions
c     --- 1 = old conformation
            iii = 1
            do j = 1,iunit
               rxuion(j,1) = rxu(self,j)
               ryuion(j,1) = ryu(self,j)
               rzuion(j,1) = rzu(self,j)
               qquion(j,1) = qqu(self,j)
            enddo

c     Begin Correction for DC-CBMC and switched beads for OLD configuration
c     --- correct the acceptance rules 
c     --- calculate the Full rcut site-site energy

c     *** excluding molecule b for first loop ***
            if (ic .eq. 1) then
               nboxi(imolb) = otherbox
            endif

c     * get total energy *

            call energy (self,imolty, v, vintra,vinter,vext,velect
     &           ,vewald,iii,box, 1,iunit,.true.,lterm,.false.,vdum
     &           ,.false.,.false.)

            if (ic .eq. 1) then
c     * return b to current box *
               nboxi(imolb) = thisbox
            endif

            if (lterm) stop 'disaster ovrlap in old conf iSWATCH'
            deleo = v - ( voldt - (voldbvib + voldbb + voldtg) ) 

            tweiold = tweiold*dexp(-beta*deleo)
            voldt     = voldt + deleo
            voldintra = vintra
            voldinter = vinter 
            voldext   = vext 
            voldelect = velect
            voldewald = vewald

c     End Correction for DC-CBMC and switched beads for OLD configuration

c     *** save the trial energies ***
            vnbox(box)   = vnbox(box)  - voldt
            vninte(box)  = vninte(box) - voldinter
            vnintr(box)  = vnintr(box) - voldintra
            vnvibb(box)  = vnvibb(box) - voldbvib 
            vntgb(box)   = vntgb(box)  - voldtg 
            vnextb(box)  = vnextb(box) - voldext
            vnbend(box)  = vnbend(box) - voldbb 
            vnelect(box) = vnelect(box) - voldelect
            vnewald(box) = vnewald(box) - voldewald

         enddo

c     *****************************
c     *** Ewald sum corrections ***
c     *****************************

c     --- Perform the Ewald sum reciprical space corrections
         if ( lewald ) then
c     --- added into tweight even though it really contains new-old

c     *** store the reciprocal space vector
            call recip(box,vdum,vdum,3)

c     *** Position 1 ***
            oldchain = imola
            newchain = imolb
            oldunit = nunit(imolta)
            newunit = nunit(imoltb)

            do j = 1,oldunit
               rxuion(j,1) = rxu(oldchain,j)
               ryuion(j,1) = ryu(oldchain,j)
               rzuion(j,1) = rzu(oldchain,j)
               qquion(j,1) = qqu(oldchain,j)
            enddo
            moltion(1) = imolta
            do j = 1,newunit
               rxuion(j,2) = rxut(2,j)
               ryuion(j,2) = ryut(2,j)
               rzuion(j,2) = rzut(2,j)
               qquion(j,2) = qqu(newchain,j)
            enddo
            moltion(2) = imoltb

            call recip(box,vrecipn,vrecipo,1)

            delen = vrecipn - vrecipo 
            tweight = tweight * dexp(-beta*delen)

            vnewald(box) = vnewald(box) + delen
            vnbox(box) = vnbox(box) + delen

c     * update the reciprocal space terms *
            call recip(box,vdum,vdum,2)

c     *** Position 2 ***
            oldchain = imolb
            newchain = imola
            oldunit = nunit(imoltb)
            newunit = nunit(imolta)

            do j = 1,oldunit
               rxuion(j,1) = rxu(oldchain,j)
               ryuion(j,1) = ryu(oldchain,j)
               rzuion(j,1) = rzu(oldchain,j)
               qquion(j,1) = qqu(oldchain,j)
            enddo
            moltion(1) = imoltb
            do j = 1,newunit
               rxuion(j,2) = rxut(1,j)
               ryuion(j,2) = ryut(1,j)
               rzuion(j,2) = rzut(1,j)
               qquion(j,2) = qqu(newchain,j)
            enddo
            moltion(2) = imolta

            call recip(box,vrecipn,vrecipo,1)

            delen = vrecipn - vrecipo 
            tweight = tweight * dexp(-beta*delen)

            vnewald(box) = vnewald(box) + delen
            vnbox(box) = vnbox(box) + delen

         endif
c     End Ewald-sum Corrections

c     ----------------------------------------------------------------------

c     *************************
c     *** Ratio Calculation ***
c     *************************

         wnlog = dlog10( tweight )
         wolog = dlog10( tweiold )
         wdlog = wnlog - wolog
         if ( wdlog .lt. -softcut ) then
c     write(iou,*) '### underflow in wratio calculation ###'
            call recip(box,vdum,vdum,4)
            return
         endif
         
         wswat = ( tweight / tweiold ) 

         if ( random() .le. wswat ) then
c     *** we can now accept !!!!! ***
            bsswat(iparty,ipairb) = bsswat(iparty,ipairb) + 1.0d0
c     write(iou,*) 'SWATCH ACCEPTED',imola,imolb

            vbox(box)     = vbox(box)    + vnbox(box)
            vinterb(box)  = vinterb(box) + vninte(box)
            vintrab(box)  = vintrab(box) + vnintr(box)
            vvibb(box)    = vvibb(box)   + vnvibb(box)
            vtgb(box)     = vtgb(box)    + vntgb(box)
            vextb(box)    = vextb(box)   + vnextb(box)
            vbendb(box)   = vbendb(box)  + vnbend(box)
            velectb(box)  = velectb(box) + vnelect(box) + vnewald(box)

c     *** update book keeping ***
c     *** assign new geometries ***
            do ic = 1,iunita
               rxu(imola,ic) = rxut(1,ic)
               ryu(imola,ic) = ryut(1,ic)
               rzu(imola,ic) = rzut(1,ic)
c     write(iou,*) 'imola:',imola
c     write(iou,*) rxu(imola,ic),ryu(imola,ic),rzu(imola,ic)
            enddo

            do ic = 1,iunitb
               rxu(imolb,ic) = rxut(2,ic)
               ryu(imolb,ic) = ryut(2,ic)
               rzu(imolb,ic) = rzut(2,ic)
c     write(iou,*) 'imolb:',imolb
c     write(iou,*) rxu(imolb,ic),ryu(imolb,ic),rzu(imolb,ic)
            enddo
            if ( lewald ) then
c     -- update reciprocal-space sum
               call recip(iboxi,vdum,vdum,2)
            endif

c     *** update center of mass
            call ctrmas(.false.,box,imolb,8)
            call ctrmas(.false.,box,imola,8)

            if (licell .and. (box .eq. boxlink)) 
     &           stop 'not yet implemented!'

c     --- call nearest neighbor list
            if ( lneigh ) then
               call updnn( imola )
               call updnn( imolb )
            endif
         else
c     *** recover the reciprocal space vectors
c     if the move is not accepted ***

            call recip(box,vdum,vdum,4)

         endif

c ****************************************************************
c * END INTRABOX SWATCH ADD *
c ****************************************************************

      else

c     --- check if the particle types are in their boxes
         if ( (ncmt(boxa,imolta) .eq. 0) .or. 
     &        (ncmt(boxb,imoltb) .eq. 0) ) then
            lempty = .true.
c     write(iou,*) 'one box out of swatch particle'
         else
c     ---get particle from box a

            iboxal = idint( dble(ncmt(boxa,imolta))*random() ) + 1
            iboxa = parbox(iboxal,boxa,imolta)
            if ( moltyp(iboxa) .ne. imolta ) write(iou,*) 'screwup'
            iboxia = nboxi(iboxa)
            if (iboxia .ne. boxa) stop 'problem in swatch'

c     ---get particle from box b

            iboxbl = idint( dble(ncmt(boxb,imoltb))*random() ) + 1
            iboxb = parbox(iboxbl,boxb,imoltb)
            if ( moltyp(iboxb) .ne. imoltb ) write(iou,*) 'screwup'
            iboxib = nboxi(iboxb)
            if (iboxib .ne. boxb) stop 'problem in swatch'

ccc--!!!JLR - for test write coordinates of a and b                                              
c            open(unit=91,file='a_init.xyz',status='unknown')                                 
c            write(91,*) nunit(imolta)                                                        
c            write(91,*)                                                                      
c            do zz = 1,nunit(imolta)                                                          
c               write(91,*) 'C ', rxu(iboxa,zz),ryu(iboxa,zz),                                
c     &              rzu(iboxa,zz)                                                            
c            enddo                                                                            
c            close(91)                                                                        
c            open(unit=92,file='b_init.xyz',status='unknown')                                 
c            write(92,*) nunit(imoltb)                                                        
c            write(92,*)                                                                      
c            do zz = 1,nunit(imoltb)                                                          
c               write(92,*) 'O ', rxu(iboxb,zz),ryu(iboxb,zz),                                
c     &              rzu(iboxb,zz)                                                            
c            enddo                                                                            
c            close(92)                                                                        
ccc--!!!JLR - end of coordinate test            

         endif

c     ---add one attempt to the count for iparty
c     write(iou,*) 'iparty:',iparty,'boxa:',boxa
         bnswat(iparty,ipairb) = bnswat(iparty,ipairb) + 1.0d0

c --- JLR 12-1-09, Count the empty attempts
c         if (lempty) return
         if (lempty) then
            bnswat_empty(iparty,ipairb) =
     &           bnswat_empty(iparty,ipairb) + 1.0d0
            return
         endif
c --- END JLR 12-1-09 ---

c$$$c     --- assign from and prev for each moltyp
c$$$         from(1) = gswatc(iparty,type_a,1) 
c$$$         prev(1) = gswatc(iparty,type_a,2) 
c$$$         from(2) = gswatc(iparty,type_b,1) 
c$$$         prev(2) = gswatc(iparty,type_b,2) 

         do zz = 1,ncut(iparty,type_a)

            from(type_a+2*(zz-1)) = gswatc(iparty,type_a,1+2*(zz-1)) 
            prev(type_a+2*(zz-1)) = gswatc(iparty,type_a,2+2*(zz-1)) 

         enddo

         do zz = 1,ncut(iparty,type_b)

            from(type_b+2*(zz-1)) = gswatc(iparty,type_b,1+2*(zz-1)) 
            prev(type_b+2*(zz-1)) = gswatc(iparty,type_b,2+2*(zz-1)) 

         enddo
c     ---store number of units in iunita and iunitb
         iunita = nunit(imolta)
         iunitb = nunit(imoltb)

c     ---store number of each type in the boxes
         orgaia = ncmt(boxa, imolta)
         orgbia = ncmt(boxa, imoltb)
         orgaib = ncmt(boxb, imolta)
         orgbib = ncmt(boxb, imoltb)
         tweight = 1.0d0
         tweiold = 1.0d0

c     --- set the trial energies to zero
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

c     --- compute the rosenbluth weights for each molecule type in each box
         do ic = 1,2
            if ( ic .eq. 1 ) then
               self = iboxa
               other = iboxb
c               s_type = imolta
c               o_type = imoltb
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
c               s_type = imoltb
c               o_type = imolta
               s_type = type_b
               o_type = type_a
               iboxnew = iboxia
               iboxold = iboxib
               imolty = imoltb
               igrow = nugrow(imoltb)
               iunit = nunit(imoltb)
            endif

c     --- store the beads that are identical
            do icbu = 1, nsampos(iparty)
               if ( ic .eq. 1 ) then
c     --- new conformation is for type a
                  new = splist(iparty,icbu,type_a)
                  old = splist(iparty,icbu,type_b)
               else
c     --- new conformation is for type b
                  new = splist(iparty,icbu,type_b)
                  old = splist(iparty,icbu,type_a)
               endif

               rxnew(new) = rxu(other,old)
               rynew(new) = ryu(other,old)
               rznew(new) = rzu(other,old)
            enddo

c        --- set up growth schedule
            ifirst = from(ic)
            iprev = prev(ic)

c --- rigid molecule add on

            if (lrigid(imolty)) then
ccc--- JLR 11-24-09
ccc---adding two if statements for rigid swaps
               if (nsampos(iparty).lt.iunit) then

                  if ( (rindex(imolty).eq.0) .or.
     &                 (ifirst.lt.riutry(imolty,1))  ) then

                     if (nsampos(iparty).ge.3) then

                        call align_planes(iparty,self,other,
     &                       s_type,o_type,rxnew,rynew,rznew)

                     else

c     --- calculate new vector from initial bead
                        do j = 1,iunit

                           rxnew(j) = rxnew(ifirst)
     &                          - (rxu(self,ifirst) - rxu(self,j))
                           rynew(j) = rynew(ifirst)
     &                          - (ryu(self,ifirst) - ryu(self,j))
                           rznew(j) = rznew(ifirst)
     &                          - (rzu(self,ifirst) - rzu(self,j))
                           enddo

                     endif      !nsampos.lt.3
                  endif         !ifirst.lt.riutry...
               endif            !nsampos.lt.iunit 
               call schedule(igrow,imolty,islen,ifirst,iprev,4)
ccc---END JLR 11-24-09
            else
               call schedule(igrow,imolty,islen,ifirst,iprev,3)
            endif !lrigid
         
c           --- Adding in multiple end regrowths

            if (ncut(iparty,s_type) .gt. 1) then
               do zz = 2,ncut(iparty,s_type)
                  ifirst = from(s_type+2*(zz-1))
                  iprev = prev(s_type+2*(zz-1))
                  
                  call schedule(igrow,imolty,islen,ifirst,iprev,5)
                  
               enddo
            endif

            waddnew = 1.0d0
c --- JLR 11-24-09 New stuff for rigid swatch
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there
                  !don't regrow anything in rosenbluth 
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1)) ) then
c                 rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
c                 rigid part is not grown, do rigrot
                  icallrose = 2
               endif
            else
c              flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            endif


            call rosenbluth(.true.,lterm,other,self,imolty,islen,
     &           iboxnew,igrow,waddnew,.false.,vdum2,icallrose)
c --- END JLR 11-24-09
c        --- termination of cbmc attempt due to walk termination ---
            if ( lterm ) return

c        --- propagate new rosenbluth weight
            tweight = tweight*weight*waddnew

c --- end rigid add on 

c     --- save the new coordinates
c     write(iou,*) 'new:',moltyp(self),iunit
            do jj = 1,igrow
               rxut(ic,jj) = rxnew(jj)
               ryut(ic,jj) = rynew(jj)
               rzut(ic,jj) = rznew(jj)
c     write(iou,*) rxnew(jj),rynew(jj),rznew(jj)
            enddo
c     write(iou,*) 'old:',moltyp(other),nunit(moltyp(other))
c     do jj = 1,nunit(moltyp(other))
c     write(iou,*) rxu(other,jj),ryu(other,jj),rzu(other,jj)
c     enddo

c     --- Corrections for switched beads, and DC-CBMC
c     --- Assign all of the grown new and old beads to rxuion
c     --- with rxuion: new = 2
            iii = 2
            do j=1,igrow
               rxuion(j,iii) = rxnew(j)
               ryuion(j,iii) = rynew(j)
               rzuion(j,iii) = rznew(j)
               qquion(j,iii) = qqu(self,j)
            enddo

            if ( iunit .ne. igrow ) then

c     --- for explicit-hydrogen model, put on the hydrogens
c     --- use phony number iins and call explct to add constrained hydrogens

               iins = nchain + 1
               moltyp(iins) = imolty
               do j=1,igrow
                  rxu(iins,j) = rxnew(j)
                  ryu(iins,j) = rynew(j)
                  rzu(iins,j) = rznew(j)
               enddo
               call explct(iins,vdum,.false.,.false.)
               do j = igrow + 1, iunit
                  rxuion(j,iii) = rxu(iins,j)
                  ryuion(j,iii) = ryu(iins,j)
                  rzuion(j,iii) = rzu(iins,j)
                  qquion(j,iii) = qqu(self,j)
                  rxut(ic,j) = rxu(iins,j)
                  ryut(ic,j) = ryu(iins,j)
                  rzut(ic,j) = rzu(iins,j)
c     write(iou,*) 'new:',rxut(ic,j),ryut(ic,j),rzut(ic,j)
               enddo
            endif

c     --- Begin DC-CBMC, explicit-hydrogen and
c     --- switched bead Corrections for NEW configuration
c     --- calculate the true site-site energy

c     ??? energy problem for cases which involve the change of the bending type
c     and torsional type for those units swatched!!!!!!

            call ctrmas(.false.,iboxnew,other,8)
            call energy (other,imolty, v, vintra,vinter,vext
     &           ,velect,vewald,iii,iboxnew,1,iunit,.true.,lterm,.false.
     &           ,vdum,.false.,.false.)
            
            if (lterm) then
c     write(iou,*) 'other ',other,' self ',self
c     stop 'interesting screwup in CBMC swatch'
               return
            endif
            delen = v - ( vnewt - (vnewbvib + vnewbb + vnewtg )) 
            tweight = tweight*dexp(-beta*delen)

            vnewt     = vnewt + delen
            vnewinter = vinter 
            vnewintra = vintra
            vnewext   = vext 
            vnewelect = velect
            vnewewald = vewald
            
c     End DC-CBMC and switched bead Corrections for NEW configuration

c     --- save the trial energies
            vnbox(iboxnew)   = vnbox(iboxnew) + vnewt
            vninte(iboxnew)  = vninte(iboxnew) + vnewinter
            vnintr(iboxnew)  = vnintr(iboxnew) + vnewintra
            vnvibb(iboxnew)  = vnvibb(iboxnew) + vnewbvib
            vntgb(iboxnew)   = vntgb(iboxnew) + vnewtg
            vnextb(iboxnew)  = vnextb(iboxnew) + vnewext
            vnbend(iboxnew)  = vnbend(iboxnew) + vnewbb
            vnelect(iboxnew) = vnelect(iboxnew) + vnewelect
            vnewald(iboxnew) = vnewald(iboxnew) + vnewewald 
c --- rigid add on

            waddold = 1.0d0

c        --- grow old chain conformation
c --- JLR 11-24-09 New stuff for rigid swatch
            if (lrigid(imolty)) then

               if(nsampos(iparty).ge.iunit) then
                  !molecule is all there
                  !don't regrow anything in rosenbluth
                  icallrose = 4
               elseif( (nsampos(iparty).ge.3)      .or.
     &              (ifirst.ge.riutry(imolty,1) )) then
c                 rigid part is grown, don't do rigrot in rosebluth
                  icallrose = 3
               else
c                 rigid part is not grown, do rigrot
                  icallrose = 2
               endif
            else
c              flexible molecule call rosenbluth in normal fashion
               icallrose = 2
            endif

            call rosenbluth(.false.,lterm,self,self,imolty,
     &           islen,iboxold,igrow,waddold,.false.,vdum2,icallrose)
c --- END JLR 11-24-09

c        --- termination of old walk due to problems generating orientations
            if ( lterm ) then
               write(iou,*) 'SWATCH: old growth rejected'
               return
            endif

c        --- propagate old rosenbluth weight
            tweiold = tweiold*weiold*waddold

c     --- store the old grown beads and explict placed beads positions
c     --- 1 = old conformation
            iii = 1
            do j = 1,iunit
               rxuion(j,1) = rxu(self,j)
               ryuion(j,1) = ryu(self,j)
               rzuion(j,1) = rzu(self,j)
               qquion(j,1) = qqu(self,j)
            enddo

c     Begin Correction for DC-CBMC and switched beads for OLD configuration
c     --- correct the acceptance rules 
c     --- calculate the Full rcut site-site energy

            call energy (self,imolty, v, vintra,vinter,vext,velect
     &           ,vewald,iii,iboxold, 1,iunit,.true.,lterm,.false.,vdum
     &           ,.false.,.false.)

            if (lterm) stop 'disaster ovrlap in old conf SWATCH'
            deleo = v - ( voldt - (voldbvib + voldbb + voldtg) ) 
            
            tweiold = tweiold*dexp(-beta*deleo)

            voldt     = voldt + deleo
            voldintra = vintra
            voldinter = vinter 
            voldext   = vext 
            voldelect = velect
            voldewald = vewald

c     End Correction for DC-CBMC and switched beads for OLD configuration

c     --- save the trial energies
            vnbox(iboxold)   = vnbox(iboxold)  - voldt
            vninte(iboxold)  = vninte(iboxold) - voldinter
            vnintr(iboxold)  = vnintr(iboxold) - voldintra
            vnvibb(iboxold)  = vnvibb(iboxold) - voldbvib 
            vntgb(iboxold)   = vntgb(iboxold)  - voldtg 
            vnextb(iboxold)  = vnextb(iboxold) - voldext
            vnbend(iboxold)  = vnbend(iboxold) - voldbb 
            vnelect(iboxold) = vnelect(iboxold) - voldelect
            vnewald(iboxold) = vnewald(iboxold) - voldewald

         enddo

c     --- Perform the Ewald sum reciprical space corrections
         if ( lewald ) then
c     --- added into tweight even though it really contains new-old
c     --- Box A
            oldchain = iboxa
            newchain = iboxb
            oldunit = nunit(imolta)
            newunit = nunit(imoltb)

            do j = 1,oldunit
               rxuion(j,1) = rxu(oldchain,j)
               ryuion(j,1) = ryu(oldchain,j)
               rzuion(j,1) = rzu(oldchain,j)
               qquion(j,1) = qqu(oldchain,j)
            enddo
            moltion(1) = imolta
            do j = 1,newunit
               rxuion(j,2) = rxut(2,j)
               ryuion(j,2) = ryut(2,j)
               rzuion(j,2) = rzut(2,j)
               qquion(j,2) = qqu(newchain,j)
            enddo
            moltion(2) = imoltb

            call recip(iboxia,vrecipn,vrecipo,1)

            delen = vrecipn - vrecipo 

            tweight = tweight * dexp(-beta*delen)

            vnewald(iboxia) = vnewald(iboxia) + delen
            vnbox(iboxia) = vnbox(iboxia) + delen

c     --- Box B
            oldchain = iboxb
            newchain = iboxa
            oldunit = nunit(imoltb)
            newunit = nunit(imolta)

            do j = 1,oldunit
               rxuion(j,1) = rxu(oldchain,j)
               ryuion(j,1) = ryu(oldchain,j)
               rzuion(j,1) = rzu(oldchain,j)
               qquion(j,1) = qqu(oldchain,j)
            enddo
            moltion(1) = imoltb
            do j = 1,newunit
               rxuion(j,2) = rxut(1,j)
               ryuion(j,2) = ryut(1,j)
               rzuion(j,2) = rzut(1,j)
               qquion(j,2) = qqu(newchain,j)
            enddo
            moltion(2) = imolta

            call recip(iboxib,vrecipn,vrecipo,1)

            delen = vrecipn - vrecipo 
            tweight = tweight * dexp(-beta*delen)

            vnewald(iboxib) = vnewald(iboxib) + delen
            vnbox(iboxib) = vnbox(iboxib) + delen

         endif
c     End Ewald-sum Corrections

C     ----------------------------------------------------------------------

         if (ltailc) then
c     ---    add tail corrections
            if (lpbcz) then
               if (lsolid(boxa) .and. .not. lrect(boxa)) then
                  vola = (hmat(boxa,1) * (hmat(boxa,5) * hmat(boxa,9) -
     +                 hmat(boxa,8)*hmat(boxa,6))+
     &                 hmat(boxa,4)*(hmat(boxa,8)
     +                 * hmat(boxa,3)-hmat(boxa,2)*
     &                 hmat(boxa,9))+hmat(boxa,7)
     +                 * (hmat(boxa,2)*hmat(boxa,6)-
     &                 hmat(boxa,5)*hmat(boxa,3)))
               else
                  vola=boxlx(boxa)*boxly(boxa)*boxlz(boxa)
               endif

               if (lsolid(boxb) .and. .not. lrect(boxb)) then
                  volb = (hmat(boxb,1) * (hmat(boxb,5) * hmat(boxb,9) -
     +                 hmat(boxb,8)*hmat(boxb,6))+
     &                 hmat(boxb,4)*(hmat(boxb,8)
     +                 * hmat(boxb,3)-hmat(boxb,2)*
     &                 hmat(boxb,9))+hmat(boxb,7)
     +                 * (hmat(boxb,2)*hmat(boxb,6)-
     &                 hmat(boxb,5)*hmat(boxb,3)))
               else
                  volb=boxlx(boxb)*boxly(boxb)*boxlz(boxb)
               endif
            else
               vola=boxlx(boxa)*boxly(boxa)
               volb=boxlx(boxb)*boxly(boxb)
            endif

c     - for new BOXINS with inserted particle
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
               endif
c --- JLR 11-24-09 don't do tail corrections for ideal box
               if (.not.lideal(ibox)) then 
                  do imt = 1, nmolty
                     do jmt = 1, nmolty
c     --- new logic for tail correction (same answer) MGM 3-25-98
                        rho = dble( ncmt(ibox,jmt) )
                        if ( jmt .eq. imolin ) rho = rho + 1.0d0
                        if ( jmt .eq. imolrm ) rho = rho - 1.0d0

                        rho = rho / dvol

                        dicount = ncmt(ibox,imt)
                        if ( imt .eq. imolin ) dicount = dicount + 1
                        if ( imt .eq. imolrm ) dicount = dicount - 1

                        dinsta = dinsta + 
     &                       dicount * coru(imt,jmt,rho,ibox)
                     enddo
                  enddo
               else
                  dinsta = 0.0d0
               endif
c --- END JLR 11-24-09              
               dinsta = dinsta - vtailb( ibox )

               tweight=tweight*dexp(-beta*dinsta)

               vntail(ibox) = dinsta 
            enddo
         else
            vntail(boxa) = 0.0d0
            vntail(boxb) = 0.0d0
         endif

         wnlog = dlog10( tweight )
         wolog = dlog10( tweiold )
         wdlog = wnlog - wolog
         if ( wdlog .lt. -softcut ) then
c     write(iou,*) '### underflow in wratio calculation ###'
            return
         endif
         wswat = ( tweight / tweiold ) * ( dble(orgaia*orgbib) /
     &        dble((orgbia+1)*(orgaib+1)) ) *
     &        dexp(beta*(eta2(boxa,imolta)
     &        +eta2(boxb,imoltb)-eta2(boxa,imoltb)
     &        -eta2(boxb,imolta)))

c     write(iou,*) 'imolta,imoltb',imolta,imoltb
c     write(iou,*) 'wswat,tweight,tweiold',wswat,tweight,tweiold
         if ( random() .le. wswat ) then
c     *** we can now accept !!!!! ***
            bsswat(iparty,ipairb) = bsswat(iparty,ipairb) + 1.0d0
c     write(iou,*) 'SWATCH ACCEPTED',iboxa,iboxb
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
            enddo

c     ---update book keeping
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
            enddo
            do ic = 1,iunitb
               rxu(iboxb,ic) = rxut(2,ic)
               ryu(iboxb,ic) = ryut(2,ic)
               rzu(iboxb,ic) = rzut(2,ic)
            enddo
            if ( lewald ) then
c     -- update reciprocal-space sum
               call recip(boxa,vdum,vdum,2)
               call recip(boxb,vdum,vdum,2)
            endif

c     ---update center of mass
            call ctrmas(.false.,boxa,iboxb,8)
            call ctrmas(.false.,boxb,iboxa,8)

            if (licell .and. (boxa .eq. boxlink .or. boxb .eq. boxlink))
     &           stop 'not yet implemented!'

c     --- call nearest neighbor list
            if ( lneigh ) then
               call updnn( iboxa )
               call updnn( iboxb )
            endif

ccc--!!!JLR - for test 2 write final coordinates of a and b                                      
c         open(unit=93,file='a_final.xyz',status='unknown')                                   
c         write(93,*) iunita                                                                  
c         write(93,*)                                                                         
c         do zz = 1,iunita                                                                    
c            write(93,*) 'C ', rxu(iboxa,zz),ryu(iboxa,zz),rzu(iboxa,zz)                      
c         enddo                                                                               
c                                                                                             
c         open(unit=93,file='b_final.xyz',status='unknown')                                   
c         write(93,*) iunitb                                                                  
c         write(93,*)                                                                         
c         do zz = 1,iunitb                                                                    
c            write(93,*) 'C ', rxu(iboxb,zz),ryu(iboxb,zz),rzu(iboxb,zz)                      
c         enddo                                                                               
c         STOP 'END OF SWATCH TEST'                                                           
ccc--!!!JLR - end of coordinate test         

         endif
c -----------------------------------------------------------------

      endif

c      write(iou,*) 'end SWATCH'

      return
      end




