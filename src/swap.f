      subroutine swap(bsswap,bnswap,bnswap_in,bnswap_out,cnt_wf1,
     &     cnt_wf2,cnt_wra1,cnt_wra2,qelect)
      
c swap
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

c    ********************************************************************
c    ** removes a molecule from one box and inserts it into the other  **
c    ** using CBMC insertion techniques.  Works for linear or branched **
c    ** molecules and for DC-CBMC and Explicit atom                    **
c    ** Rewritten from old swap and swapbr subroutines by M.G. Martin  **
c    ** 9-18-97                                                        ** 
c    ********************************************************************
 
      implicit none
      
      include 'control.inc'
      include 'coord.inc'
      include 'conver.inc'
      include 'coord2.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'cbmc.inc'
      include 'rosen.inc' 
      include 'boltzmann.inc'
      include 'external.inc'
      include 'connect.inc'
      include 'inputdata.inc'
      include 'ewaldsum.inc'
      include 'poten.inc'
      include 'fepsi.inc'
      include 'clusterbias.inc'
      include 'neigh.inc'
      include 'cell.inc'
      include 'eepar.inc'
      
      logical ovrlap,lterm,lnew,lempty,ldone,ltors,lovrh,lfavor,
     &     laccept,lswapinter,lrem_out,lins_in,lneighij,linsk_in,
     &     lremk_in,lrem_clu,lins_clu,lfixnow

     
      integer boxins,boxrem,imol,ichoi,ip,iwalk,idum,iins1,imolty1
      integer istt,iett,ncount,itype,ipair,ipairb,beg,flagon

      integer iutry,icbu,ifrom,irem,iins,glist,findex
     &     ,iii,j,ibox,iunit,ic,pointp,imolty,imt,jmt,igrow
     &     ,pointp2,jins,jmolty,neighj_num,neighk_num
     &     ,joffset,koffset,kmolty,kins,target,cnt_wf1
     &     ,cnt_wf2,neigh_old,cnt_wra1,cnt_wra2,k

      dimension glist(numax),cnt_wf1(0:6,0:6,4),cnt_wf2(0:6,0:6,4),
     &     cnt_wra1(1000,4),cnt_wra2(1000,4)

      double precision sx,sy,sz,ddum(27)

      double precision v,vintra,vinter,vext,velect,vtorold,vtornew
     &     ,vewald,vflucq,delen,deleo,rpair
      double precision vnewflucq,voldflucq,qion,ctorfo,ctorfn
      dimension qion(numax)
      double precision rxuold,ryuold,rzuold
      dimension rxuold(numax),ryuold(numax),rzuold(numax)
      double precision bsswap,bnswap,bnswap_in,bnswap_out
      double precision random,rmol,qelect,rbf,bsum
      double precision waddnew,waddold

      double precision total_NBE,vintran,velectn,vewaldn,vtgn
      double precision vbendn,vvibn 



      double precision v1insext,v1remext,v1ins,w1ins,v1rem,w1rem
     &     ,v1insint,v1remint,v1insewd,v1remewd
     &     ,wnlog,wolog,wdlog,wratio,vinsta,vremta
     &     ,volins,volrem,rho,arg,coru,v1inselc,v1remelc
      double precision rvol,x,y,z,rijsq,wbias_ins,wbias_rem,r
     &     ,xi1,xi2,xisq
      dimension bsswap(ntmax,npabmax,nbxmax*2),
     &     bnswap(ntmax,npabmax,nbxmax*2),bnswap_in(ntmax,2),
     &     bnswap_out(ntmax,2)
      dimension qelect(nntype)
      double precision vrecipn,vrecipo,vdum,whins,whrem
      double precision rxuh,ryuh,rzuh,delenh,vtrhext,vtrhintra
     &     ,vtrhinter,vtrhelect,vtrhewald,vtrhtg,bfach
     
      dimension bfach(nchmax),delenh(nchmax),vtrhinter(nchmax)
     &     ,vtrhext(nchmax),vtrhintra(nchmax),vtrhelect(nchmax)
     &     ,vtrhewald(nchmax),vtrhtg(nchmax)
     
      dimension lovrh(nchmax)
      dimension rxuh(numax,nchmax),ryuh(numax,nchmax)
     +     ,rzuh(numax,nchmax)

C --------------------------------------------------------------------

c      write(2,*) 'START SWAP'
c      write(11,*) '1:',neigh_cnt(18)

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
         endif
c       write(6,*) 'ee val', imolty, irem, pointp, boxrem, boxins
         goto 300
      endif
                                                                                
      wee_ratio = 1.0d0
   
c --- select a molecule typ with probabilities given in pmswmt
      rmol = random()
      ldone = .false.
      do imol = 1, nmolty
         if ( rmol .lt. pmswmt(imol) ) then
            if ( .not. ldone ) then
               imolty = imol
               ldone = .true.
            endif
         endif
      enddo

      if (temtyp(imolty).eq.0) return
      imolty1 = imolty
         
c ---    select a box given in pmswatyp
      if ( nswapb(imolty) .gt. 1 ) then
         rpair = random()
         do 96 ipair = 1, nswapb(imolty)
            if ( rpair .lt. pmswapb(imolty,ipair) ) then
               ipairb = ipair
               goto 97
            endif
 96      continue
      else
         ipairb = 1
      endif

 97   if (random().lt.0.5d0) then
         boxins=box1(imolty,ipairb)
         boxrem=box2(imolty,ipairb)
      else
         boxins=box2(imolty,ipairb)
         boxrem=box1(imolty,ipairb)
      endif
      if ( boxins .eq. boxrem ) then
         lswapinter = .false.
      else
         lswapinter = .true.
      endif
c      write(2,*) 'boxins:',boxins,'boxrem:',boxrem
      if ( .not. (lgibbs .or. lgrand) .and. lswapinter ) 
     &     stop 'no interbox swap if not gibbs/grand ensemble!'
         
c *** select a chain in BOXREM at random ***

      if ( ncmt(boxrem,imolty) .eq. 0 ) then
         lempty = .true.
         if ( .not. lswapinter .or. lrigid(imolty) ) return
      elseif ( lswapinter .or. lavbmc1(imolty) .and. .not.  
     &        (lavbmc2(imolty) .or. lavbmc3(imolty)) ) then

c *** for the advanced AVBMC algorithm, this particle will be selected in
c *** sub-regions defined by Vin

 197     pointp = idint( dble(ncmt(boxrem,imolty))*random() ) + 1
         if (lexpee) then
            if ((pointp.eq.eepointp).and.
     &               (boxrem.eq.box_state(mstate)).and.
     &               (ncmt(boxrem,imolty).gt.1)) then
               goto 197
            else
               return
            endif
         endif

         irem = parbox(pointp,boxrem,imolty)
         if ( moltyp(irem) .ne. imolty ) 
     &        write(2,*) 'screwup swap, irem:',irem,moltyp(irem),imolty
         ibox = nboxi(irem)
         if ( ibox .ne. boxrem ) stop 'problem in swap'
      endif

c$$$      write(2,*) 'particle ',irem,' is being removed, imolty is:',
c$$$     &     imolty,' and the box is:',boxrem


c ===>  for both gibbs and grand-canonical we have:
c --- insert a chain in box: boxins 
c --- remove one in box: boxrem
c      write(2,*) 'boxrem',boxrem,' imolty',imolty,' lempty',lempty
c     bnswap(imolty,X) decoder X = 1 is # attempts into box 1
c      X = 2 is # attempts into box 2 X=3 is success into box 1
c      X = 4 is success into box 2
c     bsswap is the same thing but keeps track of successful growths

      if (.not. lempty) bnswap(imolty,ipairb,boxins) 
     &     = bnswap(imolty,ipairb,boxins) + 1.0d0
      bsswap(imolty,ipairb,boxins) = bsswap(imolty,ipairb,boxins)+1.0d0

 300  continue
      
c     *** store number of units in iunit ***
      iunit = nunit(imolty)
      igrow = nugrow(imolty)
c     *** give i a phony number ***
      if ( lswapinter ) then
         iins = nchain + 1
         iins1 = iins
         moltyp(iins) = imolty
c     give charges to phony number 
         if ((.not.lexpee).and.(.not.leemove)) then
            if ( lempty ) then
               do icbu = 1, iunit
                  qqu(iins,icbu) = qelect(ntype(imolty,icbu))
               enddo
            else
               do icbu = 1, iunit
                  qqu(iins,icbu) = qqu(irem,icbu)
               enddo
            endif
         else
            do icbu = 1, iunit
               qqu(iins,icbu) = ee_qqu(icbu,nstate)
            enddo
         endif
      elseif ( lavbmc2(imolty) .or. lavbmc3(imolty) ) then
         iins = 0
      else
         iins = irem
         if (lexpee.and.leemove) then
            iins1 = iins
            do icbu = 1, iunit
               qqu(iins,icbu) = ee_qqu(icbu,nstate)
            enddo
         endif
      endif

c *** select a position of the first/starting unit at RANDOM ***
c *** and calculate the boltzmann weight                     ***
c *** for the chain to be INSERTED                           ***


      if (lrigid(imolty)) then
         beg = riutry(imolty,1)
      else
         beg = 1
      endif

      wbias_ins = 1.0d0
      wbias_rem = 1.0d0

      ichoi = nchoi1(imolty)

      if ( .not. lswapinter .and. lbias(imolty)) then

         if (boxins .ne. boxrem) then
            write(2,*) 'avbmc, boxins, boxrem:',boxins,boxrem
            stop
         endif
         
c *********************
c * First AVBMC Steps *
c *********************
         rmol = random()
         ldone = .false.
         do imol = 1, nmolty
            if ( rmol .lt. pmbsmt(imol) ) then
               if ( .not. ldone ) then
                  jmolty = imol
                  ldone = .true.
               endif
            endif
         enddo
         
         if ( ncmt(boxins,jmolty) .eq. 0 .or. (jmolty .eq. imolty
     &        .and. ncmt(boxins,jmolty) .eq. 1) ) then
               
c ??? what shall we do now?
            lempty = .true.
            return
            
         else
 111        pointp2=idint( dble(ncmt(boxins,jmolty))*random())+1
            jins = parbox(pointp2,boxins,jmolty)
            if ( jins .eq. iins ) goto 111
            if ( moltyp(jins) .ne. jmolty ) 
     &           write(2,*) 'screwup swap, jins:'
     &           ,jins,moltyp(jins),jmolty
            if ( nboxi(jins) .ne. boxins ) 
     &           stop 'problem in swap with jins'
         endif

         if ( lavbmc3(imolty) ) then

c *** define a second bonding region bounded by kins with kmolty type
            
            rmol = random()
            ldone = .false.
            do imol = 1, nmolty
               if ( rmol .lt. pmbsmt(imol) ) then
                  if ( .not. ldone ) then
                     kmolty = imol
                     ldone = .true.
                  endif
               endif
            enddo
            if ( ncmt(boxins,kmolty) .eq. 0 .or. (kmolty .eq. imolty
     &           .and. jmolty .eq. imolty .and. 
     &           ncmt(boxins,kmolty) .eq. 2) .or. (kmolty .eq. imolty
     &           .and. ncmt(boxins,kmolty) .eq. 1) ) then
               
c ??? what shall we do now?
               lempty = .true.
               return
            else
 112           pointp2=idint( dble(ncmt(boxins,kmolty))*random())+1
               kins = parbox(pointp2,boxins,kmolty)

c *** make sure the two regions bounded by jins and kins do not 
c *** overlap (sampling from cluster to cluster)
c *** ??? Problems: potentially infinite loop if there is only
c *** one region in the whole system

               if ( jins .eq. kins ) goto 112
c               do ip = 1,neighj_num
c                  do ic = 1,neighk_num
c                     if ( neighbor(ip,jins) .eq. kins .or.
c     &                    neighbor(ic,kins) .eq. jins .or.
c     &                    neighbor(ip,jins) .eq. neighbor(ic,kins))
c     &                    goto 112
c                  enddo
c               enddo
               x = rxu(jins,1) - rxu(kins,1)
               y = ryu(jins,1) - ryu(kins,1)
               z = rzu(jins,1) - rzu(kins,1)
               if ( lpbc ) call mimage(x,y,z,boxins)
               rijsq = x*x + y*y + z*z
               if ( rijsq .lt. (2.0d0*rbsmax)**2 ) goto 112

               if ( moltyp(kins) .ne. kmolty ) 
     &              write(2,*) 'screwup swap, kins:'
     &              ,kins,moltyp(kins),kmolty
               if ( nboxi(kins) .ne. boxins ) 
     &              stop 'problem in swap with kins'
            endif
            neighk_num = neigh_cnt(kins)
                  
         endif
            

         if ( lavbmc2(imolty) .or. lavbmc3(imolty) ) then
c            neighj_num = neigh_cnt(jins,imolty)
            neighj_num = neigh_cnt(jins)
            if (jmolty .eq. imolty) then 
               joffset = 1
            else
               joffset = 0
            endif
            if (lavbmc3(imolty) .and. kmolty .eq. imolty) then
               koffset = 1
            else
               koffset = 0 
            endif
         endif

         if ( lpbc ) call setpbc(boxins)
         
         if ( random() .lt. pmbias(imolty) ) then
            
            wbias_rem = pmbias(imolty) * (1.0d0/vol_eff)

            lins_in = .true.
            
            if (lavbmc2(imolty) .or. (lavbmc3(imolty) 
     &           .and. random() .gt. pmbias2(imolty)) ) then
              
c *** select a particle in the out region and move this particle into
c *** the in region defined by the particle jins's bonding region
c *** out -> j case

               wbias_rem=wbias_rem
     &              /dble(ncmt(boxrem,imolty)-neighj_num-joffset) 

               lrem_out = .true.
 119           pointp=idint(dble(ncmt(boxrem,imolty))*random())+1
               irem = parbox(pointp,boxrem,imolty)
               if ( irem .eq. jins ) goto 119
               x = rxu(irem,1) - rxu(jins,1)
               y = ryu(irem,1) - ryu(jins,1)
               z = rzu(irem,1) - rzu(jins,1)
               if ( lpbc ) call mimage(x,y,z,boxins)
               rijsq = x*x + y*y + z*z
               if ( rijsq .lt. rbsmax**2 .and. rijsq .gt. rbsmin**2) 
     &              goto 119

               if ( moltyp(irem) .ne. imolty ) 
     &              write(2,*) 'screwup swap1, irem:',irem,
     &              moltyp(irem),imolty
               ibox = nboxi(irem)
               if ( ibox .ne. boxrem ) stop 'problem in swap'
               iins = irem
               lremk_in = .false.
   
            elseif ( lavbmc3(imolty) ) then
               
c *** select a particle in the region bounded by kins and move this particle
c *** into the in region defined by the particle jins's bonding region
c *** k -> j case
               
               lrem_out = .false.
               lremk_in = .true.
   
               if ( neighk_num .eq. 0 .or. neighk_num .eq. 1 .and. 
     &              neighbor(1,kins) .eq. jins) then
c ??? what shall we do now?
                  lempty = .true.
                  return
               else
 113              pointp=idint(dble(neighk_num)*random())+1
c                  irem = neighbor(pointp,kins,imolty)
c                  write(2,*) 'kins,irem:',kins,irem,neighk_num
                   irem = neighbor(pointp,kins)
                   if ( irem .eq. jins ) goto 113
                   iins = irem
                endif
                wbias_rem=wbias_rem/dble(neighk_num) 


            endif               
            
c            write(2,*) 'move in'
c            write(2,*) '3:',wbias_rem,boxlx(boxins),vol_eff
            
            do icbu = 1,ichoi
c *** choose a random association distance
               rvol = random()
               r = (rbsmax*rbsmax*rbsmax*rvol + 
     &              (1.0d0-rvol)*rbsmin*rbsmin*rbsmin)**(1.0d0/3.0d0)
         
c --- calculate random vector on the unit sphere ---
 109           xi1 = ( 2.0d0 * random() ) - 1.0d0
               xi2 = ( 2.0d0 * random() ) - 1.0d0
               xisq = xi1**2 + xi2**2
               if ( xisq .lt. 1.0d0 ) then
                  x = r * 2.0d0 * xi1 * dsqrt( 1.0d0 - xisq )
                  y = r * 2.0d0 * xi2 * dsqrt( 1.0d0 - xisq )
                  z = r * ( 1.0d0 - 2.0d0 * xisq )
               else
                  goto 109
               endif
               rxp(1,icbu) = rxu(jins,1) + x
               ryp(1,icbu) = ryu(jins,1) + y
               rzp(1,icbu) = rzu(jins,1) + z
            enddo
               
         else

            lins_in = .false.
               
            if ( lavbmc2(imolty) .or. lavbmc3(imolty) ) then

c *** select a particle in the in region defined by the particle 
c *** jins's bonding region through neighbor list 
c *** and move to the out region or the region bounded by kins in AVBMC3
               if ( neighj_num .eq. 0 .or. (neighj_num .eq. 1
     &              .and. neighbor(1,jins) .eq. kins) ) then

c ??? what shall we do now?
                  lempty = .true.
                  return
               else
 114              pointp=idint(dble(neighj_num)*random())+1
c                  irem = neighbor(pointp,jins,imolty)
c                  write(2,*) 'jins,irem:',jins,irem,neighj_num
                  irem = neighbor(pointp,jins)
                  if ( irem .eq. kins ) goto 114
                  if ( moltyp(irem) .ne. imolty ) 
     &                 write(2,*) 'screwup swap2, irem:',irem,
     &                 moltyp(irem),imolty,neighj_num,pointp,jins
                  ibox = nboxi(irem)
                  if ( ibox .ne. boxrem ) stop 'problem in swap'
                  iins = irem
               endif
               lrem_out = .false.
               lremk_in = .false.

            endif
c            write(2,*) 'move out'
c            write(2,*) '4:',wbias_rem,boxlx(boxins),vol_eff
            

            if ( lavbmc3(imolty) .and. 
     &           random() .lt. pmbias2(imolty) ) then

c *** move to the region bounded by kins in AVBMC3
c *** j -> k case

               linsk_in = .true.
               wbias_rem = (1-pmbias(imolty))/vol_eff/dble(neighj_num)

               do icbu = 1,ichoi
c *** choose a random association distance
                  rvol = random()
                  r = (rbsmax*rbsmax*rbsmax*rvol + 
     &                 (1-rvol)*rbsmin*rbsmin*rbsmin)**(1.0/3.0d0)
         
c --- calculate random vector on the unit sphere ---
 139              xi1 = ( 2.0d0 * random() ) - 1.0d0
                  xi2 = ( 2.0d0 * random() ) - 1.0d0
                  xisq = xi1**2 + xi2**2
                  if ( xisq .lt. 1.0d0 ) then
                     x = r * 2.0d0 * xi1 * dsqrt( 1.0d0 - xisq )
                     y = r * 2.0d0 * xi2 * dsqrt( 1.0d0 - xisq )
                     z = r * ( 1.0d0 - 2.0d0 * xisq )
                  else
                     goto 139
                  endif
                  rxp(1,icbu) = rxu(kins,1) + x
                  ryp(1,icbu) = ryu(kins,1) + y
                  rzp(1,icbu) = rzu(kins,1) + z
               enddo

            else

               linsk_in = .false.
c *** move to the out region
c *** j -> out case

               wbias_rem = (1-pmbias(imolty))/
     &              (boxlx(boxins)*boxly(boxins)*
     &              boxlz(boxins)-vol_eff)
               if ( lavbmc2(imolty) .or. lavbmc3(imolty) ) 
     &              wbias_rem = wbias_rem/dble(neighj_num)            

               do icbu = 1,ichoi
 222              rxp(1,icbu) = boxlx(boxins) * random()
                  ryp(1,icbu) = boxly(boxins) * random()
                  if (lpbcz .or. lslit) then
                     rzp(1,icbu) = boxlz(boxins) * random()
                  else if ( lsami .or. lmuir .or. ljoe ) then
                     if ( lempty ) then
                        rzp(1,icbu) = 20*random()-10
                     else
                        rzp(1,icbu) = rzu(irem,1)
                     endif
                  else
                     rzp(1,icbu) = 0.0d0
                  endif
            
c *** determine whether it is inside the sphere with particle jins
                  x = rxp(1,icbu) - rxu(jins,1)
                  y = ryp(1,icbu) - ryu(jins,1)
                  z = rzp(1,icbu) - rzu(jins,1)
                  if ( lpbc ) call mimage(x,y,z,boxins)
                  rijsq = x*x + y*y + z*z
                  if (rijsq .lt. rbsmax**2 .and. rijsq .gt. rbsmin**2) 
     &                 goto 222
               enddo
            endif
         endif

c *************************
c * end first Avbmc calcs *
c *************************

      else

         if (lsolid(boxins) .and. .not. lrect(boxins)) then
            ibox = boxins
            do icbu = 1,ichoi

c$$$               if (lbnow.and.boxins.eq.boxsite) then
c$$$                  x = random()
c$$$                  y = random()
c$$$                  z = random()
c$$$
c$$$c     --- gaussian picking probabilities
c$$$c                  x = adev * dsqrt(- 2.0d0 * dlog(x*adev*twopiroot))
c$$$c                  y = adev * dsqrt(- 2.0d0 * dlog(y*adev*twopiroot))
c$$$c                  z = adev * dsqrt(- 2.0d0 * dlog(z*adev*twopiroot))
c$$$
c$$$c     --- old picking probabilities
c$$$                  rvol = random()
c$$$                  r = (rbsmax2**3*rvol +
c$$$     &                 (1-rvol)*rbsmin*rbsmin*rbsmin)**(1.0/3.0d0)
c$$$
c$$$c --- calculate random vector on the unit sphere ---
c$$$ 129              xi1 = ( 2.0d0 * random() ) - 1.0d0
c$$$                  xi2 = ( 2.0d0 * random() ) - 1.0d0
c$$$                  xisq = xi1**2 + xi2**2
c$$$                  if ( xisq .lt. 1.0d0 ) then
c$$$                     x = r * 2.0d0 * xi1 * dsqrt( 1.0d0 - xisq )
c$$$                     y = r * 2.0d0 * xi2 * dsqrt( 1.0d0 - xisq )
c$$$                     z = r * ( 1.0d0 - 2.0d0 * xisq )
c$$$                  else
c$$$                     goto 129
c$$$                  endif
c$$$
c$$$                  rxp(1,icbu) = rxu(jins,1) + x
c$$$                  ryp(1,icbu) = ryu(jins,1) + y
c$$$                  rzp(1,icbu) = rzu(jins,1) + z
c$$$
c$$$               else

                  sx = random()
                  sy = random()
                  sz = random()

                  rxp(1,icbu) = sx*hmat(ibox,1)+sy*hmat(ibox,4)
     &                 +sz*hmat(ibox,7)
                  ryp(1,icbu) = sx*hmat(ibox,2)+sy*hmat(ibox,5)
     &                 +sz*hmat(ibox,8)
                  rzp(1,icbu) = sx*hmat(ibox,3)+sy*hmat(ibox,6)
     &                 +sz*hmat(ibox,9)

c$$$               endif
            enddo

         else

            do icbu = 1,ichoi
               rxp(1,icbu) = boxlx(boxins) * random()
               ryp(1,icbu) = boxly(boxins) * random()
               if (lpbcz .or. lslit) then
                  rzp(1,icbu) = boxlz(boxins) * random()
               else if ( lsami .or. lmuir .or. ljoe ) then
                  if ( lempty ) then
                     rzp(1,icbu) = 20*random()-10
                  else
                     rzp(1,icbu) = rzu(irem,1)
                  endif
               else
                  rzp(1,icbu) = 0.0d0
               endif
            enddo
         endif

      endif

c     *** select starting unit ***
c     --- always using bead 1 as the starting unit
      iutry = beg
      lnew = .true.
      glist(1) = beg
c     --  insert the first atom

      call boltz(lnew,.true.,ovrlap,iins,iins,imolty,boxins,
     &     ichoi,idum,1,glist,0.0d0)

      if ( lanes .and. lflucq(imolty) ) then
c *** lfavor is used to set up the favor and favor2 to preferentially sample 
c *** the electronic degrees of freedom
         lfavor = .true.
      else
         lfavor = .false.
         if ( lswapinter )
     &        bnchem(boxins,imolty) = bnchem(boxins,imolty) + 1.0d0
      endif

      if ( ovrlap ) return

      
c     *** perform the walk according to the availibility of the choices ***
c     *** and calculate the correct weight for the trial walk           ***
      
      w1ins = 0.0d0
      do ip = 1, ichoi
         w1ins = w1ins + bfac(ip)
      enddo
      
c     --- check for termination of walk ---
      if ( w1ins .lt. softlog ) then
         write(2,*) 'caught in swap'
         return
      endif
      

c     --- select one position at random ---
      if ( ichoi .gt. 1 ) then
         rbf = w1ins * random()
         bsum = 0.0d0 
         do ip = 1, ichoi
            if ( .not. lovr(ip) ) then
               bsum = bsum + bfac(ip)
               if ( rbf .lt. bsum ) then
c     --- select ip position ---
                  iwalk = ip
                  goto 180
               endif
            endif
         enddo
         write(2,*) 'w1ins:',w1ins,'rbf:',rbf
         stop 'big time screwup -- w1ins'
      else
         iwalk = 1
      endif
      
 180  v1ins =  vtry(iwalk)  
      v1insext = vtrext(iwalk)
      v1insint = vtrinter(iwalk)
      v1inselc = vtrelect(iwalk)
      v1insewd = vtrewald(iwalk)

c      neigh_icnt = ntr_icnt(iwalk)
c      do ip = 1,neigh_icnt
c         neighi(ip) = ntr(ip,iwalk)
c      enddo
      rxnew(beg) = rxp(1,iwalk)
      rynew(beg) = ryp(1,iwalk)
      rznew(beg) = rzp(1,iwalk)

      if (lrigid(imolty)) then
c     --- calculate new vector from initial bead

         do j = beg,iunit
            rxnew(j) = rxnew(beg) 
     &           - (rxu(irem,beg) - rxu(irem,j))
            rynew(j) = rynew(beg) 
     &           - (ryu(irem,beg) - ryu(irem,j))
            rznew(j) = rznew(beg) 
     &           - (rzu(irem,beg) - rzu(irem,j))
            if ((.not.lexpee).and.(.not.leemove)) then
               qqu(iins,j) = qqu(irem,j)
            else
               qqu(iins,j) = ee_qqu(j,nstate)
            endif
         enddo
         call schedule(igrow,imolty,ifrom,iutry,0,4)
      elseif (lring(imolty)) then
         lfixnow = .true.
         call safeschedule(igrow,imolty,ifrom,iutry,findex,2)
      else
         call schedule(igrow,imolty,ifrom,iutry,0,2)
      endif
      

c     ------------------------------------------------------------------
      
      waddnew = 1.0d0
      lterm = .false.

      call rosenbluth( .true.,lterm,iins,iins,imolty,ifrom
     &     ,boxins,igrow,waddnew,lfixnow,ctorfn,2 )

c     --- termination of cbmc attempt due to walk termination ---
      if ( lterm ) return

      if ( ldual .or. lewald .or. iunit .ne. igrow 
     &     .or. ((.not. lchgall) .and. lelect(imolty)) ) then
c     --- Put on hydrogens for explicit AA model for calculation of COM
c     --- and assign all of the grown new and old beads to rxuion
c     --- with rxuion: new = 2
         iii = 2
         do j=1,igrow
            rxuion(j,iii) = rxnew(j)
            ryuion(j,iii) = rynew(j)
            rzuion(j,iii) = rznew(j)
            qquion(j,iii) = qqu(iins,j)
         enddo
         if ( igrow .ne. iunit .and. .not. lrigid(imolty) ) then
            idum = nchain+1
            do j=1,igrow
               rxu(idum,j) = rxnew(j)
               ryu(idum,j) = rynew(j)
               rzu(idum,j) = rznew(j)
c               write(2,*) rxu(idum,j),ryu(idum,j),rzu(idum,j)
            enddo
            moltyp(idum) = imolty
            call explct(idum,vtornew,.false.,.false.)
            do j=igrow+1,iunit
               rxuion(j,iii) = rxu(idum,j)
               ryuion(j,iii) = ryu(idum,j)
               rzuion(j,iii) = rzu(idum,j)
               qquion(j,iii) = qqu(iins,j)
c               write(2,*) rxu(idum,j),ryu(idum,j),rzu(idum,j)
            enddo
         endif
         moltion(iii) = imolty
         
         ibox=boxins
         nboxi(iins) = ibox
      endif
      
c     --- Begin DC-CBMC Corrections for NEW configuration
c      if (ldual .or. ((.not. lchgall) .and. lelect(imolty))) then 

      if (ldual .or. ((.not. lchgall) .and. lelect(imolty))
     &     .or. (lchgall .and. lewald .and. (.not. ldual))) then

c     calculate the true site-site energy
         istt = 1
         iett = igrow
!         write(6,*) igrow
         
         call energy (iins,imolty, v, vintra,vinter,vext
     &        ,velect,vewald,iii,ibox, istt, iett, .true.,ovrlap
     &        ,.false.,vdum,.false.,lfavor)
         
         if (ovrlap) then
            write(2,*) 'iins',iins,'irem',irem
            stop 'strange screwup in DC-CBMC swap'
         endif
c v1insewd, vnewewald and vnewintra now accounted for in v from energy
c$$$         delen = v - ( vnewinter + vnewext +vnewelect) 
c$$$     &        - (v1ins - v1insewd)
         delen = v - ( vnewinter + vnewext +vnewelect + vnewintra
     &        + vnewewald + v1ins) 
         waddnew = waddnew*dexp(-beta*delen)
         vnewt     = vnewt + delen
         vnewinter = vinter - v1insint
         vnewext   = vext - v1insext
         vnewelect = velect - v1inselc
         vnewewald = vewald - v1insewd
         vnewintra = vintra
      endif
c     End DC-CBMC Corrections for NEW configuration
      
c     Begin Explicit Atom Corrections for NEW configuration
      if ( iunit .ne. igrow ) then
c     calculate the true Lennard-Jones energy for the hydrogens
c     iii=2 new conformation
         
         istt = igrow+1
         iett = iunit
         ltors = .false.
         whins = 0.0d0
         ichoi = nchoih(imolty)
         
         do ip = 1, ichoi
c        --- place the explicit hydrogens nchoih times and compute the 
c        --- full intermolecular and intramolecular energy.  Then select 
c        --- one with the usual rosenbluth scheme

            if ( ip .eq. 1) then
c           ---  ip .eq. 1, uses the hydrogen positions generated above
               do j = igrow+1, iunit
                  rxuh(j,ip) = rxuion(j,iii)
                  ryuh(j,ip) = ryuion(j,iii)
                  rzuh(j,ip) = rzuion(j,iii)
               enddo
            else
               idum = nchain+1
c              --- generate positions for the hydrogens
               call explct(idum,vtornew,.false.,.false.)
               do j = igrow+1, iunit
                  rxuion(j,iii) = rxu(idum,j)
                  ryuion(j,iii) = ryu(idum,j)
                  rzuion(j,iii) = rzu(idum,j)
                  rxuh(j,ip) = rxu(idum,j)
                  ryuh(j,ip) = ryu(idum,j)
                  rzuh(j,ip) = rzu(idum,j)
               enddo
            endif

c ??? problem here on calculation of favor and favor2 when ichoi > 1

            call energy (iins,imolty,v, vintra,vinter,vext,velect
     &           ,vewald
     &           ,iii,ibox,istt,iett, .true.,ovrlap,ltors,vdum
     &           ,.true.,lfavor)
c            write(2,*) 'ovrlap:',ovrlap
            
c            if ( iins .eq. 118) write(2,*) 'vinter:',vinter

            if (ovrlap) then
               lovrh(ip) = .true.
            else
               delenh(ip) = v + vtornew
               
               if ( delenh(ip)*beta .gt. (3.3d0*softcut) ) then
                  lovrh(ip) = .true.
               else
c                 --- calculate the boltzman and rosenbluth weight
                  bfach(ip) = dexp(-beta*delenh(ip))
                  whins = whins + bfach(ip)
                  lovrh(ip) = .false.
                  vtrhintra(ip) = vintra
                  vtrhinter(ip) = vinter
                  vtrhext(ip) = vext
                  vtrhtg(ip) = vtornew
                  vtrhelect(ip) = velect
                  vtrhewald(ip) = vewald
               endif
            endif
         enddo

c        --- check for termination of walk
         if ( whins .lt. softlog ) return
         
         if ( ichoi .gt. 1 ) then
            rbf = whins * random()
            bsum = 0.0d0
            do ip = 1, ichoi
               if ( .not. lovrh(ip) ) then
                  bsum = bsum + bfach(ip)
                  if ( rbf .lt. bsum ) then
                     iwalk = ip
                     goto 190
                  endif
               endif
            enddo
            stop 'screw up in explicit hydrogen'
         else
            iwalk = 1
         endif
         
 190     waddnew = waddnew * whins
         vnewt     = vnewt + delenh(iwalk)
         vnewintra = vnewintra + vtrhintra(iwalk)
         vnewinter = vnewinter + vtrhinter(iwalk)
         vnewext   = vnewext + vtrhext(iwalk)
         vnewtg = vnewtg + vtrhtg(iwalk)
         vnewelect = vnewelect + vtrhelect(iwalk)
         vnewewald = vnewewald + vtrhewald(iwalk)
         do j = igrow+1, iunit
            rxuion(j,iii) = rxuh(j,iwalk)
            ryuion(j,iii) = ryuh(j,iwalk)
            rzuion(j,iii) = rzuh(j,iwalk)
         enddo
      endif
      
c     End Explicit Atom Corrections for NEW configuration

c     Begin Ewald-sum Corrections

      if ( lewald ) then
c --- reciprocal space sum
c --- prepare qquion(jj,1) etc
         moltion(1) = imolty
         if ( lswapinter ) then
            do j = 1,iunit
               qquion(j,1) = 0.0d0
            enddo
            call recip(boxins,vrecipn,vrecipo,1)
            delen = vrecipn - vrecipo 
            waddnew = waddnew*dexp(-beta*delen)
            vnewelect = vnewelect + delen
            vnewt = vnewt + delen
         endif
      endif
      
c     End Ewald-sum Corrections

      if (lpbcz) then
          if (lsolid(boxins) .and. .not. lrect(boxins)) then
             ibox = boxins
             volins = cell_vol(ibox)
         else
            volins=boxlx(boxins)*boxly(boxins)*boxlz(boxins)
         endif
      else
         volins=boxlx(boxins)*boxly(boxins)
      endif
      
c     Begin Fluctuating Charge corrections for NEW configuration

      vnewflucq = 0.0d0
      if ( lflucq(imolty) ) then
         do j = 1,iunit
            qion(j) = qqu(iins,j)
         enddo
         call charge(iins, qion, vflucq,vdum)
c        --- once we go to fully flexible will need to compute this in boltz
         vnewflucq = vflucq
         vnewt = vnewt + vflucq
         waddnew = waddnew*dexp(-beta*vflucq)
      endif

c     End Fluctuating Charge corrections for NEW configuration

c     Begin Tail corrections for BOXINS with inserted particle

      if (ltailc .and. lswapinter) then
         vinsta = 0.0d0
         do imt = 1, nmolty
            do jmt = 1, nmolty


               if ( jmt .eq. imolty ) then
                  rho = dble( ncmt(boxins,jmt) + 1 ) / volins
               else
                  rho = dble( ncmt(boxins,jmt) ) / volins
               endif
               if ( imt .eq. imolty ) then
                  vinsta = vinsta + 
     &              dble( ncmt(boxins,imt) + 1 ) * coru(imt,jmt,rho
     &                                               ,boxins)
               else
                  vinsta = vinsta + 
     &                 dble( ncmt(boxins,imt) ) * coru(imt,jmt,rho
     &                                                ,boxins)
               endif
            enddo
         enddo

c         if(LSOLPAR.and.(boxins.eq.2))then
c           vinsta = 0.0d0
c         endif

         vinsta = vinsta - vtailb( boxins )
         waddnew = waddnew*dexp(-beta*vinsta)
         vnewt = vnewt + vinsta
         vnewinter = vnewinter + vinsta
         
c     this approximate method of tail correction was used for chem. pot.
c     until 1-26-98 MGM
c      dtest = 0.0d0
c         do jmt = 1, nmolty
c            if ( jmt .eq. imolty ) then
c               rho = dble(ncmt(boxins,jmt)+1)/volins
c            else
c               rho = dble(ncmt(boxins,jmt))/volins
c            endif
c            dtest = dtest + 2.0d0*coru(imolty,jmt,rho)
c            arg=arg*dexp(-beta*2.0d0*coru(imolty,jmt,rho))
c         enddo
c         write(2,*) 'dtest',dtest
      else
         vinsta = 0.0d0
      endif

c     End Tail corrections for BOXINS with inserted particle

      if ( .not. lanes ) then
         if ( lswapinter ) then
            arg = w1ins * waddnew * weight * volins 
     &           / dble( ncmt(boxins,imolty)+1 )
            acchem(boxins,imolty) = acchem(boxins,imolty)+arg
         endif
      endif

      if (lexpee.and.leemove) then
         imolty = ee_moltyp(mstate)
         moltion(1) = imolty
         do icbu = 1, iunit
            qqu(irem,icbu) = ee_qqu(icbu,mstate)
         enddo
         goto 500
      endif

      bsswap(imolty,ipairb,boxins+nbox) = 
     &     bsswap(imolty,ipairb,boxins+nbox) + 1.0d0

C     Compute weights for the molecule to be removed from boxrem

c *** check that there is at least one molecule in BOXREM ***
      if ( lempty ) then
c         write(2,*) 'no molecule in BOXREM'
         if (lgrand) then
	    if (boxrem.eq.2) write(2,*) ' ERROR ***** array too low !'
         endif
         return
      endif

 500  continue


c *** select a position of the first/starting unit at RANDOM ***
c *** and calculate the boltzmann weight                     ***
c *** for the chain to be REMOVED                           ***

      rxp(1,1) = rxu(irem,beg)
      ryp(1,1) = ryu(irem,beg)
      rzp(1,1) = rzu(irem,beg)

      ichoi = nchoi1(imolty)

      if ( .not. lswapinter ) then

c *********************
c * second AVMBC part *
c *********************

         if (boxins .ne. boxrem) then
            write(2,*) 'avbmc, boxins, boxrem:',boxins,boxrem
            stop
         endif
         
         if ( lbias(imolty) ) then
            if (lavbmc1(imolty)) then
         
c *** determine whether it is in or out

               x = rxu(irem,1) - rxu(jins,1)
               y = ryu(irem,1) - ryu(jins,1)
               z = rzu(irem,1) - rzu(jins,1)
               if ( lpbc ) call mimage(x,y,z,boxins)
               rijsq = x*x + y*y + z*z
               if ( rijsq .gt. rbsmax*rbsmax .or. 
     &              rijsq .lt. rbsmin*rbsmin ) then
                  lrem_out = .true.
               else
                  lrem_out = .false.
               endif
            endif
            if ( lrem_out ) then

c *** out-> j case

               wbias_ins = (1-pmbias(imolty))/
     &              (boxlx(boxins)*boxly(boxins)*boxlz(boxins)-
     &              vol_eff)

               if (lavbmc2(imolty) .or. lavbmc3(imolty) ) 
     &              wbias_ins=wbias_ins/dble(neighj_num+1)

c            write(2,*) 'originally out'
c            write(2,*) '1:',wbias_ins,boxlx(boxins),vol_eff

               do icbu = 2,ichoi
 232              rxp(1,icbu) = boxlx(boxins) * random()
                  ryp(1,icbu) = boxly(boxins) * random()
                  if (lpbcz .or. lslit) then
                     rzp(1,icbu) = boxlz(boxins) * random()
                  else if ( lsami .or. lmuir .or. ljoe ) then
                     if ( lempty ) then
                        rzp(1,icbu) = 20*random()-10
                     else
                        rzp(1,icbu) = rzu(irem,1)
                     endif
                  else
                     rzp(1,icbu) = 0.0d0
                  endif
            
c *** determine whether it is inside the sphere with particle jins
                  x = rxp(1,icbu) - rxu(jins,1)
                  y = ryp(1,icbu) - ryu(jins,1)
                  z = rzp(1,icbu) - rzu(jins,1)
                  if ( lpbc ) call mimage(x,y,z,boxins)
                  rijsq = x*x + y*y + z*z
                  if ( rijsq .lt. rbsmax**2 .and. rijsq .gt. rbsmin**2 ) 
     &                 goto 232
               enddo
               
            else

               if (lavbmc3(imolty) .and. lremk_in) then

c *** selected from the region bounded by kins
c *** k -> j case 

                  target = kins
                  wbias_ins=(1-pmbias(imolty))/vol_eff
     &                 /dble(neighj_num+1)

               else

c *** selected from the region bounded by jins

                  target = jins

                  wbias_ins = pmbias(imolty)/vol_eff
               
                  if ( lavbmc2(imolty) ) 
     &                 wbias_ins=wbias_ins/dble(
     &                 ncmt(boxrem,imolty)-neighj_num-joffset+1)

                  if ( lavbmc3(imolty) ) then
                     if ( linsk_in ) then
c *** j -> k case
                        wbias_ins = wbias_ins/(neighk_num+1)
                     else
c *** j -> out case
                        wbias_ins=wbias_ins/dble(
     &                       ncmt(boxrem,imolty)-neighj_num-joffset+1)
                     endif
                  endif
               endif
              
c            write(2,*) 'originally in'
c            write(2,*) '2:',wbias_ins,boxlx(boxins),vol_eff

               do icbu = 2,ichoi
c *** choose a random association distance
                  rvol = random()
                  r = (rbsmax*rbsmax*rbsmax*rvol + 
     &                 (1-rvol)*rbsmin*rbsmin*rbsmin)**(1.0/3.0d0)
c --- calculate random vector on the unit sphere ---
 129              xi1 = ( 2.0d0 * random() ) - 1.0d0
                  xi2 = ( 2.0d0 * random() ) - 1.0d0
                  xisq = xi1**2 + xi2**2
                  if ( xisq .lt. 1.0d0 ) then
                     x = r * 2.0d0 * xi1 * dsqrt( 1.0d0 - xisq )
                     y = r * 2.0d0 * xi2 * dsqrt( 1.0d0 - xisq )
                     z = r * ( 1.0d0 - 2.0d0 * xisq )
                  else
                     goto 129
                  endif
                  rxp(1,icbu) = rxu(target,1) + x
                  ryp(1,icbu) = ryu(target,1) + y
                  rzp(1,icbu) = rzu(target,1) + z
               enddo
               
            endif
         endif

c *************************
c * end second AVMBC part *
c *************************

      else

         if (lsolid(boxrem) .and. .not. lrect(boxrem)) then

            ibox = boxrem
            do icbu = 2,ichoi
               sx = random()
               sy = random()
               sz = random()
               rxp(1,icbu) = sx*hmat(ibox,1)+sy*hmat(ibox,4)
     &              +sz*hmat(ibox,7)
               ryp(1,icbu) = sx*hmat(ibox,2)+sy*hmat(ibox,5)
     &              +sz*hmat(ibox,8)
               rzp(1,icbu) = sx*hmat(ibox,3)+sy*hmat(ibox,6)
     &              +sz*hmat(ibox,9)
            enddo

         else

            do icbu = 2,ichoi
               rxp(1,icbu) = boxlx(boxrem) * random()
               ryp(1,icbu) = boxly(boxrem) * random()
               if (lpbcz .or. lslit) then
                  rzp(1,icbu) = boxlz(boxrem) * random()
               else if ( lsami .or. lmuir .or. ljoe ) then
                  if ( lempty ) then
                     rzp(1,icbu) = 20*random()-10
                  else
                     rzp(1,icbu) = rzu(irem,1)
                  endif
               else
                  rzp(1,icbu) = 0.0d0
               endif
            enddo
         endif
         
      endif

      if ( .not. lswapinter .and. lbias(imolty) ) then
         if ( lrem_out .and. lins_in ) then
            bnswap_in(imolty,1) = bnswap_in(imolty,1) + 1.0d0
         elseif ( (.not. lrem_out) .and. (.not. lins_in ) ) then
            bnswap_out(imolty,1) = bnswap_out(imolty,1) + 1.0d0
         endif
      endif


c *** calculate the boltzmann weight of first bead          ***

      lnew = .false.

      call boltz(lnew,.true.,ovrlap,irem,irem,imolty,boxrem,ichoi,idum
     &     ,1,glist,0.0d0)

      if ( ovrlap ) then
         write(2,*) 'disaster: overlap for 1st bead in SWAP'
      endif
c *** calculate the correct weight for the  old  walk ***

      w1rem = 0.0d0
      do ip = 1, ichoi
         w1rem = w1rem + bfac(ip)
      enddo

c --- check for termination of walk ---
      if ( w1rem .lt. softlog ) then 
         write(2,*) ' run problem : soft overlap in old'
      endif

      v1rem = vtry(1)
      v1remint = vtrinter(1)
      v1remext = vtrext(1)
      v1remelc = vtrelect(1)
      v1remewd = vtrewald(1)

      waddold = 1.0d0

c     --- call rosenbluth for old conformation

      call rosenbluth(.false.,lterm,irem,irem,imolty,ifrom
     &     ,boxrem,igrow,waddold,lfixnow,ctorfo,2 )

      if ( lterm ) then 
         write(2,*) 'SWAP: rosenbluth old rejected'
         return
      endif

      if ( ldual .or. lewald .or. igrow .ne. iunit 
     &     .or. ((.not. lchgall) .and. lelect(imolty)) ) then
c     --- store the old grown beads and explict placed beads positions
c     --- 1 = old conformation
         iii = 1
         do j = 1,iunit
            rxuion(j,1) = rxu(irem,j)
            ryuion(j,1) = ryu(irem,j)
            rzuion(j,1) = rzu(irem,j)
            qquion(j,1) = qqu(irem,j)
         enddo
      endif

c     Begin Correction for DC-CBMC for OLD configuration
c      if (ldual .or. ((.not. lchgall) .and. lelect(imolty)) ) then 


      if (ldual .or. ((.not. lchgall) .and. lelect(imolty))
     &     .or. (lchgall .and. lewald .and. (.not. ldual))) then


c     --- correct the acceptance rules 
c     --- calculate the Full rcut site-site energy
         istt=1
         iett = igrow

         call energy (irem,imolty, v, vintra,vinter,vext,velect
     &        ,vewald,iii, boxrem, istt, iett, .true.,ovrlap
     &        ,.false.,vtorold,.false.,lfavor)
            
         if (ovrlap) stop 'disaster ovrlap in old conf SWAP'
c v now includes vnewintra,v1remewd and voldewald, take out
c$$$         deleo = v - ( voldinter + voldext +voldelect) 
c$$$     &        - (v1rem - v1remewd)
         deleo = v - ( voldinter + voldext +voldelect + voldintra 
     &        + voldewald + v1rem) 
         waddold = waddold*dexp(-beta*deleo)
         voldt     = voldt + deleo
         voldinter = vinter - v1remint
         voldext   = vext - v1remext
         voldelect = velect - v1remelc
         voldewald = vewald - v1remewd
         voldintra = vintra
      endif
c     End Correction for DC-CBMC for OLD configuration

c     Begin Correction for Explicit Atom for OLD configuration
      if ( iunit .ne. igrow ) then
c        --- calculate the true Lennard-Jones energy for the hydrogens
c        --- Calculate the energy of the non-backbone beads 
c        --- iii=1 old conformation

         ibox = boxrem
         ltors = .true.
         istt = igrow + 1
         iett = iunit
         
         call energy (irem,imolty,v, vintra,vinter,vext,velect
     &        ,vewald
     &        ,iii,ibox,istt,iett,.true.,ovrlap,ltors,vtorold
     &        ,.true.,lfavor)

c         if (irem .eq. 118)  write(2,*) 'for old',vinter
         if (ovrlap) stop 'disaster ovrlap in old conf SWAP'
         deleo = v + vtorold
         voldt     = voldt + deleo
         voldintra = voldintra + vintra
         voldinter = voldinter + vinter
         voldext   = voldext + vext
         voldtg    = voldtg + vtorold
         voldelect = voldelect + velect
         voldewald = voldewald + vewald
         
         whrem = dexp(-beta*deleo)
         ichoi = nchoih(imolty)

         if ( ichoi .ne. 1 ) then
c        --- torsion is calculated in explicit for the placed atoms
            ltors = .false.
c        store the old coordinates for the hydrogens
            do j = igrow+1, iunit
               rxuold(j) = rxuion(j,iii)
               ryuold(j) = ryuion(j,iii)
               rzuold(j) = rzuion(j,iii)
            enddo

            do ip = 2,ichoi
c        --- rosenbluth weight for multiple placement of explicit hydrogens 
               call explct(irem,vtorold,.false.,.false.)
               do j = igrow + 1, iunit
                  rxuion(j,iii) = rxu(irem,j)
                  ryuion(j,iii) = ryu(irem,j)
                  rzuion(j,iii) = rzu(irem,j)
               enddo

c ??? problem here on calculation of favor and favor2               

               call energy (irem,imolty,v, vintra,vinter,vext
     &              ,velect
     &              ,vewald,iii,ibox,istt,iett, .true.,ovrlap,ltors
     &              ,vdum,.true.,lfavor)
               
               deleo = v + vtorold
               if ( .not. ovrlap ) then
                  whrem = whrem + dexp(-beta*deleo)
               endif
            enddo
            if ( ichoi .gt. 1 ) then
c              --- restore the old hydrogen positions
               do j = igrow+1, iunit
                  rxuion(j,iii) = rxuold(j)
                  ryuion(j,iii) = ryuold(j)
                  rzuion(j,iii) = rzuold(j)
                  rxu(irem,j) = rxuold(j)
                  ryu(irem,j) = ryuold(j)
                  rzu(irem,j) = rzuold(j)
               enddo
            endif
         endif

         waddold = waddold*whrem
      endif
      
c     End Correction for Explicit Atom for OLD configuration

c     Begin Ewald-sum Corrections for OLD configuration

      if ( lewald ) then
c --- reciprocal space sum on r*uion
c --- prepare qquion(jj,1) etc
         if ( lswapinter ) then
            do j = 1,iunit
               qquion(j,2) = 0.0d0
            enddo
         endif
         call recip(boxrem,vrecipn,vrecipo,1)
         deleo = vrecipo - vrecipn
         voldt = voldt + deleo 
         voldelect = voldelect + deleo
         waddold = waddold * dexp(-beta*deleo)
      endif

c     End Ewald-sum Corrections for OLD configuration
      
c     Begin Fluctuating Charge corrections for OLD configuration

      voldflucq = 0.0d0
      if ( lflucq(imolty) ) then
c        --- this is a temporary fix - eventually should implement fully 
c        --- flexible fluctuating charge calculation into boltz
c        --- right now the old flucq energy is the same as the new flucq
         voldflucq = vnewflucq
         voldt = voldt + voldflucq
         waddold = waddold*dexp(-beta*voldflucq)
      endif

      
c     End Fluctuating Charge corrections for OLD configuration

      if (lpbcz) then
         if (lsolid(boxrem) .and. .not. lrect(boxrem)) then
            ibox = boxrem
            volrem = cell_vol(ibox)
         else
            volrem=boxlx(boxrem)*boxly(boxrem)*boxlz(boxrem)
         endif
      else
         volrem=boxlx(boxrem)*boxly(boxrem)
      endif
      
c     Start of intermolecular tail correction for boxrem

      if (ltailc .and. lswapinter) then
c     --- BOXREM without removed particle
         vremta = 0.0d0
         do imt = 1, nmolty
            do jmt = 1, nmolty
                  
               if ( jmt .eq. imolty ) then
                  rho = dble( ncmt(boxrem,jmt) - 1 ) / volrem
               else
                  rho = dble( ncmt(boxrem,jmt) ) / volrem
               endif
               if ( imt .eq. imolty ) then
                  vremta = vremta + 
     &              dble( ncmt(boxrem,imt) - 1 ) * coru(imt,jmt,rho,
     &                                                boxrem)
               else
                  vremta = vremta + 
     &              dble( ncmt(boxrem,imt) ) * coru(imt,jmt,rho
     &                                              ,boxrem )
               endif
            enddo
         enddo
c         if(LSOLPAR.and.(boxrem.eq.2)) then
c           vremta = 0.0d0
c         endif
   
         vremta = - vremta + vtailb( boxrem )
         waddold=waddold*dexp(-beta*vremta) 
         voldt = voldt + vremta
         voldinter = voldinter + vremta
      else
         vremta = 0.0d0
      endif
c     End of intermolecular tail correction for boxrem

c     --- Add contributions of the first bead and additional beads:

      vnewt     = vnewt  + v1ins
      vnewinter = vnewinter + v1insint
      vnewext   = vnewext + v1insext
      vnewelect = vnewelect + v1inselc
      vnewewald = vnewewald + v1insewd
      
      voldt     = voldt  + v1rem
      voldinter = voldinter+(v1remint)
      voldext   = voldext+v1remext
      voldelect = voldelect + v1remelc
      voldewald = voldewald + v1remewd

      weight= w1ins * waddnew * weight
      weiold= w1rem * waddold * weiold

      wnlog = dlog10 ( weight )
      wolog = dlog10 ( weiold )
      wdlog = wnlog - wolog
      
      if ( wdlog .lt. -softcut ) then
c         write(2,*) '### underflow in wratio calculation ###'
         return
      endif

c *** For ANES algorithm, do the Fluctuating charge moves.
      
      if ( lanes ) then
         if ( lswapinter ) then
            nboxi(irem) = boxins         
            parbox(ncmt(boxins,imolty)+1,boxins,imolty)= irem
            parbox(pointp,boxrem,imolty)=
     &           parbox(ncmt(boxrem,imolty),boxrem,imolty)
            parbox(ncmt(boxrem,imolty),boxrem,imolty)=0         
            nchbox(boxins) = nchbox(boxins) + 1
            nchbox(boxrem) = nchbox(boxrem) - 1
            ncmt(boxins,imolty) = ncmt(boxins,imolty) + 1
            ncmt(boxrem,imolty) = ncmt(boxrem,imolty) - 1         
         endif
         
         favor(irem) = 1.0d0
         favor2(irem) = 1.0d0

         call anes(irem,boxins,boxrem,3,laccept,vdum,vdum,vdum,
     &        vdum,vdum,vdum,vdum,vdum,vdum,vinsta,vremta,
     &        vnewflucq,voldflucq,lswapinter)

         if ( lswapinter ) then
            arg = weight * volins 
     &           / dble( ncmt(boxins,imolty) )
            acchem(boxins,imolty) = acchem(boxins,imolty)+arg
            bnchem(boxins,imolty) = bnchem(boxins,imolty) + 1.0d0
         endif

         if ( laccept ) then
            bnswap(imolty,ipairb,boxins+nbox)=
     &           bnswap(imolty,ipairb,boxins+nbox)+1.0d0

         elseif ( lswapinter ) then
            nboxi(irem) = boxrem            
            parbox(ncmt(boxrem,imolty)+1,boxrem,imolty)=irem
            parbox(ncmt(boxins,imolty),boxins,imolty)=0
            nchbox(boxins) = nchbox(boxins) - 1
            nchbox(boxrem) = nchbox(boxrem) + 1
            ncmt(boxins,imolty) = ncmt(boxins,imolty) - 1
            ncmt(boxrem,imolty) = ncmt(boxrem,imolty) + 1         
         endif    
         return         
      endif

      if ( lswapinter ) then
         if (lgibbs.and.((.not.leemove).and.(.not.lexpee))) then
c     --- Note: acceptance based on only molecules of type imolty
            wratio = ( weight / weiold ) * wee_ratio *
     +           ( volins * dble( ncmt(boxrem,imolty) ) /
     +           ( volrem * dble( ncmt(boxins,imolty1) + 1 ) ) )
     +           * dexp(beta*(eta2(boxrem,imolty)-
     +           eta2(boxins,imolty1)))
         else if (lgibbs.and.(leemove.and.lexpee)) then
            wratio = ( weight / weiold ) * wee_ratio *
     +            volins/volrem * dexp(beta*(eta2(boxrem,imolty)-
     +           eta2(boxins,imolty1)))
                                                                                
         else if (lgrand) then
            if (boxins.eq.1) then
c              --- molecule added to box 1
               wratio = (weight /  weiold ) *
     +              volins * B(imolty) / (ncmt(boxins,imolty)+1)
            else
c              --- molecule removed from box 1
               wratio = (weight /weiold)*
     +              dble(ncmt(boxrem,imolty))/(volrem*B(imolty))
                                                                                
            endif
         endif
      else
         wratio = (weight*wbias_ins)/(weiold*wbias_rem)*wee_ratio
      endif

      if (.not. lswapinter) then

         neigh_old = neigh_cnt(irem)
         if ( neigh_old .le. 6 .and. neigh_icnt .le. 6 ) then
            if ( lavbmc1(imolty) ) then
               if ( lrem_out ) then
                  if ( lins_in ) then
                     ip = 1
                  else
                     ip = 2
                  endif
               else
                  if ( lins_in ) then
                     ip = 3
                  else
                     ip = 4
                  endif
               endif
            elseif ( lavbmc2(imolty) ) then
               if ( lins_in ) then
                  ip = 1
               else
                  ip = 2
               endif
            else
               if ( lins_in ) then
                  if ( lremk_in ) then
                     ip = 1
                  else
                     ip = 2
                  endif
               elseif ( linsk_in ) then
                  ip =3
               else
                  ip =4
               endif
            endif
c         write(2,*) neigh_old,neigh_icnt,ip
c            if ( .not. lins_in .and. neigh_old .eq. 0 ) then
c               write(2,*) '####error:',irem,jins
c            endif
            cnt_wf1(neigh_old,neigh_icnt,ip)
     &           = cnt_wf1(neigh_old,neigh_icnt,ip)+1
            if (neigh_old .eq. 0 .and. neigh_icnt .eq. 1) then
               wdlog = dlog10 (wratio)
               ic = dint((wdlog+95.0d0)/0.1d0)+1
               if ( ic .lt. 1 ) ic = 1
               if ( ic .gt. 1000 ) ic = 1000
               cnt_wra1(ic,ip) = cnt_wra1(ic,ip) + 1
            elseif (neigh_old .eq. 1 .and. neigh_icnt .eq. 0) then
               wdlog = dlog10 (wratio)
               ic = dint((wdlog+95.0d0)/0.1d0) + 1
               if ( ic .lt. 1 ) ic = 1
               if ( ic .gt. 1000 ) ic = 1000
               cnt_wra2(ic,ip) = cnt_wra2(ic,ip) + 1
            endif
         endif

      endif

      if (lfixnow) then
         wratio = wratio * ctorfo / ctorfn
      endif
       
c         wratio = 1.0   

      if ( random() .le. wratio ) then
c         write(2,*) 'SWAP MOVE ACCEPTED',irem
c *** we can now accept !!!!! ***
         if ((.not.leemove).and.(.not.lexpee)) then
            bnswap(imolty,ipairb,boxins+nbox) = 
     &        bnswap(imolty,ipairb,boxins+nbox) + 1.0d0
         endif
         if ( .not. lswapinter .and. lbias(imolty) ) then
            if ( lrem_out .and. lins_in ) then
               bnswap_in(imolty,2) = bnswap_in(imolty,2) + 1.0d0
            elseif ( (.not. lrem_out) .and. (.not. lins_in ) ) then
               bnswap_out(imolty,2) = bnswap_out(imolty,2) + 1.0d0
            endif
         endif

c---update the position, it will be used to get the bonded energy
         do ic = 1,igrow
            rxu(irem,ic) = rxnew(ic)
            ryu(irem,ic) = rynew(ic)
            rzu(irem,ic) = rznew(ic)
            if (leemove.and.lexpee) qqu(irem,ic) = qqu(iins1,ic)
         enddo
         do ic = igrow+1,iunit
            rxu(irem,ic) = rxuion(ic,2)
            ryu(irem,ic) = ryuion(ic,2)
            rzu(irem,ic) = rzuion(ic,2)
         enddo

         vvibn  = 0.0d0
         vbendn = 0.0d0
         vtgn   = 0.0d0
         
         if (lrigid(imolty).and.(pm_atom_tra.gt.0.000001)) then 
           call Intra_energy(irem,imolty, vdum ,vintran, vdum,vdum,
     &       velectn,vewaldn,flagon, boxins, 1, iunit,.true.,ovrlap,
     &       .false.
     &       ,vdum,.false.,.false.,vvibn,vbendn,vtgn) 
         endif
c         total_NBE = vintran+velectn+vewaldn+vtgn+vbendn+vvibn 
          total_NBE = vtgn+vbendn+vvibn
!         write(6,*) vintran,velectn,vewaldn       

!         write(6,*) 'irem', irem  

c ---    update energies:

         vbox(boxrem)     = vbox(boxrem)     - voldt - total_NBE
	 vinterb(boxrem)  = vinterb(boxrem)  - voldinter
         vtailb(boxrem)   = vtailb(boxrem)   - vremta
	 vintrab(boxrem)  = vintrab(boxrem)  - voldintra
	 vvibb(boxrem)    = vvibb(boxrem)    - voldbvib - vvibn
	 vtgb(boxrem)     = vtgb(boxrem)     - voldtg - vtgn
	 vextb(boxrem)    = vextb(boxrem)    - voldext
	 vbendb(boxrem)   = vbendb(boxrem)   - voldbb - vbendn
         velectb(boxrem)  = velectb(boxrem)  - (voldelect+voldewald)
         vflucqb(boxrem)  = vflucqb(boxrem)  - voldflucq
	
 
         vbox(boxins)     = vbox(boxins)     + vnewt + total_NBE
	 vinterb(boxins)  = vinterb(boxins)  + vnewinter
	 vtailb(boxins)   = vtailb(boxins)   + vinsta
	 vintrab(boxins)  = vintrab(boxins)  + vnewintra
	 vvibb(boxins)    =  vvibb(boxins)   + vnewbvib + vvibn
	 vtgb(boxins)     = vtgb(boxins)     + vnewtg + vtgn
	 vextb(boxins)    = vextb(boxins)    + vnewext
	 vbendb(boxins)   = vbendb(boxins)   + vnewbb + vbendn
         velectb(boxins)  = velectb(boxins)  + (vnewelect+vnewewald)
         vflucqb(boxins)  = vflucqb(boxins)  + vnewflucq

c ---    update book keeping

         if ( lswapinter ) then
            nboxi(irem) = boxins
         
            if ((.not.leemove).and.(.not.lexpee)) then
               parbox(ncmt(boxins,imolty)+1,boxins,imolty)= irem
                parbox(pointp,boxrem,imolty)=
     &           parbox(ncmt(boxrem,imolty),boxrem,imolty)
               parbox(ncmt(boxrem,imolty),boxrem,imolty)=0
            
               nchbox(boxins) = nchbox(boxins) + 1
               nchbox(boxrem) = nchbox(boxrem) - 1
               ncmt(boxins,imolty) = ncmt(boxins,imolty) + 1
               ncmt(boxrem,imolty) = ncmt(boxrem,imolty) - 1
            else
               parbox(ncmt(boxins,imolty1)+1,boxins,imolty1)= irem
                parbox(pointp,boxrem,imolty)=
     &           parbox(ncmt(boxrem,imolty),boxrem,imolty)
               parbox(ncmt(boxrem,imolty),boxrem,imolty)=0
            
               nchbox(boxins) = nchbox(boxins) + 1
               nchbox(boxrem) = nchbox(boxrem) - 1
               ncmt(boxins,imolty1) = ncmt(boxins,imolty1) + 1
               ncmt(boxrem,imolty) = ncmt(boxrem,imolty) - 1
               eepointp = ncmt(boxins,imolty1)
            endif

            if ( lexpand(imolty) ) then
               itype = eetype(imolty)
               ncmt2(boxins,imolty,itype) = 
     &              ncmt2(boxins,imolty,itype) + 1
               ncmt2(boxrem,imolty,itype) = 
     &              ncmt2(boxrem,imolty,itype) + 1
            endif
         endif

!         do ic = 1,igrow
!            rxu(irem,ic) = rxnew(ic)
!            ryu(irem,ic) = rynew(ic)
!            rzu(irem,ic) = rznew(ic)
!         enddo
!         do ic = igrow+1,iunit
!            rxu(irem,ic) = rxuion(ic,2)
!            ryu(irem,ic) = ryuion(ic,2)
!            rzu(irem,ic) = rzuion(ic,2)
!         enddo
         
         if ( lewald ) then
c *** update reciprocal-space sum
            if ( lswapinter ) then
               call recip(boxins,vdum,vdum,2)
               call recip(boxrem,vdum,vdum,2)
            else
               call recip(boxins,vdum,vdum,2)
            endif
         endif

c ---    update center of mass
         call ctrmas(.false.,boxins,irem,3)
c *** update linkcell, if applicable
         if ( licell .and. ((boxins .eq. boxlink) .or. (boxrem .eq.  
     &        boxlink))) then
            call linkcell(2,irem,vdum,vdum,vdum,ddum)
         endif
         
c         write(2,*) 'lneighbor:',lneighbor

         if ( lneigh ) call updnn( irem )
         if ( lneighbor ) then
            neigh_old = neigh_cnt(irem)
            if ( neigh_old .le. 6 .and. neigh_icnt .le. 6 ) then
               cnt_wf2(neigh_old,neigh_icnt,ip)
     &              = cnt_wf2(neigh_old,neigh_icnt,ip)+1
            endif

            do 10 ic = 1, neigh_old
               j = neighbor(ic,irem)
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. irem ) then
                     neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                     neigh_cnt(j) = neigh_cnt(j)-1
                     goto 10
                  endif
               enddo
 10         continue
            neigh_cnt(irem) = neigh_icnt
c            write(2,*) 'irem:',irem,neigh_icnt
            do ic = 1,neigh_icnt
               j = neighi(ic)
               neighbor(ic,irem)=j
               lneighij = .false.
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. irem ) then
                     lneighij = .true.
                  endif
               enddo
               if ( .not. lneighij ) then
                  neigh_cnt(j) = neigh_cnt(j)+1
                  neighbor(neigh_cnt(j),j) = irem
               endif
            enddo
         endif

c         write(2,*) irem,'end SWAP'

      endif
 
c -----------------------------------------------------------------
      return
      end

