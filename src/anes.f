      subroutine anes(i,ibox,boxrem,mtype,laccept,deltv,vintern,vintran,
     &     vextn,velectn,vintero,vintrao,vexto,velecto,vinsta,vremta,
     &     vnewflucq,voldflucq,lswapinter)

c anes 
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

c    *********************************************************************
c    ** optimize the electronic configuration for trans, rot, and swap  **
c    ** (config, swatch in the future) moves and accept/reject the      **
c    ** combined move.                                                  **
c    ** written on June 25/99 by Bin Chen.                              **
c    *********************************************************************

      implicit none
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ensemble.inc'
      include 'rosen.inc'
      include 'inputdata.inc'
      include 'system.inc'

      logical::laccept,lswapinter
      integer::i,ibox,boxrem,mtype,imolty,iunit,ichoiq,ip,ibox2,j
      real(8)::deltv,vintern,vintran,vextn,velectn,vintero,
     &     vintrao,vexto,velecto,vinsta,vremta,vnewflucq,voldflucq,
     &     vboxo(nbxmax),vinterbo(nbxmax),vintrabo(nbxmax),
     &     vextbo(nbxmax),velectbo(nbxmax),vflucqbo(nbxmax),
     &     vtailbo(nbxmax),vvibbo(nbxmax),vtgbo(nbxmax),
     &     vbendbo(nbxmax),rxuo(numax),ryuo(numax),rzuo(numax),
     &     xcmo,ycmo,zcmo,vdum,wratio,volins,volrem,deltvb,random,
     &     vnewt2,voldt2
      real(8)::qquo(nmax,numax) 

c      write(iou,*) 'START the optimization of the charge configuration'

      imolty = moltyp(i)
      iunit = nunit(imolty)

c *** store the old energy, old coordinates and ewald sum
      do ibox2 = 1,nbox
         vboxo(ibox2) = vbox(ibox2)
         vinterbo(ibox2) = vinterb(ibox2)
         vintrabo(ibox2) = vintrab(ibox2)
         vextbo(ibox2) = vextb(ibox2)
         vflucqbo(ibox2) = vflucqb(ibox2)
         velectbo(ibox2) = velectb(ibox2)

         if ( mtype .eq. 3 ) then
            vtailbo(ibox2) = vtailb(ibox2)
            vvibbo(ibox2) = vvibb(ibox2)
            vtgbo(ibox2) = vtgb(ibox2)
            vbendbo(ibox2) = vbendb(ibox2)
         endif
         
      enddo
      do j = 1,iunit
         rxuo(j) = rxu(i,j)
         ryuo(j) = ryu(i,j)
         rzuo(j) = rzu(i,j)
      enddo
      xcmo = xcm(i)
      ycmo = ycm(i)
      zcmo = zcm(i)
c *** store the old charges
      do ip = 1,nchain
         do j = 1,nunit(moltyp(ip))
            qquo(ip,j) = qqu(ip,j)
         enddo
      enddo
      if (lewald) then
c *** store the reciprocal-space sum
         do ibox2 = 1,nbox
            call recip(ibox2,vdum,vdum,3)
         enddo
         if ( ldielect ) then
            call dipole(ibox,2)
         endif
      endif

c *** on the new coordinates, continue to use the fluctuating charge
c *** algorithm to optimize the charge configurations, update the
c *** energy, coordinates and the ewald sum
      
      if ( mtype .eq. 3 ) then
c *** for swap move
         vbox(ibox)     = vbox(ibox) + vnewt
         vinterb(ibox)  = vinterb(ibox) + vnewinter
         vintrab(ibox)  = vintrab(ibox) + vnewintra
         vextb(ibox)    = vextb(ibox)   + vnewext
         velectb(ibox)   = velectb(ibox)  + vnewelect+vnewewald
         vtailb(ibox)   = vtailb(ibox)   + vinsta
         vvibb(ibox)    =  vvibb(ibox)   + vnewbvib
         vtgb(ibox)     = vtgb(ibox)     + vnewtg
	 vbendb(ibox)   = vbendb(ibox)   + vnewbb
         vflucqb(ibox)  = vflucqb(ibox)  + vnewflucq

         vbox(boxrem)     = vbox(boxrem)     - voldt 
         vinterb(boxrem)  = vinterb(boxrem)  - voldinter
         vtailb(boxrem)   = vtailb(boxrem)   - vremta
         vintrab(boxrem)  = vintrab(boxrem)  - voldintra
         vvibb(boxrem)    = vvibb(boxrem)    - voldbvib
         vtgb(boxrem)     = vtgb(boxrem)     - voldtg
         vextb(boxrem)    = vextb(boxrem)    - voldext
         vbendb(boxrem)   = vbendb(boxrem)   - voldbb
         velectb(boxrem)  = velectb(boxrem)  - 
     &        (voldelect+voldewald)
         vflucqb(boxrem)  = vflucqb(boxrem)  - voldflucq
      else

         vbox(ibox)     = vbox(ibox) + deltv
         vinterb(ibox)  = vinterb(ibox) + (vintern - vintero)
         vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
         vextb(ibox)    = vextb(ibox)   + (vextn   - vexto)
         velectb(ibox)   = velectb(ibox)  + (velectn - velecto)
      endif
      
      do j = 1,iunit
         rxu(i,j) = rxuion(j,2)
         ryu(i,j) = ryuion(j,2)
         rzu(i,j) = rzuion(j,2)
      enddo
c *** update chain center of mass
      call ctrmas(.false.,ibox,i,mtype)

      if (lewald) then
c *** update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
         if (mtype .eq. 3) call recip(boxrem,vdum,vdum,2)
         if ( ldielect ) then
c *** update the dipole term
            call dipole(ibox,1)
         endif
      endif            

c *** begin to optimize the charge configuration

      if ( mtype .eq. 3 ) then
         do ichoiq = 1,500
            call flucq(-1,0)
         enddo
c         do ichoiq = 1,0
c            call flucq(-2,0)
c         enddo
         do ichoiq = 1,500
            call flucq(2,0) 
         enddo
         deltv = vbox(ibox) - vboxo(ibox)
         weight = dexp(-deltv*beta)
         deltv = vboxo(boxrem) - vbox(boxrem)
         weiold = dexp(-deltv*beta)
         volins=boxlx(ibox)*boxly(ibox)*boxlz(ibox)
         volrem=boxlx(boxrem)*boxly(boxrem)*boxlz(boxrem)

         if ( lswapinter ) then
            if (lgibbs) then
c     --- Note: acceptance based on only molecules of type imolty
               wratio = ( weight / weiold ) *
     +              ( volins * dble( ncmt(boxrem,imolty)+1 ) / 
     +              ( volrem * dble( ncmt(ibox,imolty) ) ) )
            elseif (lgrand) then
               if (ibox.eq.1) then
c           --- molecule added to box 1
                  wratio = (weight /  weiold ) * 
     +                 volins * B(imolty) / (ncmt(ibox,imolty)) 
               else
c            --- molecule removed from box 1
                  wratio = (weight /  weiold ) * 
     +                 (ncmt(boxrem,imolty)+1)/ (B(imolty)*volrem) 
               endif
            endif
         else
               wratio = weight / weiold
         endif
         if ( wratio .gt. random() ) then
            laccept = .true.
         else
            laccept = .false.
         endif

      else

         do ichoiq = 1,nchoiq(ibox)
            call flucq(0,ibox)
         enddo
         vnewt2 = 0.0d0
         voldt2 = 0.0d0
         do ibox2 = 1, nbox
            vnewt2 = vnewt2 + vbox(ibox2)
            voldt2 = voldt2 + vboxo(ibox2)
         enddo
         deltv = vnewt2 - voldt2

         deltvb = beta * deltv
         if ( deltvb .gt. (2.3d0*softcut) ) then
            laccept = .false.
         elseif ( deltv .le. 0.0d0 ) then
            laccept = .true.
         elseif ( dexp(-deltvb) .gt. random() ) then
            laccept = .true.
         else
            laccept = .false.
         endif

      endif

      if ( laccept ) then

c *** combined move can be accepted now !!!

      else
c *** restore the old energy and old coordinates and ewald sum
         do ibox2 = 1,nbox
            vbox(ibox2) = vboxo(ibox2)
            vinterb(ibox2) = vinterbo(ibox2)
            vintrab(ibox2) = vintrabo(ibox2)
            vextb(ibox2) = vextbo(ibox2)
            velectb(ibox2) = velectbo(ibox2)
            vflucqb(ibox2) = vflucqbo(ibox2)
            if ( mtype .eq. 3 ) then
               vtailb(ibox2) = vtailbo(ibox2)
               vvibb(ibox2) = vvibbo(ibox2)
               vtgb(ibox2) = vtgbo(ibox2)
               vbendb(ibox2) = vbendbo(ibox2)
            endif
         enddo
         do j = 1,iunit
            rxu(i,j) = rxuo(j)
            ryu(i,j) = ryuo(j)
            rzu(i,j) = rzuo(j)
         enddo
         do ip = 1,nchain
            do j = 1,nunit(moltyp(ip))
               qqu(ip,j) = qquo(ip,j)
            enddo
         enddo
         xcm(i) = xcmo
         ycm(i) = ycmo
         zcm(i) = zcmo
         if (lewald) then
c *** restore the reciprocal-space sum
            do ibox2 = 1,nbox
               call recip(ibox2,vdum,vdum,4)
            enddo
            if ( ldielect ) then
c *** restore old dipole moment
               call dipole(ibox,3)
            endif
         endif
      endif    

c      write(iou,*) 'END the optimization of the charge configuration'

      return         
      end





