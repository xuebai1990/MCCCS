      subroutine rotxyz (lx,ly,lz )

c rotxyz
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

c    *******************************************************************
c    ** makes a rotational movement around "x" space-fixed axis.      **
c    ** the maximum displacement is controlled by rmrotx and the      **
c    ** number of successful rotation is given by bsrotx.             **
c    **                                                               **
c    ** rotxyz chooses one of the three space-fixed axes at random    **
c    ** and rotates the molecule around this axis by dgamma radians.  **
c    ** the maximum angular displacement is dgamax.                   **
c    *******************************************************************
      implicit none

c *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ensemble.inc'
      include 'neigh2.inc'
      include 'system.inc' 
      include 'inputdata.inc'
      include 'bnbsma.inc'
      include 'neigh.inc'
      include 'ipswpar.inc'
      include 'eepar.inc'

      logical::lx,ly,lz,ovrlap,lneighij,lclu_cmp,lexclude(nmax)
      integer::i,ibox,flagon,iunit,j,imolty,iuroty,icbu,ic,ip,k
      real(8)::rx,ry,rz,dchain,rchain,random,vnew,vold
     +                ,vintrao,dgamma,rxorig,ryorig,rzorig
     +                ,rxnew2,rynew2,rznew2,vintran,disvsq,deltv
     +                ,deltvb,vintern,vintero,vextn,vexto,vdum
     &                ,velectn,velecto
     & ,velectn_intra,velectn_inter,velecto_intra,velecto_inter
c *** further variable definitions
      real(8)::cosdg, sindg, rmrot
      real(8)::vrecipn,vrecipo


      logical::laccept

C --------------------------------------------------------------------

c      write(iou,*) 'start ROTXYZ'
      ovrlap = .false.
      if (lgrand) then
c ---    select a chain at random in box 1!
c         (in box 2 is an ideal gas!)
         ibox = 1
         rchain  = random()
         do icbu = 1,nmolty
            if ( rchain .lt. pmromt(icbu) ) then
               imolty = icbu
               rchain = 2.0d0
            endif
         enddo
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) bnrotx(imolty,ibox) = bnrotx(imolty,ibox) + 1.0d0
            if(ly) bnroty(imolty,ibox) = bnroty(imolty,ibox) + 1.0d0
            if(lz) bnrotz(imolty,ibox) = bnrotz(imolty,ibox) + 1.0d0
            return
         endif

 10      dchain = dble(temtyp(imolty))
         i = int( dchain*random() + 1 )
         i = parall(imolty,i)
         if (nboxi(i) .ne. 1) goto 10
         iuroty = iurot(imolty)
      else 
c ***    select a chain type at random ***
         rchain  = random()
         do icbu = 1,nmolty
            if ( rchain .lt. pmromt(icbu) ) then
               imolty = icbu
               rchain = 2.0d0
            endif
         enddo

         if ((lexpee).and.(imolty.ge.nmolty1))
     &      imolty = ee_moltyp(mstate)
                                                                                
         if (temtyp(imolty).eq.0) return
         
         dchain = dble(temtyp(imolty))
         i = int( dchain*random() + 1 )
         i = parall(imolty,i)

         ibox = nboxi(i)
         iuroty = iurot(imolty)
      endif

ckea 6/4/09 -- for multiple rotation centers
      if(iuroty.lt.0) then
         if(nrotbd(imolty).gt.1) then
            rchain = random()
            do icbu = 1,nrotbd(imolty)
               if( rchain .lt. pmrotbd(icbu,imolty)) then
                  iuroty = irotbd(icbu,imolty)
               endif
            enddo
         else
            iuroty = irotbd(1,imolty)
         endif
      endif

c *** store number of units of i in iunit ***
      iunit = nunit(imolty)

c *** store current positions in old-new array #1# ***
      do j = 1, iunit
         rxuion(j,1) = rxu(i,j)
         ryuion(j,1) = ryu(i,j)
         rzuion(j,1) = rzu(i,j)
         qquion(j,1) = qqu(i,j)
      enddo
      moltion(1) = imolty
 
c *** choose a random angular displacement ***
 
      if (lx) then
         rmrot = rmrotx(imolty,ibox)
         bnrotx(imolty,ibox) = bnrotx(imolty,ibox) + 1.0d0
      else if (ly) then
         rmrot = rmroty(imolty,ibox)
         bnroty(imolty,ibox) = bnroty(imolty,ibox) + 1.0d0
      else if (lz) then
         rmrot = rmrotz(imolty,ibox)
         bnrotz(imolty,ibox) = bnrotz(imolty,ibox) + 1.0d0
      endif
      dgamma = ( 2.0d0*random() - 1.0d0 ) * rmrot
 
c *** set up the rotation marix ***
 
      cosdg = dcos( dgamma )
      sindg = dsin( dgamma )
 
c *** Determine the rotation coordinates ***
      if (iuroty .eq. 0) then
c *** Use the center of mass for rotation
         rxorig = xcm(i)
         ryorig = ycm(i)
         rzorig = zcm(i)
      else
c *** Use iurot for rotation
         rxorig = rxuion(iuroty,1)
         ryorig = ryuion(iuroty,1)
         rzorig = rzuion(iuroty,1)
      endif

c      write(iou,*) 'before rotating'
c      write(iou,*) xcm(i),ycm(i),zcm(i)
c      write(iou,*) rxu(i,1),ryu(i,1),rzu(i,1)


      if (lx) then 
 
C ***    ROTATE UNITS OF I AROUND X-AXIS ***
 
         do  j = 1, iunit
            ry = ryuion(j,1) - ryorig
            rz = rzuion(j,1) - rzorig
            rynew2 = cosdg * ry + sindg * rz
            rznew2 = cosdg * rz - sindg * ry
               
            rxuion(j,2) = rxuion(j,1)
            ryuion(j,2) = ryorig + rynew2
            rzuion(j,2) = rzorig + rznew2
            qquion(j,2) = qquion(j,1) 
         enddo
      endif
      if (ly) then 
 
C ***    ROTATE UNITS OF I AROUND y-AXIS ***
 
         do  j = 1, iunit
            rx = rxuion(j,1) - rxorig
            rz = rzuion(j,1) - rzorig
            rxnew2 = cosdg * rx - sindg * rz
            rznew2 = cosdg * rz + sindg * rx

            rxuion(j,2) = rxorig + rxnew2
            ryuion(j,2) = ryuion(j,1)
            rzuion(j,2) = rzorig + rznew2
            qquion(j,2) = qquion(j,1)
         enddo
      endif
      if (lz) then 
 
C ***    ROTATE UNITS OF I AROUND z-AXIS ***
 
         do  j = 1, iunit
            rx = rxuion(j,1) - rxorig
            ry = ryuion(j,1) - ryorig
            rxnew2 = cosdg * rx + sindg * ry
            rynew2 = cosdg * ry - sindg * rx

            rxuion(j,2) = rxorig + rxnew2
            ryuion(j,2) = ryorig + rynew2
            rzuion(j,2) = rzuion(j,1)
            qquion(j,2) = qquion(j,1)
         enddo
      endif 
      moltion(2) = imolty
 
c  *** calculate the energy of i in the new configuration ***

      flagon = 2
      call energy(i,imolty, vnew,vintran,vintern,vextn,velectn,vdum
     &     ,flagon, ibox,1,iunit,.false.,ovrlap,.false.,vdum,
     +     .false.,.false.)
      if (ovrlap) return
 
c *** calculate the energy of i in the old configuration ***
      flagon = 1
      call energy(i,imolty, vold,vintrao,vintero,vexto,velecto,vdum
     &     ,flagon,ibox, 1, iunit,.false.,ovrlap,.false.,vdum,
     &     .false.,.false.)

      if (ovrlap) stop 'disaster- overlap for old conf in ROTXYZ'
      if ( lewald .and. lelect(imolty) ) then
         call recip(ibox,vrecipn,vrecipo,1)
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
           vrecipn  =  (1.0d0-(1.0d0-etais)*lambdais)*vrecipn
           vrecipo  =  (1.0d0-(1.0d0-etais)*lambdais)*vrecipo
         elseif (lstageb) then
           vrecipn  =  etais*vrecipn
           vrecipo  =  etais*vrecipo
         elseif (lstagec) then
           vrecipn  =  (etais+(1.0d0-etais)*lambdais)*vrecipn
           vrecipo  =  (etais+(1.0d0-etais)*lambdais)*vrecipo
         endif
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      endif

c *** check for acceptance ***
 
      deltv  = vnew - vold
      deltvb = beta * deltv

c *** For ANES algorithm, do the Fluctuating charge moves.

      if ( lanes ) then
         call anes(i,ibox,ibox,2,laccept,deltv,vintern,vintran,vextn,
     &        velectn,vintero,vintrao,vexto,velecto,vdum,vdum,vdum,
     &        vdum,.false.)
         if (laccept) then
            if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0d0
            if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0d0
            if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0d0
         endif
         return         
      endif

      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
c        --- move accepted
      elseif ( dexp(-deltvb) .gt. random() ) then
c        --- move accepted
      else
c        --- move rejected
         return
      endif

c      write(iou,*) 'ROTXYZ accepted',i
      vbox(ibox) = vbox(ibox) + deltv
      vinterb(ibox)  = vinterb(ibox) + (vintern -  vintero)
      vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
      vextb(ibox)    = vextb(ibox)   + (vextn-vexto)
      velectb(ibox)  = velectb(ibox) + (velectn-velecto)
      vipswb(ibox) = vipswb(ibox) + (vipswn-vipswo)
      vwellipswb(ibox) = vwellipswb(ibox) + (vwellipswn-vwellipswo)
      vipsw = vipswb(ibox)
      vwellipsw = vwellipswb(ibox)
 
      do j = 1, iunit
         rxu(i,j) = rxuion(j,2)
         ryu(i,j) = ryuion(j,2)
         rzu(i,j) = rzuion(j,2)
      enddo

      if (lewald .and. lelect(imolty)) then
c *** update reciprocal-space sum
         call recip(ibox,vdum,vdum,2)
      endif
 
      if ( ldielect ) then
           call dipole(ibox,1)
      endif
 

c *** update center of mass
      call ctrmas(.false.,ibox,i,2)

      if (lx) bsrotx(imolty,ibox) = bsrotx(imolty,ibox) + 1.0d0
      if (ly) bsroty(imolty,ibox) = bsroty(imolty,ibox) + 1.0d0
      if (lz) bsrotz(imolty,ibox) = bsrotz(imolty,ibox) + 1.0d0

      if ( lneigh ) then
c *** check for update of near neighbour bitmap ***
c *** check for headgroup ***
         disvec(1,i,1) = disvec(1,i,1) + rxuion(1,2) - rxuion(1,1)
         disvec(1,i,2) = disvec(1,i,2) + ryuion(1,2) - ryuion(1,1)
         disvec(1,i,3) = disvec(1,i,3) + rzuion(1,2) - rzuion(1,1)
         disvsq = disvec(1,i,1) * disvec(1,i,1) +
     +        disvec(1,i,2) * disvec(1,i,2) +
     +        disvec(1,i,3) * disvec(1,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
c *** check for last unit ***
         disvec(2,i,1) = disvec(2,i,1)+rxuion(iunit,2)-rxuion(iunit,1)
         disvec(2,i,2) = disvec(2,i,2)+ryuion(iunit,2)-ryuion(iunit,1)
         disvec(2,i,3) = disvec(2,i,3)+rzuion(iunit,2)-rzuion(iunit,1)
         disvsq = disvec(2,i,1) * disvec(2,i,1) +
     +        disvec(2,i,2) * disvec(2,i,2) +
     +        disvec(2,i,3) * disvec(2,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
      endif

      if ( lneighbor ) then
c         write(iou,*) 'in rotxyz:',i,neigh_cnt(i)
         do 11 ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                  neigh_cnt(j) = neigh_cnt(j)-1
                  goto 11
               endif
            enddo
 11      continue
         neigh_cnt(i) = neigh_icnt
         do ic = 1,neigh_icnt
            j = neighi(ic)
            neighbor(ic,i)=j
            lneighij = .false.
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  lneighij = .true.
               endif
            enddo
            if ( .not. lneighij ) then
               neigh_cnt(j) = neigh_cnt(j)+1
               neighbor(neigh_cnt(j),j) = i
            endif
         enddo
      endif
 
c      write(iou,*) 'end ROTXYZ'

      return
      end

