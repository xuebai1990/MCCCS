      subroutine flucq (ichoice,boxi)

c flucq
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
c    ** select one (or two in charge transfer) molecule and displace  **
c    ** the charge magnitude of charge sites on selected molecule(s)  **
c    ** according to charge nutrality and preferential strategy.      **
c    ** rewritten by Bin Chen at 6-25-99.                             **
c    *******************************************************************
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ensemble.inc'
      include 'system.inc' 
      include 'inputdata.inc'
      include 'bnbsma.inc'
      include 'ewaldsum.inc'
      include 'fepsi.inc'
      include 'poten.inc'
      
      logical linterqt,ovrlap
      integer i,ibox,iunit,j,imolty,icbu,mainunit,ic,ii
     &     ,ncount,ichoiq,ichoice,qunit,boxi,flagon
      double precision dchain,random,vnew,vold,deltv,deltvb
     &  ,velectn,velecto,vflucqn,vflucqo
     &  ,dispbig,displit,vintern,vintero,vewaldn,vewaldo
      double precision qion
      double precision velectn_intra,velectn_inter,velecto_intra,
     &  velecto_inter 
      dimension qion(numax)
      double precision vrecipn,vrecipo,sumr,sumi,arg
      dimension sumr(2),sumi(2)
      
      integer maini,mainj,jchain
      double precision qionj(numax),vinterjo,vflucqjo,velectjo,
     &     voldj,vinterjn,vflucqjn,velectjn,vnewj,corr,rij,erfunc,
     &     rxuij,ryuij,rzuij,vdum,vewaldjn,vewaldjo
      double precision qoldj2,vnewi,velectni,vinterni,voldi,
     &     velectoi,vinteroi

C --------------------------------------------------------------------

c      write(2,*) 'start FLUCQ'
c ***    select a chain at random ***
      dchain  = random()
      do icbu = 1,nmolty
         if ( dchain .lt. pmfqmt(icbu) ) then
            imolty = icbu
            dchain = 2.0d0
         endif
      enddo
      
      if (lgrand) then
c ---    select a chain at random in box 1!
c         (in box 2 is an ideal gas!)
         ibox = 1
         if (nchbox(ibox).eq.0) then
            bnflcq(imolty,ibox) = bnflcq(imolty,ibox) + 1.0d0
            return
         endif
         i = dint( dble(ncmt(1,imolty))*random() ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(2,*) 'screwup'

      elseif ( lanes ) then

         if ( ichoice .eq. -1 ) then
c *** equal probability charge moves in two boxes (for swap)
c *** preferential charge moves according to favor
            dchain = dble(temtyp(imolty))
 77         i = int( dchain*random() ) + 1
            i = parall(imolty,i)
            if ( random() .gt. favor(i)) goto 77
            ibox = nboxi(i)
         elseif ( ichoice .eq. -2 ) then
c *** equal probability charge moves in two boxes (for swap)
c *** preferential charge moves according to favor2
            dchain = dble(temtyp(imolty))
 66         i = int( dchain*random() ) + 1
            i = parall(imolty,i)
            if ( random() .gt. favor2(i)) goto 66
            ibox = nboxi(i)
         elseif ( ichoice .eq. 0 ) then
c *** equal probability charge moves in separate box (for trans,rot)
 88         dchain = dble(temtyp(imolty))
            i = int( dchain*random() ) + 1
            i = parall(imolty,i)
            ibox = nboxi(i)
            if ( ibox .ne. boxi ) goto 88
         elseif (ichoice .eq. 2) then
            dchain = dble(temtyp(imolty))
            i = int( dchain*random() ) + 1
            i = parall(imolty,i)
            ibox = nboxi(i)
         endif
      else
         dchain = dble(temtyp(imolty))
         i = int( dchain*random() ) + 1
         i = parall(imolty,i)
         ibox = nboxi(i)
         
      endif

c *** For charge transfer case, select a second molecule to perform the
c *** fluctuating-charge move. The second molecule should be in the same
c *** box as chain i.

      If ( lqtrans(imolty) ) then
         if (lgrand) then
c ---    select a chain at random in box 1!
c     (in box 2 is an ideal gas!)
            ibox = 1
            if (nchbox(ibox).eq.0) then
               bnflcq(imolty,ibox) = bnflcq(imolty,ibox) + 1.0d0
               return
            endif
            jchain = dint( dble(ncmt(1,imolty))*random() ) + 1
            jchain = parbox(jchain,1,imolty)
            if ( moltyp(jchain) .ne. imolty ) write(2,*) 'screwup'

         elseif ( lanes ) then

            if ( ichoice .eq. -1 ) then
c *** equal probability charge moves in two boxes (for swap)
c *** preferential charge moves according to favor
               dchain = dble(temtyp(imolty))
 70            jchain = int( dchain*random() ) + 1
               jchain = parall(imolty,jchain)
               if ( random() .gt. favor(jchain)) goto 70
               if ( nboxi(jchain) .ne. ibox ) goto 70 
            elseif ( ichoice .eq. -2 ) then
c *** equal probability charge moves in two boxes (for swap)
c *** preferential charge moves according to favor
               dchain = dble(temtyp(imolty))
 60            jchain = int( dchain*random() ) + 1
               jchain = parall(imolty,jchain)
               if ( random() .gt. favor2(jchain)) goto 60
               if ( nboxi(jchain) .ne. ibox ) goto 60
            elseif ( ichoice .eq. 0 ) then
c *** equal probability charge moves in separate box (for trans,rot)
 80            dchain = dble(temtyp(imolty))
               jchain = int( dchain*random() ) + 1
               jchain = parall(imolty,jchain)
               if ( nboxi(jchain) .ne. boxi ) goto 80
            elseif (ichoice .eq. 2) then
 90            dchain = dble(temtyp(imolty))
               jchain = int( dchain*random() ) + 1
               jchain = parall(imolty,jchain)
               if ( nboxi(jchain) .ne. ibox ) goto 90
            endif
         else
 100        dchain = dble(temtyp(imolty))
            jchain = int( dchain*random() ) + 1
            jchain = parall(imolty,jchain)
            if ( nboxi(jchain) .ne. ibox ) goto 100
         endif
      endif

      if ( lqtrans(imolty) .and. (jchain .ne. i ) ) then
         linterqt = .true.
      else
         linterqt = .false.
      endif

c *** store number of units of i in iunit ***

      iunit = nunit(imolty)
c *** count the number of charge sites
      qunit = 0
      do j = 1, iunit
         if ( lqchg(ntype(imolty,j)) ) qunit = qunit + 1
      enddo

c     --- store the charges in qion
      do  j = 1, iunit
         qion(j)   = qqu(i,j)         
      enddo

c *** calculate the polariztion energy of i in the old configuration ***
 
      call charge(i, qion, vflucqo, vewaldo)
      
      If ( linterqt ) then
         do j = 1,iunit
            qionj(j) = qqu(jchain,j)
         enddo

c *** calculate the polarization energy of jchain in the old configuration ***

         call charge(jchain, qionj, vflucqjo, vewaldjo)

         vflucqo = vflucqo + vflucqjo
         vewaldo = vewaldo + vewaldjo
      endif

c --- Choose one of the units as the main charge transfer site

 30   mainunit = int( dble(iunit)*random() ) + 1
c *** for unit which is not a charge site 
      if ( .not. lqchg(ntype(imolty,mainunit)) ) goto 30
      bnflcq(imolty,ibox) = bnflcq(imolty,ibox) + 1.0d0
      dispbig = ( 2.0d0*random() - 1.0d0 )*rmflcq(imolty,ibox)

      if ( linterqt ) then
c *** For charge transfer case, i molecule increases by dispbig and
c *** jchain molecule decreases by dispbig. 
c *** correction for the reptition of the calculation of the 
c *** coulombic real space term between maini and mainj
         maini = mainunit
 32      mainunit = int( dble(iunit)*random() ) + 1
c *** for unit which is not a charge site
         if ( .not. lqchg(ntype(imolty,mainunit)) ) goto 32
         mainj = mainunit
         corr = (qion(maini)-qionj(mainj))*dispbig 
     &        - qion(maini)*qionj(mainj)
         qion(maini) = qion(maini) + dispbig
         qionj(mainj) = qionj(mainj) - dispbig
         corr = corr + qion(maini)*qionj(mainj)
         rxuij = rxu(i,maini)-rxu(jchain,mainj)
         ryuij = ryu(i,maini)-ryu(jchain,mainj)
         rzuij = rzu(i,maini)-rzu(jchain,mainj)
         if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
         rij = dsqrt( rxuij*rxuij + ryuij*ryuij + rzuij*rzuij )
         corr = qqfact*erfunc(calp(ibox)*rij)*corr/rij

c *** ??? problems here for charge transfer

      else

         displit = (-dispbig) / (dble(qunit)-1.0d0)

c --- displace the charges on the molecule
      
         do j = 1, iunit
            if ( j .eq. mainunit ) then
               qion(j) = qion(j) + dispbig
            elseif ( lqchg(ntype(imolty,j)) ) then
               qion(j) = qion(j) + displit
            endif
            rxuion(j,1) = rxu(i,j)
            ryuion(j,1) = ryu(i,j)
            rzuion(j,1) = rzu(i,j)
            qquion(j,1) = qqu(i,j)
            rxuion(j,2) = rxuion(j,1)
            ryuion(j,2) = ryuion(j,1)
            rzuion(j,2) = rzuion(j,1)
            qquion(j,2) = qion(j)
         enddo
      endif
      moltion(1) = imolty
      moltion(2) = imolty

c  *** calculate the polarization energy of i in the new configuration ***
 
      call charge(i, qion, vflucqn, vewaldn)

      if ( linterqt ) then
         
c  *** calculate the polarization energy of jchain in the new configuration ***

         call charge(jchain, qionj, vflucqjn, vewaldjn)
         vflucqn = vflucqn + vflucqjn
         vewaldn = vewaldn + vewaldjn

c *** calculate the energy of i in the new configuration ***
         flagon = 2 
         rxuion(maini,flagon) = rxu(i,maini)
         ryuion(maini,flagon) = ryu(i,maini)
         rzuion(maini,flagon) = rzu(i,maini)
         qquion(maini,flagon) = qion(maini)
         qoldj2 = qqu(jchain,mainj)
         qqu(jchain,mainj) = qionj(mainj)

         call energy(i,imolty, vnew,vdum, vintern,vdum,velectn
     &     ,vdum,flagon, ibox,maini, maini,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)
         vnewi = vnew + vflucqn + vewaldn
         velectni = velectn + vewaldn
         vinterni = vintern

c *** calculate the energy of i in the old configuration ***
         flagon = 1
         rxuion(maini,flagon) = rxu(i,maini)
         ryuion(maini,flagon) = ryu(i,maini)
         rzuion(maini,flagon) = rzu(i,maini)
         qquion(maini,flagon) = qqu(i,maini)

         call energy(i,imolty,vold,vdum,vintero,vdum,velecto
     &        ,vdum,flagon,ibox,maini, maini,.false.,ovrlap,.false.
     &        ,vdum,.false.,.false.)
         
         voldi = vold + vflucqo + vewaldo
         velectoi = velecto + vewaldo
         vinteroi = vintero

c  *** restore the old charges

         qqu(jchain,mainj) = qoldj2

c *** calculate the energy of jchain in the new configuration ***

         flagon = 2 
         rxuion(mainj,flagon) = rxu(jchain,mainj)
         ryuion(mainj,flagon) = ryu(jchain,mainj)
         rzuion(mainj,flagon) = rzu(jchain,mainj)
         qquion(mainj,flagon) = qionj(mainj)
         call energy(jchain,imolty, vnew,vdum, vintern,vdum,velectn
     &     ,vdum,flagon, ibox,mainj, mainj,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)

         vnew = vnew + vnewi
         velectn = velectn + velectni
         vintern = vintern + vinterni

c *** calculate the energy of jchain in the old configuration ***

         flagon = 1
         rxuion(mainj,flagon) = rxu(jchain,mainj)
         ryuion(mainj,flagon) = ryu(jchain,mainj)
         rzuion(mainj,flagon) = rzu(jchain,mainj)
         qquion(mainj,flagon) = qqu(jchain,mainj)

         call energy(jchain,imolty, vold,vdum, vintero,vdum,velecto
     &     ,vdum,flagon, ibox,mainj, mainj,.false.,ovrlap,.false.
     &     ,vdum,.false.,.false.)

         vold = vold + voldi
         velecto = velecto + velectoi
         vintero = vintero + vinteroi
  
c *** prepare the rxuion etc for ewald sum and dielectric constant
c *** and for energy calculation

         if ( linterqt ) then
            rxuion(1,1) = rxu(i,maini)
            ryuion(1,1) = ryu(i,maini)
            rzuion(1,1) = rzu(i,maini)
            qquion(1,1) = qqu(i,maini)
            rxuion(2,1) = rxu(jchain,mainj)
            ryuion(2,1) = ryu(jchain,mainj)
            rzuion(2,1) = rzu(jchain,mainj)
            qquion(2,1) = qqu(jchain,mainj)
            rxuion(1,2) = rxu(i,maini)
            ryuion(1,2) = ryu(i,maini)
            rzuion(1,2) = rzu(i,maini)
            qquion(1,2) = qion(maini)
            rxuion(2,2) = rxu(jchain,mainj)
            ryuion(2,2) = ryu(jchain,mainj)
            rzuion(2,2) = rzu(jchain,mainj)
            qquion(2,2) = qionj(mainj)
            
            do j = 3, iunit
               rxuion(j,1) = rxuion(j,2)
               ryuion(j,1) = ryuion(j,2)
               rzuion(j,1) = rzuion(j,2)
               qquion(j,1) = qquion(j,2)
            enddo
         endif

      else

c *** calculate the energy of i in the new configuration ***
         flagon = 2 
         call energy(i,imolty, vnew,vdum, vintern,vdum,velectn
     &        ,vdum,flagon, ibox,1, iunit,.false.,ovrlap,.false.
     &        ,vdum,.false.,.false.)
         if (ovrlap) return
         vnew = vnew + vflucqn + vewaldn
         velectn = velectn + vewaldn

c *** calculate the energy of i in the old configuration ***
         flagon = 1
         call energy(i,imolty,vold,vdum,vintero,vdum,velecto
     &        ,vdum,flagon,ibox,1, iunit,.false.,ovrlap,.false.
     &        ,vdum,.false.,.false.)
         vold = vold + vflucqo + vewaldo
         velecto = velecto + vewaldo

      endif

c *** Begin Ewald-sum correction
      if ( lewald ) then
         call recip(ibox,vrecipn,vrecipo,1)
         vnew = vnew + vrecipn
         vold = vold + vrecipo
         velectn = velectn + vrecipn
         velecto = velecto + vrecipo
      endif
         

c *** check for acceptance ***
 
      deltv  = vnew - vold
c --- use the thermostat temperature instead of real temp
      deltvb = fqbeta * deltv

c      if ( deltv .lt. -100.0d0) then
c         write(2,*) i,favor(i),deltv
c      endif 
      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
c        accept move
      elseif ( dexp(-deltvb) .gt. random() ) then
c        accept move
      else
         return
      endif

      vbox(ibox)     = vbox(ibox) + deltv
      velectb(ibox)  = velectb(ibox)  + (velectn - velecto)
      vflucqb(ibox)  = vflucqb(ibox) + (vflucqn - vflucqo)
      vinterb(ibox) = vinterb(ibox) + (vintern - vintero)
c      write(2,*) 'this move has been accepted!!!'
      do j = 1,iunit
         qqu(i,j) = qion(j)
         if ( linterqt ) 
     &        qqu(jchain,j) = qionj(j)
      enddo
c --- update the reciprocal-space sum
      if ( lewald ) then
         call recip(ibox,vdum,vdum,2)
      endif

      if ( ldielect ) then
          call dipole(ibox,1)
      endif
 


      bsflcq(imolty,ibox) = bsflcq(imolty,ibox) + 1.0d0

c      write(2,*) 'end FLUCQ'

      return
      end


