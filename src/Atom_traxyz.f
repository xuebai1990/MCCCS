      subroutine Atom_traxyz (lx,ly,lz )

c traxyz
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
c    ** makes a translational movement in x,y,or z-direction.         **
c    ** the maximum displacement is controlled by rAtrax(yz) and the  **
c    ** number of successful trial moves is stored in Abstrax(yz).     **
c    ** The attempts are stored in Abntrax(yz)                         **
c    *******************************************************************
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ensemble.inc'
      include 'neigh2.inc'
      include 'system.inc' 
      include 'inputdata.inc'
      include 'bnbsma.inc'
      include 'neigh.inc'

      logical lx,ly,lz,ovrlap,idum,ddum

      logical lneighij,lclu_cmp,lexclude(nmax)

      integer i,ibox,flagon,iunit,j,imolty,icbu,ncount,ic,ip,k
      integer pick_unit, pick_chain
      double precision rx,ry,rz,dchain,ddx,ddy,ddz,random,vnew,vold
     +                 ,vintran,vintrao,deltv,deltvb,disvsq
     +                 ,vintern,vintero,vextn,vexto,rchain
     &                 ,velectn,velecto,vdum
     &                 ,vrecipo,vrecipn   
     & ,velectn_intra,velectn_inter,velecto_intra,velecto_inter
 
      double precision vvibn,vbendn,vtgn,vvibo,vbendo,vtgo

      double precision vewaldn, vewaldo   

      dimension ddum(27)

      logical laccept

C --------------------------------------------------------------------

c      write(2,*) 'start TRAXYZ'
      ovrlap = .false.
c     ***    select a chain at random ***
      rchain  = random()
      do icbu = 1,nmolty
         if ( rchain .lt. pmtrmt(icbu) ) then
            imolty = icbu
            rchain = 2.0d0
         endif
      enddo
      
      if (lgrand) then
c ---    select a chain at random in box 1!
c         (in box 2 is an ideal gas!)
         ibox = 1
         if (ncmt(ibox,imolty).eq.0) then
            if(lx) Abntrax = Abntrax + 1.0d0
            if(ly) Abntray = Abntray + 1.0d0
            if(lz) Abntraz = Abntraz + 1.0d0
            return
         endif
         pick_chain = idint( dble(ncmt(1,imolty))*random() ) + 1
         pick_chain = parbox(pick_chain,1,imolty)
         if ( moltyp(pick_chain) .ne. imolty ) 
     &           write(2,*) 'screwup traxyz'


      else

         dchain = dble(temtyp(imolty))
         pick_chain = int( dchain*random() + 1 )
         pick_chain = parall(imolty,pick_chain)
         ibox = nboxi(pick_chain)

      endif

c *** store number of units of i in iunit ***

      iunit = nunit(imolty)

      pick_unit = int(dble(iunit*random()) + 1 )

c      write(2,*) pick_unit, imolty, pick_chain

      i = pick_chain 

      do j = 1,iunit  
        rxuion(j,1) = rxu(i,j)
        ryuion(j,1) = ryu(i,j)
        rzuion(j,1) = rzu(i,j)
        qquion(j,1) = qqu(i,j)
      enddo

      moltion(1) = imolty

c *** move i ***
      if (lx) then
         rx =  ( 2.0*random() - 1.0d0 ) * Armtrax
         Abntrax = Abntrax + 1.0d0
      else
         rx=0
      endif
      if (ly) then
         ry =  ( 2.0*random() - 1.0d0 ) * Armtray
         Abntray = Abntray + 1.0d0
      else
         ry=0
      endif
      if (lz) then
         rz =  ( 2.0*random() - 1.0d0 ) * Armtraz
         Abntraz = Abntraz + 1.0d0
      else
         rz=0
      endif

      do j = 1,iunit
        if (j .eq. pick_unit) then
           rxuion(j,2) = rxuion(j,1) + rx 
           ryuion(j,2) = ryuion(j,1) + ry 
           rzuion(j,2) = rzuion(j,1) + rz 
           qquion(j,2) = qquion(j,1)
        else
           rxuion(j,2) = rxuion(j,1) 
           ryuion(j,2) = ryuion(j,1) 
           rzuion(j,2) = rzuion(j,1) 
           qquion(j,2) = qquion(j,1)
        endif
      enddo

      moltion(2) = imolty
      
c  *** calculate the energy of i in the new configuration ***
      flagon = 2 
       call Atom_energy(i,imolty, vnew,vintran, vintern,vextn,velectn
     &     ,vewaldn,flagon, ibox,pick_unit, pick_unit,.true.,ovrlap,
     &      .false.
     &     ,vdum,.false.,.false.,vvibn,vbendn,vtgn)
      if (ovrlap) return

c *** calculate the energy of i in the old configuration ***
      flagon = 1
      call Atom_energy(i,imolty,vold,vintrao,vintero,vexto,velecto
     &     ,vewaldo,flagon,ibox,pick_unit, pick_unit,.true.,ovrlap,
     &       .false.
     &     ,vdum,.false.,.false.,vvibo,vbendo,vtgo)

      if (ovrlap) stop 'disaster ovrlap in old conf of TRAXYZ'
      
      if ( lewald .and. lelect(imolty) ) then
         call recip_atom(ibox,vrecipn,vrecipo,1,pick_unit)
         velectn = velectn + vrecipn + vewaldn
         velecto = velecto + vrecipo + vewaldo
         vnew = vnew + vrecipn
         vold = vold + vrecipo
      endif
c *** check for acceptance ***
 
      deltv  = vnew - vold
      deltvb = beta * deltv

c *** For ANES algorithm, do the Fluctuating charge moves.
c *** For time being it will not work for atom disp [Neeraj]****
      if ( lanes ) then
         call anes(i,ibox,ibox,1,laccept,deltv,vintern,vintran,vextn,
     &        velectn,vintero,vintrao,vexto,velecto,vdum,vdum,
     &        vdum,vdum,.false.)
         if ( laccept ) then
            if (lx) Abstrax = Abstrax + 1.0d0
            if (ly) Abstray = Abstray + 1.0d0
            if (lz) Abstraz = Abstraz + 1.0d0
         endif
         return         
      endif
 
      if ( deltvb .gt. (2.3d0*softcut) ) return

      if ( deltv .le. 0.0d0 ) then
c        --- accept move
      elseif ( dexp(-deltvb) .gt. random() ) then
c        --- accept move
      else
c        --- move rejected
         return
      endif

c      write(2,*) 'TRAXYZ accepted i',i
      vbox(ibox)     = vbox(ibox) + deltv
      vinterb(ibox)  = vinterb(ibox) + (vintern - vintero)
      vintrab(ibox)  = vintrab(ibox) + (vintran - vintrao)
      vextb(ibox)    = vextb(ibox)   + (vextn   - vexto)
      velectb(ibox)   = velectb(ibox)  + (velectn - velecto)

      vtgb(ibox) = vtgb(ibox) + (vtgn-vtgo)
      vbendb(ibox) = vbendb(ibox) + (vbendn-vbendo)
      vvibb(ibox) = vvibb(ibox) + (vvibn-vvibo)

      
      if (lx) rxu(i,pick_unit) = rxuion(pick_unit,2)
      if (ly) ryu(i,pick_unit) = ryuion(pick_unit,2)
      if (lz) rzu(i,pick_unit) = rzuion(pick_unit,2)

      if (lewald .and. lelect(imolty)) then
c *** update reciprocal-space sum
         call recip_atom(ibox,vdum,vdum,2,pick_unit)
      endif

      if ( ldielect ) then
c *** update the dipole term
         call dipole(ibox,1)
      endif


c *** update chain center of mass

      call ctrmas(.false.,ibox,i,10)

      if (lx) Abstrax = Abstrax + 1.0d0
      if (ly) Abstray = Abstray + 1.0d0
      if (lz) Abstraz = Abstraz + 1.0d0

      if ( licell .and. (ibox.eq.boxlink)) then
c     --- update linkcell list
         call linkcell(2,i,vdum,vdum,vdum,ddum)
      endif

      if ( lneigh ) then
c *** check for update of near neighbour bitmap ***
c *** check for headgroup ***
         disvec(1,i,1) = disvec(1,i,1) + rx
         disvsq = disvec(1,i,1) * disvec(1,i,1) +
     +        disvec(1,i,2) * disvec(1,i,2) +
     +        disvec(1,i,3) * disvec(1,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
c *** check for last unit ***
         disvec(2,i,1) = disvec(2,i,1) + rx
         disvsq = disvec(2,i,1) * disvec(2,i,1) +
     +        disvec(2,i,2) * disvec(2,i,2) +
     +        disvec(2,i,3) * disvec(2,i,3)
         if (disvsq .gt. upnnsq) call updnn( i )
      endif

      if ( lneighbor ) then
         
         do 10 ic = 1, neigh_cnt(i)
            j = neighbor(ic,i)
c            write(2,*) ic,i,'j:',j
            do ip = 1,neigh_cnt(j)
               if ( neighbor(ip,j) .eq. i ) then
                  neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                  neigh_cnt(j) = neigh_cnt(j)-1
                  goto 10
               endif
            enddo
 10      continue
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

c      write(2,*) 'end TRAXYZ',i

      return
      end
