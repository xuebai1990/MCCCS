      subroutine prvolume

c prvolume
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

c    *************************************************************
c    ** makes an isotropic volume change under const. pressure  **
c    ** the maximum change is controlled by rmtrax and the      **
c    ** number of successful trial moves is stored in bsvol.    **
c    *************************************************************
c
c --- perform change of the volume: random walk in ln(vol)
c
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'ensemble.inc'
      include 'system.inc' 
      include 'inputdata.inc'
      include 'bnbsma.inc'
      include 'ewaldsum.inc'
      include 'neigh.inc'
      include 'cell.inc'
ckea
      include 'garofalini.inc'

      logical::ovrlap,lvol,lx,ly,lz,la,lb,lc

      integer::i,j,imolty,ic,ncount,jhmat,numvecto
      real(8)::bxo,byo,bzo,hbxo,hbyo,hbzo
      real(8)::volo,vboxo,dfac,
     +                 vintero,vtailo,vexto,velecto,vflucqo
      real(8)::kxo(vectormax),kyo(vectormax),
     +                 kzo(vectormax),prefacto(vectormax)
      real(8)::voln,vboxn,
     +                 vintern,vtailn,vextn,velectn,vflucqn
      real(8)::rxuo(nmax,numax),ryuo(nmax,numax)
     &                ,rzuo(nmax,numax),qquo(nmax,numax)
      real(8)::random,df,dx0,dy0,dz0,dx,dy,dz,v,dele,vdum
     &                 ,vinter,vtail,vext,vminim,velect,vflucq
     &                 ,velect_intra,velect_inter
      real(8)::xcmo(nmax),ycmo(nmax),zcmo(nmax),calpo
     &                 ,xcmoi,ycmoi,zcmoi,xcmi,ycmi,zcmi
     &                 ,xcmj,ycmj,zcmj,rxuij,ryuij,rzuij
      real(8)::rbcut,rbox
      real(8)::w(3),vx,vy,vz
      integer::boxvch,ichoiq,ibox
      integer::neigho_cnt(nmax),neigho(100,nmax)
      real(8)::hmato(9),hmatio(9)

      integer::idum
      real(8)::lddum,lddum2(27),min_boxl
c KEA
      real(8)::v3o,v3n

C --------------------------------------------------------------------

c      write(iou,*) 'start VOLUME of LNPTGIBBS'

c     Select a box at  random to change the volume of box


      lx = .false.
      ly = .false.
      lz = .false.


      rbox = random()
      do ibox = 1,nbox
         if (rbox .lt. pmvlmt(ibox) ) then
            boxvch=ibox
            goto 99
         endif
      enddo

      
 99   bnvol(boxvch) = bnvol(boxvch) + 1.0d0

      if ( lsolid(boxvch) ) then
c *** volume move independently in x, y, z directions
         rbox = random()
         if ( rbox .le. pmvolx ) then
            lx = .true.
            ly = .false.
            lz = .false.
         elseif ( rbox .le. pmvoly ) then
            lx = .false.
            ly = .true.
            lz = .false.
         else
            lx = .false.
            ly = .false.
            lz = .true.
         endif
         if ( .not. lrect(boxvch) ) then
            do i = 1,9
               hmato(i) = hmat(boxvch,i)
               hmatio(i) = hmati(boxvch,i)
            enddo
         endif
      endif

c --- store old box lengths, energy, configuration etc

c---- in the box of volume change ----
      bxo    = boxlx(boxvch)
      byo    = boxly(boxvch)
      if ( lpbcz ) bzo = boxlz(boxvch)

      if (lsolid(boxvch) .and. .not. lrect(boxvch)) then
         volo = cell_vol(boxvch)
      else
         if ( lpbcz ) then
            volo   = bxo*byo*bzo
         else
            volo   = bxo*byo
         endif
      endif

      vboxo    = vbox(boxvch)
      vintero  = vinterb(boxvch)
      vtailo   = vtailb(boxvch)
      vexto    = vextb(boxvch)
      velecto  = velectb(boxvch)
      vflucqo = vflucqb(boxvch)
ckea
      v3o = v3garob(boxvch)
c --- store neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_o(i) = neigh_cnt(i)
            do j=1,neigh_o(i)
               neighboro(j,i) = neighbor(j,i)
               ndijo(j,i) = ndij(j,i)
               nxijo(j,i) = nxij(j,i)
               nyijo(j,i) = nyij(j,i)
               nzijo(j,i) = nzij(j,i)
            enddo
         enddo
      endif

c --- store old k vectors and reciprocal sum
      if ( lewald ) then
         calpo = calp(boxvch)
         ncount = numvect(boxvch)
         numvecto = numvect(boxvch)
         call recip(boxvch,vdum,vdum,3)
         do ic = 1,ncount
            kxo(ic) = kx(ic,boxvch)
            kyo(ic) = ky(ic,boxvch)
            kzo(ic) = kz(ic,boxvch)
            prefacto(ic) = prefact(ic,boxvch)
         enddo
      endif
      do i = 1, nchain
c ----- Check if the chain i is in the correct box
         if (nboxi(i). eq. boxvch) then

            imolty = moltyp(i)
            xcmo(i) = xcm(i)
            ycmo(i) = ycm(i)
            if (lpbcz) zcmo(i) = zcm(i)
            do j = 1, nunit(imolty)
               rxuo(i,j) = rxu(i,j)
               ryuo(i,j) = ryu(i,j)
               if ( lpbcz ) rzuo(i,j) = rzu(i,j)
               qquo(i,j) = qqu(i,j)
            enddo
            if (lneighbor) then
               neigho_cnt(i) = neigh_cnt(i)
               do j = 1,neigho_cnt(i)
                  neigho(j,i)=neighbor(j,i)
               enddo
            endif
         endif
      enddo

      if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then
c *** select one of the cell edge
         rbox = 3.0d0*random()
         if ( lx ) then
            if ( rbox .le. 1.0d0 ) then
               la = .true. 
               lb = .false.
               lc = .false.
            elseif (rbox .le. 2.0d0 ) then
               la = .false.
               lb = .true.
               lc = .false.
            else
               la = .false.
               lb = .false.
               lc = .true.
            endif
         elseif ( ly ) then
            la = .false.
            if ( rbox .le. 1.5d0 ) then
               lb = .true.
               lc = .false.
            else
               lb = .false.
               lc = .true.
            endif
         else
            la = .false.
            lb = .false.
            lc = .true.
         endif

         if ( la ) then
            if ( lx ) jhmat = 1
            if ( ly ) jhmat = 2
            if ( lz ) jhmat = 3
         elseif ( lb ) then
            if ( lx ) jhmat = 4
            if ( ly ) jhmat = 5
            if ( lz ) jhmat = 6
         elseif ( lc ) then
            if ( lx ) jhmat = 7
            if ( ly ) jhmat = 8
            if ( lz ) jhmat = 9
         endif

         hmat(boxvch,jhmat) = hmat(boxvch,jhmat) + rmhmat(boxvch,jhmat)*
     +        ( 2.0d0*random() - 1.0d0 )
         bnhmat(boxvch,jhmat) = bnhmat(boxvch,jhmat) + 1.0d0
          
         call matops(boxvch)
            
         voln = cell_vol(boxvch)

         w(1) = min_width(boxvch,1)
         w(2) = min_width(boxvch,2)
         w(3) = min_width(boxvch,3)

         rbcut = rcut(boxvch)

         if (rbcut/w(1) .gt. 0.5d0 .or.
     &        rbcut/w(2) .gt. 0.5d0 .or.
     &        rbcut/w(3) .gt. 0.5d0) then
            write(iou,*) 'Problem with line 275 prvolume.f'
            write(iou,*) 'non-rectangular prvolume move rejected-',
     & ' box width below cutoff size'
            write(iou,*) 'w1:',w(1),'w2:',w(2),'w3:',w(3)
            hmat(boxvch,jhmat) = hmato(jhmat)
            call dump
            call cleanup('')
c            goto 500
         endif


c *** determine the displacement of the COM
         do i = 1,nchain
            imolty = moltyp(i)
            if (nboxi(i) .eq. boxvch) then
               if ( lx ) then
                  dx = sxcm(i)*(hmat(boxvch,1)-hmato(1))+
     +                 sycm(i)*(hmat(boxvch,4)-hmato(4))+
     +                 szcm(i)*(hmat(boxvch,7)-hmato(7))
                  xcm(i) = xcm(i) + dx
                  do j = 1, nunit(imolty)
                     rxu(i,j) = rxu(i,j) + dx
                  enddo
               elseif ( ly ) then
                  dy = sxcm(i)*(hmat(boxvch,2)-hmato(2))+
     +                 sycm(i)*(hmat(boxvch,5)-hmato(5))+
     +                 szcm(i)*(hmat(boxvch,8)-hmato(8))
                  ycm(i) = ycm(i) + dy
                  do j = 1, nunit(imolty)
                     ryu(i,j) = ryu(i,j) + dy
                  enddo
               else
                  dz = sxcm(i)*(hmat(boxvch,3)-hmato(3))+
     +                 sycm(i)*(hmat(boxvch,6)-hmato(6))+
     +                 szcm(i)*(hmat(boxvch,9)-hmato(9))
                  zcm(i) = zcm(i) + dz
                  do j = 1, nunit(imolty)
                     rzu(i,j) = rzu(i,j) + dz
                  enddo
               endif
            endif
         enddo

      else

c --- calculate new volume         
         voln = volo +
     +        rmvol(boxvch) * ( 2.0d0*random() - 1.0d0 )
         rbcut = rcut(boxvch)
         if ( lpbcz ) then
            if (lsolid(boxvch).and.lrect(boxvch)) then
c *** volume move independently in x, y, z directions
               dfac=voln/volo
            else
               dfac = (voln/volo)**(1.0d0/3.0d0)
            endif
c            vminim = rbcut**(3.0d0)
         else
            dfac= dsqrt(voln/volo)
c            vminim = rbcut**(2.0d0)
         endif

c         if ( voln .lt. vminim ) then
c            write(iou,*) 'prvolume move rejected - below cut-off size' 
c            return
c         endif

         if ( lsolid(boxvch).and.lrect(boxvch) ) then
            if (lx) boxlx(boxvch) = boxlx(boxvch) * dfac
            if (ly) boxly(boxvch) = boxly(boxvch) * dfac
            if (lz) boxlz(boxvch) = boxlz(boxvch) * dfac
         else
            boxlx(boxvch) = boxlx(boxvch) * dfac
            boxly(boxvch) = boxly(boxvch) * dfac
            if ( lpbcz ) then
               boxlz(boxvch) = boxlz(boxvch) * dfac
            endif
         endif


         rbcut = 2.0d0*rbcut

         if (boxlx(boxvch) .lt. rbcut .or.
     &        boxly(boxvch) .lt. rbcut .or.
     &        (lpbcz .and. boxlz(boxvch) .lt. rbcut) ) then
            boxlx(boxvch) = bxo
            boxly(boxvch) = byo
            if ( lpbcz ) then
               boxlz(boxvch) = bzo
            endif
            write(iou,*) 'boxvch',boxvch
            write(iou,*) 'Problem in line 381 of subroutine prvolume.f'
            write(iou,*) 'A move was attempted that would lead to a 
     & boxlength less than twice rcut'
            call dump
            call cleanup('')
            return
         endif

c *** determine new positions of the molecules
c *** calculate centre of mass and its displacement

c - WARNING
         if ( .not. lfold ) 
     &        call cleanup('volume move only correct with folded coordinates')

         df = dfac - 1.0d0

         do i = 1, nchain
            imolty = moltyp(i)

c ----- Check if the chain i is in the correct box
            if (nboxi(i) .eq. boxvch) then
            
               if (lsolid(boxvch).and.lrect(boxvch)) then
                  if ( lx ) then
                     dx = xcm(i) * df
                     xcm(i) = xcm(i) + dx
                     do j = 1, nunit(imolty)
                        rxu(i,j) = rxu(i,j) + dx
                     enddo
                  endif
                  if ( ly ) then
                     dy = ycm(i) * df
                     ycm(i) = ycm(i) + dy
                     do j = 1, nunit(imolty)
                        ryu(i,j) = ryu(i,j) + dy
                     enddo
                  endif
                  if ( lz ) then
                     dz = zcm(i) * df
                     zcm(i) = zcm(i) + dz
                     do j = 1, nunit(imolty)
                        rzu(i,j) = rzu(i,j) + dz
                     enddo
                  endif
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
                  enddo
               endif
            endif
         enddo
         
      endif

      lvol = .true.
      if ( lchgall ) then
           if (lsolid(boxvch).and.(.not.lrect(boxvch))) then
             min_boxl = min(min_width(boxvch,1),min_width(boxvch,2),
     &                    min_width(boxvch,3))
           else
              min_boxl = min(boxlx(boxvch),boxly(boxvch),boxlz(boxvch))
           endif
           calp(boxvch) = kalp(boxvch)/boxlx(boxvch)
      endif    
      call sumup( ovrlap, v, vinter, vtail, vdum,vdum,
     +     vdum,vdum,vext,velect,vdum, boxvch, lvol)
      if ( ovrlap ) then
c         write(iou,*) 'move rejected due to overlap in PRVOLUME'
         goto 500
      endif

      vintern  = vinter
      vtailn   = vtail
      vextn    = vext  
      velectn  = velect
      v3n = v3garo
      vboxn    = vboxo + (vintern-vintero) + (vextn-vexto) 
     &     + (velectn-velecto) + (v3n-v3o)
c      write (6,*) 'new  energy',  vboxn

      if ( lanes ) then
c *** for ANES algorithm, optimize the charge configuration         
c *** on the new coordinates, continue to use the fluctuating charge
c *** algorithm to optimize the charge configurations, update the
c *** energy, coordinates and the ewald sum

         vbox(boxvch)    = vbox(boxvch) + (vboxn - vboxo)
         vinterb(boxvch) = vinterb(boxvch) + (vintern-vintero)
         vtailb(boxvch)  = vtailb(boxvch) + (vtailn-vtailo)
         vextb(boxvch)   = vextb(boxvch) + (vextn-vexto)  
         velectb(boxvch) = velectb(boxvch) + (velectn-velecto)
         do ichoiq = 1,nchoiq(boxvch)
            call flucq(0,boxvch)
         enddo
         
         dele = (vbox(boxvch) - vboxo)+ express(boxvch)*(voln-volo) 
     +        - ((nchbox(boxvch)+ghost_particles(boxvch))
     +              * dlog(voln/volo) / beta )

      else
         
         dele = ( vboxn - vboxo ) +  express(boxvch)*(voln-volo) 
     +        - ((nchbox(boxvch)+ghost_particles(boxvch)) 
     +              * dlog(voln/volo)/beta )
      endif

c --- acceptance test

c--- check problem--
c      write (6,*) 'length of box',      boxlx(boxvch)
c      write (6,*) 'change of  vol',     voln - volo
c      write (6,*) 'change of  energy',  vboxn, vboxo
c      write (6,*) 'dele',     dele


      if (random() .lt. dexp(-(beta*dele)) ) then
c --      accepted
          if ( lsolid(boxvch) .and. .not. lrect(boxvch) ) then
             bshmat(boxvch,jhmat) = bshmat(boxvch,jhmat) + 1.0d0
          else
             bsvol(boxvch) = bsvol(boxvch) + 1.0d0
          endif
          if ( .not. lanes ) then
             vbox(boxvch)    = vbox(boxvch) + (vboxn - vboxo)
             vinterb(boxvch) = vinterb(boxvch) + (vintern-vintero)
             vtailb(boxvch)  = vtailb(boxvch) + (vtailn-vtailo)
             vextb(boxvch)   = vextb(boxvch) + (vextn-vexto)  
             velectb(boxvch) = velectb(boxvch) + (velectn-velecto)
             v3garob(boxvch) = v3garob(boxvch) + (v3n-v3o)
          endif
c ---     update centers of mass
          call ctrmas(.true.,boxvch,0,5)
c *** update linkcell, if applicable
          if (licell .and. (boxvch .eq. boxlink)) then
             call linkcell(1,1,lddum,lddum,lddum,lddum2)
          endif
          return
      endif

c ---     rejected
c --- re store old box lengths, energy, configuration etc
 500  continue
      if ( lsolid(boxvch) .and. .not. lrect( boxvch)) then
         do i = 1,9
            hmat(boxvch,i) = hmato(i)
         enddo
         call matops(boxvch)  
      else
         boxlx(boxvch)   = bxo
         boxly(boxvch)   = byo
         if ( lpbcz ) boxlz(boxvch)   = bzo
      endif
      vbox(boxvch) = vboxo
      vinterb(boxvch) = vintero
      vtailb(boxvch) = vtailo
      vextb(boxvch) = vexto
      velectb(boxvch) = velecto
      vflucqb(boxvch) = vflucqo
      v3garob(boxvch) = v3o

c --- restore old neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_cnt(i) = neigh_o(i)
            do j=1,neigh_cnt(i)
               neighbor(j,i) = neighboro(j,i)
               ndij(j,i) = ndijo(j,i)
               nxij(j,i) = nxijo(j,i)
               nyij(j,i) = nyijo(j,i)
               nzij(j,i) = nzijo(j,i)
            enddo
         enddo
      endif

      if ( lewald ) then
c --- restore old k vectors and reciprocal sum and calp
         calp(boxvch) = calpo
c         if (numvect(boxvch) .ne. numvecto) then
c            write(iou,*) '### numvect changing in prvol'
c         endif
         numvect(boxvch) = numvecto
         ncount = numvect(boxvch)
         call recip(boxvch,vdum,vdum,4)
         do ic = 1,ncount
            kx(ic,boxvch) = kxo(ic)
            ky(ic,boxvch) = kyo(ic)
            kz(ic,boxvch) = kzo(ic)
            prefact(ic,boxvch) = prefacto(ic)
         enddo
      endif
      do i = 1, nchain
         imolty = moltyp(i)
c ----- Check if the chain i is in the correct box
         if (nboxi(i) .eq. boxvch) then
            xcm(i) = xcmo(i)
            ycm(i) = ycmo(i)
            if ( lpbcz ) zcm(i) = zcmo(i)
            do j = 1, nunit(imolty)
               rxu(i,j) = rxuo(i,j)
               ryu(i,j) = ryuo(i,j)
               if ( lpbcz ) rzu(i,j) = rzuo(i,j)
               qqu(i,j) = qquo(i,j)
            enddo
            if (lneighbor) then
               neigh_cnt(i) = neigho_cnt(i)
               do j = 1,neigh_cnt(i)
                  neighbor(j,i)=neigho(j,i)
               enddo
            endif
         endif
      enddo

c      write(iou,*) 'end VOLUME of LNPTGIBBS'

      return
      end



