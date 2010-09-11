      subroutine volume

! volume
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
! John Stubbs, and Collin Wick and Ilja Siepmann  
!                     
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to 
!
! Free Software Foundation, Inc. 
! 59 Temple Place - Suite 330
! Boston, MA  02111-1307, USA.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!    *************************************************************
!    ** makes an isotropic volume change                        **
!    ** the maximum change is controlled by rmtrax and the      **
!    ** number of successful trial moves is stored in bsvol.    **
!    *************************************************************
!
! --- perform change of the volume: random walk in ln(vol)
!
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'ensemble.inc'
      include 'system.inc' 
      include 'inputdata.inc'
      include 'bnbsma.inc'
      include 'ewaldsum.inc'
      include 'cell.inc'
!kea
      include 'garofalini.inc'
      include 'neigh.inc'

!      logical::ovrlap,lvol,lx(nbxmax),ly(nbxmax),lz(nbxmax),lncubic
      logical::ovrlap,lvol,lx,ly,lz,lncubic
      logical::la,lb,lc

      integer::i,j,ibox,imolty,ic,ncount,boxa,boxb,ipair,ipairb
      real(8)::bxo(nbxmax),byo(nbxmax),bzo(nbxmax)
      real(8)::volo(nbxmax),vboxo(nbxmax),dfac(nbxmax),
     &                 voln(nbxmax),vflucqo(nbxmax),vintero(nbxmax),
     &                 vtailo(nbxmax),vexto(nbxmax),velecto(nbxmax)
      real(8)::kxo(vectormax,nbxmax),kyo(vectormax,nbxmax),
     &                 kzo(vectormax,nbxmax),prefacto(vectormax,nbxmax)
      real(8)::vboxn(nbxmax),vintern(nbxmax),
     &                 vtailn(nbxmax),vextn(nbxmax),velectn(nbxmax)
      real(8)::rxuo(nmax,numax),ryuo(nmax,numax)
     &                ,rzuo(nmax,numax),qquo(nmax,numax)
      real(8)::volt,expdv,random,df,dx,dy,dz,v,dele,
     &                 vinter,vtail,vext,boxlen,vdum,velect
     &                 ,velect_intra,velect_inter
      real(8)::vminim(nbxmax)
      real(8)::xcmo,ycmo,zcmo,calpo(nbxmax),numvecto(nbxmax)
      real(8)::rbcut(nbxmax),rbcuta,rbcutb,rpair,rm

      real(8)::w(3),vx,vy,vz,min_boxl

      dimension xcmo(nmax),ycmo(nmax),zcmo(nmax)

      integer::ichoiq
      integer::hbox,jbox,jhmat
      real(8)::rbox,hmato(9),hmatio(9)

      integer::idum
      real(8)::lddum,lddum2(27)
!kea
      real(8)::v3n(nbxmax),v3o(nbxmax)

! --------------------------------------------------------------------

!      write(iou,*) 'start VOLUME'

! *** select pair of boxes to do the volume move
      if ( nvolb .gt. 1 ) then
         rpair = random()
         do 96 ipair = 1, nvolb
            if ( rpair .lt. pmvolb(ipair) ) then
               ipairb = ipair
               goto 97
            end if
 96      continue
      else
         ipairb = 1
      end if
      
 97   boxa = box5(ipairb)
      boxb = box6(ipairb)

      bnvol(ipairb) = bnvol(ipairb) + 1.0d0

      lx = .false.
      ly = .false.
      lz = .false.

      if ( lsolid(boxa) ) then
! *** volume move independently in x, y, z directions
! * changing to only move in z direction:
         rm = random()
         if ( rm .le. pmvolx ) then
!            lx(boxa) = .true.
!            ly(boxa) = .false.
!            lz(boxa) = .false.
            lx = .true.
            ly = .false.
            lz = .false.
         elseif ( rm .le. pmvoly ) then
!            lx(boxa) = .false.
!            ly(boxa) = .true.
!            lz(boxa) = .false.
            lx = .false.
            ly = .true.
            lz = .false.
         else
!            lx(boxa) = .false.
!            ly(boxa) = .false.
!            lz(boxa) = .true.
            lx = .false.
            ly = .false.
            lz = .true.
         end if
         if (.not. lrect(boxa)) then
            do i = 1,9
               hmato(i) = hmat(boxa,i)
               hmatio(i) = hmati(boxa,i)
            end do
         end if
      end if
      if ( lsolid(boxb) ) then
! *** volume move independently in x, y, z directions
! * changing to only move in z direction:
         rm = random()
         if ( rm .le. pmvolx ) then
!            lx(boxb) = .true.
!            ly(boxb) = .false.
!            lz(boxb) = .false.
            lx = .true.
            ly = .false.
            lz = .false.
         elseif ( rm .le. pmvoly ) then
!            lx(boxb) = .false.
!            ly(boxb) = .true.
!            lz(boxb) = .false.
            lx = .false.
            ly = .true.
            lz = .false.
         else
!            lx(boxb) = .false.
!            ly(boxb) = .false.
!            lz(boxb) = .true.
            lx = .false.
            ly = .false.
            lz = .true.
         end if
         if ( .not. lrect(boxb) ) then
            do i = 1,9
               hmato(i) = hmat(boxb,i)
               hmatio(i) = hmati(boxb,i)
            end do
         end if
      end if

      if (lsolid(boxa) .and. .not. lrect(boxa)) then
         if (lsolid(boxb) .and. .not. lrect(boxb)) then
            write(iou,*) 'can not perform volume move between',
     &           ' two non-rectangular boxes'
            call cleanup('')
         end if
      end if

! --- store old box lengths, energy, configuration etc
      lncubic = .false.

      do ibox = 1, 2
         if (ibox .eq. 1) i = boxa
         if (ibox .eq. 2) i = boxb

         bxo(i) = boxlx(i)
         byo(i) = boxly(i)

         if ( lpbcz ) then
            bzo(i) = boxlz(i)
         end if

         if (lsolid(i) .and. .not. lrect(i)) then
             volo(i) = cell_vol(i)
             lncubic = .true.
         else
            if ( lpbcz ) then
               volo(i)   = bxo(i)*byo(i)*bzo(i)
            else
               volo(i)   = bxo(i)*byo(i)
            end if
         end if
         vboxo(i)    = vbox(i)
         vintero(i)  = vinterb(i)
         vtailo(i)   = vtailb(i)
         vexto(i)    = vextb(i)  
         velecto(i)  = velectb(i)
         vflucqo(i) = vflucqb(i)
         v3o(i)      = v3garob(i)

! --- store neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_o(i) = neigh_cnt(i)
            do j=1,neigh_o(i)
               neighboro(j,i) = neighbor(j,i)
               ndijo(j,i) = ndij(j,i)
               nxijo(j,i) = nxij(j,i)
               nyijo(j,i) = nyij(j,i)
               nzijo(j,i) = nzij(j,i)
            end do
         end do
      end if

! --- store old k vectors and reciprocal sum
         if ( lewald ) then
            calpo(i) = calp(i)
            numvecto(i) = numvect(i)
            ncount = numvect(i)
            call recip(i,vdum,vdum,3)
            do ic = 1,ncount
               kxo(ic,i) = kx(ic,i)
               kyo(ic,i) = ky(ic,i)
               kzo(ic,i) = kz(ic,i)
               prefacto(ic,i) = prefact(ic,i)
            end do
         end if
      end do

      do i = 1, nchain
         ibox = nboxi(i)
         if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
            imolty = moltyp(i)
               xcmo(i) = xcm(i)
               ycmo(i) = ycm(i)
               if (lpbcz) zcmo(i) = zcm(i)
               do j = 1, nunit(imolty)
                  rxuo(i,j) = rxu(i,j)
                  ryuo(i,j) = ryu(i,j)
                  if ( lpbcz ) rzuo(i,j) = rzu(i,j)
                  qquo(i,j) = qqu(i,j)
               end do
         end if

      end do


!      write(iou,*) 'before lncubic', lx,ly,lz

! --- calculate total volume
      volt = volo(boxa) + volo(boxb)

      if ( lncubic ) then

! *** select one of the cell edge
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
            end if
         elseif ( ly ) then
            la = .false.
            if ( rbox .le. 1.5d0 ) then
               lb = .true.
               lc = .false.
            else
               lb = .false.
               lc = .true.
            end if
         else
            la = .false.
            lb = .false.
            lc = .true.
         end if

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
         end if

         if ( lsolid(boxa) .and. .not. lsolid(boxb) ) then
            hbox = boxa
            jbox = boxb

            hmat(boxa,jhmat) = hmat(boxa,jhmat) + rmhmat(boxa,jhmat)*
     &           ( 2.0d0*random() - 1.0d0 )
            bnhmat(boxa,jhmat) = bnhmat(boxa,jhmat) + 1.0d0

         else
            hbox = boxb
            jbox = boxa

            hmat(boxb,jhmat) = hmat(boxb,jhmat) + rmhmat(boxb,jhmat)*
     &           ( 2.0d0*random() - 1.0d0 )
            bnhmat(boxb,jhmat) = bnhmat(boxb,jhmat) + 1.0d0

         end if

         call matops(hbox)

         voln(hbox) = cell_vol(hbox)

         rbcut(hbox) = rcut(hbox)

         w(1) = min_width(hbox,1) 
         w(2) = min_width(hbox,2)
         w(3) = min_width(hbox,3)


         if (rbcut(hbox)/w(1) .gt. 0.5d0 .or.
     &        rbcut(hbox)/w(2) .gt. 0.5d0 .or.
     &        rbcut(hbox)/w(3) .gt. 0.5d0) then
            write(iou,*) 'Problem with line 381 in volume.f'
            write(iou,*) 'non-rectangular volume move rejected-',
     &                 ' box width below cutoff size'
            write(iou,*) 'w1:',w(1),'w2:',w(2),'w3:',w(3)
            hmat(hbox,jhmat) = hmato(jhmat)
            call dump
            call cleanup('')
!            goto 500
         end if

         voln(jbox) = volt-voln(hbox)
!         if ( la ) boxlx(hbox) = hmat(hbox,1)
!         if ( lb ) boxly(hbox) = dsqrt(hmat(hbox,4)*hmat(hbox,4)+
!     &                 hmat(hbox,5)*hmat(hbox,5))
!         if ( lc ) boxlz(hbox) = dsqrt(hmat(hbox,7)*hmat(hbox,7)
!     &        +hmat(hbox,8)*hmat(hbox,8) +
!     &        hmat(hbox,9)*hmat(hbox,9))

! *** determine the displacement of the COM
        
         df=(voln(jbox)/volo(jbox))**(1.0d0/3.0d0)
         boxlx(jbox) = boxlx(jbox)*df
         boxly(jbox) = boxly(jbox)*df
         boxlz(jbox) = boxlz(jbox)*df
         df = df - 1.0d0

         do i = 1,nchain
            imolty = moltyp(i)
            if (nboxi(i) .eq. hbox) then
               if ( lx ) then
                  dx = sxcm(i)*(hmat(hbox,1)-hmato(1))+
     &                 sycm(i)*(hmat(hbox,4)-hmato(4))+
     &                 szcm(i)*(hmat(hbox,7)-hmato(7))
                  xcm(i) = xcm(i) + dx
                  do j = 1, nunit(imolty)
                     rxu(i,j) = rxu(i,j) + dx
                  end do
               elseif ( ly ) then
                  dy = sxcm(i)*(hmat(hbox,2)-hmato(2))+
     &                 sycm(i)*(hmat(hbox,5)-hmato(5))+
     &                 szcm(i)*(hmat(hbox,8)-hmato(8))
                  ycm(i) = ycm(i) + dy
                  do j = 1, nunit(imolty)
                     ryu(i,j) = ryu(i,j) + dy
                  end do
               else
                  dz = sxcm(i)*(hmat(hbox,3)-hmato(3))+
     &                 sycm(i)*(hmat(hbox,6)-hmato(6))+
     &                 szcm(i)*(hmat(hbox,9)-hmato(9))
                  zcm(i) = zcm(i) + dz
                  do j = 1, nunit(imolty)
                     rzu(i,j) = rzu(i,j) + dz
                  end do
               end if
            elseif (nboxi(i) .eq. jbox) then
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
               end do
            end if

         end do

      else
      
        
! --- calculate new volume
         expdv = dexp(dlog(volo(boxa)/volo(boxb))
     &        + rmvol(ipairb)*(2.0d0*random()-1.0d0))
         voln(boxa)= expdv*volt/(1+expdv)
         voln(boxb)= volt-voln(boxa)
         rbcut(boxa) = rcut(boxa)
         rbcut(boxb) = rcut(boxb)
         
         if ( lpbcz ) then

            if ( lsolid(boxa).and.lrect(boxa) ) then
! *** volume move independently in x, y, z directions
               dfac(boxa)=voln(boxa)/volo(boxa)
            else
               dfac(boxa)= (voln(boxa)/volo(boxa))**(1.0d0/3.0d0)
            end if
         

            if ( lsolid(boxb).and.lrect(boxb) ) then
! *** volume move independently in x, y, z directions
               dfac(boxb)=voln(boxb)/volo(boxb)
            else
               dfac(boxb)= (voln(boxb)/volo(boxb))**(1.0d0/3.0d0)
            end if
            vminim(boxa) = (2.0d0*rbcut(boxa))**(3.0d0)
         else
            dfac(boxa)= dsqrt(voln(boxa)/volo(boxa))
            dfac(boxb)= dsqrt(voln(boxb)/volo(boxb))
            vminim(boxb) = (2.0d0*rbcut(boxb))**(2.0d0)
         end if
         
         if ( lsolid(boxa).and. lrect(boxa) ) then
            if (lx) boxlx(boxa) = boxlx(boxa) * dfac(boxa)
            if (ly) boxly(boxa) = boxly(boxa) * dfac(boxa)
            if (lz) boxlz(boxa) = boxlz(boxa) * dfac(boxa)
         else
            
            boxlx(boxa) = boxlx(boxa) * dfac(boxa)
            boxly(boxa) = boxly(boxa) * dfac(boxa)
            if ( lpbcz ) then
               boxlz(boxa) = boxlz(boxa) * dfac(boxa)
            end if
         end if

         if ( lsolid(boxb).and.lrect(boxb) ) then
            if (lx) boxlx(boxb) = boxlx(boxb) * dfac(boxb)
            if (ly) boxly(boxb) = boxly(boxb) * dfac(boxb)
            if (lz) boxlz(boxb) = boxlz(boxb) * dfac(boxb)
         else
            boxlx(boxb) = boxlx(boxb) * dfac(boxb)
            boxly(boxb) = boxly(boxb) * dfac(boxb)
            if ( lpbcz ) then
               boxlz(boxb) = boxlz(boxb) * dfac(boxb)
            end if
         end if

         rbcuta = 2.0d0*rbcut(boxa)
         rbcutb = 2.0d0*rbcut(boxb)
         if ( boxlx(boxa) .lt. rbcuta .or. 
     &        boxly(boxa) .lt. rbcuta .or. 
     &        (lpbcz .and. boxlz(boxa) .lt. rbcuta) .or.
     &        boxlx(boxb) .lt. rbcutb .or. 
     &        boxly(boxb) .lt. rbcutb .or. 
     &        (lpbcz .and. boxlz(boxb) .lt. rbcutb) ) then
                        
            write(iou,*) 'Problem in line 552 of subroutine volume.f'
            write(iou,*) 'A move was attempted that would lead to a 
     & boxlength less than twice rcut'

            boxlx(boxa) = bxo(boxa)
            boxlx(boxb) = bxo(boxb)
            boxly(boxa) = byo(boxa)
            boxly(boxb) = byo(boxb)
            if ( lpbcz ) then 
               boxlz(boxa) = bzo(boxa)
               boxlz(boxb) = bzo(boxb)
            end if
            call dump 
            call cleanup('')
            return
         end if


! *** determine new positions of the molecules
! *** calculate centre of mass and its displacement

! - WARNING
         if ( .not. lfold ) 
     &        call cleanup('volume move only correct with folded coordinates')

         do i = 1, nchain

            ibox = nboxi(i)

            if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
               
               imolty = moltyp(i)
               df = dfac(ibox) - 1.0d0
               
!     if ( lsolid(ibox) ) then
               if ( lsolid(ibox).and.lrect(ibox) ) then
!     if ( lx(ibox) ) then
                  if ( lx ) then
                     dx = xcm(i) * df
                     xcm(i) = xcm(i) + dx
                     do j = 1, nunit(imolty)
                        rxu(i,j) = rxu(i,j) + dx
                     end do
                  end if
!     if ( ly(ibox) ) then
                  if ( ly ) then
                     dy = ycm(i) * df
                     ycm(i) = ycm(i) + dy
                     do j = 1, nunit(imolty)
                        ryu(i,j) = ryu(i,j) + dy
                     end do
                  end if
!     if ( lz(ibox) ) then
                  if ( lz ) then
                     dz = zcm(i) * df
                     zcm(i) = zcm(i) + dz
                     do j = 1, nunit(imolty)
                        rzu(i,j) = rzu(i,j) + dz
                     end do
                  end if
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
                  end do
               end if
            end if
         end do

      end if

      lvol = .true.
      if ( lchgall ) then
         if (lsolid(boxa).and.(.not.lrect(boxa))) then
             min_boxl = min(min_width(boxa,1),min_width(boxa,2),
     &                    min_width(boxa,3))
         else
              min_boxl = min(boxlx(boxa),boxly(boxa),boxlz(boxa))
         end if
         calp(boxa) = kalp(boxa)/min_boxl
         if (lsolid(boxb).and.(.not.lrect(boxb))) then
             min_boxl = min(min_width(boxb,1),min_width(boxb,2),
     &                    min_width(boxb,3))
         else
              min_boxl = min(boxlx(boxb),boxly(boxb),boxlz(boxb))
         end if
         calp(boxb) = kalp(boxb)/min_boxl
      end if

      do i = 1,2
         if ( i .eq. 1 ) ibox = boxa
         if ( i .eq. 2 ) ibox = boxb
         call sumup( ovrlap, v, vinter,vtail, vdum,vdum,
     &                  vdum,vdum,vext,velect,vdum, ibox, lvol)
         if ( ovrlap ) goto 500
         vintern(ibox) = vinter
         vtailn(ibox)  = vtail
         vextn(ibox)   = vext  
         velectn(ibox) = velect
         v3n(ibox) = v3garo
         vboxn(ibox)   = vboxo(ibox) + (vintern(ibox)-vintero(ibox))
     &       + (vextn(ibox)-vexto(ibox)) + (velectn(ibox)-velecto(ibox))
     &        + (v3n(ibox)-v3o(ibox))
!kea
      end do

      if ( lanes ) then
! *** for ANES algorithm, optimize the charge configuration
         do i = 1,2
            if ( i .eq. 1 ) ibox = boxa
            if ( i .eq. 2 ) ibox = boxb
! *** on the new coordinates, continue to use the fluctuating charge
! *** algorithm to optimize the charge configurations, update the
! *** energy, coordinates and the ewald sum

            vbox(ibox) = vbox(ibox) + (vboxn(ibox) - vboxo(ibox)) 
            vinterb(ibox)  = vinterb(ibox) + 
     &           (vintern(ibox) - vintero(ibox))
            vtailb(ibox) = vtailb(ibox) + (vtailn(ibox) - vtailo(ibox))
            vextb(ibox) = vextb(ibox) + (vextn(ibox) - vexto(ibox))
            velectb(ibox) = velectb(ibox) + 
     &           (velectn(ibox) - velecto(ibox))
            do ichoiq = 1,nchoiq(ibox)
               call flucq(0,ibox)
            end do
         end do
         
         dele = (vbox(boxa) - vboxo(boxa))+( vbox(boxb)- vboxo(boxb))
     &        - ((nchbox(boxa)+1+ghost_particles(boxa))
     &        *dlog(voln(boxa)/volo(boxa))/beta)
     &        - ((nchbox(boxb)+1+ghost_particles(boxb))
     &        *dlog(voln(boxb)/volo(boxb))/beta)

      elseif (lncubic) then

         dele = (vboxn(boxa)-vboxo(boxa)) + (vboxn(boxb)-vboxo(boxb))
     &        - ((nchbox(boxa)+ghost_particles(boxa))
     &        *dlog(voln(boxa)/volo(boxa))/beta)
     &        - ((nchbox(boxb)+ghost_particles(boxb))
     &        *dlog(voln(boxb)/volo(boxb))/beta)

      else
         
         dele = (vboxn(boxa)-vboxo(boxa)) + (vboxn(boxb)-vboxo(boxb))
     &        - ((nchbox(boxa)+1+ghost_particles(boxa))
     &        *dlog(voln(boxa)/volo(boxa))/beta)
     &        - ((nchbox(boxb)+1+ghost_particles(boxb))
     &        *dlog(voln(boxb)/volo(boxb))/beta)
      end if

! --- acceptance test

!      write(iou,*) 'dele',dele
      if (random() .lt. dexp(-beta*dele) ) then
! --      accepted
         if ( lncubic ) then
            bshmat(hbox,jhmat) = bshmat(hbox,jhmat) + 1.0d0
         else
            bsvol(ipairb) = bsvol(ipairb) + 1.0d0
         end if

          if ( .not. lanes ) then
             do ibox = 1,2
                if ( ibox .eq. 1 ) i = boxa
                if ( ibox .eq. 2 ) i = boxb
                vbox(i)    = vbox(i) + (vboxn(i) - vboxo(i)) 
                vinterb(i) = vinterb(i) + (vintern(i) - vintero(i))
                vtailb(i)  = vtailb(i) + (vtailn(i) - vtailo(i))
                vextb(i)   = vextb(i) + (vextn(i) - vexto(i))  
                velectb(i) = velectb(i) + (velectn(i) - velecto(i))
                v3garob(i)     = v3garob(i) + (v3n(i)-v3o(i))
             end do
          end if

! --      update centers of mass
          call ctrmas(.true.,boxa,0,5)
          call ctrmas(.true.,boxb,0,5)
! *** update linkcell, if applicable
          if (licell .and. (boxa .eq. boxlink .or. boxb .eq. boxlink)) 
     &         then
             call linkcell(1,1,lddum,lddum,lddum,lddum2)
          end if
          return
      end if
! ---     rejected
! --- restore old energy, box lengths
 500  do ibox = 1, 2
         if ( ibox .eq. 1 ) i = boxa
         if ( ibox .eq. 2 ) i = boxb
         vbox(i)    = vboxo(i)
         vinterb(i)  = vintero(i)
         vtailb(i)   = vtailo(i)
         vextb(i)    = vexto(i)  
         velectb(i)  = velecto(i)
         vflucqb(i) = vflucqo(i)
         v3garob(i)      = v3o(i)

         if (lsolid(i) .and. .not. lrect(i)) then
            do j = 1,9
               hmat(ibox,j) = hmato(j)
            end do
            call matops(i)
         end if
        
         boxlx(i)   = bxo(i)
         boxly(i)   = byo(i)
         if ( lpbcz ) boxlz(i)   = bzo(i)

         if ( lewald ) then
! --- restore old k vectors and reciprocal sum and calp
            calp(i) = calpo(i)
!           ncount = numvect(i)
            ncount = numvecto(i)
            numvect(i) = numvecto(i)
            call recip(i,vdum,vdum,4)
            do ic = 1,ncount
               kx(ic,i) = kxo(ic,i)
               ky(ic,i) = kyo(ic,i)
               kz(ic,i) = kzo(ic,i)
               prefact(ic,i) = prefacto(ic,i)
            end do
         end if
      end do

! --- restore old neighbor list for garofalini --- KEA
      if (lgaro) then
         do i=1,nchain
            neigh_cnt(i) = neigh_o(i)
            do j=1,neigh_cnt(i)
               neighbor(j,i) = neighboro(j,i)
               ndij(j,i) = ndijo(j,i)
               nxij(j,i) = nxijo(j,i)
               nyij(j,i) = nyijo(j,i)
               nzij(j,i) = nzijo(j,i)
            end do
         end do
      end if

      do i = 1, nchain
         ibox = nboxi(i) 
         if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
            imolty = moltyp(i)
            xcm(i) = xcmo(i)
            ycm(i) = ycmo(i)
            if ( lpbcz ) zcm(i) = zcmo(i)
            do j = 1, nunit(imolty)
               rxu(i,j) = rxuo(i,j)
               ryu(i,j) = ryuo(i,j)
               if ( lpbcz ) rzu(i,j) = rzuo(i,j)
               qqu(i,j) = qquo(i,j)
            end do
         end if
      end do

!      write(iou,*) 'end VOLUME'
!      call dump

      return
      end

