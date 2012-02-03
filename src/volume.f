      subroutine volume

c volume
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
c    ** makes an isotropic volume change                        **
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
      include 'cell.inc'

c      logical ovrlap,lvol,lx(nbxmax),ly(nbxmax),lz(nbxmax),lncubic
      logical ovrlap,lvol,lx,ly,lz,lncubic
      logical la,lb,lc

      integer i,j,ibox,imolty,ic,ncount,boxa,boxb,ipair,ipairb
      double precision bxo(nbxmax),byo(nbxmax),bzo(nbxmax)
      double precision volo(nbxmax),vboxo(nbxmax),dfac(nbxmax),
     +                 voln(nbxmax),vflucqo(nbxmax),vintero(nbxmax),
     +                 vtailo(nbxmax),vexto(nbxmax),velecto(nbxmax)
      double precision kxo(vectormax,nbxmax),kyo(vectormax,nbxmax),
     +                 kzo(vectormax,nbxmax),prefacto(vectormax,nbxmax)
      double precision vboxn(nbxmax),vintern(nbxmax),
     +                 vtailn(nbxmax),vextn(nbxmax),velectn(nbxmax)
      double precision rxuo(nmax,numax),ryuo(nmax,numax)
     &                ,rzuo(nmax,numax),qquo(nmax,numax)
      double precision volt,expdv,random,df,dx,dy,dz,v,dele,
     &                 vinter,vtail,vext,boxlen,vdum,velect
      double precision vminim(nbxmax)
      double precision xcmo,ycmo,zcmo,calpo(nbxmax),numvecto(nbxmax)
      double precision rbcut(nbxmax),rbcuta,rbcutb,rpair,rm

      double precision w(3),vx,vy,vz

      dimension xcmo(nmax),ycmo(nmax),zcmo(nmax)

      integer ichoiq
      integer hbox,jbox,jhmat
      double precision rbox,hmato(9),hmatio(9)

      integer idum
      double precision lddum,lddum2(27)

C --------------------------------------------------------------------

c      write(6,*) 'start VOLUME'

c *** select pair of boxes to do the volume move
      if ( nvolb .gt. 1 ) then
         rpair = random()
         do 96 ipair = 1, nvolb
            if ( rpair .lt. pmvolb(ipair) ) then
               ipairb = ipair
               goto 97
            endif
 96      continue
      else
         ipairb = 1
      endif
      
 97   boxa = box5(ipairb)
      boxb = box6(ipairb)

      bnvol(ipairb) = bnvol(ipairb) + 1.0d0

      lx = .false.
      ly = .false.
      lz = .false.

      if ( lsolid(boxa) ) then
c *** volume move independently in x, y, z directions
c * changing to only move in z direction:
         rm = random()
         if ( rm .le. pmvolx ) then
c            lx(boxa) = .true.
c            ly(boxa) = .false.
c            lz(boxa) = .false.
            lx = .true.
            ly = .false.
            lz = .false.
         elseif ( rm .le. pmvoly ) then
c            lx(boxa) = .false.
c            ly(boxa) = .true.
c            lz(boxa) = .false.
            lx = .false.
            ly = .true.
            lz = .false.
         else
c            lx(boxa) = .false.
c            ly(boxa) = .false.
c            lz(boxa) = .true.
            lx = .false.
            ly = .false.
            lz = .true.
         endif
         if (.not. lrect(boxa)) then
            do i = 1,9
               hmato(i) = hmat(boxa,i)
               hmatio(i) = hmati(boxa,i)
            enddo
         endif
      endif
      if ( lsolid(boxb) ) then
c *** volume move independently in x, y, z directions
c * changing to only move in z direction:
         rm = random()
         if ( rm .le. pmvolx ) then
c            lx(boxb) = .true.
c            ly(boxb) = .false.
c            lz(boxb) = .false.
            lx = .true.
            ly = .false.
            lz = .false.
         elseif ( rm .le. pmvoly ) then
c            lx(boxb) = .false.
c            ly(boxb) = .true.
c            lz(boxb) = .false.
            lx = .false.
            ly = .true.
            lz = .false.
         else
c            lx(boxb) = .false.
c            ly(boxb) = .false.
c            lz(boxb) = .true.
            lx = .false.
            ly = .false.
            lz = .true.
         endif
         if ( .not. lrect(boxb) ) then
            do i = 1,9
               hmato(i) = hmat(boxb,i)
               hmatio(i) = hmati(boxb,i)
            enddo
         endif
      endif

      if (lsolid(boxa) .and. .not. lrect(boxa)) then
         if (lsolid(boxb) .and. .not. lrect(boxb)) then
            write(6,*) 'can not perform volume move between',
     &           ' two non-rectangular boxes'
            stop
         endif
      endif


c --- store old box lengths, energy, configuration etc
      lncubic = .false.

      do ibox = 1, 2
         if (ibox .eq. 1) i = boxa
         if (ibox .eq. 2) i = boxb

         bxo(i) = boxlx(i)
         byo(i) = boxly(i)
         if ( lpbcz ) then
            bzo(i) = boxlz(i)
         endif

         if (lsolid(i) .and. .not. lrect(i)) then
            volo(i) = (hmat(i,1) * (hmat(i,5) * hmat(i,9) -
     +           hmat(i,8) * hmat(i,6)) + hmat(i,4) * (hmat(i,8)
     +           * hmat(i,3) - hmat(i,2) * hmat(i,9)) + hmat(i,7)
     +           * (hmat(i,2) * hmat(i,6) - hmat(i,5)*hmat(i,3)))
            lncubic = .true.
         else
            if ( lpbcz ) then
               volo(i)   = bxo(i)*byo(i)*bzo(i)
            else
               volo(i)   = bxo(i)*byo(i)
            endif
         endif
         vboxo(i)    = vbox(i)
         vintero(i)  = vinterb(i)
         vtailo(i)   = vtailb(i)
         vexto(i)    = vextb(i)  
         velecto(i)  = velectb(i)
         vflucqo(i) = vflucqb(i)

c --- store old k vectors and reciprocal sum
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
            enddo
         endif
      enddo

      do i = 1, nchain
         ibox = nboxi(i)
         if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
            imolty = moltyp(i)
c$$$            if ( lsolid(ibox) ) then
c$$$               if ( lx(ibox) ) then
c$$$                  xcmo(i) = xcm(i)
c$$$                  do j = 1, nunit(imolty)
c$$$                     rxuo(i,j) = rxu(i,j)
c$$$                     qquo(i,j) = qqu(i,j)
c$$$                  enddo
c$$$               endif
c$$$               if ( ly(ibox) ) then
c$$$                  ycmo(i) = ycm(i)
c$$$                  do j = 1, nunit(imolty)
c$$$                     ryuo(i,j) = ryu(i,j)
c$$$                     qquo(i,j) = qqu(i,j)
c$$$                  enddo
c$$$               endif
c$$$               if ( lz(ibox) ) then
c$$$                  zcmo(i) = zcm(i)
c$$$                  do j = 1, nunit(imolty)
c$$$                     rzuo(i,j) = rzu(i,j)
c$$$                     qquo(i,j) = qqu(i,j)
c$$$                  enddo
c$$$               endif
c$$$            else
               xcmo(i) = xcm(i)
               ycmo(i) = ycm(i)
               if (lpbcz) zcmo(i) = zcm(i)
               do j = 1, nunit(imolty)
                  rxuo(i,j) = rxu(i,j)
                  ryuo(i,j) = ryu(i,j)
                  if ( lpbcz ) rzuo(i,j) = rzu(i,j)
                  qquo(i,j) = qqu(i,j)
               enddo
c$$$            endif
         endif

      enddo


c      write(6,*) 'before lncubic', lx,ly,lz

c --- calculate total volume
      volt = volo(boxa) + volo(boxb)

      if ( lncubic ) then

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

         if ( lsolid(boxa) .and. .not. lsolid(boxb) ) then
            hbox = boxa
            jbox = boxb

            hmat(boxa,jhmat) = hmat(boxa,jhmat) + rmhmat(boxa,jhmat)*
     +           ( 2.0d0*random() - 1.0d0 )
            bnhmat(boxa,jhmat) = bnhmat(boxa,jhmat) + 1.0d0

         else
            hbox = boxb
            jbox = boxa

            hmat(boxb,jhmat) = hmat(boxb,jhmat) + rmhmat(boxb,jhmat)*
     +           ( 2.0d0*random() - 1.0d0 )
            bnhmat(boxb,jhmat) = bnhmat(boxb,jhmat) + 1.0d0

         endif
         voln(hbox) = (hmat(hbox,1) * (hmat(hbox,5) * hmat(hbox,9) -
     +        hmat(hbox,8) * hmat(hbox,6)) + hmat(hbox,4)*(hmat(hbox,8)
     +        * hmat(hbox,3)-hmat(hbox,2) * hmat(hbox,9))+hmat(hbox,7)
     +        * (hmat(hbox,2)*hmat(hbox,6)-hmat(hbox,5)*hmat(hbox,3)))

c * calculate new boxwidths, compare to cutoff
c * v = u2 x u3                     
         vx = hmat(hbox,5)*hmat(hbox,9) - 
     &        hmat(hbox,6)*hmat(hbox,8)
         vy = hmat(hbox,6)*hmat(hbox,7) - 
     &        hmat(hbox,4)*hmat(hbox,9)
         vz = hmat(hbox,4)*hmat(hbox,8) - 
     &        hmat(hbox,5)*hmat(hbox,7)
         
         w(1) = voln(hbox) / dsqrt(vx**2 + vy**2 + vz**2)
         
c * v = u3 x u1                     
         vx = hmat(hbox,8)*hmat(hbox,3) - 
     &        hmat(hbox,9)*hmat(hbox,2)
         vy = hmat(hbox,9)*hmat(hbox,1) - 
     &        hmat(hbox,7)*hmat(hbox,3)
         vz = hmat(hbox,7)*hmat(hbox,2) - 
     &        hmat(hbox,8)*hmat(hbox,1)
         
         w(2) = voln(hbox) / dsqrt(vx**2 + vy**2 + vz**2)

c * v = u1 x u2
         vx = hmat(hbox,2)*hmat(hbox,6) - 
     &        hmat(hbox,3)*hmat(hbox,5)
         vy = hmat(hbox,3)*hmat(hbox,4) - 
     &        hmat(hbox,1)*hmat(hbox,6)
         vz = hmat(hbox,1)*hmat(hbox,5) - 
     &        hmat(hbox,2)*hmat(hbox,4)
         
         w(3) = voln(hbox) / dsqrt(vx**2 + vy**2 + vz**2)

         if ( rcut .gt. rcutchg(hbox) .or. lchgall ) then
            rbcut(hbox) = rcut
         else
            rbcut(hbox) = rcutchg(hbox)
         endif

         if (rbcut(hbox)/w(1) .gt. 0.5d0 .or.
     &        rbcut(hbox)/w(2) .gt. 0.5d0 .or.
     &        rbcut(hbox)/w(3) .gt. 0.5d0) then
            write(6,*) 'Problem with line 381 in volume.f'
            write(6,*) 'non-rectangular volume move rejected-',
     &' box width below cutoff size'
            write(6,*) 'w1:',w(1),'w2:',w(2),'w3:',w(3)
            hmat(hbox,jhmat) = hmato(jhmat)
            call dump
            stop
c            goto 500
         endif

         voln(jbox) = volt-voln(hbox)
         if ( la ) boxlx(hbox) = hmat(hbox,1)
         if ( lb ) boxly(hbox) = dsqrt(hmat(hbox,4)*hmat(hbox,4)+
     &                 hmat(hbox,5)*hmat(hbox,5))
         if ( lc ) boxlz(hbox) = dsqrt(hmat(hbox,7)*hmat(hbox,7)
     &        +hmat(hbox,8)*hmat(hbox,8) +
     &        hmat(hbox,9)*hmat(hbox,9))

         hmati(hbox,1) = (hmat(hbox,5)*hmat(hbox,9)-hmat(hbox,8)
     &        *hmat(hbox,6))/voln(hbox)
         hmati(hbox,5) = (hmat(hbox,1)*hmat(hbox,9)-hmat(hbox,7)
     &        *hmat(hbox,3))/voln(hbox)
         hmati(hbox,9) = (hmat(hbox,1)*hmat(hbox,5)-hmat(hbox,4)
     &        *hmat(hbox,2))/voln(hbox)
         hmati(hbox,4) = (hmat(hbox,7)*hmat(hbox,6)-hmat(hbox,4)
     &        *hmat(hbox,9))/voln(hbox)
         hmati(hbox,2) = (hmat(hbox,3)*hmat(hbox,8)-hmat(hbox,2)
     &        *hmat(hbox,9))/voln(hbox)
         hmati(hbox,7) = (hmat(hbox,4)*hmat(hbox,8)-hmat(hbox,7)
     &        *hmat(hbox,5))/voln(hbox)
         hmati(hbox,3) = (hmat(hbox,2)*hmat(hbox,6)-hmat(hbox,3)
     &        *hmat(hbox,5))/voln(hbox)
         hmati(hbox,8) = (hmat(hbox,7)*hmat(hbox,2)-hmat(hbox,8)
     &        *hmat(hbox,1))/voln(hbox)
         hmati(hbox,6) = (hmat(hbox,3)*hmat(hbox,4)-hmat(hbox,6)
     &        *hmat(hbox,1))/voln(hbox)
c *** determine the displacement of the COM
        
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
     +                 sycm(i)*(hmat(hbox,4)-hmato(4))+
     +                 szcm(i)*(hmat(hbox,7)-hmato(7))
                  xcm(i) = xcm(i) + dx
                  do j = 1, nunit(imolty)
                     rxu(i,j) = rxu(i,j) + dx
                  enddo
               elseif ( ly ) then
                  dy = sxcm(i)*(hmat(hbox,2)-hmato(2))+
     +                 sycm(i)*(hmat(hbox,5)-hmato(5))+
     +                 szcm(i)*(hmat(hbox,8)-hmato(8))
                  ycm(i) = ycm(i) + dy
                  do j = 1, nunit(imolty)
                     ryu(i,j) = ryu(i,j) + dy
                  enddo
               else
                  dz = sxcm(i)*(hmat(hbox,3)-hmato(3))+
     +                 sycm(i)*(hmat(hbox,6)-hmato(6))+
     +                 szcm(i)*(hmat(hbox,9)-hmato(9))
                  zcm(i) = zcm(i) + dz
                  do j = 1, nunit(imolty)
                     rzu(i,j) = rzu(i,j) + dz
                  enddo
               endif
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
               enddo
            endif

         enddo

      else
      
        
c --- calculate new volume
         expdv = dexp(dlog(volo(boxa)/volo(boxb))
     +        + rmvol(ipairb)*(2.0d0*random()-1.0d0))
         voln(boxa)= expdv*volt/(1+expdv)
         voln(boxb)= volt-voln(boxa)
         if ( rcut .gt. rcutchg(boxa) .or. lchgall ) then
            rbcut(boxa) = rcut
         else
            rbcut(boxa) = rcutchg(boxa)
         endif
         if ( rcut .gt. rcutchg(boxb) .or. lchgall ) then
            rbcut(boxb) = rcut
         else
            rbcut(boxb) = rcutchg(boxb)
         endif
         
         if ( lpbcz ) then

            if ( lsolid(boxa).and.lrect(boxa) ) then
c *** volume move independently in x, y, z directions
               dfac(boxa)=voln(boxa)/volo(boxa)
            else
               dfac(boxa)= (voln(boxa)/volo(boxa))**(1.0d0/3.0d0)
            endif
         

            if ( lsolid(boxb).and.lrect(boxb) ) then
c *** volume move independently in x, y, z directions
               dfac(boxb)=voln(boxb)/volo(boxb)
            else
               dfac(boxb)= (voln(boxb)/volo(boxb))**(1.0d0/3.0d0)
            endif
            vminim(boxa) = (2.0d0*rbcut(boxa))**(3.0d0)
         else
            dfac(boxa)= dsqrt(voln(boxa)/volo(boxa))
            dfac(boxb)= dsqrt(voln(boxb)/volo(boxb))
            vminim(boxb) = (2.0d0*rbcut(boxb))**(2.0d0)
         endif
         
         if ( lsolid(boxa).and. lrect(boxa) ) then
            if (lx) boxlx(boxa) = boxlx(boxa) * dfac(boxa)
            if (ly) boxly(boxa) = boxly(boxa) * dfac(boxa)
            if (lz) boxlz(boxa) = boxlz(boxa) * dfac(boxa)
         else
            
            boxlx(boxa) = boxlx(boxa) * dfac(boxa)
            boxly(boxa) = boxly(boxa) * dfac(boxa)
            if ( lpbcz ) then
               boxlz(boxa) = boxlz(boxa) * dfac(boxa)
            endif
         endif

         if ( lsolid(boxb).and.lrect(boxb) ) then
            if (lx) boxlx(boxb) = boxlx(boxb) * dfac(boxb)
            if (ly) boxly(boxb) = boxly(boxb) * dfac(boxb)
            if (lz) boxlz(boxb) = boxlz(boxb) * dfac(boxb)
         else
            boxlx(boxb) = boxlx(boxb) * dfac(boxb)
            boxly(boxb) = boxly(boxb) * dfac(boxb)
            if ( lpbcz ) then
               boxlz(boxb) = boxlz(boxb) * dfac(boxb)
            endif
         endif

         rbcuta = 2.0d0*rbcut(boxa)
         rbcutb = 2.0d0*rbcut(boxb)
         if ( boxlx(boxa) .lt. rbcuta .or. 
     &        boxly(boxa) .lt. rbcuta .or. 
     &        (lpbcz .and. boxlz(boxa) .lt. rbcuta) .or.
     &        boxlx(boxb) .lt. rbcutb .or. 
     &        boxly(boxb) .lt. rbcutb .or. 
     &        (lpbcz .and. boxlz(boxb) .lt. rbcutb) ) then
                        
            write(6,*) 'Problem in line 552 of subroutine volume.f'
            write(6,*) 'A move was attempted that would lead to a 
     & boxlength less than twice rcut'

            boxlx(boxa) = bxo(boxa)
            boxlx(boxb) = bxo(boxb)
            boxly(boxa) = byo(boxa)
            boxly(boxb) = byo(boxb)
            if ( lpbcz ) then 
               boxlz(boxa) = bzo(boxa)
               boxlz(boxb) = bzo(boxb)
            endif
            call dump 
            stop
            return
         endif


c *** determine new positions of the molecules
c *** calculate centre of mass and its displacement

c - WARNING
         if ( .not. lfold ) 
     &        stop 'volume move only correct with folded coordinates'

         do i = 1, nchain

            ibox = nboxi(i)

            if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
               
               imolty = moltyp(i)
               df = dfac(ibox) - 1.0d0
               
c     if ( lsolid(ibox) ) then
               if ( lsolid(ibox).and.lrect(ibox) ) then
c     if ( lx(ibox) ) then
                  if ( lx ) then
                     dx = xcm(i) * df
                     xcm(i) = xcm(i) + dx
                     do j = 1, nunit(imolty)
                        rxu(i,j) = rxu(i,j) + dx
                     enddo
                  endif
c     if ( ly(ibox) ) then
                  if ( ly ) then
                     dy = ycm(i) * df
                     ycm(i) = ycm(i) + dy
                     do j = 1, nunit(imolty)
                        ryu(i,j) = ryu(i,j) + dy
                     enddo
                  endif
c     if ( lz(ibox) ) then
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
         calp(boxa) = kalp(boxa)/boxlx(boxa)
         calp(boxb) = kalp(boxb)/boxlx(boxb)
      endif
      do i = 1,2
         if ( i .eq. 1 ) ibox = boxa
         if ( i .eq. 2 ) ibox = boxb
         call sumup( ovrlap, v, vinter,vtail, vdum,vdum,
     +                  vdum,vdum,vext,velect,vdum, ibox, lvol )
         if ( ovrlap ) goto 500
         vintern(ibox) = vinter
         vtailn(ibox)  = vtail
         vextn(ibox)   = vext  
         velectn(ibox) = velect
         vboxn(ibox)   = vboxo(ibox) + (vintern(ibox)-vintero(ibox))
     &   + (vextn(ibox)-vexto(ibox)) + (velectn(ibox)-velecto(ibox))
      enddo

      if ( lanes ) then
c *** for ANES algorithm, optimize the charge configuration
         do i = 1,2
            if ( i .eq. 1 ) ibox = boxa
            if ( i .eq. 2 ) ibox = boxb
c *** on the new coordinates, continue to use the fluctuating charge
c *** algorithm to optimize the charge configurations, update the
c *** energy, coordinates and the ewald sum

            vbox(ibox) = vbox(ibox) + (vboxn(ibox) - vboxo(ibox)) 
            vinterb(ibox)  = vinterb(ibox) + 
     &           (vintern(ibox) - vintero(ibox))
            vtailb(ibox) = vtailb(ibox) + (vtailn(ibox) - vtailo(ibox))
            vextb(ibox) = vextb(ibox) + (vextn(ibox) - vexto(ibox))
            velectb(ibox) = velectb(ibox) + 
     &           (velectn(ibox) - velecto(ibox))
            do ichoiq = 1,nchoiq(ibox)
               call flucq(0,ibox)
            enddo
         enddo
         
         dele = (vbox(boxa) - vboxo(boxa))+( vbox(boxb)- vboxo(boxb))
     +        - ((nchbox(boxa)+1)*dlog(voln(boxa)/volo(boxa))/beta)
     +        - ((nchbox(boxb)+1)*dlog(voln(boxb)/volo(boxb))/beta)

      elseif (lncubic) then

         dele = (vboxn(boxa)-vboxo(boxa)) + (vboxn(boxb)-vboxo(boxb))
     +        - (nchbox(boxa)*dlog(voln(boxa)/volo(boxa))/beta)
     +        - (nchbox(boxb)*dlog(voln(boxb)/volo(boxb))/beta)

      else
         
         dele = (vboxn(boxa)-vboxo(boxa)) + (vboxn(boxb)-vboxo(boxb))
     +        - ((nchbox(boxa)+1)*dlog(voln(boxa)/volo(boxa))/beta)
     +        - ((nchbox(boxb)+1)*dlog(voln(boxb)/volo(boxb))/beta)
      endif

c --- acceptance test

c      write(6,*) 'dele',dele
      if (random() .lt. dexp(-beta*dele) ) then
c --      accepted
         if ( lncubic ) then
            bshmat(hbox,jhmat) = bshmat(hbox,jhmat) + 1.0d0
         else
            bsvol(ipairb) = bsvol(ipairb) + 1.0d0
         endif

          if ( .not. lanes ) then
             do ibox = 1,2
                if ( ibox .eq. 1 ) i = boxa
                if ( ibox .eq. 2 ) i = boxb
                vbox(i)    = vbox(i) + (vboxn(i) - vboxo(i)) 
                vinterb(i) = vinterb(i) + (vintern(i) - vintero(i))
                vtailb(i)  = vtailb(i) + (vtailn(i) - vtailo(i))
                vextb(i)   = vextb(i) + (vextn(i) - vexto(i))  
                velectb(i) = velectb(i) + (velectn(i) - velecto(i))
             enddo
          endif

c --      update centers of mass
          call ctrmas(.true.,boxa,0,5)
          call ctrmas(.true.,boxb,0,5)
c *** update linkcell, if applicable
          if (licell .and. (boxa .eq. boxlink .or. boxb .eq. boxlink)) 
     &         then
             call linkcell(1,1,lddum,lddum,lddum,lddum2)
          endif
          return
      endif
c ---     rejected
c --- restore old energy, box lengths
 500  do ibox = 1, 2
         if ( ibox .eq. 1 ) i = boxa
         if ( ibox .eq. 2 ) i = boxb
         vbox(i)    = vboxo(i)
         vinterb(i)  = vintero(i)
         vtailb(i)   = vtailo(i)
         vextb(i)    = vexto(i)  
         velectb(i)  = velecto(i)
         vflucqb(i) = vflucqo(i)

         if (lsolid(i) .and. .not. lrect(i)) then
            do j = 1,9
               hmat(ibox,j) = hmato(j)
               hmati(ibox,j) = hmatio(j)
            enddo
         endif

         boxlx(i)   = bxo(i)
         boxly(i)   = byo(i)
         if ( lpbcz ) boxlz(i)   = bzo(i)

         if ( lewald ) then
c --- restore old k vectors and reciprocal sum and calp
            calp(i) = calpo(i)
c            ncount = numvect(i)
            ncount = numvecto(i)
            numvect(i) = numvecto(i)
            call recip(i,vdum,vdum,4)
            do ic = 1,ncount
               kx(ic,i) = kxo(ic,i)
               ky(ic,i) = kyo(ic,i)
               kz(ic,i) = kzo(ic,i)
               prefact(ic,i) = prefacto(ic,i)
            enddo
         endif
      enddo

      do i = 1, nchain
         ibox = nboxi(i) 
         if ( ibox .eq. boxa .or. ibox .eq. boxb ) then
            imolty = moltyp(i)
c$$$            if ( lsolid(ibox) ) then
c$$$               if ( lx(ibox) ) then
c$$$                  xcm(i) = xcmo(i)
c$$$                  do j = 1, nunit(imolty)
c$$$                     rxu(i,j) = rxuo(i,j)
c$$$                     qqu(i,j) = qquo(i,j)
c$$$                  enddo
c$$$               endif
c$$$               if ( ly(ibox) ) then
c$$$                  ycm(i) = ycmo(i)
c$$$                  do j = 1, nunit(imolty)
c$$$                     ryu(i,j) = ryuo(i,j)
c$$$                     qqu(i,j) = qquo(i,j)
c$$$                  enddo
c$$$               endif
c$$$               if ( lz(ibox) ) then
c$$$                  zcm(i) = zcmo(i)
c$$$                  do j = 1, nunit(imolty)
c$$$                     rzu(i,j) = rzuo(i,j)
c$$$                     qqu(i,j) = qquo(i,j)
c$$$                  enddo
c$$$               endif
c$$$            else
            xcm(i) = xcmo(i)
            ycm(i) = ycmo(i)
            if ( lpbcz ) zcm(i) = zcmo(i)
            do j = 1, nunit(imolty)
               rxu(i,j) = rxuo(i,j)
               ryu(i,j) = ryuo(i,j)
               if ( lpbcz ) rzu(i,j) = rzuo(i,j)
               qqu(i,j) = qquo(i,j)
            enddo
c$$$            endif
         endif
      enddo

c      write(6,*) 'end VOLUME'
      call dump

      return
      end

