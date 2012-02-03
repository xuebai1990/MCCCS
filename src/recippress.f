      subroutine recippress(ibox,repress,pxx,pyy,pzz,pxy,pyx,pxz,pzx,
     &                      pyz,pzy)
c recippress
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

      implicit none
c    ********************************************************************
c    ** calculates the reciprocal space contribution to pressure using **
c    ** thermodynamic definition. See J. Chem. Phys. Vol. 109 P2791.   **
c    ** written in 1998 by Bin Chen.                                   **
c    ** modified to calculate surface tension, 11/24/03 JMS            **
c    ********************************************************************

      include 'control.inc'
      include 'ewaldsum.inc'
      include 'coord.inc'
      include 'poten.inc'
      integer ncount,ibox,i,ii,imolty
      double precision factor,repress,repressx,repressy,repressz
     &     ,recipintra,piix,piiy,piiz,xcmi,ycmi,zcmi,arg

      double precision pxx,pyy,pzz,intraxx,intrayy,intrazz,intraxy
     &     ,intraxz,intrazy,intrayz,intrayx,intrazx,pxy,pyx,pyz,pzy
     &     ,pxz,pzx

      repress  = 0.0d0
      repressx = 0.0d0
      repressy = 0.0d0
      repressz = 0.0d0
      recipintra = 0.0d0
      pxy = 0.0d0
      pxz = 0.0d0
      pyx = 0.0d0
      pyz = 0.0d0
      pzx = 0.0d0
      pzy = 0.0d0

      intraxx = 0.0d0
      intrayy = 0.0d0
      intrazz = 0.0d0
      intraxy = 0.0d0
      intrazy = 0.0d0
      intraxz = 0.0d0
      intrazx = 0.0d0
      intrayz = 0.0d0
      intrayx = 0.0d0

      do ncount = 1,numvect(ibox)
         factor = prefact(ncount,ibox)*(ssumr(ncount,ibox)*
     &        ssumr(ncount,ibox) + ssumi(ncount,ibox)*
     &        ssumi(ncount,ibox))
         repressx = repressx + factor*(1.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kx(ncount,ibox)*kx(ncount,ibox))
         repressy = repressy + factor*(1.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*ky(ncount,ibox)*ky(ncount,ibox))
         repressz = repressz + factor*(1.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kz(ncount,ibox)*kz(ncount,ibox))
         pxy = pxy + factor*(0.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kx(ncount,ibox)*ky(ncount,ibox))
         pxz = pxz + factor*(0.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kx(ncount,ibox)*kz(ncount,ibox))
         pyz = pyz + factor*(0.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*ky(ncount,ibox)*kz(ncount,ibox))
      enddo
      repress = repressx + repressy + repressz
c * keep x,y,z separate for surface tension calculation
      pxx = repressx
      pyy = repressy
      pzz = repressz
      pyx = pxy
      pzx = pxz
      pzy = pyz

c --- the intramolecular part should be substracted
      do 10 i = 1, nchain
c ### check if i is in relevant box ###
         if ( nboxi(i) .eq. ibox ) then

            imolty = moltyp(i)
            if ( .not. lelect(imolty) ) goto 10
            xcmi = xcm(i)
            ycmi = ycm(i)
            zcmi = zcm(i)  

c --- loop over all beads ii of chain i 
            do ii = 1, nunit(imolty)
               
c --- compute the vector of the bead to the COM (p)
               
               piix = rxu(i,ii) - xcmi
               piiy = ryu(i,ii) - ycmi
               piiz = rzu(i,ii) - zcmi
             
               do ncount = 1,numvect(ibox)

c --- compute the dot product of k and r
                  
                  arg = kx(ncount,ibox)*rxu(i,ii) +
     &                 ky(ncount,ibox)*ryu(i,ii) +
     &                 kz(ncount,ibox)*rzu(i,ii)
                  factor = prefact(ncount,ibox)*2.0d0*
     &                 (-ssumr(ncount,ibox)*dsin(arg)
     &                 +ssumi(ncount,ibox)*dcos(arg))*qqu(i,ii)
                  recipintra = recipintra + factor*
     &                 (kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy
     &                 +kz(ncount,ibox)*piiz)
c * keep x,y and z separate for surface tension calculation
                  intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                  intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                  intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                  intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                  intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                  intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                  intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                  intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                  intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)

               enddo
            enddo
         endif
 10   continue
      repress = (repress + recipintra)*qqfact
      
      pxx = (pxx + intraxx)*qqfact
      pyy = (pyy + intrayy)*qqfact
      pzz = (pzz + intrazz)*qqfact

      pxy = pxy + intraxy
      pyx = pyx + intrayx
      pxz = pxz + intraxz
      pzx = pzx + intrazx
      pyz = pyz + intrayz
      pzy = pzy + intrazy

c      write(2,*) 'internal part:',intraxx,intrayy,intrazz
      
      return
      end



