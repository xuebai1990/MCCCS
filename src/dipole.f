      subroutine dipole(ibox,mtype)

c dipole
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
      integer ibox,mtype,i,imolty,zz,ii
      double precision dipox(2),dipoy(2),dipoz(2)
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ewaldsum.inc'
      
      if ( mtype .eq. 0 ) then

c in sumup, initiliaze dipole to be 0 and then sum up all the dipoles

         dipolex(ibox) = 0.0d0
         dipoley(ibox) = 0.0d0
         dipolez(ibox) = 0.0d0
         
         do i = 1, nchain
            if ( nboxi(i) .eq. ibox ) then
               imolty = moltyp(i)
               do ii = 1,nunit(imolty)
                  dipolex(ibox) = dipolex(ibox) 
     &                 + qqu(i,ii)*rxu(i,ii)
                  dipoley(ibox) = dipoley(ibox) 
     &                 + qqu(i,ii)*ryu(i,ii)
                  dipolez(ibox) = dipolez(ibox) 
     &                 + qqu(i,ii)*rzu(i,ii)
               enddo
            endif
         enddo

      elseif(mtype .eq. 1 ) then
         
c *** calculate the dipole moment after the traslation, rotation and 
c *** charge move

         do zz = 1,2
            dipox(zz) = 0.0d0
            dipoy(zz) = 0.0d0
            dipoz(zz) = 0.0d0
            imolty = moltion(zz)
            do i = 1,nunit(imolty)
               dipox(zz) = dipox(zz) + 
     &              qquion(i,zz)*rxuion(i,zz)
               dipoy(zz) = dipoy(zz) +
     &              qquion(i,zz)*ryuion(i,zz)
               dipoz(zz) = dipoz(zz) +
     &              qquion(i,zz)*rzuion(i,zz)
            enddo
         enddo
         dipolex(ibox) = dipolex(ibox) - dipox(1) + dipox(2) 
         dipoley(ibox) = dipoley(ibox) - dipoy(1) + dipoy(2)
         dipolez(ibox) = dipolez(ibox) - dipoz(1) + dipoz(2)

      elseif(mtype .eq. 2) then

c *** store the old dipole moment

         dipolexo = dipolex(ibox) 
         dipoleyo = dipoley(ibox)
         dipolezo = dipolez(ibox)
         
      elseif(mtype .eq. 3) then

c *** restore the old dipole moment

         dipolex(ibox) = dipolexo 
         dipoley(ibox) = dipoleyo
         dipolez(ibox) = dipolezo

      endif
      
       
      return
      end
