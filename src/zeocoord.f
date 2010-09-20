      subroutine zeocoord

! zeocoord
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

      implicit none
      include 'zeolite.inc'
      include 'zeopoten.inc'
      include 'mpi.inc'
      integer::count,frac,izeo,bonding(8),atomtype
      real(8)::wzeo,charge,alpha,beta,gamma
      character::atom*4

      open (unit = 47, file = 'zeolite.cssr')

      znum=0

      if (myid.eq.0) write(16,100)
      read(47,*)    zeorx,zeory,zeorz
      if (myid.eq.0) write(16,101) zeorx,zeory,zeorz,zeorx*zeory*zeorz
      read(47,*)    alpha,beta,gamma
      if (myid.eq.0) write(16,102) alpha,beta,gamma
      if (alpha.ne.90.or.beta.ne.90.or.gamma.ne.90) 
     &     call cleanup('** zeocoord: not cubic **')

      read(47,*)   nzeo,frac,nx,ny,nz
      

      if (myid.eq.0) write(16,103) nzeo
      if (nzeo.gt.nzeomax)
     &     call cleanup('** zeocoord: nzeo gt nzeomax **')

!     === Calculate zeolite density from assumption no of Si = 0.5* no O

      wzeo = (nzeo*16.00 + 0.5*nzeo*28.09)/(6.023e23)      
      if (myid.eq.0) write(16,104) wzeo,1000.0/(wzeo*6.023e23)

!     === Converting to absolute coordinates within [0,ri>

      do izeo = 1,nzeo
         read(47,99) count,atom,zeox(izeo),zeoy(izeo),
     &        zeoz(izeo)
         if (frac.eq.1) then
            zeox(izeo) = mod(zeox(izeo)+1.0d0,1.0d0)*zeorx
            zeoy(izeo) = mod(zeoy(izeo)+1.0d0,1.0d0)*zeory
            zeoz(izeo) = mod(zeoz(izeo)+1.0d0,1.0d0)*zeorz
         else
            zeox(izeo) = mod(zeox(izeo)+zeorx,zeorx)
            zeoy(izeo) = mod(zeoy(izeo)+zeory,zeory)
            zeoz(izeo) = mod(zeoz(izeo)+zeorz,zeorz)
         end if 
         idzeo(izeo)=atomtype(zntype,atom)
!         if (myid.eq.0) write(16,*) zeox(izeo),zeoy(izeo),
!     &        zeoz(izeo),idzeo(izeo)
      end do

      if (myid.eq.0) then
         write(16,105) nx,ny,nz
         write(16,*) 'number of Si:',znum(1),'number of O:',znum(2)
      end if

      if (myid.eq.0) close(16)

! 99   format(i4,1x,a4,2x,3(f9.5,1x),8i4,1x,f7.3)
 99   format(i4,1x,a4,2x,3(f9.5,1x))
 100  format(/,' READING ZEOLITE LATTICE FROM FILE zeolite.cssr:',/,
     &     ' --------------------------------------------------')
 101  format(  ' box dimensions                    = ',3f10.3,
     &     ' Angstrom',/,
     &     ' simulation volume                 = ', f10.1,
     &     20x,' Angst**3')
 102  format(  ' box angles                        = ',3f10.3,
     &     ' degrees')
 103  format(  ' number of zeolite atoms           = ',i10)
 104  format(  ' mass zeolite                      = ',e10.3,
     &     ' grams',/,
     &     ' one adsorbed molecule in sim box  = ',f10.3 ,
     &     ' mmol/gram')
 105  format(  ' number of zeolite cells           = ',3i5,/,
     &     ' --------------------------------------------------')

      return
      end


