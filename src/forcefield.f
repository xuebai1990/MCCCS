      subroutine forcefield(rczeo)

c forcefield
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
   
      include 'zeolite.inc'
      include 'zeopoten.inc'
      include 'control.inc'

      double precision epsilon,sigma,rczeo
      integer          itype,jtype,idz

C     In atomtype.f atoms O of the zeolite have been given id=1
C     Hopefully in some other place the guestmolecules CH4 have been
C     assigned id=2


c      if (zntype.ne.2) stop '** forcefield: ntype ne 2 not allowed **'

C     Force field data:

C     The interaction parameters are in this case calculated from 
C     atomic values

C     Fill the matrices

C Needed are epsilon and sigma such that in Angstroms and kcal/mol:
C Energy = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

      read(25,*)  
      read(25,*)
      idz=4
      do itype = 1,zntype
	 read(25,*) epsilon,sigma
         zeps(itype,idz) = epsilon
         zsig2(itype,idz) = sigma*sigma
      enddo

C Write epsilons and sigma's:

      write(16,100)
      write( 6,100)
      jtype=idz
      do itype = 1,zntype
         write(16,1001) itype,jtype,zeps(itype,jtype),
     +                  itype,jtype,dsqrt(zsig2(itype,jtype))
         write( 6,1001) itype,jtype,zeps(itype,jtype),
     +                  itype,jtype,dsqrt(zsig2(itype,jtype))
      enddo
 
C Force field cutoff radius rczeo for interactions between guests
C and zeolite atoms:

      write(16,101) rczeo
      write( 6,101) rczeo
      do itype = 1,zntype
        zrc2(itype,idz)=rczeo**2
      enddo
c
c === calculate cutt-off of the potential
c
      if (lshift) then
        write(6,102)
        do itype = 1,zntype
          zencut(itype,idz)=4.*zeps(itype,idz)*
     +         ( (zsig2(itype,idz)/zrc2(itype,idz))**6
     +          -(zsig2(itype,idz)/zrc2(itype,idz))**3)
          write(6,103) itype,zencut(itype,idz)
        enddo
      endif     

  100 format(/,' FORCE FIELD DEFINITION:',/,
     +         ' ------------------------------------------------')
  101 format(/,' Force field cutoffs: rczeo = ',f6.1,/,
     +         '                      rcads = ',f6.1,' Angstrom',/,
     +         ' ------------------------------------------------',//)
 1000 format(/,' Lennard-Jones Parameters (K, Angstrom):')
 1001 format(  ' epsilon(',i2,',',i2,') = ',f6.1,
     +         ' K &  sigma(',i2,',',i2,') = ',f6.1,' Angstrom')
 102  format(' Value at cut-off distance ' )
 103  format('    interaction ',i3, '   : ',f8.2,'[K]')        
      return
      end

