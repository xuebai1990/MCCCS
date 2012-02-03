      subroutine charge ( i, qion, vflucq, vewald )    

c charge
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
c    ** calculates the intramolcular polarization energies            **
c    ** for fluctuating charge moves                                  **
c    ** written by Marcus G. Martin 9-16-97                           **
c    ** rewritten by Bin Chen 6-25-99                                 **
c    *******************************************************************
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'poten.inc'
      include 'conver.inc'
      include 'ewaldsum.inc'

      integer i,imolty,iunit,ii,jj,ntii,ntjj,ntij,ibox
      double precision vflucq,qion(numax),qqii,vewald,rxui,ryui,rzui,
     &     rxuij,ryuij,rzuij,rij,erfunc
      vflucq = 0.0d0
      vewald = 0.0d0
      imolty = moltyp(i)
      iunit = nunit(imolty)
      ibox = nboxi(i)

C *************************************
C *** INTRACHAIN FLUCQ INTERACTIONS ***
C *************************************

c --- calculate intramolecular flucq energy for chain i      
      do ii = 1, iunit
         ntii = ntype(imolty,ii)
         qqii = qion(ii)
         
         do jj = ii, iunit
            if ( ii .eq. jj ) then
               vflucq = vflucq + xiq(ntii)*qqii
     &              + jayself(ntii)*qqii*qqii
            else
               ntjj = ntype(imolty,jj)
               ntij = (ntii-1)*nntype + ntjj
               
               vflucq = vflucq + jayq(ntij)*qqii*qion(jj)
            endif
         enddo
      enddo
c     --- remove the ground state gas phase energy
      vflucq = vflucq - fqegp(imolty)

      if ( lewald ) then
         do ii = 1, iunit
            rxui = rxu(i,ii)
            ryui = ryu(i,ii)
            rzui = rzu(i,ii)
            qqii = qion(ii)
            
            do jj = ii+1, iunit
c              --- correction term in ewald sum
               rxuij = rxui - rxu(i,jj)
               ryuij = ryui - ryu(i,jj)
               rzuij = rzui - rzu(i,jj)
               rij = dsqrt(rxuij*rxuij + ryuij*ryuij
     &              + rzuij*rzuij)
               vewald = vewald + qqii*qion(jj)*
     &              (erfunc(calp(ibox)*rij)-1.0d0)/rij
            enddo
c --- self term in ewald sum ---
            vewald = vewald - qqii*qqii*calp(ibox)/dsqrt(onepi)
         enddo
      endif
      vewald = qqfact*vewald
      return
      end






