      subroutine charge ( i, qion, vflucq, vewald )    

! charge
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
 
!    *******************************************************************
!    ** calculates the intramolcular polarization energies            **
!    ** for fluctuating charge moves                                  **
!    ** written by Marcus G. Martin 9-16-97                           **
!    ** rewritten by Bin Chen 6-25-99                                 **
!    *******************************************************************
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'poten.inc'
      include 'conver.inc'
      include 'ewaldsum.inc'
      include 'connect.inc'      

      integer::i,imolty,iunit,ii,jj,ntii,ntjj,ntij,ibox
      real(8)::vflucq,qion(numax),qqii,vewald,rxui,ryui,rzui,
     &     rxuij,ryuij,rzuij,rij,erfunc
      vflucq = 0.0d0
      vewald = 0.0d0
      imolty = moltyp(i)
      iunit = nunit(imolty)
      ibox = nboxi(i)

! *************************************
! *** INTRACHAIN FLUCQ INTERACTIONS ***
! *************************************

! --- calculate intramolecular flucq energy for chain i      
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
            end if
         end do
      end do
!     --- remove the ground state gas phase energy
      vflucq = vflucq - fqegp(imolty)

      if ( lewald ) then
         do ii = 1, iunit
            rxui = rxu(i,ii)
            ryui = ryu(i,ii)
            rzui = rzu(i,ii)
            qqii = qion(ii)
            
            do jj = ii+1, iunit
!              --- correction term in ewald sum
               rxuij = rxui - rxu(i,jj)
               ryuij = ryui - ryu(i,jj)
               rzuij = rzui - rzu(i,jj)
               rij = dsqrt(rxuij*rxuij + ryuij*ryuij
     &              + rzuij*rzuij)
               if (.not.lqinclu(imolty,ii,jj)) then
                   vewald = vewald + qqii*qion(jj)*
     &              (erfunc(calp(ibox)*rij)-1.0d0)/rij
               else
                   vewald = vewald + (1.0d0-qscale2(imolty,ii,jj))*qqii*
     &                       qion(jj)*
     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
               end if 
            end do
! --- self term in ewald sum ---
            vewald = vewald - qqii*qqii*calp(ibox)/dsqrt(onepi)
         end do
      end if
      vewald = qqfact*vewald
      return
      end






