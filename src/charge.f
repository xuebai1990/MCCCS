      subroutine charge ( i, qion, vflucq, vewald )    

!    *******************************************************************
!    ** calculates the intramolcular polarization energies            **
!    ** for fluctuating charge moves                                  **
!    ** written by Marcus G. Martin 9-16-97                           **
!    ** rewritten by Bin Chen 6-25-99                                 **
!    *******************************************************************
 
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'poten.inc'
!$$$      include 'conver.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'connect.inc'      

      integer(KIND=normal_int)::i,imolty,iunit,ii,jj,ntii,ntjj,ntij,ibox
      real(KIND=double_precision)::vflucq,qion(numax),qqii,vewald,rxui ,ryui,rzui,rxuij,ryuij,rzuij,rij
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
               vflucq = vflucq + xiq(ntii)*qqii + jayself(ntii)*qqii*qqii
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
               rij = dsqrt(rxuij*rxuij + ryuij*ryuij + rzuij*rzuij)
               if (.not.lqinclu(imolty,ii,jj)) then
                   vewald = vewald + qqii*qion(jj)* (erfunc(calp(ibox)*rij)-1.0d0)/rij
               else
                   vewald = vewald + (1.0d0-qscale2(imolty,ii,jj))*qqii* qion(jj)* (erfunc(calp(ibox)*rij)-1.0d0)/rij
               end if 
            end do
! --- self term in ewald sum ---
            vewald = vewald - qqii*qqii*calp(ibox)/dsqrt(onepi)
         end do
      end if
      vewald = qqfact*vewald
      return
      end






