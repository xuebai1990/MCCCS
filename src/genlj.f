      function genlj (rijsq,sr2,epsilon2)

!  ********************************************************************
!  ***  calculates the energy of the Generalized Lennard Jones
!       potential                                               ***
!  ***      parameters defined  poten.inc           ***
!  *******************************************************************
!   Also you should make sure that the rcut you choose
!  is .gt. sigma*(2**(2/n0)) (because currently tail corrections
!  are evaluated with this assumption which is rarely incorrect.)
!  ********************************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'control.inc'
!$$$      include 'poten.inc'
      real(KIND=double_precision)::rijsq,rij,srij,sr2,epsilon2,genlj
     
      rij=(dsqrt(rijsq))
      srij=dsqrt(sr2)

      if ( (rij) .le.(rij*srij)*2.0d0**(2.0d0/n0) ) then
         genlj = 4.0d0*epsilon2*(((srij)**n0)-((srij)**(n0/2.0d0)))
      else
         genlj =epsilon2*(((2.0d0**((4.0d0*n1/n0)))*((srij)**
     &  (2.0d0*n1)))-((2.0d0**((2.0d0*(n1/n0))+1.0d0))*((srij)**(n1))))
      end if

      
! In reduced units      
!      xij=(1.0d0/dsqrt(sr2))
!      
!      if ( (xij) .le. 2.0d0**(2.0d0/n0) ) then
!         genlj = 4.0d0 * (((1/xij)**n0)-((1/xij)**(n0/2.0d0)))
!      else
!         genlj =(( (2.0d0**((4.0d0*n1/n0)))*(1/xij)**(2.0d0*n1))-
!     & (2.0d0**((2.0d0*(n1/n0))+1.0d0)*((1/xij)**(n1))))
!      end if
!      vinter=vinter+e
       

!      not usually used in MC    
!      if (lshift) then
!         genlj = genlj - shiftgenlj()
!      end if
      
      return
      end

      
