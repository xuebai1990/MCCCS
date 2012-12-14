      subroutine bondlength(vibtype,requil,kvib,betaT,length,vvib)

!     **************************************************************
!     *** computes the bond length for a given vibration type    ***
!     *** vibtype is the vibration type                          ***
!     *** requil is the equilibrium bond length                  ***
!     *** kvib is the force constant for the bond length         ***
!     *** length is the bond length returned by this subroutine  ***
!     *** betaT is 1/kT                                           ***
!     *** vvib is the vibration energy for this bond             ***
!     *** M.G. Martin  2-4-98                                    ***
!     **************************************************************

      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'tabulated.inc'

      integer(KIND=normal_int)::vibtype
      real(KIND=double_precision)::length,bond,random,bf,vvib,betaT ,kvib,requil,tabulated_vib

      vvib = 0.0d0

      if ( kvib .gt. 0.1d0 ) then
!     --- random bond length from Boltzmann distribution ---
 107     bond = (0.2d0*random()+0.9d0)
         
!     --- correct for jacobian by dividing by the max^3
!     --- 1.331  = (1.1)^3, the equilibrium bond length cancels out
         bf = bond*bond*bond/(1.331d0)
         if ( random() .ge. bf ) goto 107
         bond = bond * requil
         
!     --- correct for the bond energy
         vvib = kvib * (bond-requil )**2
         bf = dexp ( -(vvib * betaT) )
         if ( random() .ge. bf ) goto 107
         length = bond

!  --- tabulated bond stretching potential
!  --- added 12/02/08 by KM
      else if (L_vib_table) then
!     --- random bond length from Boltzmann distribution ---
!     ---  +/- 25% of equilibrium bond length
 108     bond = (0.5d0*random()+0.75d0)
         
!     --- correct for jacobian by dividing by the max^3
!     --- 1.953125  = (1.25)^3, the equilibrium bond length cancels out
         bf = bond*bond*bond/(1.953125d0)
         if ( random() .ge. bf ) goto 108
         
         bond = bond * requil
         call lininter_vib(bond, tabulated_vib, vibtype)
         vvib=tabulated_vib

         bf = dexp ( -(vvib * betaT) )
         
         if ( random() .ge. bf ) goto 108
         length = bond

      else
!     --- fixed bond length
         length = requil
      end if
      
      return

      end
