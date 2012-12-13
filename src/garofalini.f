      function garofalini(rijsq,ntij,qa,qb,aa,bb)

!     **************************************************************
!     ***  calculates the energy using the garofalini (SiO2/H2O) ***
!     ***  exp-6 potential modified using articles in suijtab    ***
!     ***  parameters are defined in suijtab.f  KE ANDERSON      ***
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

!$$$      include 'control.inc'
!$$$      include 'garofalini.inc'
!$$$      include 'poten.inc'

      real(KIND=double_precision)::rijsq,rij,hterm,coul,qa,qb
      real(KIND=double_precision)::hterma,garofalini
      integer(KIND=normal_int)::ntij,aa,bb,i

!      write(iou,*) 'input',rijsq,ntij,qa,qb
      rij = dsqrt(rijsq)
      hterm = 0.0d0
      coul = 0.0d0
      garofalini = 0.0d0


!      write(iou,*) aa,bb,ntij,qa,qb,gbeta(ntij),galpha(ntij),grho(ntij)
!     H term
      do i=1,3
         hterma =   ga(ntij,i)/(1+dexp(gb(ntij,i)*(rij-gc(ntij,i))))
!         write(iou,*) i,hterma,' (',ntij,')'
         hterm = hterm + hterma
      end do

      coul = qa*qb*erfunc(rij/gbeta(ntij))/rij
!      write(iou,*) 'erfunc',coul*rij
      coul = coul * qqfact

      garofalini = galpha(ntij)*dexp(-rij/grho(ntij)) + hterm  + coul

!      write(iou,*) 'i,j,v2',aa,bb,'H term:',hterm,' Coul term:'
!     &     ,coul,' the rest:',garofalini-hterm-coul,
!     &     ' Total:',garofalini
!      write(iou,*) '                     ',rij

!      write(iou,*) 'leaving garofalini'

      return
      end

      
