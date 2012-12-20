MODULE energy_external
  use const_math,only:twopi,fourpi
  use util_math,only:mbessel
  use sim_system,only:sig2ij,epsij,qqu,lcorreg
  implicit none
  private
  save
  public::slitpore,exgrph,v_elect_field,ntsubst

! EXTERNAL.INC
  real::a1,delta,rsol
  integer::ntsubst
! -- Slitpore
! -- a1, delta are in Angstroms, rsol [=] 1/A^3
  parameter (a1 = 2.460d0)
  parameter (delta = 3.40d0)
  parameter (ntsubst = 190)
  parameter (rsol = 0.114d0)
! Joe Hautman's parameters are read from standard input, to be used with ljoe = .true.
! AT PRESENT no parameters for polymeric surfactants (see LJPSUR)

contains
! -- calculates the surface energy of a bend with a featureless
! -- graphite surface
  function slitpore(z,ntij)
!$$$	include 'control.inc'
!$$$	include 'external.inc'
!$$$	include 'poten.inc'
	integer::ntij
	real::vgs,z
	real::slitpore
	real::sig
!	real::coef1,coef2,coef3,coef4

	sig = sqrt(sig2ij(ntij))

!	coef1 = twopi*rsol*epsij(ntij)*sig2ij(ntij)*delta
!	coef2 = 2.0/5.0*(sig/z)**10
!	coef3 = (sig/z)**4
!	coef4 = sig2ij(ntij)**2/(3*delta*(z+0.61*delta)**3)

	vgs = twopi*rsol*epsij(ntij)*sig2ij(ntij)*delta*(2.0/5.0*(sig/z)**10 		-(sig/z)**4	-(sig2ij(ntij)**2/(3*delta*(z+0.61*delta)**3)))

    slitpore = vgs
  end function slitpore

! - calculates the energy of a bead with a graphite surface
  function exgrph(x,y,z,ntij)
!$$$        include 'control.inc'
!$$$        include 'external.inc'
!$$$        include 'poten.inc'
        real::aa,aa2
        real::a1sq
        real::e0,e1
        real::exgrph
        real::fxy
        real::x,y,z
        real::bb,cc,dd
        real::k2,k5
        real::sz2
        real::zzz
        integer::ntij

        exgrph = 0.0d0
        e0 = 0.0d0
        e1 = 0.0d0
        fxy = 0.0d0

!       write(81,*) ntij,sqrt(sig2ij(ntij))

        sz2 = sig2ij(ntij)/(z**2)

        aa = twopi * rsol * delta * sig2ij(ntij)

        e0 = aa*epsij(ntij)*((2.0d0/5.0d0)*(sz2**5) - (sz2**2) - (sig2ij(ntij)**2/(3.0d0*delta*(0.61*delta+z)**3)))

!       write(82,*) e0,aa,delta,z
        if ( lcorreg ) then
                a1sq = a1**2

                aa2 = (sig2ij(ntij)/a1sq)**3

                bb = aa2*fourpi*epsij(ntij)/sqrt(3.0d0)

!               bb = fourpi*epsij(ntij)*sig2ij(ntij)**3/
!     +                 (sqrt(3.0d0)*a1**6)

                cc = aa2/(30.0d0*(twopi/sqrt(3.0d0))**5)

!               cc = sig2ij(ntij)**6/
!     +         (30.0d0*a1**6*(twopi/sqrt(3.0d0)**5))

                dd = 2.0d0*(twopi/sqrt(3.0d0))**2

                zzz = fourpi*z/(sqrt(3.0d0)*a1)

                k2 = mbessel(zzz,2.0d0)

                k5 = mbessel(zzz,5.0d0)
!               write(84,*) zzz,k2,k5
                e1 = bb*(cc * k5 * (a1/z)**5 - dd * k2 * (a1/z)**2)
!               write(82,*) bb,cc,dd,e1,k2,k5
                fxy = -2.0d0*(cos(twopi*(x/a1 + y/sqrt(3.0d0)/a1)) + cos(twopi*(x/a1 - y/sqrt(3.0d0)/a1)) + cos(fourpi*y/sqrt(3.0d0)/a1))

!       write(82,'(6g12.5)') x,y,z,fxy,e1,e0
                exgrph = e0 + e1*fxy
        else
!       write(83,'(6g12.5)') x,y,z,e0
                exgrph = e0
        end if
  end function exgrph

! **********************************************************************
! **  calculates interaction of molecule i with an external field E  ***
! **  added 06/24/07 by KM                                           ***
! **********************************************************************
  function v_elect_field(i, j, rzfield,E)
!$$$      include 'control.inc'
!$$$      include 'external.inc'
!$$$      include 'coord.inc'

      real::v_elect_field,rzfield, E
      integer::i,j


! ********************************************
! **  units
! **  E in V/A, q in e, rz in A
! **  E*q*rz = V*e
! **  1 V*e = 11600 K
! ********************************************

      v_elect_field = -E*rzfield*qqu(i,j)

!      write(io_output,*) 'E ', E, ' exfield ', exfield

      return
  end function v_elect_field
end MODULE energy_external
