	function slitpore(z,ntij)
! -- calculates the surface energy of a bend with a featureless
! -- graphite surface

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
	
!$$$	include 'control.inc'
!$$$	include 'external.inc'
!$$$	include 'poten.inc'	
	integer(KIND=normal_int)::twopi
	integer(KIND=normal_int)::ntij
	parameter (twopi=6.283185307179586d0)
	real(KIND=double_precision)::vgs,z
	real(KIND=double_precision)::slitpore
	real(KIND=double_precision)::sig
!	real(KIND=double_precision)::coef1,coef2,coef3,coef4
	
	sig = sqrt(sig2ij(ntij))
	
!	coef1 = twopi*rsol*epsij(ntij)*sig2ij(ntij)*delta
!	coef2 = 2.0/5.0*(sig/z)**10
!	coef3 = (sig/z)**4
!	coef4 = sig2ij(ntij)**2/(3*delta*(z+0.61*delta)**3)
	
	vgs = twopi*rsol*epsij(ntij)*sig2ij(ntij)*delta
     &		*(2.0/5.0*(sig/z)**10
     &		-(sig/z)**4
     &		-(sig2ij(ntij)**2/(3*delta*(z+0.61*delta)**3)))
     
     	slitpore = vgs
	end
