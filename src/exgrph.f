        function exgrph(x,y,z,ntij)

! - calculates the energy of a bead with a graphite surface

        use sim_system
        use var_type
        use const_math
        use util_math
        implicit none
!$$$        include 'control.inc'
!$$$        include 'external.inc'
!$$$        include 'poten.inc'
        real(KIND=double_precision)::aa,aa2
        real(KIND=double_precision)::a1sq
        real(KIND=double_precision)::e0,e1
        real(KIND=double_precision)::exgrph
        real(KIND=double_precision)::fxy
        real(KIND=double_precision)::x,y,z
        real(KIND=double_precision)::bb,cc,dd
        real(KIND=double_precision)::k2,k5
        real(KIND=double_precision)::sz2
        real(KIND=double_precision)::zzz
        integer(KIND=normal_int)::ntij
        
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


