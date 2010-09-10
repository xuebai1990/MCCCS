        function exgrph (x,y,z,ntij)

c exgrph
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

c - calculates the energy of a bead with a graphite surface


        implicit none
        
        include 'control.inc'
        include 'external.inc'
        include 'poten.inc'
        real(8)::aa,aa2
        real(8)::a1sq
        real(8)::e0,e1
        real(8)::exgrph
        real(8)::fxy
        real(8)::x,y,z
        real(8)::pi2,pi4
        real(8)::bb,cc,dd
        real(8)::k2,k5
        real(8)::mbessel
        real(8)::sz2
        real(8)::zz
        integer::ntij
        
        parameter (pi2 = 6.28318530718d0, pi4 = 12.5663706144d0)

        exgrph = 0.0d0
        e0 = 0.0d0
        e1 = 0.0d0
        fxy = 0.0d0

c       write(81,*) ntij,sqrt(sig2ij(ntij))
        
        sz2 = sig2ij(ntij)/(z**2)
        
        aa = pi2 * rsol * delta * sig2ij(ntij)

        e0 = aa*epsij(ntij)*((2.0d0/5.0d0)*(sz2**5) - (sz2**2) -
     +          (sig2ij(ntij)**2/(3.0d0*delta*(0.61*delta+z)**3)))

c       write(82,*) e0,aa,delta,z       
        if ( lcorreg ) then
                a1sq = a1**2
        
                aa2 = (sig2ij(ntij)/a1sq)**3
        
                bb = aa2*pi4*epsij(ntij)/sqrt(3.0d0)
                
c               bb = pi4*epsij(ntij)*sig2ij(ntij)**3/
c     +                 (sqrt(3.0d0)*a1**6)     

                cc = aa2/(30.0d0*(pi2/sqrt(3.0d0))**5)
                
c               cc = sig2ij(ntij)**6/
c     +         (30.0d0*a1**6*(pi2/sqrt(3.0d0)**5))

                dd = 2.0d0*(pi2/sqrt(3.0d0))**2
                
                zz = pi4*z/(sqrt(3.0d0)*a1)
                
                k2 = mbessel(zz,2.0d0)
                
                k5 = mbessel(zz,5.0d0)
c               write(84,*) zz,k2,k5            
                e1 = bb*(cc * k5 * (a1/z)**5 - dd * k2 * (a1/z)**2)
c               write(82,*) bb,cc,dd,e1,k2,k5
                fxy = -2.0d0*(cos(pi2*(x/a1 + y/sqrt(3.0d0)/a1)) +
     +                  cos(pi2*(x/a1 - y/sqrt(3.0d0)/a1)) +
     +                  cos(pi4*y/sqrt(3.0d0)/a1))

c       write(82,'(6g12.5)') x,y,z,fxy,e1,e0
                exgrph = e0 + e1*fxy
        else
c       write(83,'(6g12.5)') x,y,z,e0
                exgrph = e0     
        endif
        end function exgrph


