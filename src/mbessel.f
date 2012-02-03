      function mbessel(z,nu)

c mbessel
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

        double precision z,nu
        double precision mbessel
        double precision pi
        parameter (pi = 3.14159265359d0)
        
c       mbessel = sqrt(pi/(2.0d0*z))*exp(-z)*
c     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
c     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))

c -- simple form     
        mbessel = sqrt(pi/(2.0d0*z))*exp(-z)
c     +         (1.0d0 + (4.0d0*nu**2-1)/(8.0d0*z) +
c     +         (4.0d0*nu**2-1)*(4.0d0*nu**2-9.0d0)/(2.0d0*64.0d0*z**2))
        end function mbessel

