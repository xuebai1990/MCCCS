      subroutine calctor(iu1,iu2,iu3,iu4,jttor,vtor)

c calctor
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

      include 'control.inc'
      include 'fix.inc'

      integer iu1,iu2,iu3,iu4,jttor
      
      double precision thetac,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2
     &     ,daa1,da1a2,dot,vtor,vtorso
      
c     --- calculate cross products d_a x d_a-1 
      xaa1 = yvec(iu2,iu1) * zvec(iu3,iu2)
     &     + zvec(iu2,iu1) * yvec(iu2,iu3)
      yaa1 = zvec(iu2,iu1) * xvec(iu3,iu2) 
     &     + xvec(iu2,iu1) * zvec(iu2,iu3)
      zaa1 = xvec(iu2,iu1) * yvec(iu3,iu2) 
     &     + yvec(iu2,iu1) * xvec(iu2,iu3)
      
c     --- calculate cross products d_a-1 x d_a-2
      xa1a2 = yvec(iu2,iu3) * zvec(iu3,iu4) -
     +     zvec(iu2,iu3) * yvec(iu3,iu4)
      ya1a2 = zvec(iu2,iu3) * xvec(iu3,iu4) -
     +     xvec(iu2,iu3) * zvec(iu3,iu4)
      za1a2 = xvec(iu2,iu3) * yvec(iu3,iu4) -
     +     yvec(iu2,iu3) * xvec(iu3,iu4)
   
c     --- calculate lengths of cross products ***
      daa1 = dsqrt ( xaa1**2 + yaa1**2 + zaa1**2 )
      da1a2 = dsqrt ( xa1a2**2 + ya1a2**2 
     &     + za1a2**2 )
      
c     --- calculate dot product of cross products ***
      dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
      thetac = - (dot / ( daa1 * da1a2 ))

      vtor = vtorso(thetac,jttor)

      return
      end
