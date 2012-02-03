      subroutine Intra_energy ( i,imolty, v, vintra, vinter,vext
     &     ,velect,vewald,flagon,ibox, istart, iend,lljii,ovrlap
     &     ,ltors,vtors,lcharge_table,lfavor,vvib,vbend,vtg)

c energy
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
 
c    *******************************************************************
c    ** calculates the total potential energy for a configuration.    **
c    *******************************************************************
 
      implicit none

c *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'neigh.inc'
      include 'poten.inc'
      include 'coord2.inc' 
      include 'external.inc'
      include 'connect.inc'
      include 'ewaldsum.inc'
      include 'fepsi.inc'
      include 'qqlist.inc'
      include 'clusterbias.inc'
      include 'nsix.inc'
      include 'peboco.inc'
      include 'cell.inc'

      logical lqimol,lqjmol,lexplt,lcoulo,lfavor,lij2,liji,lqchgi
      logical lljii,ovrlap,ltors,lcharge_table,lt,lfound

      integer growii,growjj,k,cellinc,jcell,ic,nmole
      integer i,ibox, istart, iend,ii,ntii,flagon,jjj,iii
     +       ,j,jj,ntjj,ntij,ntj,imolty,jmolty,ncell
      integer iivib,jjtor,ip1,ip2,ip3,it,nchp2,acellinc

      integer jjvib,jjben  

      double precision vvib,vbend,vtg,theta

      double precision ljsami,ljpsur,ljmuir,v,vintra, vinter,vext 
     +                ,rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq
     +                ,ryuij,rzuij,sr2,sr6,rij,rijsq,dzui,dz3,dz12
     +                ,exgrph,exsami,exmuir,exzeo,vtors,exsix,velect
     +                ,vewald,mmff,rbcut,rvdwsq,rchgsq,ninesix
      double precision erfunc,qave
      double precision xvec,yvec,zvec,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2
     &     ,daa1,da1a2,dot,thetac,vtorso
      double precision xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,epsilon2,sigma2
      double precision sx,sy,sz
      double precision slitpore	
      double precision distij(numax,numax)
      double precision xcc,ycc,zcc,tcc,spltor

      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)
      dimension lcoulo(numax,numax),cellinc(cmax),jcell(nmax)
      dimension acellinc(numax,27)

C --------------------------------------------------------------------

c      write(2,*) 'start ENERGY'
      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut = rcut(ibox)
      if (ldual) rcinsq = rcutin*rcutin
      rminsq = rmin * rmin

      v = 0.0d0
      vvib = 0.0d0
      vbend = 0.0d0
      vtg = 0.0d0
 
c - branched and linear molecules with connectivity table -
c - go through entire chain -
c - calculate all bonds vectors and lengths
c - calculate all stretching, bending, and torsional potentials
c - that have an end-bead with an index smaller than the current bead
             do ii = 1, nunit(imolty)
                rxui=rxu(i,ii)
                ryui=ryu(i,ii)
                rzui=rzu(i,ii)
                do iivib = 1, invib(imolty,ii)
                   jj = ijvib(imolty,ii,iivib)
c                   xvec(ii,jj) = rxu(i,jj) - rxui
c                   yvec(ii,jj) = ryu(i,jj) - ryui
c                   zvec(ii,jj) = rzu(i,jj) - rzui
                   xvec(ii,jj) = rxu(i,jj) - rxui
                   yvec(ii,jj) = ryu(i,jj) - ryui
                   zvec(ii,jj) = rzu(i,jj) - rzui 
                   distij(ii,jj) = dsqrt( xvec(ii,jj)**2
     +                 + yvec(ii,jj)**2 + zvec(ii,jj)**2 )

                   if ( nunit(imolty) .ne. nugrow(imolty) )then
c                  --- account for explct atoms in opposite direction
                      xvec(jj,ii)   = -xvec(ii,jj)
                      yvec(jj,ii)   = -yvec(ii,jj)
                      zvec(jj,ii)   = -zvec(ii,jj)
                      distij(jj,ii) = distij(ii,jj)
                   endif
                enddo
             enddo

c - stretching -
c             if ( brvibk(1) .gt. 0.01d0 .or. lninesix) then
                do j = 2, nunit(imolty)
                   do jjvib = 1, invib(imolty,j)
                      ip1 = ijvib(imolty,j,jjvib)
                      it  = itvib(imolty,j,jjvib)
                      if ( ip1 .lt. j ) vvib = vvib +
     +             brvibk(it) * ( distij(ip1,j) - brvib(it) )**2
                   enddo
                enddo
c             endif


c - bending -
c ### molecule with bond bending
             do j = 2, nunit(imolty)
                do jjben = 1, inben(imolty,j)
                   ip2 = ijben3(imolty,j,jjben)
                   if ( ip2 .lt. j ) then
                      ip1 = ijben2(imolty,j,jjben)
                      it  = itben(imolty,j,jjben)
                      thetac = ( xvec(ip1,j)*xvec(ip1,ip2) +
     +                     yvec(ip1,j)*yvec(ip1,ip2) +
     +                     zvec(ip1,j)*zvec(ip1,ip2) ) /
     +                     ( distij(ip1,j)*distij(ip1,ip2) )
                      if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
                      if ( thetac .le. -1.0d0 ) thetac = -1.0d0

                      theta = dacos(thetac)
                      vbend = vbend +
     +                     brbenk(it) * (theta-brben(it))**2

c                      write(2,*) 'ip2,ip1,j',ip2,ip1,j
c                      write(2,*) 'bend energy, theta '
c     &                     ,brbenk(it) * (theta-brben(it))**2,theta
                   endif
                enddo
             enddo

c - torsions -
c ### molecule with dihedral potenials ###
             do j = 2, nunit(imolty)
                do jjtor = 1, intor(imolty,j)
                   ip3 = ijtor4(imolty,j,jjtor)
                   if ( ip3 .lt. j ) then
                      ip1 = ijtor2(imolty,j,jjtor)
                      ip2 = ijtor3(imolty,j,jjtor)
                      it  = ittor(imolty,j,jjtor)
c*** calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                      xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     +                     zvec(ip1,j) * yvec(ip1,ip2)
                      yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     +                     xvec(ip1,j) * zvec(ip1,ip2)
                      zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     +                     yvec(ip1,j) * xvec(ip1,ip2)
                      xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     +                     zvec(ip1,ip2) * yvec(ip3,ip2)
                      ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     +                     xvec(ip1,ip2) * zvec(ip3,ip2)
                      za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     +                     yvec(ip1,ip2) * xvec(ip3,ip2)
c *** calculate lengths of cross products ***
                      daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                      da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
c *** calculate dot product of cross products ***
                      dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                      thetac = - dot / ( daa1 * da1a2 )
c     KEA -- added for extending range to +/- 180
c     additional definitions for torsions
                     if (L_tor_table) then
c     *** calculate cross product of cross products ***
                       xcc = yaa1*za1a2 - zaa1*ya1a2
                       ycc = zaa1*xa1a2 - xaa1*za1a2
                       zcc = xaa1*ya1a2 - yaa1*xa1a2
c     *** calculate scalar triple product ***
                       tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2)
     &                    + zcc*zvec(ip1,ip2)
                       theta = dacos(thetac)
                       if (tcc .lt. 0.0d0) theta = -theta
                       if (L_spline) then
                          call splint(theta,spltor,it)
                       elseif(L_linear) then
                          call lininter(theta,spltor,it)
                       endif
                       vtg = vtg + spltor
                     else
                      vtg = vtg + vtorso( thetac, it )
!                       write(17,*) j,it,vtg,
!     &                           vtorso( thetac, it )
                     endif
                   endif
                enddo
             enddo

C----------------------------------------
     
            velect = velect*qqfact
            vewald = vewald*qqfact

 
c      velect = velect * qqfact
c      vewald = vewald * qqfact

c     note that vintra is only computed when the flag lljii is true
      v = vinter + vext + vintra + velect + vewald + vvib + vbend + vtg 
C  NEERAJ: Debugging start

c      write(2,*) 'vinter:',vinter,'vext:',vext,'vintra:',vintra,'velect'
c     & ,velect,'vewald:'vewald

C  NEERAJ: DEbugging end

c      write(2,*) 'end ENERGY'

      return
      end





