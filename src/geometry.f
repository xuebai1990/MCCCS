      subroutine geometry(lnew,iw,i,imolty,angstart,iuprev,glist
     &     ,bondlen,bendang,phi,vvibtr,vbbtr, maxlen, wei_bend )

c geometry
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

c     ***********************************************************************
c     ** determines the new geometry of the bond lengths and angles to be  **
c     ** rotated on the cone for rosenb.f                                  **
c     ** for old computes the rosenbluth weight for bending and determines **
c     ** the old bond lenghts, angles, and phi for growth                  **
c     ** bondlen(count) is the bondlengths from the grow bead to count     **
c     ** bendang(count) is the bond angle between iuprev,iufrom, and count **
c     ** phi(count) is the angle around the cone between count and 1       **
c     ** written by M.G. Martin 7-10-98 from geomnew and geomold           **
c     ***********************************************************************

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'connect.inc'
      include 'conver.inc'
      include 'rosen.inc'
      include 'cbmc.inc'

c     --- variables passed to/from the subroutine

      logical lnew
      integer iw,imolty,angstart,iuprev,glist,i
      double precision bondlen,bendang,phi,vvibtr,vbbtr,maxlen
      double precision wei_bend

      dimension glist(numax)
      dimension bondlen(numax),bendang(numax),phi(numax)

c     --- local variables
      
      integer count,ntogrow,iugrow,iufrom,iv,juvib
     &     ,jtvib,iu2back,ib,iulast,type,aaa,iuone

      double precision equil,kforce,thetaone,thetatwo
     &     ,vvib,length,angle,phione,phitwo,vangle,random,vphi

c     --- new variables
      integer ibend,nchben_a,nchben_b,start,start_ang
      double precision bsum_try,rsint,ang_trial,bfactor,rbf,bs

      dimension ang_trial(nchbn_max),bfactor(nchbn_max)

c     --- variables from geomold
      double precision rxui,ryui,rzui,rxuij,ryuij,rzuij
     &     ,xvecprev,yvecprev,zvecprev,distprev,xvecgrow
     &     ,yvecgrow,zvecgrow,distgrow,anglec
      double precision xub,yub,zub,dum,ux,uy,uz,alpha,gamma

c     --- assign nchben_a and ncben_b
      nchben_a = nchbna(imolty)
      nchben_b = nchbnb(imolty)

c      write(2,*) 'START GEOMETRY'

c     --- initialize trial energies
      vvibtr = 0.0d0
      vbbtr = 0.0d0

c     --- assign grownum and growfrom to local variables
      ntogrow = grownum(iw)
      iufrom = growfrom(iw)

      if ( .not. lnew ) then
c        --- OLD store r*ui positions for unit iufrom
         rxui = rxu(i,iufrom)
         ryui = ryu(i,iufrom)
         rzui = rzu(i,iufrom)
      endif

c     --- Begin Bond length selection based on Boltzmann rejection
c     --- determine the bond lengths of the beads to be grown
      maxlen = 0.0d0
      do count = 1,ntogrow
         iugrow = growlist(iw,count)
         glist(count) = iugrow
            
c        --- determine the vibration (bond) type as jtvib
         do iv = 1, invib(imolty,iugrow)
            juvib = ijvib(imolty,iugrow,iv)
            if ( juvib .eq. iufrom ) then
               jtvib = itvib(imolty,iugrow,iv)
            endif
         enddo

         if ( lnew ) then
c           --- compute bond length
            equil = brvib(jtvib)
            kforce = brvibk(jtvib)
            call bondlength( jtvib,equil,kforce,beta,length,vvib )

         else
c           --- compute bond length
            rxuij = rxu(i,iugrow) - rxui
            ryuij = ryu(i,iugrow) - ryui
            rzuij = rzu(i,iugrow) - rzui

c           --- do not need mimage for intramolecular
            length = dsqrt(rxuij*rxuij + ryuij*ryuij + rzuij*rzuij)
         
c           --- compute vibration energy
            equil = brvib(jtvib)
            kforce = brvibk(jtvib)
            vvib = kforce * (length-equil )**2
         endif
            
c        --- adjust maximum bond length of those being grown
         if ( length .gt. maxlen ) maxlen = length

c        --- assign bondlength and add up vibrational energy
         bondlen(count) = length
         vvibtr = vvibtr + vvib

      enddo

c     --- Finished Bond length selection
      
c     --- determine iuprev
      if ( growprev(iw) .eq. 0 ) then
c        --- need to choose one bead to grow on unit sphere
         iuprev = growlist(iw,1)
         angstart = 2
      else
         iuprev = growprev(iw)
         angstart = 1
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- Begin Bond angle biased selection ---                 c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if ( .not. lnew ) then
c        --- compute the vector from iufrom to iuprev
         xvecprev = rxu(i,iuprev) - rxui
         yvecprev = ryu(i,iuprev) - ryui
         zvecprev = rzu(i,iuprev) - rzui
         distprev = dsqrt( xvecprev*xvecprev + yvecprev*yvecprev 
     &        + zvecprev*zvecprev )
      endif

c     --- initialize wei_bend
      wei_bend = 1.0d0

c     --- determine the iugrow-iufrom-iuprev angles 
      do count = angstart,ntogrow
         iugrow = growlist(iw,count)
         do ib = 1, inben(imolty,iugrow)
            iulast = ijben2(imolty,iugrow,ib)
            if ( iulast .eq. iufrom ) then
               iu2back = ijben3(imolty,iugrow,ib)
               if ( iu2back .eq. iuprev ) then
                  type = itben(imolty,iugrow,ib)
                  equil = brben(type)
                  kforce = brbenk(type)
               endif
            endif
         enddo

         if ( kforce .gt. 0.1d0 ) then
c           --- flexible bond angle
c           --- initialize bsum_try
            bsum_try = 0.0d0

            if ( .not. lnew ) then
c              --- first ibend is the old conformation
c              --- compute vector from iufrom to iugrow
               xvecgrow = rxu(i,iugrow) - rxui
               yvecgrow = ryu(i,iugrow) - ryui
               zvecgrow = rzu(i,iugrow) - rzui
               distgrow = bondlen(count)

c              --- dot product divided by lengths gives cos(angle)
               anglec = ( xvecprev*xvecgrow + yvecprev*yvecgrow 
     &              + zvecprev*zvecgrow ) / (distprev*distgrow)
               angle = dacos(anglec)

c              --- compute the energy of this angle
               vangle = kforce * (angle - equil)**2
               ang_trial(1) = angle
               bfactor(1) = dexp( -beta*vangle )
               bsum_try = bsum_try + bfactor(1)

c              --- skip first ibend in next loop
               start = 2

            else
c              --- new conformation start at 1
               start = 1
            endif

c           --- compute trial angles and energies
            do ibend = start,nchben_a
c               --- choose the angle uniformly on sin(theta)
               rsint = 2.0d0*random() - 1.0d0
               angle = dacos(rsint)
               ang_trial(ibend) = angle

c              --- calculate the bond angle energy
               vangle = kforce * (angle - equil)**2
               bfactor(ibend) = dexp(-beta*vangle)
               bsum_try = bsum_try + bfactor(ibend)

            enddo

            if ( lnew ) then
c              --- select one of the trial sites via bias
               rbf = random()*bsum_try
               bs = 0.0d0
               do ibend = 1,nchben_a
                  bs = bs + bfactor(ibend)
                  if ( rbf .lt. bs ) then
                     angle = ang_trial(ibend)
                     vangle = dlog(bfactor(ibend))/(-beta)
                     goto 10
                  endif
               enddo
 10            continue
            else
c              --- select the old conformation
               angle = ang_trial(1)
               vangle = dlog(bfactor(1))/(-beta)
            endif

c           --- propagate the rosenbluth weight
            wei_bend = wei_bend * bsum_try/dble(nchben_a)

         else
c           --- fixed bond angle
            angle = equil
            vangle = 0.0d0
         endif
         
         bendang(count) = angle
         vbbtr = vbbtr + vangle
         
      enddo

      if ( lnew ) then
c        --- assign phi(angstart) to 0.0
         phi(angstart) = 0.0d0
c        --- skip angstart in the loop below
         start_ang = angstart+1
      else
c        --- set up the cone using iuprev
         xub = -xvecprev/distprev
         yub = -yvecprev/distprev
         zub = -zvecprev/distprev
         
         call cone(1,xub,yub,zub,dum,dum,dum,dum,dum)
         
c        --- need to determine angstart
         start_ang = angstart
      endif
         
c     --- determine the angles of the grown beads relative to anglestart
      do count = start_ang,ntogrow
         iugrow = growlist(iw,count)

c        --- initialize bsum_try
         bsum_try = 0.0d0

         if ( .not. lnew ) then
c           --- compute vector from iufrom to iugrow
            xvecgrow = rxu(i,iugrow) - rxui
            yvecgrow = ryu(i,iugrow) - ryui
            zvecgrow = rzu(i,iugrow) - rzui
            distgrow = bondlen(count)

c           --- turn this into a unit vector
            ux = xvecgrow/distgrow
            uy = yvecgrow/distgrow
            uz = zvecgrow/distgrow
            alpha = bendang(count)

c           --- compute gamma (phi) for these unit vectors
            call cone(3,dum,dum,dum,alpha,gamma,ux,uy,uz)

            phitwo = gamma

            call cone(2,dum,dum,dum,alpha,gamma,ux,uy,uz)

            if ( angstart .eq. count ) then
c              --- this is the first phi, no need to compute bias
               vphi = 0.0d0
               goto 20
            else
               vphi = 0.0d0
               do aaa = angstart,count-1
                  iuone = growlist(iw,aaa)
                  phione = phi(aaa)
                  thetaone = bendang(aaa)
                  thetatwo = bendang(count)
                  call coneangle(thetaone,phione,thetatwo,phitwo
     &                 ,angle)
            
                  do ib = 1, inben(imolty,iugrow)
                     iulast = ijben2(imolty,iugrow,ib)
                     if ( iulast .eq. iufrom ) then
                        iu2back = ijben3(imolty,iugrow,ib)
                        if ( iu2back .eq. iuone ) then
                           type = itben(imolty,iugrow,ib)
c                          --- calculate the bond angle energy
                           vphi = vphi + brbenk(type) 
     &                          * (angle - brben(type))**2
                        endif
                     endif
                  enddo
               enddo

               ang_trial(1) = phitwo
               bfactor(1) = dexp( -beta * vphi )
               bsum_try = bsum_try + bfactor(1)

            endif
c           --- skip first ibend for OLD
            start = 2
         else
c           --- perform all ibend for NEW
            start = 1
         endif

c        --- compute trial energies and weights
         do ibend = start,nchben_b

c           --- determine a value of phitwo
            phitwo = random()*twopi
            vphi = 0.0d0
            do aaa = angstart,count-1
               iuone = growlist(iw,aaa)
               phione = phi(aaa)
               thetaone = bendang(aaa)
               thetatwo = bendang(count)
               call coneangle(thetaone,phione,thetatwo,phitwo
     &              ,angle)
            
               do ib = 1, inben(imolty,iugrow)
                  iulast = ijben2(imolty,iugrow,ib)
                  if ( iulast .eq. iufrom ) then
                     iu2back = ijben3(imolty,iugrow,ib)
                     if ( iu2back .eq. iuone ) then
                        type = itben(imolty,iugrow,ib)
c                       --- calculate the bond angle energy
                        vphi = vphi + brbenk(type) 
     &                       * (angle - brben(type))**2
                     endif
                  endif
               enddo
            enddo

c           --- store the boltzmann factors and phi
            bfactor(ibend) = dexp(-beta*vphi)
            ang_trial(ibend) = phitwo
            bsum_try = bsum_try + bfactor(ibend)

         enddo
   
         if ( lnew ) then
c           --- select a value of phitwo in a biased fashion
            rbf = random()*bsum_try
            bs = 0.0d0
            do ibend = 1,nchben_b
               bs = bs + bfactor(ibend)
               if ( rbf .lt. bs ) then
                  phitwo = ang_trial(ibend)
                  vphi = dlog( bfactor(ibend) )/(-beta)
                  goto 15
               endif
            enddo
 15         continue
         else
c           --- select the OLD value of phitwo
            phitwo = ang_trial(1)
            vphi = dlog( bfactor(1) )/(-beta)
         endif

c        --- propagate angle weight
         wei_bend = wei_bend * (bsum_try/dble(nchben_b))

 20      continue
c        --- store the angle for thetatwo
         phi(count) = phitwo
         vbbtr = vbbtr + vphi
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     --- End Bond angle biased selection ---           c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      write(2,*) 'FINISH GEOMETRY'

      return
      end
      subroutine bendangle(equil,kforce,beta,angle,vangle )

c bendagle 
c Copyright (C) 2000 Marcus Martin, Bin Chen, Collin Wick, and Ilja Siepmann


c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.

c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.

c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.



c     *************************************************************
c     *** choose a bond angle                                   ***
c     *** equil is equilibrium bond angle                       ***
c     *** kforce is the force constant                          ***
c     *** beta is 1/kT                                          ***
c     *** angle is the angle returned by this subroutine        ***
c     *** vangle is the energy of the angle                     ***
c     *** M.G. Martin                                           ***
c     *************************************************************

      implicit none

      double precision equil, kforce, beta, angle, vangle, rr, v1, v2
     &     ,random

c      write(2,*) 'start BENDANGLE'

      if ( kforce .gt. 0.1d0 ) then
c        --- find a vector inside the unit sphere
 80      v1 = 2.0d0*random()-1.0d0
         v2 = 2.0d0*random()-1.0d0
         rr = v1*v1 + v2*v2
         if (rr .ge. 1.0d0 ) goto 80
         
c        --- select angle from a gaussian distribution
         angle = equil + v1*dsqrt( (-dlog(rr))
     &        /( kforce*beta*rr) )
         
         if (angle .le. 0.0d0 .or. angle .ge. 3.1415 ) then
            write(2,*) 'chose angle outside of 0,Pi in bendangle'
            goto 80
         endif

c        --- correct for the phase space of the angle
         if ( random() .gt. dsin(angle) ) goto 80
         vangle = kforce * (angle - equil)**2
      else
c        --- fixed bond angle
         angle = equil
         vangle = 0.0d0
      endif

c      write(2,*) 'end BENDANGLE'

      return
      end
