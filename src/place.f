      subroutine place(lnew,lterm,i,imolty,ibox,index,wplace)

c place
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

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'cbmc.inc'
      include 'connect.inc'
      include 'fix.inc'
      include 'system.inc'
      include 'rosen.inc'

      logical lnew,lterm,ovrlap

      integer i,j,imolty,count,counta,iu,ju,ku,jtvib,start,iv,dir
     &     ,index,ivib,nchvib,ibend,ib,type,site,ip,ichoi,niplace
     &     ,iw,iufrom,it,jut2,jut3,jut4,jttor,ja,ibox,glist,iwalk
     &     ,iuprev,list,nchben_a, nchben_b,iuback2,iuone,max,iu2back
 
      parameter(max=10)

      double precision wplace,equil,kforce,bsum_try,mincb,delcb
     &     ,ux,uy,uz,r,vvib,bfactor,third,random,length,bs,rbf
     &     ,vvibtr,wei_vib,bendang,vangle,vbbtr,angle,vphi,phione
     &     ,thetac,rx,ry,rz,rsint,dist
     &     ,wei_bend,ang_trial,angles,vctor,vdha
     &     ,vbend,vtorsion,bsum,alpha,gamma,twopi,dum,phi,thetaone
     &     ,thetatwo,phitwo


      dimension r(nchbn_max),bfactor(nchbn_max)
     &     ,bendang(numax,numax),ang_trial(nchbn_max),dist(max)
     &     ,niplace(numax),vbend(nchmax),vtorsion(nchmax),phi(max)
     &     ,list(max),glist(max)

c     ***********************************************************
c     **  Places hydrogens after the growth of the backbone of **
c     **  a molecule for linear, branched or cylic molecules.  **
c     **           -- Uses CDCBMC to grow them --              **
c     ***********************************************************


c      write(2,*) 'START PLACE'


      nchvib = nchbna(imolty)
      nchben_a = nchvib
      nchben_b = nchbnb(imolty)
      third = 1.0d0 / 3.0d0
      twopi = dacos(-1.0d0) * 2.0d0
      wplace = 1.0d0
     
      ichoi = nchoi(imolty)

      do j = 1, nunit(imolty)
         niplace(j) = 0
      enddo

      do iw = 1, nplace
         do count = 1, pnum(iw)
            iu = iplace(iw,count)
            niplace(iu) = iw     
         enddo
      enddo
      
      do iw = 1, nplace
         vvibtr = 0.0d0
         wei_vib = 1.0d0
         iufrom = pfrom(iw)
         do count = 1, pnum(iw)

            iu = iplace(iw,count)

            if (invib(imolty,iu).gt.1) then
               write(2,*) 'iu,invib',iu,invib(imolty,iu)
               stop 'invib can no be larger than one for hydrogen'
            endif

c     --- determine bond lengths
            iv = 1
            
            ju = ijvib(imolty,iu,iv)

            if (iufrom.ne.ju) then
               write(2,*) 'iu,ju,iufrom',iu,ju,iufrom
               stop 'ju not equal to iufrom'
            endif

            jtvib = itvib(imolty,iu,iv)
            
            equil = brvib(jtvib)
            kforce = brvibk(jtvib)

            if (kforce.gt.0.001) then
c     --- we will use flexible bond lengths
               bsum_try = 0.0d0
               mincb = brvibmin(jtvib)**3
               delcb = brvibmax(jtvib)**3 - mincb
               
               if (.not. lnew) then
                  ux = rxu(i,ju) - rxu(i,iu)
                  uy = ryu(i,ju) - ryu(i,iu)
                  uz = rzu(i,ju) - rzu(i,iu)
                  r(1) = dsqrt(ux**2 + uy**2 + uz**2)
                  vvib = kforce * (r(1) - equil)**2
                  bfactor(1) = dexp(-beta*vvib)
                  bsum_try = bsum_try + bfactor(1)
                  start = 2
               else
                  start = 1
               endif
               
               do ivib = start, nchvib
                  r(ivib) = (mincb + random()*delcb)**third
                  vvib = kforce * ( r(ivib) - equil )**2
                  bfactor(ivib) = dexp(-beta*vvib)
                  bsum_try = bsum_try + bfactor(ivib)
               enddo
               
               wei_vib = wei_vib * bsum_try
               
               if (lnew) then
c     --- select one of the trial sites via bias
                  rbf = random()*bsum_try
                  bs = 0.0d0
                  do ivib = 1, nchvib
                     bs = bs + bfactor(ivib)
                     if (rbf .lt. bs ) then
                        length = r(ivib)
                        vvib = dlog(bfactor(ivib))/(-beta)
                        goto 5
                     endif
                  enddo
 5                continue
               else
c     --- select old conformation
                  length = r(1)
                  vvib = dlog(bfactor(1))/(-beta)
               endif
               vvibtr = vvibtr + vvib

            else
               
c     --- our bond lengths are fixed
               length = equil
               
            endif
            
            distij(iu,ju) = length
            distij(ju,iu) = length
         enddo
         if (lnew) then
            vnewt = vnewt + vvibtr
            vnewbvib  = vnewbvib  + vvibtr
         else
            voldt = voldt + vvibtr
            voldbvib = voldbvib + vvibtr
         endif

      enddo

      do iw = 1, nplace

         iufrom = pfrom(iw)
         iuprev = pprev(iw)

c     --- first, set up cone
         dist(2) = distij(iuprev,iufrom)

         rx = xvec(iuprev,iufrom) / dist(2) 
         ry = yvec(iuprev,iufrom) / dist(2)
         rz = zvec(iuprev,iufrom) / dist(2)
         
         call cone(1,rx,ry,rz,dum,dum,dum,dum,dum)
               
c     --- now that we set up cone, we must determine other beads grown
c     --- from iufrom
         counta = 2

         do iv = 1, invib(imolty,iufrom)
            ku = ijvib(imolty,iufrom,iv)

c     --- make sure that ku is not equal to a site we are growing or iuprev
            if (ku.eq.iuprev) goto 95
            do count = 1, pnum(iw)
               iu = iplace(iw,count)
               if (iu.eq.ku) goto 95
            enddo
            
c     --- we must determine the angle associated with this one
            counta = counta + 1
            dist(counta) = distij(iufrom,ku)
            ux = xvec(iufrom,ku) / dist(counta)
            uy = yvec(iufrom,ku) / dist(counta)
            uz = zvec(iufrom,ku) / dist(counta)
           
c     --- determine angle with iuprev
            thetac = -(ux*rx + uy*ry + uz*rz)

            bendang(ku,iuprev) = dacos(thetac)
            
            alpha = bendang(ku,iuprev)

            call cone(3,dum,dum,dum,alpha,gamma,ux,uy,uz)

            phi(counta) = gamma
            list(counta) = ku

 95         continue
         enddo

       
         do ip = 1, ichoi

            wei_bend = 1.0d0
            vbbtr = 0.0d0
            vdha = 0.0d0
            
            do count = 1, pnum(iw)
               
               iu = iplace(iw,count)
      
               ju = ijvib(imolty,iu,1)
            
               do ib = 1, inben(imolty,iu)
                              
                  ku = ijben3(imolty,iu,ib)
               
                  if (ku.eq.pprev(iw)) then

                     if (ju.ne.ijben2(imolty,iu,ib)) then
                        write(2,*) 'ju,ijben2',ju,ijben2(imolty,iu,ib)
                        stop 'ju not equal to ijben2 in place'
                     endif

                     type = itben(imolty,iu,ib)
                     equil = brben(type)
                     kforce = brbenk(type)
               
c     --- initialize bsum_try
                     bsum_try = 0
                  
                     if (.not.lnew.and.ip.eq.1) then
c     --- first ibend is the old conformation
c     --- compute vector from iufrom to iugrow
                        ux = rxu(i,iu) - rxu(i,ju)
                        uy = ryu(i,iu) - ryu(i,ju)
                        uz = rzu(i,iu) - rzu(i,ju)
                        dist(1) = dsqrt(ux**2 + uy**2 + uz**2)
                        
c     --- dot product divided by lengths gives cos(angle)
                        thetac = -( ux*rx + uy*ry 
     &                       + uz*rz ) 
     &                       / (dist(1))
                        angle = dacos(thetac)
                        
c     --- compute the energy of this angle
                        vangle = kforce * (angle - equil)**2
                        ang_trial(1) = angle
                        bfactor(1) = dexp( -beta*vangle )
                        bsum_try = bsum_try + bfactor(1)
                        
c     --- skip first ibend in next loop
                        start = 2
                        
                     else
c     --- new conformation start at 1
                        start = 1
                     endif
                     
c     --- compute trial angles and energies
                     do ibend = start,nchben_a
c     --- choose the angle uniformly on sin(theta)
                        rsint = 2.0d0*random() - 1.0d0
                        angle = dacos(rsint)
                        ang_trial(ibend) = angle
                        
c     --- calculate the bond angle energy
                        vangle = kforce * (angle - equil)**2
                        bfactor(ibend) = dexp(-beta*vangle)
                        bsum_try = bsum_try + bfactor(ibend)
                     enddo
                  
                     if ( lnew.or.ip.ne.1 ) then
c     --- select one of the trial sites via bias
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
 10                     continue
                     else
c     --- select the old conformation
                        angle = ang_trial(1)
                        vangle = dlog(bfactor(1))/(-beta)
                     endif
                  
c     --- propagate the rosenbluth weight
                     wei_bend = wei_bend * bsum_try/dble(nchben_a)
                     
                     bendang(iu,ku) = angle

                     vbbtr = vbbtr + vangle
                  endif
               
               enddo

            enddo
            
c     --- now we must determine the second posible position for our sites
                        
            do count = 1, pnum(iw)

               iu = iplace(iw,count)
               
c     --- initialize bsum_try
               bsum_try = 0.0d0

               if ( .not. lnew .and.ip.eq.1) then
                  ux = rxu(i,iu) - rxu(i,iufrom)
                  uy = ryu(i,iu) - ryu(i,iufrom)
                  uz = rzu(i,iu) - rzu(i,iufrom)
                  
                  dist(1) = dsqrt(ux**2 + uy**2 + uz**2)

                  ux = ux / dist(1)
                  uy = uy / dist(1)
                  uz = uz / dist(1)
                  alpha = bendang(iu,iuprev)
                  thetatwo = alpha

                  call cone(3,dum,dum,dum,alpha,gamma,ux,uy,uz)

                  phitwo = gamma

                  vphi = 0.0d0

                  do j = 3, counta + count - 1
                     ku = list(j)

                     phione = phi(j)
                     thetaone = bendang(ku,iuprev)
                     

                     call coneangle(thetaone,phione,thetatwo,phitwo
     &                    ,angle)
                     
                     do ib = 1, inben(imolty,iu)
                        iu2back = ijben3(imolty,iu,ib)
                        
                        if (iu2back.eq.ku) then
                           type = itben(imolty,iu,ib)

                           vphi = vphi + brbenk(type)
     &                          * (angle - brben(type))**2

                        endif
                     enddo
                  enddo
                  
                  ang_trial(1) = phitwo
                  bfactor(1) = dexp( -beta * vphi )
                  bsum_try = bsum_try + bfactor(1)

                  start = 2
               else
                  start = 1
                  thetatwo = bendang(iu,iuprev)
               endif
                  
               do ibend = start, nchben_b
                  phitwo = random() * twopi
                  vphi = 0.0d0
                  do j = 3, counta + count - 1
                     ku = list(j)
                     phione = phi(j)
                     thetaone = bendang(ku,iuprev)
                                          
                     call coneangle(thetaone,phione,thetatwo,phitwo
     &                    ,angle)

                     do ib = 1, inben(imolty,iu)
                        iu2back = ijben3(imolty,iu,ib)
                        if (iu2back.eq.ku) then
                           type = itben(imolty,iu,ib)
                           vphi = vphi + brbenk(type) 
     &                          * (angle - brben(type))**2
                        endif
                     enddo
                  enddo
c           --- store the boltzmann factors and phi
                  bfactor(ibend) = dexp(-beta*vphi)
                  ang_trial(ibend) = phitwo
                  bsum_try = bsum_try + bfactor(ibend)
               enddo

               if ( lnew.or.ip.ne.1 ) then
                  rbf = random() * bsum_try
                  bs = 0.0d0
                  do ibend = 1, nchben_b
                     bs = bs + bfactor(ibend)
                     if (rbf .lt. bs) then
                        phitwo = ang_trial(ibend)
                        vphi = dlog( bfactor(ibend)) / (-beta)
                        goto 15
                     endif
                  enddo
               else
                  phitwo = ang_trial(1)
                  vphi = dlog( bfactor(1) )/ (-beta)
               endif

 15            continue

               wei_bend = wei_bend * (bsum_try/dble(nchben_b))

               vbbtr = vbbtr + vphi

               phi(counta+count) = phitwo
               list(counta+count) = iu

c     --- determine vectors associated with this
               
               call cone(2,dum,dum,dum,thetatwo,phitwo,ux,uy,uz)
               
               xvec(iufrom,iu) = ux * distij(iufrom,iu)
               yvec(iufrom,iu) = uy * distij(iufrom,iu)
               zvec(iufrom,iu) = uz * distij(iufrom,iu)

               xvec(iu,iufrom) = -ux 
               yvec(iu,iufrom) = -uy
               zvec(iu,iufrom) = -uz
                        
c     --- now to calculate all torsions

               do it = 1, intor(imolty,iu)

                  jut2 = ijtor2(imolty,iu,it)
                  jut3 = ijtor3(imolty,iu,it)
                  jut4 = ijtor4(imolty,iu,it)

                  if (niplace(jut4).lt.niplace(iu)) then
                        
                     if (.not. lexist(jut4)) then
                        write(2,*) 'jut4,jut3,jut2,iu',
     &                       jut4,jut3,jut2,iu
                        stop 'trouble jut4 in place'
                     endif
                     
                     jttor = ittor(imolty,iu,it)

                     call calctor(iu,jut2,jut3,jut4,jttor,vctor)
                     vdha = vdha + vctor
                  endif
               enddo
            enddo
                    
            do count = 1, pnum(iw)
               iu = iplace(iw,count)
               if (lnew) then
                  rxp(count,ip) = xvec(iufrom,iu) + rxnew(iufrom)
                  ryp(count,ip) = yvec(iufrom,iu) + rynew(iufrom)
                  rzp(count,ip) = zvec(iufrom,iu) + rznew(iufrom)
               else
                  rxp(count,ip) = xvec(iufrom,iu) + rxu(i,iufrom)
                  ryp(count,ip) = yvec(iufrom,iu) + ryu(i,iufrom)
                  rzp(count,ip) = zvec(iufrom,iu) + rzu(i,iufrom)
               endif
               glist(count) = iu
            enddo
            vtorsion(ip) = vdha
            vbend(ip) = vbbtr
            bsum_tor(ip) = dexp(-beta * vdha) * wei_bend
         enddo

c     --- now we calculate the nonbonded interactions
         call boltz( lnew, .false., ovrlap,i,i,imolty,ibox,ichoi
     &        ,iufrom,pnum(iw),glist)

         if (ovrlap) then
            lterm = .true.
            return
         endif

         bsum = 0.0d0

         do ip = 1, ichoi
            if (.not. lovr(ip)) then
               bsum = bsum + bfac(ip) * bsum_tor(ip)
            endif
         enddo

c     --- update new rosenbluth weight + vibrations
         wplace = wplace * bsum * wei_vib
         
         if (wplace .lt. softlog) then
            lterm = .true.
            return
         endif
         
         if (lnew) then
            rbf = bsum * random()
            bs = 0.0d0
            
            do ip = 1, ichoi
               if (.not. lovr(ip)) then
                  bs = bs + bfac(ip) * bsum_tor(ip)
                  if (rbf .lt. bs) then
                     iwalk = ip
                     goto 20
                  endif
               endif
            enddo

            stop 'BIG TIME SCREWUP IN PLACE'

         endif
 20      continue

         if (lnew) then

            vnewt = vnewt + vbend(iwalk) + vtorsion(iwalk)
     &           + vtrintra(iwalk)
            vnewbb    = vnewbb    + vbend(iwalk)
            vnewtg    = vnewtg    + vtorsion(iwalk)
            vnewext   = vnewext   + vtrext(iwalk)
            vnewintra = vnewintra + vtrintra(iwalk)
            vnewinter = vnewinter + vtrinter(iwalk)
            vnewelect = vnewelect + vtrelect(iwalk)
            vnewewald = vnewewald + vtrewald(iwalk)
         else
            voldt = voldt + vbend(1) + vtorsion(1)
     &           + vtrintra(1)
            voldbb    = voldbb    + vbend(1)
            voldtg    = voldtg    + vtorsion(1)
            voldext   = voldext   + vtrext(1)
            voldintra = voldintra + vtrintra(1)
            voldinter = voldinter + vtrinter(1)
            voldelect = voldelect + vtrelect(1)

            voldewald = voldewald + vtrewald(1)
         endif

         do count = 1, pnum(iw)
            iu = iplace(iw,count)
         
            if (lnew) then
               rxnew(iu) = rxp(count,iwalk)
               rynew(iu) = ryp(count,iwalk)
               rznew(iu) = rzp(count,iwalk)
            endif

            lexist(iu) = .true.

            if (lnew) then
               xvec(iu,iufrom) = rxnew(iufrom) - rxnew(iu)
               yvec(iu,iufrom) = rynew(iufrom) - rynew(iu)
               zvec(iu,iufrom) = rznew(iufrom) - rznew(iu)
            else
               xvec(iu,iufrom) = rxu(i,iufrom) - rxu(i,iu)
               yvec(iu,iufrom) = ryu(i,iufrom) - ryu(i,iu)
               zvec(iu,iufrom) = rzu(i,iufrom) - rzu(i,iu)
            endif
            distij(iu,iufrom) = dsqrt( xvec(iu,iufrom)**2
     +           + yvec(iu,iufrom)**2 + zvec(iu,iufrom)**2 )

            xvec(iufrom,iu) = - xvec(iu,iufrom)
            yvec(iufrom,iu) = - yvec(iu,iufrom)
            zvec(iufrom,iu) = - zvec(iu,iufrom)
            distij(iufrom,iu) = distij(iu,iufrom)
         enddo
      enddo

c      write(2,*) 'END PLACE'
      
      return
      end




