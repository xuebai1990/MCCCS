      subroutine safecbmc(iinit,lnew,i,iw,igrow,imolty,count
     &     ,ux,uy,uz,vphi,vtor,wei_bv,lterm,movetype)

c safecbmc
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

c     **************************************************************
c     **   Finshes the last two steps for Fixed Endpoint CBMC     **
c     **************************************************************
c     **     Originally completed by Collin Wick on 1-1-2000      **
c     **************************************************************
c     **     --- SEE safeschedule.f FOR MORE INFORMATION ---      **
c     **************************************************************

c     --- lshit is used for diagnistics ---

c     iinit = 1  initial setup for two beads to go
c     iinit = 2  calculates closing probabilities
c     iinit = 3  does final crankshaft move


      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'connect.inc'
      include 'fix.inc'
      include 'rosen.inc'
      include 'cbmc.inc'

      logical::lnew,lshit,lterm,ldo,lreturn

      integer::igrow,imolty,count,counta,j,ja,ivib,iufrom,iuprev
     +     ,iinit,iu,ju,ku,i,iv,juvib,jtvib,type,iu2,ib,iw,ntogrow
     +     ,itor,ip,ichoi,ichtor,countb,bin,max,nu,iu1,dir,diracc
     +     ,start,nchben_a,nchben_b,ibend,iopen,last,iclose,nchvib

      integer::jttor,it,jut2,jut3,jut4,movetype,lu,k,opencount

      real(8)::vdha,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2,dot
     +     ,daa1,da1a2,phicrank,bf_tor,vtorsion,vbend,rbf
     +     ,vtorso,ran_tor,bs,ang_bend,bfactor,bsum_bend,wei_bv
     +     ,bsum_try,third,vibtr,lengthc

c     *** to conserve memory, max is the maximum number of endpoints
c     *** possible in place of numax
      parameter(max=10)

      real(8)::flength,x,y,z,equil,kforce,length
     +     ,vvib,equilb,kforceb,ux,uy,uz,hdist,lengtha,lengthb
     +     ,vtor,vphi,thetac,angle,equila,kforcea,ovphi
     +     ,alpha,twopi,phidisp,dum,random,rxt,ryt,rzt
     +     ,phiacc,rxa,rya,rza,angles,bangles,vctor
     +     ,r,mincb,delcb,vkforce,vequil,vvibration,ovvib

      dimension flength(numax,numax),alpha(max,numax)
     +     ,equilb(numax,numax),kforceb(numax,numax),equila(max)
     +     ,rxa(max,max),rya(max,max),rza(max,max)
     +     ,rxt(max),ryt(max),rzt(max),vbend(2*nchtor_max)
     +     ,phicrank(2*nchtor_max,max),phiacc(max)
     +     ,vtorsion(2*nchtor_max),bf_tor(2*nchtor_max)
     +     ,dir(2*nchtor_max,max),diracc(max),kforcea(max)
     +     ,bfactor(nchbn_max),ang_bend(nchbn_max),angles(3)
     +     ,bangles(max,3),iopen(2),r(nchbn_max)
     +     ,vkforce(numax,numax),vequil(numax,numax)
     +     ,vvibration(2*nchtor_max)

      save kforceb,equilb,flength,vequil,vkforce

c     ---------------------------------------------------------------
    
c      write(iou,*) 'START SAFECMBC ',iinit,'!'
c      print*,'iinit',iinit,'iw',iw,'igrow',igrow,'count',count
      vphi = 0.0d0
      ovphi = 0.0d0
      wei_bv = 1.0d0
      
      lreturn = .false.
      
 400  continue
      
      if (iinit.eq.1.or.(.not.lreturn.and.lcrank)) then
         
         third = 1.0d0 / 3.0d0
         nchvib = nchbna(imolty)
         ntogrow = grownum(iw)
         
         if (lcrank) then
            ntogrow = 1
         else
            ntogrow = grownum(iw)
         endif
         
c     *** lets first determine our bond distances ***

c     --- find vibrations for iu - ibef
         do 155 count = 1, ntogrow
            
            ja = 0
            
            if (lcrank) then
               wbefnum = 1
               iwbef(1) = growfrom(iw)
               iu = growfrom(iw)
            else
               iu = growlist(iw,count)
            endif
            do j = 1, wbefnum
               if (iu.eq.iwbef(j)) then
c     --- ja defines iwbef in counter
                  ja = j
               endif
            enddo
            
c     --- we do not close with this bead
            if (ja.eq.0) goto 155
            
            do iv = 1, invib(imolty,iu)
               juvib = ijvib(imolty,iu,iv)
               do counta = 1, befnum(ja)
                  ju = ibef(ja,counta)
                  if (juvib.eq.ju) then
                     jtvib = itvib(imolty,iu,iv)
                     equil = brvib(jtvib)
                     kforce = brvibk(jtvib)
                     if (kforce.gt.0.1d0) then
c     we will use flexible bond lengths
                        bsum_try = 0.0d0
                        mincb = brvibmin(jtvib)**3
                        delcb = brvibmax(jtvib)**3 - mincb
                        
                        if (.not.lnew) then
                           x = rxu(i,iu) - rxu(i,ju)
                           y = ryu(i,iu) - ryu(i,ju)
                           z = rzu(i,iu) - rzu(i,ju)
                           r(1) = dsqrt(x**2 + y**2 + z**2)
                           vvib = kforce * ( r(1) - equil )**2
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
                        
                        if (lnew) then
c     --- select one of the trial sites via vias
                           rbf = random()*bsum_try
                           bs = 0.0d0
                           do ivib = 1, nchvib
                              bs = bs + bfactor(ivib)
                              if (rbf .lt. bs ) then
                                 flength(iu,ju) = r(ivib)
                                 vvib = dlog(bfactor(ivib))/(-beta)
                                 goto 4
                              endif
                           enddo
 4                         continue
                        else
c     --- select old conformation
                           flength(iu,ju) = r(1)
                           vvib = dlog(bfactor(1))/(-beta)
                        endif
                        
c     --- add up vibrational energy
                        vphi = vphi + vvib
                        
c     --- propogate rosenbluth weight
                        wei_bv = wei_bv * bsum_try/dble(nchvib)
                        
                     else
                        if (lnew) then
c     --- compute new bond length
                           call bondlength(jtvib,equil,kforce,beta,
     $                          length,vvib)
                           flength(iu,ju) = length

                        else
c     --- compute old bond length
                           x = rxu(i,ju) - rxu(i,iu)
                           y = ryu(i,ju) - ryu(i,iu)
                           z = rzu(i,ju) - rzu(i,iu)
                           flength(iu,ju) = dsqrt(x**2+y**2+z**2)
                        endif
                     endif 
                  endif
               enddo
            enddo
            
c     --- find vibrations for ibef - iend
            do counta = 1, befnum(ja)
               ju = ibef(ja,counta)
               do iv = 1, invib(imolty,ju)
                  juvib = ijvib(imolty,ju,iv)
                  do j = 1, fcount(ju)
                     ku = fclose(ju,j)
                     if (juvib.eq.ku) then
                        jtvib = itvib(imolty,ju,iv)
                        equil = brvib(jtvib)
                        kforce = brvibk(jtvib)
                        if (kforce.gt.0.1d0) then
c     --- we have flexible bond lengths
                              
                           bsum_try = 0.0d0
                           mincb = brvibmin(jtvib)**3
                           delcb = brvibmax(jtvib)**3 - mincb
                           
                           if (j.gt.1) then
                              vequil(ju,ku) = equil
                              vkforce(ju,ku) = kforce
                              goto 112
                           endif
                              
                           if (.not.lnew) then
                              x = rxu(i,ku) - rxu(i,ju)
                              y = ryu(i,ku) - ryu(i,ju)
                              z = rzu(i,ku) - rzu(i,ju)
                              r(1) = dsqrt(x**2 + y**2 + z**2)
                              vvib = kforce * ( r(1) - equil )**2
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
                           
                           if (lnew) then
c     --- select one of the trial sites via vias
                              rbf = random()*bsum_try
                              bs = 0.0d0
                              do ivib = 1, nchvib
                                 bs = bs + bfactor(ivib)
                                 if (rbf .lt. bs ) then
                                    flength(ju,ku) = r(ivib)
                                    vvib = dlog(bfactor(ivib)) /(-beta)
                                    goto 6
                                 endif
                              enddo
 6                            continue
                           else
c     --- select old conformation
                              flength(ju,ku) = r(1)
                              vvib = dlog(bfactor(1))/(-beta)
                              
                           endif
                           
c     --- add up vibrational energy
                           vphi = vphi + vvib
                           
c     --- propogate rosenbluth weight
                           wei_bv = wei_bv * bsum_try/dble(nchvib)
                        
 112                       continue
                           
                        else
                           
                           if (lnew) then
c     --- compute new bond length
                              call bondlength(jtvib,equil,kforce,beta
     &                             ,length,vvib)
                              flength(ju,ku) = length
                           else
c     --- compute old bond length
                              x = rxu(i,ku) - rxu(i,ju)
                              y = ryu(i,ku) - ryu(i,ju)
                              z = rzu(i,ku) - rzu(i,ju)
                              
                              flength(ju,ku) = dsqrt(x**2 + y**2 + z**2)
                              
                           endif
                        endif
                     endif
                  enddo
               enddo
            enddo
            
c        --- determine angles for iwbef-ibef-iend

            if (lcrank) then
               do ib = 1, inben(imolty,iu)
                  iu2 = ijben3(imolty,iu,ib)
                  type = itben(imolty,iu,ib)
                  do counta = 1, grownum(iw)
                     ju = growlist(iw,counta)
                     if (fcount(ju).ne.0) then
                        do j = 1, fcount(ju)
                           ku = fclose(ju,j)
                           if (ku.eq.iu2) then
                              equilb(iu,ku) = brben(type)
                              kforceb(iu,ku) = brbenk(type)
                           endif
                        enddo
                     endif
                  enddo
               enddo
            else
               do ib = 1, inben(imolty,iu)
                  iu2 = ijben3(imolty,iu,ib)
                  type = itben(imolty,iu,ib)
                  do j = 1, fcount(iu)
                     ju = fclose(iu,j)
                     if (ju.eq.iu2) then
                        equilb(iu,ju) = brben(type)
                        kforceb(iu,ju) = brbenk(type)
                     endif
                  enddo
               enddo
            endif
 155     continue
      
         if (lcrank) then
c     --- we need to calculate the new angle
            iufrom = growfrom(iw)
            do count = 1, grownum(iw)
               iu = growlist(iw,count)
               if (fcount(iu).ne.0) then
                  j = 1
                  ju = fclose(iu,j)
                  lengtha = flength(iufrom,iu)
                  lengthb = flength(iu,ju)
                  
                  x = rxu(i,ju) - rxnew(iufrom)
                  y = ryu(i,ju) - rynew(iufrom)
                  z = rzu(i,ju) - rznew(iufrom)
                  
                  hdist = dsqrt( x**2 + y**2 + z**2 )
c     --- use law of cosines to calculate bond angle
                  thetac = (lengtha**2 + lengthb**2 - hdist**2) / (2.0d0
     $                 * lengtha * lengthb)
                  
c     --- check to make sure this will give a number
                  
                  if (abs(thetac).gt.1.0d0) then
                     vtor = 0
                     vphi = 0
                     lterm = .true.
                     return
                  endif
                  angle = dacos(thetac)
                  
                  ovphi = ovphi + kforceb(iufrom,ju) * (angle -
     $                 equilb(iufrom,ju))**2
             
c                  write(iou,*) iufrom,iu,ju,kforceb(iufrom,ju)*(angle
c     &                 - equilb(iufrom,ju))**2
     
                  wei_bv = wei_bv * dexp( - beta * ovphi )
               endif
            enddo
            vibtr= vphi
            
            lreturn = .true.
            goto 400
         endif

c     ********************************************************************
      elseif (iinit.eq.2.or.iinit.eq.4) then
c     --- lets determine closing energy for this bead alone
         lshit = .false.

c         if (iinit.eq.4) then
c            lshit = .true.
c         endif

         vtor = 1.0d0

c     --- this case iw and count is sent in from rosenbluth
         
         iu = growlist(iw,count)

c     --- determine iwbef count
         do j = 1, wbefnum
            if (iu.eq.iwbef(j)) then
               ja = j
               goto 100
            endif
         enddo
 100     continue

         do j = 1, fcount(iu)
            ku = fclose(iu,j)

            if (movetype.eq.2.and.lnew) then
               x = rxnew(ku) - ux
               y = rynew(ku) - uy
               z = rznew(ku) - uz
            else
               x = rxu(i,ku) - ux
               y = ryu(i,ku) - uy
               z = rzu(i,ku) - uz
            endif
               
            hdist = dsqrt( x**2 + y**2 + z**2 )

            if (j.gt.1) then
c     --- we will use a phoney probability since we don't know our bond
c     --- distances yet
               bin = anint( hdist * 10.0d0 )
               vtor = vtor * probf(iu,ku,bin)
            else
c     --- we can calculate an angle here
               do counta = 1, befnum(ja)
                  ju = ibef(ja,counta)
                  
                  do countb = 1, fcount(ju)
                     nu = fclose(ju,countb)
                     if (nu.eq.ku) goto 150
                  enddo
                  goto 175
                  
 150              continue
                  
                  lengtha = flength(iu,ju)
                  lengthb = flength(ju,ku)
                  
c     --- use law of cosines to calculate bond angle

                  thetac = (lengtha**2 + lengthb**2 - hdist**2) / (2
     $                 .0d0 * lengtha * lengthb)
                  
c     --- check to make sure this will give a number

                  if (abs(thetac).gt.1.0d0) then
                     vtor = 0
                     vphi = 0
                     return
                  endif
                  angle = dacos( thetac )
                  
                  vphi = vphi + kforceb(iu,ku) * ( angle - equilb(iu,ku)
     $                 )**2
                  
 175              continue

               enddo
            endif
         
c     --- determine torsion interaction with growpast if it exists

            if (pastnum(ku).ne.0) then
               
               do counta = 1, pastnum(ku)
                  nu = ipast(ku,counta)
                  if (.not.lplace(imolty,nu)) then
                     if (movetype.eq.2.and.lnew) then
                        x = rxnew(nu) - ux
                        y = rynew(nu) - uy
                        z = rznew(nu) - uz
                     else
                        x = rxu(i,nu) - ux
                        y = ryu(i,nu) - uy
                        z = rzu(i,nu) - uz
                     endif
                     hdist = dsqrt( x**2 + y**2 + z**2 )
                     bin = anint( hdist * 10.0d0 )
                     vtor = vtor * probf(iu,nu,bin)
                     
                     if (nextnum(nu).ne.0) then
                        do k = 1, nextnum(nu)
                           lu = inext(nu,k)
                           if (.not.lplace(imolty,lu)) then
                              if (movetype.eq.2.and.lnew) then
                                 x = rxnew(lu) - ux
                                 y = rynew(lu) - uy
                                 z = rznew(lu) - uz 
                              else
                                 x = rxu(i,lu) - ux
                                 y = ryu(i,lu) - uy
                                 z = rzu(i,lu) - uz
                              endif
                              hdist = dsqrt( x**2 + y**2 + z**2 )
                              bin = anint( hdist * 10.0d0 )
                              vtor = vtor * probf(iu,lu,bin)
                           endif
                        enddo
                     endif
                  endif
               enddo
            endif
         enddo

c *********************************************************************
      else
 
c     --- CRANKSHAFT MOVE
         ntogrow = grownum(iw)
         iufrom = growfrom(iw)
         iuprev = growprev(iw)
         iopen(1) = 0
         iopen(2) = 0
         iclose = 0
         ovvib = 0
c     --- determine bond distances for iuprev - iufrom

         opencount = 0

         do count = 1, ntogrow
            iu = growlist(iw,count)
            if (fcount(iu).ne.0) then

               do j = 1, fcount(iu)
                  ju = fclose(iu,j)
                  if (pastnum(ju).ne.0) then
                     do counta = 1, pastnum(ju)
                        ku = ipast(ju,counta)
                        
c     --- determine angles for iu - iend - ipast if they exist
                     
                        do ib = 1, inben(imolty,iu)
                           iu2 = ijben3(imolty,iu,ib)
                           type = itben(imolty,iu,ib)
                           if (iu2.eq.ku) then
                              equilb(iu,ku) = brben(type)
                              kforceb(iu,ku) = brbenk(type)
                              goto 125
                           endif
                        enddo
                        write(iou,*) 'iu,ju,ku',iu,ju,ku
                        call cleanup('no bond angle for these')
 125                    continue
                        
                        if (pastnum(ju).gt.1) then
                           do countb = counta+1, pastnum(ju)
                              lu = ipast(ju,countb)
                              
                              do ib = 1, inben(imolty,ku)
                                 iu2 = ijben3(imolty,ku,ib)
                                 type = itben(imolty,ku,ib)
                                 if (iu2.eq.lu) then
                                    equilb(ku,lu) = brben(type)
                                    kforceb(ku,lu) = brbenk(type)
                                 endif
                              enddo
                           enddo
                        endif

                     enddo
                  endif
                
                  if (j.gt.1) goto 126
  
c     --- calculate distances from iufrom to iend
c     *** count is iu count, and ju is iend bead
                  if (lnew) then
                     if (movetype.eq.2) then
                        xvec(iufrom,ju) = rxnew(ju) - rxnew(iufrom)
                        yvec(iufrom,ju) = rynew(ju) - rynew(iufrom)
                        zvec(iufrom,ju) = rznew(ju) - rznew(iufrom)
                     else
                        xvec(iufrom,ju) = rxu(i,ju) - rxnew(iufrom)
                        yvec(iufrom,ju) = ryu(i,ju) - rynew(iufrom)
                        zvec(iufrom,ju) = rzu(i,ju) - rznew(iufrom)
                     endif
                  else
                     xvec(iufrom,ju) = rxu(i,ju) - rxu(i,iufrom)
                     yvec(iufrom,ju) = ryu(i,ju) - ryu(i,iufrom)
                     zvec(iufrom,ju) = rzu(i,ju) - rzu(i,iufrom)
                  endif
                  hdist = dsqrt( xvec(iufrom,ju)**2 + yvec(iufrom,ju)**2
     $                 + zvec(iufrom,ju)**2 )
                  
c     --- normalize these distances to one for cone
                  xvec(iufrom,ju) = xvec(iufrom,ju) / hdist
                  yvec(iufrom,ju) = yvec(iufrom,ju) / hdist
                  zvec(iufrom,ju) = zvec(iufrom,ju) / hdist
                  
c     --- calculate alpha
c     *** count is iu count, and ju is iend bead
                  lengtha = flength(iufrom,iu) 
                  lengthb = flength(iu,ju)
                  
                  thetac = (lengtha**2 + hdist**2 - lengthb**2) / (2.0d0
     $                 * lengtha * hdist)
                                    
                  if (abs(thetac).gt.1.0d0) then
                     vphi = 0
                     vtor = 0
                     lterm = .true.
                     return
                  endif
                  alpha(count,ju) = dacos(-1.0d0) - dacos(thetac)
 126              continue
               enddo
            endif

c     --- determine angle for    - iu - iufrom - 
            do ib = 1, inben(imolty,iu)
               iu1 = ijben2(imolty,iu,ib)
               type = itben(imolty,iu,ib)
               
               if (iu1.eq.iufrom) then
                  iu2 = ijben3(imolty,iu,ib)
                  if (iu2.eq.iuprev) then
                     equila(count) = brben(type)
                     kforcea(count) = brbenk(type)
                  else
                     equilb(iu,iu2) = brben(type)
                     kforceb(iu,iu2) = brbenk(type)
                  endif
               endif
            enddo

            if (fcount(iu).gt.1) then
               do j = 1, fcount(iu) - 1
                  ju = fclose(iu,j)
                  do counta = j+1, fcount(iu)
                     ku = fclose(iu,counta)
                     do ib = 1, inben(imolty,ku)
                        iu2 = ijben3(imolty,ku,ib)
                        type = itben(imolty,ku,ib)
                     
                        if (iu2.eq.ju) then
                           equilb(ju,ku) = brben(type)
                           kforceb(ju,ku) = brbenk(type)
                        endif
                     enddo
                  enddo
               enddo
            endif
               


            


c ---------------------------------------------------------------
c     begin part that is only for closures with an open bead

            if (fcount(iu).eq.0) then
              
               opencount = opencount + 1

c     --- we first have to find the bond length here
               do iv = 1, invib(imolty,iufrom)
                  juvib = ijvib(imolty,iufrom,iv)
                  if (juvib.eq.iu)  then
                     jtvib = itvib(imolty,iufrom,iv)
                     equil = brvib(jtvib)
                     kforce = brvibk(jtvib)
                     if (kforce.gt.0.1d0) then
c     we will use flexible bond lengths
                        bsum_try = 0.0d0
                        mincb = brvibmin(jtvib)**3
                        delcb = brvibmax(jtvib)**3 - mincb

                        if (.not.lnew) then
                           x = rxu(i,iufrom) - rxu(i,iu)
                           y = ryu(i,iufrom) - ryu(i,iu)
                           z = rzu(i,iufrom) - rzu(i,iu)
                           r(1) = dsqrt(x**2 + y**2 + z**2)
                           vvib = kforce * ( r(1) - equil )**2
                           bfactor(1) = dexp(-beta*vvib)
                           bsum_try = bsum_try + bfactor(1)
                           start = 2
                        else
                           start = 1
                        endif
                      
                        do ivib = start, nchvib
                           
                           r(ivib) = (mincb 
     &                          + random()*delcb)**third
                           vvib = kforce * 
     &                          ( r(ivib) - equil )**2
                           bfactor(ivib) = dexp(-beta*vvib)
                           bsum_try = bsum_try + bfactor(ivib)
                        enddo
                         
                        if (lnew) then
c     --- select one of the trial sites via vias
                           rbf = random()*bsum_try
                           bs = 0.0d0
                           do ivib = 1, nchvib
                              bs = bs + bfactor(ivib)
                              if (rbf .lt. bs ) then
                                 flength(iufrom,iu) = r(ivib)
                                 vvib = dlog(bfactor(ivib))/(-beta)
                                 goto 61
                              endif
                           enddo
 61                        continue
                        else
c     --- select old conformation
                           flength(iufrom,iu) = r(1)
                           vvib = dlog(bfactor(1))/(-beta)
                        endif
                        
c     --- add up vibrational energy
                        ovvib = ovvib + vvib
                                                
c     --- propogate rosenbluth weight
                        wei_bv = wei_bv * bsum_try/dble(nchvib)
                        
c     *************************************
                     else
                        
                        if (lnew) then
c     --- compute new bond length
                           call bondlength(jtvib,equil,kforce,beta,
     &                          length,vvib)
                           flength(iufrom,iu) = length
                        else
c     --- compute old bond length
                           x = rxu(i,iufrom) - rxu(i,iu)
                           y = ryu(i,iufrom) - ryu(i,iu)
                           z = rzu(i,iufrom) - rzu(i,iu)
                           flength(iufrom,iu) = dsqrt(x**2+y**2+z**2)
                           
                        endif
                     endif
                  endif
               enddo
                  
c     --- find bond angles

               bsum_bend = 0
               if (.not. lnew) then               
                 
                  lengtha = flength(iufrom,iu)
                  lengthb = distij(iufrom,iuprev)
                
                  thetac = -( (rxu(i,iu) - rxu(i,iufrom))
     &                 * xvec(iuprev,iufrom) 
     &                 + (ryu(i,iu) - ryu(i,iufrom))
     &                 * yvec(iuprev,iufrom)
     &                 + (rzu(i,iu) - rzu(i,iufrom))
     &                 * zvec(iuprev,iufrom)) 
     &                 / (lengtha*lengthb)
                  
                  angle = dacos(thetac)
                  vphi =  kforcea(count) * (angle-equila(count))**2
         
c                  write(iou,*) 'b',iu,iufrom,iuprev,vphi


                  ang_bend(1) = angle
                  bfactor(1) = dexp( -beta*vphi)
                  bsum_bend = bsum_bend + bfactor(1)
                  start = 2
               else
                  start = 1
               endif
               
               nchben_a = nchbna(imolty)

               do ibend = start, nchben_a
c     --- choose angle uniformly on sin(angle)
                  thetac = 2.0d0*random() - 1.0d0
                  angle = dacos(thetac)
                  ang_bend(ibend) = angle

c     --- find bend energy
                  vphi = kforcea(count) * (angle-equila(count))**2
                  bfactor(ibend) = dexp(-beta*vphi)
                  bsum_bend = bsum_bend + bfactor(ibend)
               enddo

               if (lnew) then
c     --- select one of the trial sites at random
                  rbf = random()*bsum_bend
                  bs = 0.0d0
                  do ibend = 1, nchben_a
                     bs = bs + bfactor(ibend)
                     if (rbf.lt.bs) then
                        bangles(count,1) = ang_bend(ibend)
                        ovphi = ovphi + dlog(bfactor(ibend))/(-beta)


c                        write(iou,*) 'c',iu,iufrom,iuprev
c     &                       ,dlog(bfactor(ibend))/(-beta)

                        goto 10
                     endif
                  enddo
 10               continue
               else
                  bangles(count,1) = ang_bend(1)
                  ovphi = ovphi + dlog(bfactor(1))/(-beta)
               endif

               wei_bv = wei_bv * bsum_bend/dble(nchben_a)

c     --- now we have to determine the angle with a crankshaft bead
c     --- for the opencount = 1, or the other free bead for opencount > 1


               bsum_bend = 0

               do counta = 1, ntogrow

                  if (counta.ne.count) then
                     ju = growlist(iw,counta)

c     --- check to see if the conditions just stated are true
                     if (opencount.gt.1) then
c     --- we only want the angle with a free bead
                        if (fcount(ju).ne.0) then
                           goto 25
                        endif
                     else
c     --- we want the angle with the closing bead
                        if (fcount(ju).eq.0) then
                           goto 25
                        elseif (iclose.ne.0) then
                           goto 25
                        else
                           iclose = counta
                        endif
                     endif
                     
                     if (.not. lnew) then
                        
                        lengthb = flength(iufrom,ju)
                        
                        thetac = ( (rxu(i,ju) - rxu(i,iufrom))
     &                       * (rxu(i,iu) - rxu(i,iufrom))
     &                       + (ryu(i,ju) - ryu(i,iufrom))
     &                       * (ryu(i,iu) - ryu(i,iufrom))
     &                       + (rzu(i,ju) - rzu(i,iufrom))
     &                       * (rzu(i,iu) - rzu(i,iufrom)))
     &                       / (lengtha * lengthb)

                        angle = dacos(thetac)
                        
                        vphi = kforceb(iu,ju) * (angle-equilb(iu,ju))**2
                        
                        
c                        write(iou,*) 'd',iu,iufrom,ju,vphi

                        ang_bend(1) = angle
                        bfactor(1) = dexp( -beta*vphi )
                        bsum_bend = bsum_bend + bfactor(1)

                        start = 2
                     else
                        start = 1
                     endif

                     nchben_b = nchbnb(imolty)
                     do ibend = start, nchben_b
c     --- choose angle uniformly on sin(angle)
                        thetac = 2.0d0*random() - 1.0d0
                        angle = dacos(thetac)
                        ang_bend(ibend) = angle
                        
c     --- find bend energy
                        vphi = kforceb(iu,ju) * (angle-equilb(iu,ju))**2
                        
                        bfactor(ibend) = dexp(-beta*vphi)
                        bsum_bend = bsum_bend + bfactor(ibend)
                     enddo

                     if (lnew) then
c     --- select one of the trial sites at random
                        rbf = random()*bsum_bend
                        bs = 0.0d0
                        do ibend = 1, nchben_b
                           bs = bs + bfactor(ibend)
                           if (rbf.lt.bs) then
                              bangles(count,2) = ang_bend(ibend)
                              ovphi = ovphi 
     &                             + dlog(bfactor(ibend))/(-beta)

c                              write(iou,*) 'd',iu,iufrom,ju
c     &                             ,dlog(bfactor(ibend))/(-beta)

                              goto 20
                           endif
                        enddo
 20                     continue
                     else
                        bangles(count,2) = ang_bend(1)
                        ovphi = ovphi + dlog(bfactor(1))/(-beta)
                        
                     endif
                     
                     wei_bv = wei_bv * bsum_bend/dble(nchben_b)
                   endif
 25               continue
               enddo
            endif
         enddo

c     end part for the case with an open closing bead


c -----------------------------------------------------------------          
c     --- loop over all choices
         ichoi = nchoi(imolty)
c     --- double ichtor to give extra help for this move
         ichtor = nchtor(imolty) * 2
         twopi = 2.0d0 * dacos(-1.0d0)
         
         do ip = 1, ichoi
            
            bsum_tor(ip) = 0
            do itor = 1, ichtor

               vvib = 0
               countb = 0
               lshit = .false.

               if (lnew.and.ip.eq.19
     &              .and.itor.eq.122) then
                  lshit = .true.
               endif

               if (.not.lnew
     &              .and.ip.eq.1.and.itor.eq.1) then
                  lshit = .true.
               endif

               vdha = 0
               vphi = 0
                              
               ldo = .false.
               start = 1
               last = ntogrow
 30            continue
               do count = start, last

                  dir(itor,count) = 0

                  iu = growlist(iw,count)

                  if (fcount(iu).gt.0) then
c     --- set up cone     

                     ju = fclose(iu,1)
                     call cone(1, xvec(iufrom,ju),yvec(iufrom,ju)
     &                    ,zvec(iufrom,ju),dum,dum,dum,dum,dum)

c     --- determine phidisp
                     if (.not.lnew.and.ip.eq.1.and
     &                    .itor.eq.1) then
c     --- give old unit vector for connection
                        xx(count) = rxu(i,iu) - rxu(i,iufrom)
                        yy(count) = ryu(i,iu) - ryu(i,iufrom)
                        zz(count) = rzu(i,iu) - rzu(i,iufrom)

                        ux = rxu(i,iu)
                        uy = ryu(i,iu)
                        uz = rzu(i,iu)

                     else
                        phidisp = twopi*random()
                           
                        call cone(2,dum,dum,dum,alpha(count,ju)
     &                       ,phidisp,x,y,z)
                        
                        phicrank(itor,count) = phidisp  
                        
                        xx(count) = x * flength(iufrom,iu)
                        yy(count) = y * flength(iufrom,iu)
                        zz(count) = z * flength(iufrom,iu)
     
                        if (lnew) then
                           ux = xx(count) + rxnew(iufrom)
                           uy = yy(count) + rynew(iufrom)
                           uz = zz(count) + rznew(iufrom)
                        else
                           ux = xx(count) + rxu(i,iufrom)
                           uy = yy(count) + ryu(i,iufrom)
                           uz = zz(count) + rzu(i,iufrom)
                        endif                      
                     endif

                     if (fcount(iu).gt.1) then
c     --- calculate distances to the other endpoints

                        do j = 2, fcount(iu)

                           ju = fclose(iu,j)

                           equil = vequil(iu,ju)
                           kforce = vkforce(iu,ju)


                           if (movetype.eq.2.and.lnew) then
                              x = rxnew(ju) - ux
                              y = rynew(ju) - uy
                              z = rznew(ju) - uz
                           else
                              x = rxu(i,ju) - ux
                              y = ryu(i,ju) - uy
                              z = rzu(i,ju) - uz
                           endif
                           
                           xvec(iu,ju) = x
                           yvec(iu,ju) = y
                           zvec(iu,ju) = z
                           
                           xvec(ju,iu) = -x
                           yvec(ju,iu) = -y
                           zvec(ju,iu) = -z
                           
                           length = dsqrt(x**2 + y**2 + z**2)
                     
                           flength(iu,ju) = length
                           
                           vvib = vvib + kforce * (length - equil)**2

c     --- we need to calculate the new angle

                           lengtha = flength(iufrom,iu)
                           lengthb = length
                           
                           if (lnew) then
                              if (movetype.eq.2) then
                                 x = rxnew(ju) - rxnew(iufrom)
                                 y = rynew(ju) - rynew(iufrom)
                                 z = rznew(ju) - rznew(iufrom)
                              else
                                 x = rxu(i,ju) - rxnew(iufrom)
                                 y = ryu(i,ju) - rynew(iufrom)
                                 z = rzu(i,ju) - rznew(iufrom)
                              endif
                           else
                              x = rxu(i,ju) - rxu(i,iufrom)
                              y = ryu(i,ju) - ryu(i,iufrom)
                              z = rzu(i,ju) - rzu(i,iufrom)
                           endif

                              
                           hdist = dsqrt( x**2 + y**2 + z**2 )
c     --- use law of cosines to calculate bond angle
                           thetac = (lengtha**2 + lengthb**2 
     &                          - hdist**2)
     &                          / (2.0d0 * lengtha * lengthb)
                           
c     --- check to make sure this will give a number
                           
                           if (abs(thetac).gt.1.0d0) then
                              vtor = 0
                              vphi = 0
                              vvib = 0
                              bf_tor(itor) = 0
                              goto 190
                           endif
                           angle = dacos(thetac)
                         
                           vphi = vphi + kforceb(iufrom,ju) * (angle
     &                          - equilb(iufrom,ju))**2
                           
c                           if (lshit) then
c                              write(iou,*) iufrom,iu,ju
c     &                             ,kforceb(iufrom,ju) 
c     &                             * (angle
c     &                          - equilb(iufrom,ju))**2
c                           endif

                        enddo
                           
                     endif
                     
                  elseif (fcount(iu).eq.0) then
                     
c     --- we use the angle with the closed bond
                     if (.not.ldo) then
                        ldo = .true.
                        countb = countb + 1
                        iopen(countb) = count
                        goto 40
                     elseif (countb.eq.1.and.iopen(1).ne.count) then
                        countb = countb + 1
                        iopen(countb) = count
                        goto 40
                     elseif (countb.eq.2) then
                        countb = countb + 1
                     else
                        ldo = .false.
                        countb = countb + 1
                     endif

                     if (.not.lnew.and.ip.eq.1.and.itor.eq.1) then
c     --- give old unit vector for connection
                        xx(count) = rxu(i,iu) - rxu(i,iufrom)
                        yy(count) = ryu(i,iu) - ryu(i,iufrom)
                        zz(count) = rzu(i,iu) - rzu(i,iufrom)
                     else
                        do counta = 1, ntogrow
                           if (counta.ne.count) then

                              if (opencount.gt.1) then
                                 if (countb.eq.3) then
c     --- we only want the angle with the closing bead
                                    if (iopen(2).eq.counta) then
                                       goto 35
                                    endif
                                 elseif (countb.eq.4) then
c     --- we want the angle with the open bead
                                    if (iopen(1).ne.counta) then
                                       goto 35
                                    endif
                                 endif
                              elseif(iclose.ne.counta) then
                                 goto 35
                              endif                                 
                              ju = growlist(iw,counta)

c     --- determine position by angle with this bond

                              lengtha = distij(iufrom,iuprev)
                              lengthb = flength(iufrom,ju)

                              rxt(1) = xvec(iufrom,iuprev)/lengtha
                              ryt(1) = yvec(iufrom,iuprev)/lengtha
                              rzt(1) = zvec(iufrom,iuprev)/lengtha

                              rxt(2) = xx(counta)/lengthb
                              ryt(2) = yy(counta)/lengthb
                              rzt(2) = zz(counta)/lengthb

                              angles(1) = bangles(count,1)
                              angles(2) = bangles(count,2)

                              call close(3,rxt,ryt,rzt,dum
     &                             ,angles,lterm)

                              if (lterm) then
                                 lterm = .false.
                                 vphi = 0
                                 vdha = 0
                                 vvib = 0
                                 bf_tor(itor) = 0
                                 goto 190
                              endif

c     --- since there are two possibilities, choose one at random
                              lengtha = flength(iufrom,iu)
                              if (random().lt.0.5d0) then
                                 xx(count) = rxt(3) * lengtha
                                 yy(count) = ryt(3) * lengtha
                                 zz(count) = rzt(3) * lengtha
                                 dir(itor,count) = 1
                              else
                                 xx(count) = rxt(4) * lengtha
                                 yy(count) = ryt(4) * lengtha
                                 zz(count) = rzt(4) * lengtha
                                 dir(itor,count) = 2
                              endif     
                           
c     ------------------------------------------------------------------
 35                           continue
                           endif
                        enddo
                     endif
                                             
                  else
c     --- we have to close differently

                     goto 132


c     **************************************************************

                     if (ip.eq.1.and.itor.eq.1) then
                        if (lnew) then
                           rxt(1) = rxnew(iufrom)
                           ryt(1) = rynew(iufrom)
                           rzt(1) = rznew(iufrom)
                        else
                           rxt(1) = rxu(i,iufrom)
                           ryt(1) = ryu(i,iufrom)
                           rzt(1) = rzu(i,iufrom)
                        endif
                        
                        do j = 1, fcount(iu)
                           ju = fclose(iu,j)
                           rxt(j+1) = rxu(i,ju)
                           ryt(j+1) = ryu(i,ju)
                           rzt(j+1) = rzu(i,ju)
                        enddo
                                                      
                        call close(1,rxt,ryt,rzt,flength(iufrom,iu)
     &                       ,angles,lterm)

                        if (lterm) then
                           return
                        endif

                        if (lnew) then
                           x = rxt(1) - rxnew(iufrom)
                           y = ryt(1) - rynew(iufrom)
                           z = rzt(1) - rznew(iufrom)
                        else
                           x = rxt(1) - rxu(i,iufrom)
                           y = ryt(1) - ryu(i,iufrom)
                           z = rzt(1) - rzu(i,iufrom)
                        endif
                        
                        length = dsqrt(x**2+y**2+z**2)
                        
                        if (lnew) then
                           x = rxt(2) - rxnew(iufrom)
                           y = ryt(2) - rynew(iufrom)
                           z = rzt(2) - rznew(iufrom)
                        else
                           x = rxt(2) - rxu(i,iufrom)
                           y = ryt(2) - ryu(i,iufrom)
                           z = rzt(2) - rzu(i,iufrom)
                        endif
                           
                        length = dsqrt(x**2+y**2+z**2)
                        
                        rxa(count,1) = rxt(1)
                        rya(count,1) = ryt(1)
                        rza(count,1) = rzt(1)

                        rxa(count,2) = rxt(2)
                        rya(count,2) = ryt(2)
                        rza(count,2) = rzt(2)
                     endif

c     --- determine position at random
                        
                     j = int(2.0d0 * random()) + 1

                     if (.not.lnew.and.ip.eq.1.and
     &                    .itor.eq.1) then
c     --- give old unit vector for connection
                        xx(count) = rxu(i,iu) - rxu(i,iufrom)
                        yy(count) = ryu(i,iu) - ryu(i,iufrom)
                        zz(count) = rzu(i,iu) - rzu(i,iufrom)
                     else
                        dir(itor,count) = j
                        if (lnew) then
                           xx(count) = rxa(count,j) - rxnew(iufrom) 
                           yy(count) = rya(count,j) - rynew(iufrom) 
                           zz(count) = rza(count,j) - rznew(iufrom)
                        else
                           xx(count) = rxa(count,j) - rxu(i,iufrom) 
                           yy(count) = rya(count,j) - ryu(i,iufrom) 
                           zz(count) = rza(count,j) - rzu(i,iufrom) 
                        endif
                     endif
c     *****************************************************************


 132                 continue

                  endif
                     
                  xvec(iufrom,iu) = xx(count)
                  xvec(iu,iufrom) = -xx(count)
                  yvec(iufrom,iu) = yy(count)
                  yvec(iu,iufrom) = -yy(count)
                  zvec(iufrom,iu) = zz(count)
                  zvec(iu,iufrom) = -zz(count)

c     --- determine position of trial spot
                  if (lnew) then
                     rxt(count) = rxnew(iufrom) + xx(count)
                     ryt(count) = rynew(iufrom) + yy(count)
                     rzt(count) = rznew(iufrom) + zz(count)
                  else
                     rxt(count) = rxu(i,iufrom) + xx(count) 
                     ryt(count) = ryu(i,iufrom) + yy(count) 
                     rzt(count) = rzu(i,iufrom) + zz(count) 
                  endif
               
                  lengtha = flength(iufrom,iu)
                                    
c     *** now that we have the trial spot, determine the bending
c     *** energies that are applicable                 

                  if (iopen(1).eq.count.or.iopen(2).eq.count) then
c     --- we already calculated all necessary bending energies for these
                     goto 40
                  endif

c     --- first determine bending energy for iuprev-iufrom-iu
                  if (iuprev.ne.0) then
                     length = distij(iufrom,iuprev)
                     thetac = -(xx(count)*xvec(iuprev,iufrom) 
     &                    + yy(count)*yvec(iuprev,iufrom)
     &                    + zz(count)*zvec(iuprev,iufrom)) 
     &                    / (lengtha*length)
                     
                     angle = dacos(thetac)
                     
                     vphi = vphi + kforcea(count) * (angle-equila(count)
     $                    )**2

c                     if (lshit) then
c                        write(iou,*) iu,iufrom,iuprev,kforcea(count)
c     &                       * (angle-equila(count))**2
c                     endif
                  endif

                  do j = 1, fcount(iu)
c     --- determine vectors from trial spot to iend
                     
                     ju = fclose(iu,j)
                     if (movetype.eq.2.and.lnew) then
                        xvec(iu,ju) = rxnew(ju) - rxt(count)
                        yvec(iu,ju) = rynew(ju) - ryt(count)
                        zvec(iu,ju) = rznew(ju) - rzt(count)
                     else
                        xvec(iu,ju) = rxu(i,ju) - rxt(count)
                        yvec(iu,ju) = ryu(i,ju) - ryt(count)
                        zvec(iu,ju) = rzu(i,ju) - rzt(count)
                     endif
                     xvec(ju,iu) = - xvec(iu,ju)
                     yvec(ju,iu) = - yvec(iu,ju)
                     zvec(ju,iu) = - zvec(iu,ju)

                     lengthb = flength(iu,ju)

c     --- now determine bending energy for iu-iend-ipast if it exists    
                     if (pastnum(ju).ne.0) then
                        do counta = 1, pastnum(ju)
                           ku = ipast(ju,counta)
                           length = distij(ju,ku)
                           thetac = - (xvec(iu,ju) * xvec(ju,ku)
     &                          + yvec(iu,ju) * yvec(ju,ku)
     &                          + zvec(iu,ju) * zvec(ju,ku))
     &                          / (length * lengthb)
                           angle = dacos( thetac )
                          
                           vphi = vphi + kforceb(iu,ku) * (angle
     $                          -equilb(iu,ku))**2
                           
c                           if (lshit) then
c                              write(iou,*) iu,ju,ku
c     &                             ,kforceb(iu,ku)
c     &                          * (angle-equilb(iu,ku))**2
c                           endif
                        enddo
                     endif
                  enddo
 40               continue
               enddo
               if (ldo) then
                 
                  if (countb.eq.3) then
                     start = iopen(2)
                     last = iopen(2)
                  else
                     start = iopen(1)
                     last = iopen(1)
                  endif
                  goto 30
               endif
             
c     *** calculate four sets of torsion energies 
c     ***   iu to iend to ipast to inext

c ------------------------------------------------------------------------
c     --- first calculate torsion for iu-iufrom-
               do count = 1, ntogrow
                  iu = growlist(iw,count)
                  do it = 1, intor(imolty,iu)
                     jut2 = ijtor2(imolty,iu,it)
                     
                     if (jut2.eq.iufrom) then
                        jut3 = ijtor3(imolty,iu,it)
                        jut4 = ijtor4(imolty,iu,it)

                        if (lpnow(jut4)) goto 41

c     --- jut4 must already exist or we made a big mistake
                        if (.not. lexist(jut4)) then
                           if (jut4.gt.iring(imolty)) then
                              goto 41
                           endif
                           write(iou,*) 'jut4,jut3,jut2,iu',
     &                          jut4,jut3,jut2,iu
                           call cleanup('trouble jut4 in crankshaft')
                        endif
                        jttor = ittor(imolty,iu,it)

                        call calctor(iu,jut2,jut3,jut4,jttor,vctor)

                        vdha = vdha + vctor

 41                     continue

                     endif
                  enddo

c     --- now lets calculate torsions for iend-iu-iufrom-iuprev

                  do j = 1, fcount(iu)
                     ju = fclose(iu,j)
                     do it = 1, intor(imolty,ju)
                        jut2 = ijtor2(imolty,ju,it)
                        jut3 = ijtor3(imolty,ju,it)
                        jut4 = ijtor4(imolty,ju,it)

                        if (jut2.eq.iu.and.jut4.eq.iuprev)
     &                       then

                           jttor = ittor(imolty,ju,it)

                           call calctor(ju,jut2,jut3,jut4,jttor,vctor)
                           
c                       --- add torsion energy to vdha
                           vdha = vdha + vctor
                        endif
                     enddo

c     --- calculate torsions for ipast-iend-iu-

                  if (pastnum(ju).ne.0) then      
                  do counta = 1, pastnum(ju)
                     ku = ipast(ju,counta)
                     do it = 1, intor(imolty,ku)
                        jut2 = ijtor2(imolty,ku,it)
                        jut3 = ijtor3(imolty,ku,it)
                        
                        if (jut2.eq.ju.and.jut3.eq.iu) then
                           jut4 = ijtor4(imolty,ku,it)
                           
                           if (lpnow(jut4)) goto 42

c     --- jut4 must already exist or we made a big mistake
                           if (.not. lexist(jut4)) then
                              write(iou,*) 'jut4,jut3,jut2,iu',
     &                             jut4,jut3,jut2,iu
                              call cleanup('trouble jut4 in crankshaft')
                           endif

                           jttor = ittor(imolty,ku,it)


                           call calctor(ku,jut2,jut3,jut4,jttor,vctor)

c                       --- add torsion energy to vdha
                           vdha = vdha + vctor
 42                        continue
                        endif
                     enddo

c     --- calculate torsions for inext-ipast-iend-iu

                  if (nextnum(ku).ne.0) then
                     do ja = 1, nextnum(ku)
                     nu = inext(ku,ja)

                     do it = 1, intor(imolty,nu)
                        jut2 = ijtor2(imolty,nu,it)
                        jut3 = ijtor3(imolty,nu,it)
                        jut4 = ijtor4(imolty,nu,it)

                        if (jut2.eq.ku.and.jut3.eq.ju.
     &                       and.jut4.eq.iu) then
                        
                           jttor = ittor(imolty,nu,it)

                           call calctor(nu,jut2,jut3,jut4,jttor,vctor)
                           
c                       --- add torsion energy to vdha
                           vdha = vdha + vctor

                        endif
                     enddo
                  enddo
                  endif
c     --- leave these unidented to save space
                  enddo
                  endif
                  enddo
               enddo
c     --- done determining torsions
c ------------------------------------------------------------------------
             

c     --- determine angles for iend - iu - iend
               do count = 1, ntogrow
                  iu = growlist(iw,count)

                  if (fcount(iu).gt.1) then

                     do j = 1, fcount(iu) - 1
                        ju = fclose(iu,j)
                        lengtha = flength(iu,ju)
                        
                        do counta = j + 1, fcount(iu)

                           ku = fclose(iu,counta)
                           
                           lengthb = flength(iu,ku)
                           
                           thetac = (xvec(iu,ju)*xvec(iu,ku)
     &                          + yvec(iu,ju)*yvec(iu,ku)
     &                          + zvec(iu,ju)*zvec(iu,ku))
     &                          / (lengtha*lengthb)

                           if (abs(thetac).gt.1) then
                              write(iou,*) '*********************'
     &                             ,'****************************'
                              write(iou,*) iu,ku,xvec(iu,ku)
     &                             ,yvec(iu,ku)
     &                             ,zvec(iu,ku),lengthb
                              call cleanup('shitfuck')
                           endif

                           angle = dacos(thetac)

                           vphi = vphi + kforceb(ju,ku)
     &                          * (angle-equilb(ju,ku))**2
       
c                           if (lshit) then
c                              write(iou,*) ju,iu,ku,kforceb(ju,ku)
c     &                          * (angle-equilb(ju,ku))**2
c                           endif
                           

                        enddo
                     enddo
                  endif
               enddo

c ------------------------------------------------------------------------
c     --- now lets figure out the rest of the bending for ntogrow > 1

               if (ntogrow.gt.1) then
                  do count = 1, ntogrow - 1
                     iu = growlist(iw,count)
                     do counta = count+1, ntogrow
                        ju = growlist(iw,counta)
                        if (.not.((iopen(1).eq.count
     &                       .or.iopen(1).eq.counta)
     &                       .and.((iclose.eq.count.or
     &                       .iclose.eq.counta)
     &                       .or.iopen(2).eq.count
     &                       .or.iopen(2).eq.counta))) then
c     --- we already calculated these
                           
                           lengtha = flength(iufrom,iu)
                           lengthb = flength(iufrom,ju)

                           thetac = (xx(count) * xx(counta)
     &                          + yy(count) * yy(counta)
     &                          + zz(count) * zz(counta)) / 
     &                          (lengtha * lengthb)

                           angle = dacos(thetac)
                           
                           vphi = vphi + kforceb(iu,ju) * (angle
     $                          -equilb(iu,ju))**2
                           
c                           if (lshit) then
c                              write(iou,*) iu,iufrom,ju,kforceb(iu,ju)
c     &                          * (angle-equilb(iu,ju))**2.0d0
c                           endif
                        endif
                     enddo
                  enddo
               endif

               bf_tor(itor) = dexp(-(vvib + vphi + vdha)*beta)
 190           continue
               
               vvibration(itor) = vvib 
               vtorsion(itor) = vdha
               vbend(itor) = vphi 
               bsum_tor(ip) = bsum_tor(ip) + bf_tor(itor)
            enddo

            if (lnew.or.ip.ne.1) then
c     --- choose one of the trial sites in a biased fashion
               ran_tor = random() * bsum_tor(ip)
               bs = 0
               do itor = 1, ichtor
                  bs = bs + bf_tor(itor)
                  if (ran_tor .lt. bs ) then

c     --- save torsion energy of this trial position
                     vtgtr(ip) = vtorsion(itor)
                     vtbend(ip) = vbend(itor) + ovphi
                     vtvib(ip) = vvibration(itor) + ovvib
c     --- assign the phidisp of this trial position
                     do count = 1, ntogrow
                        phiacc(count) = phicrank(itor,count)
                        diracc(count) = dir(itor,count)
                     enddo
c     --- exit the loop
                     goto 200
                  endif
               enddo

 200           continue
            else
c     --- select old conformation
               vtgtr(ip) = vtorsion(1)
               vtbend(ip) = vbend(1) + ovphi
               vtvib(ip) = vvibration(1) + ovvib
            endif

c     --- divide bsum by ichtor 
            bsum_tor(ip) = bsum_tor(ip) / dble(ichtor)

            start = 1
            last = ntogrow
            ldo = .false.
            countb = 0            
 60         continue

c     --- for accepted phidisp set up the vectors
            do count = 1,ntogrow

               iu = growlist(iw,count)

               if (.not.lnew.and.ip.eq.1) then
                  length = flength(iufrom,iu)   
                  xx(count) = (rxu(i,iu) - rxu(i,iufrom))
     &                 / length
                  yy(count) = (ryu(i,iu) - ryu(i,iufrom))
     &                 / length
                  zz(count) = (rzu(i,iu) - rzu(i,iufrom))
     &                 / length
               else
                  if (fcount(iu).gt.0) then
                     ju = fclose(iu,1)
                     call cone(1, xvec(iufrom,ju),yvec(iufrom,ju)
     &                    ,zvec(iufrom,ju),dum,dum,dum,dum,dum)
                     phidisp = phiacc(count)
                     
                     call cone(2,dum,dum,dum,alpha(count,ju)
     &                    ,phidisp,x,y,z)

c     --- store the unit vectors in xx, yy, zz
                     xx(count) = x
                     yy(count) = y
                     zz(count) = z

c                  elseif (fcount(iu).gt.1) then

c                     j = diracc(count)

c                     if (ip.gt.1.and..not.lnew) call cleanup('')

c                     if (lnew) then
c                        xx(count) = rxa(count,j) - rxnew(iufrom)  
c                        yy(count) = rya(count,j) - rynew(iufrom)
c                        zz(count) = rza(count,j) - rznew(iufrom)
c                     else
c                        xx(count) = rxa(count,j) - rxu(i,iufrom)
c                        yy(count) = rya(count,j) - ryu(i,iufrom)
c                        zz(count) = rza(count,j) - rzu(i,iufrom)
c                     endif

c                     xx(count) = xx(count) / length
c                     yy(count) = yy(count) / length
c                     zz(count) = zz(count) / length

                  else

                     if (.not.ldo) then
                        ldo = .true.
                        countb = countb + 1
                        iopen(countb) = count
                        goto 70
                     elseif (countb.eq.1.and.iopen(1).ne.count) then
                        countb = countb + 1
                        iopen(countb) = count
                        goto 70
                     elseif (countb.eq.2) then
                        countb = countb + 1
                     else
                        ldo = .false.
                        countb = countb + 1
                     endif
                                       
                     do counta = 1, ntogrow

                        if (counta.ne.count) then
                           if (opencount.gt.1) then
                              if (countb.eq.3) then
c     --- we only want the angle with the closing bead

                                 if (iopen(2).eq.counta) then
                                    goto 65
                                 endif
                              elseif (countb.eq.4) then
c     --- we want the angle with the open bead

                                 if (iopen(1).ne.counta) then
                                    goto 65
                                 endif
                              endif
                           elseif(iclose.ne.counta) then
                              goto 65
                           endif

                           lengtha = distij(iufrom,iuprev)

                           rxt(1) = xvec(iufrom,iuprev)/lengtha
                           ryt(1) = yvec(iufrom,iuprev)/lengtha
                           rzt(1) = zvec(iufrom,iuprev)/lengtha

                           rxt(2) = xx(counta)
                           ryt(2) = yy(counta)
                           rzt(2) = zz(counta)

                           angles(1) = bangles(count,1)
                           angles(2) = bangles(count,2)

                           call close(3,rxt,ryt,rzt,dum
     &                          ,angles,lterm)


                           if (lterm) then
                              return
                           endif

c     --- choose the given possibility
                           if (diracc(count).eq.1) then
                              lengtha = flength(iufrom,iu)
                              xx(count) = rxt(3) 
                              yy(count) = ryt(3)
                              zz(count) = rzt(3) 
                           else
                              lengtha = flength(iufrom,iu)
                              xx(count) = rxt(4) 
                              yy(count) = ryt(4) 
                              zz(count) = rzt(4) 
                           endif
                        endif

 65                     continue

                     enddo
                  endif
               endif
 70            continue
            enddo

            if (ldo) then
               if (countb.eq.3) then
                  start = iopen(2)
                  last = iopen(2)
               else
                  start = iopen(1)
                  last = iopen(1)
               endif
               goto 60
            endif
            
c     --- accepted coordinates, save them in r*p(trial)
            do count = 1,ntogrow
               iu = growlist(iw,count)
               length = flength(iufrom,iu)
               if ( lnew ) then
c     --- use new positions
                  rxp(count,ip) = rxnew(iufrom) + xx(count)*length
                  ryp(count,ip) = rynew(iufrom) + yy(count)*length
                  rzp(count,ip) = rznew(iufrom) + zz(count)*length
               else
c     --- use old coordinates
                  rxp(count,ip) = rxu(i,iufrom) + xx(count)*length
                  ryp(count,ip) = ryu(i,iufrom) + yy(count)*length
                  rzp(count,ip) = rzu(i,iufrom) + zz(count)*length
               endif
            enddo
         enddo

         if (lcrank) then
            vphi = vibtr
         endif
      endif
     
c      write(iou,*) 'END SAFECMBC ',iinit,'!'
      return 

      end
      












