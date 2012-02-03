      subroutine rosenbluth ( lnew,lterm,i,icharge,imolty,ifrom,ibox
     &     ,igrow,wadd,lfixnow,cwtorf,movetype )

c rosenbluth
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
c    **   performs a configurational bias move for branched molecules **
c    *******************************************************************

      implicit none
      
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'poten.inc'
      include 'conver.inc'
      include 'cbmc.inc'
      include 'rosen.inc' 
      include 'connect.inc'
      include 'fix.inc'
      include 'ipswpar.inc'

c     --- variables passed to the subroutine
      logical lnew,lterm,lwbef
      integer i,j,ja,icharge,imolty,ifrom,ibox,igrow,tac

c     --- local variables
      
      logical ovrlap,ltorsion,lfixnow,lfixed,lreturn

      integer glist,iuprev,iufrom,ichoi,ntogrow,count
      integer iu,iv,iw,ju,ip,ichtor
     &       ,it,jut2,jut3,jut4,jttor,iwalk,ivect
      integer angstart,toracc
      
      dimension glist(numax)

      double precision dum,xub,yub,zub,length,lengtha,lengthb,wadd

      double precision vdha,x,y,z,maxlen,vorient,vtorf
     &                ,xaa1,yaa1,zaa1,daa1,xa1a2,ya1a2,za1a2,da1a2
     &                ,thetac,dot,rbf,bsum,bs,random,xcc,ycc,zcc,tcc
      double precision vbbtr,vvibtr,vtorso,wei_vib,wbendv,dist
      double precision bondlen,bendang,phi,phidisp,phinew,thetanew
      double precision cwtorf,vfbbtr,vphi,theta,spltor,rm

      dimension bondlen(numax),bendang(numax),phi(numax)
     
c     -- new stuff
      integer itor,bin,counta,movetype,ku
      double precision bf_tor,vtorsion,phitors,ran_tor
     &     ,wei_bend,jacobian,ctorf
      dimension bf_tor(nchtor_max),vtorsion(nchtor_max)
     &     ,phitors(nchtor_max),ctorf(nchmax,nchtor_max)
      dimension toracc(nchmax)
     &     ,vfbbtr(nchmax,nchtor_max)


c ------------------------------------------------------------------

c      write(2,*) 'start ROSENBLUTH'
      lterm = .false.
      cwtorf = 1.0d0
      wei_vib = 1.0d0

c *******************************************
c * Rosenbluth weight of trial conformation *
c *******************************************

c        --- initialize conformation energies and weight
      if ( lnew ) then
c        --- set the initial weight to unity ***
         weight = 1.0d0
c        --- set total energy of trial configuration to zero ***
         vnewt     = 0.0d0
         vnewtg    = 0.0d0
         vnewbb    = 0.0d0
         vnewbvib  = 0.0d0
         vnewext   = 0.0d0
         vnewintra = 0.0d0
         vnewinter = 0.0d0
         vnewelect = 0.0d0
         vnewewald = 0.0d0
         vipswn = 0.0d0
         vwellipswn = 0.0d0
      else
c        --- old conformation
c        --- set the initial weight of the old configuration to unity ***
         weiold = 1.0d0
c        --- set total energy of trial configuration to zero ***
         voldt     = 0.0d0
         voldtg    = 0.0d0
         voldbb    = 0.0d0
         voldbvib  = 0.0d0
         voldext   = 0.0d0
         voldintra = 0.0d0
         voldinter = 0.0d0
         voldelect = 0.0d0
         voldewald = 0.0d0
         vipswo = 0.0d0
         vwellipswo = 0.0d0
      endif

c     --- for rigid molecules
      if (lrigid(imolty).and.movetype.ne.1) then
         call rigrot( lnew,lterm,i,icharge,imolty,ibox,wadd )
         
         if (rindex(imolty).eq.0) then
            return
         endif
         if (lterm) then
            return
         endif
      endif

c     --- set lexist to lexshed
      do iu = 1,igrow
         lexist(iu) = lexshed(iu)
      enddo

c     --- calculate all bond vectors for lexist 
      do iu = 1, igrow
         do iv = 1, invib(imolty,iu)
            ju = ijvib(imolty,iu,iv)
            if ( lexist(iu) .and. lexist(ju) ) then
               if ( lnew ) then
c                 --- use new coordinates
                  xvec(iu,ju) = rxnew(ju) - rxnew(iu)
                  yvec(iu,ju) = rynew(ju) - rynew(iu)
                  zvec(iu,ju) = rznew(ju) - rznew(iu)
               else
c                 --- use old coordinates
                  xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
                  yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
                  zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
               endif
               distij(iu,ju) = dsqrt( xvec(iu,ju)**2
     +              + yvec(iu,ju)**2 + zvec(iu,ju)**2 )
            endif
         enddo
      enddo 

c *************************
c * loop over trial units *
c *************************
      do 200 iw = 1, ifrom
c        --- set vibration and bending energies for this growth to 0.0

         if (llrig.and.lsave(iw)) goto 200

         iufrom = growfrom(iw)
         ntogrow = grownum(iw)
         lfixed = .false.
         lwbef = .false.
         if (lfixnow) then
            do count = 1, ntogrow
               iu = growlist(iw,count)
               do j = 1, wbefnum
                  do ja = 1, befnum(j)
                     if (iu.eq.ibef(j,ja)) then
c     --- time to do final crankshaft move
                        call safecbmc(3,lnew,i,iw,igrow,imolty
     &                       ,count,x,y,z,vphi,vtorf,wbendv
     &                       ,lterm,movetype)

                        if (lterm) then
                           return
                        endif
                        lfixed = .true.
                        wei_bend = wbendv
                        if (lcrank) then
                           vvibtr = vphi
                        else
                           vvibtr = 0.0d0
                        endif
                        do counta = 1, ntogrow
                           glist(counta) = growlist(iw,counta)
                        enddo
                        ichoi = nchoi(imolty)
c     --- sometimes this loop makes the code skip geometry, which sets this
c     --- maxlen (the maximum bond length that CBMC will try to grow)
                        maxlen=2.0d0
                        goto 250
                     endif
                  enddo
               enddo
            enddo
         endif      
         
c        --- perform the biased selection of bond angles and get lengths
         call geometry(lnew,iw,i,imolty,angstart,iuprev,glist
     &        ,bondlen,bendang,phi,vvibtr,vbbtr,maxlen,wei_bend )

c         write(2,*) 'lnew, wei_bend',lnew,wei_bend

c     --- for lfixnow check if there are two beads to go

         if (lfixnow) then
            do count = 1, ntogrow
               iu = growlist(iw,count)
               do j = 1, wbefnum
                  if (iu.eq.iwbef(j)) then
c     --- lets setup for two beads to go
                     
                     call safecbmc(1,lnew,i,iw,igrow,imolty
     &                    ,count,x,y,z,vphi,vtorf,wbendv
     &                    ,lterm,movetype)
                     wei_vib = wei_vib * wbendv
                     vvibtr = vvibtr + vphi
                     lwbef = .true.
                     goto 150
                  endif
               enddo
            enddo
         endif
 150     continue


c        --- we now have the bond lengths and angles for the grown beads
c        -- select nchoi trial positions based only on torsions
         ichoi = nchoi(imolty)
         ichtor = nchtor(imolty)
         do ip = 1,ichoi
            ivect = 0

            lreturn = .false.
 205        continue
            
c           --- set up the cone based on iuprev (could be grown if no prev)
            if ( .not.lreturn.and.growprev(iw) .eq. 0 ) then
c              --- calculate random vector on the unit sphere for the first bead
               count = 1
               if ( (.not. lnew) .and. ip .eq. 1 ) then
c                 --- use old molecule position
                  iu = growlist(iw,count)
                  length = bondlen(count)
c                 --- compute unit vector to be used in cone and torsion
                  x = ( rxu(i,iu) - rxu(i,iufrom) )/length 
                  y = ( ryu(i,iu) - ryu(i,iufrom) )/length 
                  z = ( rzu(i,iu) - rzu(i,iufrom) )/length 
c                 --- store this in xx yy zz
                  xx(count) = x
                  yy(count) = y
                  zz(count) = z
               else
c                 --- choose randomly on the unit sphere
                  call sphere(x,y,z)
                  xx(count) = x
                  yy(count) = y
                  zz(count) = z
               endif

               if ( ntogrow .gt. 1 ) then
c                 ---set up the cone 
                  xub = -x 
                  yub = -y 
                  zub = -z 
                  call cone (1,xub,yub,zub,dum,dum,dum,dum,dum )
               endif

               if (lrigid(imolty)) then
                  growprev(iw)=riutry(imolty,iw)+1
                  lreturn = .true.
                  goto 205
               endif

               ltorsion = .false.
               
            else
c              --- set up the cone based on iuprev and iufrom
               length = distij(iuprev,iufrom)
               xub = xvec(iuprev,iufrom) / length 
               yub = yvec(iuprev,iufrom) / length
               zub = zvec(iuprev,iufrom) / length

               call cone (1,xub,yub,zub,dum,dum,dum,dum,dum )
               if (movetype.eq.2.and.lring(imolty)
     &              .and.iw.eq.1) then
                  ltorsion = .false.
               else
                  ltorsion = .true.
               endif
            endif

c           --- Begin loop to determine torsional angle

            if ( ltorsion ) then

c              --- initialize bsum_tor
               bsum_tor(ip) = 0.0d0
               
               do itor = 1,ichtor

                  if ( (.not. lnew) .and. ip .eq. 1 
     &                 .and. itor .eq. 1) then
c                    --- old conformation - set phidisp to 0.0d0
                     phidisp = 0.0d0
                  else
c                    --- choose a random displacement angle from anglestart
c                    --- assign the positions based on angles and lengths above
                     phidisp = twopi*random()
                  endif

                  do count = angstart,ntogrow
                     phinew = phi(count) + phidisp
                     thetanew = bendang(count)
                     
                     call cone(2,dum,dum,dum,thetanew,phinew,x,y,z)
c     --- store the unit vectors in xx, yy, zz
                     xx(count) = x
                     yy(count) = y
                     zz(count) = z
                  enddo

c              --- set energies of trial position to zero ---
                  vdha = 0.0d0
c              --- compute torsion energy for given trial conformation
                  do count = 1,ntogrow
                     iu = growlist(iw,count)

                     if (movetype.eq.2.and
     &                    .lring(imolty).and.iw.lt.3) then
                        bf_tor(itor) = 1.0d0
                        goto 300
                     endif

                     do 299 it = 1, intor(imolty,iu)
                        jut2 = ijtor2(imolty,iu,it)
                        jut3 = ijtor3(imolty,iu,it)
                        if ( jut2 .eq. iufrom .and. 
     &                       jut3 .eq. iuprev) then
                           jut4 = ijtor4(imolty,iu,it)
                           
c                       --- jut4 must already exist or we made a big mistake
                           if ( .not. lexist(jut4) )  then
c * allow regrowth where one torsion may already exist and one may not
                              goto 299
c                              write(2,*) 'jut4,jut3,jut2,iu',
c     &                             jut4,jut3,jut2,iu
c                              stop 'trouble jut4'
                           endif
                           jttor = ittor(imolty,iu,it)
                           
c                       --- calculate cross products d_a x d_a-1 
                           xaa1 = yy(count) * zvec(jut3,jut2) 
     &                          + zz(count) * yvec(jut2,jut3)
                           yaa1 = zz(count) * xvec(jut3,jut2) 
     &                          + xx(count) * zvec(jut2,jut3)
                           zaa1 = xx(count) * yvec(jut3,jut2) 
     &                          + yy(count) * xvec(jut2,jut3)
                           
c                       --- calculate cross products d_a-1 x d_a-2
                           xa1a2 = yvec(jut2,jut3) * zvec(jut3,jut4) -
     +                          zvec(jut2,jut3) * yvec(jut3,jut4)
                           ya1a2 = zvec(jut2,jut3) * xvec(jut3,jut4) -
     +                          xvec(jut2,jut3) * zvec(jut3,jut4)
                           za1a2 = xvec(jut2,jut3) * yvec(jut3,jut4) -
     +                          yvec(jut2,jut3) * xvec(jut3,jut4)

c                       --- calculate lengths of cross products ***
                           daa1 = dsqrt ( xaa1**2 + yaa1**2 + zaa1**2 )
                           da1a2 = dsqrt ( xa1a2**2 + ya1a2**2 
     &                          + za1a2**2 )
                           
c                       --- calculate dot product of cross products ***
                           dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                           thetac = - (dot / ( daa1 * da1a2 ))
                           if (thetac.gt.1.0d0) thetac=1.0d0
                           if (thetac.lt.-1.0d0) thetac=-1.0d0

c     KEA -- added for extending range to +/- 180
                           if (L_tor_table) then
c     *** calculate cross product of cross products ***
                              xcc = yaa1*za1a2 - zaa1*ya1a2
                              ycc = zaa1*xa1a2 - xaa1*za1a2
                              zcc = xaa1*ya1a2 - yaa1*xa1a2
c     *** calculate scalar triple product ***
                              tcc = xcc*xvec(jut2,jut3)
     &                             + ycc*yvec(jut2,jut3)
     &                             + zcc*zvec(jut2,jut3)
                              theta = dacos(thetac)
                              if (tcc .lt. 0.0d0) theta = -theta
                              if (L_spline) then
                                 call splint(theta,spltor,jttor)
                              elseif(L_linear) then
                                 call lininter(theta,spltor,jttor)
                              endif

c                       --- add torsion energy to vdha
                              vdha = vdha + spltor
                           else
c                       --- add torsion energy to vdha
                              vdha = vdha + vtorso( thetac, jttor )
                           endif
                        endif

c                     enddo
 299                 continue

                  enddo

c                 --- compute boltzmann factor and add it to bsum_tor
                  bf_tor(itor) = dexp ( -vdha * beta )

 300              continue

c                 --- store vtorsion and phidisp for this trial
                  vtorsion(itor) = vdha
                  phitors(itor) = phidisp
                  
c                 --- for safecbmc add extra weight to assure closure
                  if (lfixnow) then
                     
                     ctorf(ip,itor) = 1.0d0
                     
                     vfbbtr(ip,itor) = 0

                     do count = 1, ntogrow
                        length = bondlen(count)

                        if (lnew) then
                           x = rxnew(iufrom) + xx(count)*length
                           y = rynew(iufrom) + yy(count)*length
                           z = rznew(iufrom) + zz(count)*length
                        else
                           x = rxu(i,iufrom) + xx(count)*length
                           y = ryu(i,iufrom) + yy(count)*length
                           z = rzu(i,iufrom) + zz(count)*length
                        endif
                        
                        iu = growlist(iw,count)
                       
                        if (movetype.eq.2.and.lnew) then
                        if (lwbef) then

c                       --- determine special closing energies
                           
                           call safecbmc(2,lnew,i,iw,igrow,imolty
     &                          ,count,x,y,z,vphi,vtorf,wbendv
     &                          ,lterm,movetype)
                       
                           bf_tor(itor) = bf_tor(itor) * vtorf 
     &                          * dexp( - beta * vphi )
                           ctorf(ip,itor) = ctorf(ip,itor) 
     &                          * vtorf
                           
                           vfbbtr(ip,itor) =  vfbbtr(ip,itor)
     &                          + vphi
                        else
                           do j = 1, fcount(iu)
                              ju = fclose(iu,j)
                              dist = dsqrt((x-rxnew(ju))**2
     &                             + (y-rynew(ju))**2
     &                             + (z-rznew(ju))**2)
                              bin = anint(dist*10.0d0)
                              
                              bf_tor(itor) = bf_tor(itor)
     &                             * probf(iu,ju,bin)
                              ctorf(ip,itor) = ctorf(ip,itor) * 
     &                             probf(iu,ju,bin)

                              if (iw.gt.2) then
                                 do counta = 1, pastnum(ju)
                                    ku = ipast(ju,counta)
                                    if (.not.lplace(imolty,ku)) then
                                       dist = dsqrt((x-rxnew(ku))**2
     &                                      + (y-rynew(ku))**2
     &                                      + (z-rznew(ku))**2)
                                       bin = anint(dist*10.0d0)
                                    
                                       bf_tor(itor) = bf_tor(itor)
     &                                      * probf(iu,ku,bin)
                                       ctorf(ip,itor) = ctorf(ip,itor) * 
     &                                      probf(iu,ku,bin)   
                                    endif
                                 enddo
                              endif


                           enddo
                        endif

                        else
                        if (lwbef) then

c                       --- determine special closing energies
                           
                           call safecbmc(2,lnew,i,iw,igrow,imolty
     &                          ,count,x,y,z,vphi,vtorf
     &                          ,wbendv,lterm,movetype)
                           
                           bf_tor(itor) = bf_tor(itor) * vtorf 
     &                          * dexp( - beta * vphi )
                           ctorf(ip,itor) = ctorf(ip,itor) 
     &                          * vtorf
                           
                           vfbbtr(ip,itor) =  vfbbtr(ip,itor)
     &                          + vphi
                        else
                        if (fcount(iu).gt.0) then
                           do j = 1, fcount(iu)
                              ju = fclose(iu,j)
                              dist = dsqrt((x-rxu(i,ju))**2
     &                             + (y-ryu(i,ju))**2
     &                             + (z-rzu(i,ju))**2)
                              bin = anint(dist*10.0d0)
                              
                              bf_tor(itor) = bf_tor(itor)
     &                             * probf(ju,iu,bin)
                              ctorf(ip,itor) = ctorf(ip,itor) * 
     &                             probf(iu,ju,bin)
               
                              if (pastnum(ju).ne.0) then
                                 do counta = 1, pastnum(ju)
                                    ku = ipast(ju,counta)
                                    if (.not.lplace(imolty,ku)) then
                                       dist = dsqrt((x-rxu(i,ku))**2
     &                                      + (y-ryu(i,ku))**2
     &                                      + (z-rzu(i,ku))**2)
                                       bin = anint(dist*10.0d0)
                                    
                                       bf_tor(itor) = bf_tor(itor)
     &                                      * probf(iu,ku,bin)
                                       ctorf(ip,itor) = ctorf(ip,itor) * 
     &                                      probf(iu,ku,bin)   
                                    endif

                                 enddo
                              endif
                           enddo
                        endif
                        endif
                        endif

                     enddo
                  endif
                  
                  bsum_tor(ip) = bsum_tor(ip) + bf_tor(itor)


               enddo

               if ( lnew .or. ip .ne. 1 ) then
c                 --- choose one of the trial sites in a biased fashion
                  ran_tor = random()*bsum_tor(ip)
                  bs = 0.0d0
                  do itor = 1,ichtor
                     bs = bs + bf_tor(itor)
                     if ( ran_tor .lt. bs ) then
c                       --- save torsion energy of this trial position
                        vtgtr(ip) = vtorsion(itor)

                        toracc(ip) = itor

c                       --- assign the phidisp of this trial postion
                        phidisp = phitors(itor)

c                       --- exit the loop
                        goto 100
                     endif
                  enddo
 100              continue
               else
c                 --- select the old conformation
                  vtgtr(ip) = vtorsion(1)
                  phidisp = phitors(1)
               endif

c              --- divide bsum by ichtor
               bsum_tor(ip) = bsum_tor(ip) / dble(ichtor)
               
            else

c              --- no torsion energy, choose phidisp at random (except old)
               if ( (.not. lnew) .and. ip .eq. 1 ) then
c                 --- old conformation - set phidisp to 0.0d0
                  phidisp = 0.0d0
               else

c                 --- choose a random displacement angle from anglestart
c                 --- assign the positions based on angles and lengths above
                  phidisp = twopi*random()
               endif

c              --- assign the torsional energy a value of 0.0
               vtgtr(ip) = 0.0d0

c              --- set bsum_tor to 1.0d0
               bsum_tor(ip) = 1.0d0

            endif

c           --- for accepted phidisp set up the vectors
            do count = angstart,ntogrow
               phinew = phi(count) + phidisp
               thetanew = bendang(count)
               
               call cone(2,dum,dum,dum,thetanew,phinew,x,y,z)
c              --- store the unit vectors in xx, yy, zz
               xx(count) = x
               yy(count) = y
               zz(count) = z

            enddo

c           --- accepted coordinates, save them in r*p(trial)
            do count = 1,ntogrow
               length = bondlen(count)
               if ( lnew ) then
c                 --- use new positions
                  rxp(count,ip) = rxnew(iufrom) + xx(count)*length
                  ryp(count,ip) = rynew(iufrom) + yy(count)*length
                  rzp(count,ip) = rznew(iufrom) + zz(count)*length
               else
c                 --- use old coordinates
                  rxp(count,ip) = rxu(i,iufrom) + xx(count)*length
                  ryp(count,ip) = ryu(i,iufrom) + yy(count)*length
                  rzp(count,ip) = rzu(i,iufrom) + zz(count)*length
               endif
               
            enddo

         enddo

c       --- now that we have the trial site need to compute non-bonded energy

 250     continue

         call boltz ( lnew,.false.,ovrlap,i,icharge,imolty,ibox,ichoi
     &        ,iufrom ,ntogrow, glist, maxlen)

         if ( ovrlap ) then
            lterm = .true.
            return
         endif

c ---------------------------------------------------------------------
 
c *** perform the walk according to the availibility of the choices ***
c *** and calculate the correct weight for the trial walk           ***

         bsum = 0.0d0
         do ip = 1, ichoi
c           --- include both the torsional and the LJ/qq
            bsum = bsum + bfac(ip)*bsum_tor(ip)
         enddo

         if ( lnew ) then
c           --- update new rosenbluth weight - include bending weight
            weight = weight * bsum * wei_bend * wei_vib

            if ( weight .lt. softlog ) then
               lterm=.true.
               return
            endif

c           --- select one position at random ---
            rbf = bsum * random()
            bs = 0.0d0 
            do 170 ip = 1, ichoi
               if ( .not. lovr(ip) ) then
                  bs = bs + bfac(ip)*bsum_tor(ip)
                  if ( rbf .lt. bs ) then
c                    --- select ip position ---
                     iwalk = ip
                     goto 180
                  endif
               endif
 170        continue
         else
c           --- old conformation, update weiold - include wei_bend
            weiold = weiold * bsum * wei_bend * wei_vib
            if (weiold .lt. softlog) write(2,*) '###old weight too low'
         endif

 180     continue

         if (lfixed) then
c     --- determine jacobian contribution for crankshaft
            jacobian = 1.0d0
            do count = 1, ntogrow
               iu = growlist(iw,count)
               if (fcount(iu).gt.0) then
                  counta = 1
                  ju = fclose(iu,counta)
                  if (lnew) then
                     if (movetype.eq.2) then
                        x = rxnew(ju) - rxnew(iufrom)
                        y = rynew(ju) - rynew(iufrom)
                        z = rznew(ju) - rznew(iufrom)
                        length = dsqrt( x**2 + y**2 + z**2 )
                        
                        x = rxnew(ju) - rxp(count,iwalk)
                        y = rynew(ju) - ryp(count,iwalk)
                        z = rznew(ju) - rzp(count,iwalk)
                        lengtha = dsqrt( x**2 + y**2 + z**2 )
                     else
                        x = rxu(i,ju) - rxnew(iufrom)
                        y = ryu(i,ju) - rynew(iufrom)
                        z = rzu(i,ju) - rznew(iufrom)
                        length = dsqrt( x**2 + y**2 + z**2 )

                        x = rxu(i,ju) - rxp(count,iwalk)
                        y = ryu(i,ju) - ryp(count,iwalk)
                        z = rzu(i,ju) - rzp(count,iwalk)
                        lengtha = dsqrt( x**2 + y**2 + z**2 )
                     endif
                                  
                     x = rxp(count,iwalk) - rxnew(iufrom)
                     y = ryp(count,iwalk) - rynew(iufrom)
                     z = rzp(count,iwalk) - rznew(iufrom)
                     lengthb = dsqrt( x**2 + y**2 + z**2 )
                  else
                     x = rxu(i,ju) - rxu(i,iufrom)
                     y = ryu(i,ju) - ryu(i,iufrom)
                     z = rzu(i,ju) - rzu(i,iufrom)
                     length = dsqrt( x**2 + y**2 + z**2 )

                     x = rxu(i,ju) - rxu(i,iu)
                     y = ryu(i,ju) - ryu(i,iu)
                     z = rzu(i,ju) - rzu(i,iu)
                     lengtha = dsqrt( x**2 + y**2 + z**2 )

                     x = rxu(i,iu) - rxu(i,iufrom)
                     y = ryu(i,iu) - ryu(i,iufrom)
                     z = rzu(i,iu) - rzu(i,iufrom)
                     lengthb = dsqrt( x**2 + y**2 + z**2 )
                  endif
                  jacobian = jacobian / (length*lengtha*lengthb)
               endif
            enddo
            bsum = bsum * jacobian
         endif

         if ( lnew ) then
            if (lfixnow) then
               if (lwbef) then
                  tac = toracc(iwalk)
                  vbbtr = vbbtr + vfbbtr(iwalk,tac)
               endif
               if (lfixed) then
                  vbbtr = vtbend(iwalk)
                  vvibtr = vvibtr + vtvib(iwalk)
               else

                  if (.not.(movetype.eq.2.and.lring(imolty).and 
     &                 .iw.eq.1)) then

                     tac = toracc(iwalk)
                     cwtorf = cwtorf * ctorf(iwalk,tac)

                  endif
               endif
            endif

c           --- update new trial energies
            vnewt     = vnewt     + vtry(iwalk)  
     +           + vtgtr(iwalk) + vvibtr + vbbtr
            vnewbvib  = vnewbvib  + vvibtr
            vnewbb    = vnewbb    + vbbtr
            vnewtg    = vnewtg    + vtgtr(iwalk)
            vnewext   = vnewext   + vtrext(iwalk)
            vnewintra = vnewintra + vtrintra(iwalk)
            vnewinter = vnewinter + vtrinter(iwalk)
            vnewelect = vnewelect + vtrelect(iwalk)
            vnewewald = vnewewald + vtrewald(iwalk)
            vipswn = vipswn+vipswnt(iwalk)
            vwellipswn = vwellipswn+vwellipswnt(iwalk)
         else
            if (lfixnow) then
               if (lfixed) then
                  vbbtr = vtbend(1)
                  vvibtr = vvibtr + vtvib(1)
               else
                                    
                  if (.not.(movetype.eq.2.and.lring(imolty).and 
     &                 .iw.eq.1)) then
                     
                     cwtorf = cwtorf * ctorf(1,1)

                  endif
               endif
               if (lwbef) then
                  vbbtr = vbbtr + vfbbtr(1,1)
               endif
            endif

c            --- update old trail energies
            voldt     = voldt     + vtry(1)  
     +           + vtgtr(1) + vvibtr + vbbtr
            voldbvib  = voldbvib  + vvibtr
            voldbb    = voldbb    + vbbtr
            voldtg    = voldtg    + vtgtr(1)
            voldext   = voldext   + vtrext(1)
            voldintra = voldintra + vtrintra(1)
            voldinter = voldinter + vtrinter(1)
            voldelect = voldelect + vtrelect(1)
            voldewald = voldewald + vtrewald(1)
            vipswo = vipswo+vipswot(1)
            vwellipswo = vwellipswo+vwellipswot(1)
         endif

         do count = 1,ntogrow
            iu = growlist(iw,count)
            if ( lnew ) then
c              --- assign new positions to r*new
               rxnew(iu) = rxp(count,iwalk)
               rynew(iu) = ryp(count,iwalk)
               rznew(iu) = rzp(count,iwalk)
            endif

c           --- set lexist(iu) to true so ewald sum computed properly
            lexist(iu) = .true.

c           --- store new existing vectors between beads and iufrom
            if ( lnew ) then
c              --- use r*new positions
               xvec(iu,iufrom) = rxnew(iufrom) - rxnew(iu)
               yvec(iu,iufrom) = rynew(iufrom) - rynew(iu)
               zvec(iu,iufrom) = rznew(iufrom) - rznew(iu)
            else
c              --- use r*u positions
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

            if (lfixnow) then
c     --- we must store new vectors with endpoints
               if (fcount(iu).gt.0) then
                  do j = 1, fcount(iu)
                     ju = fclose(iu,j)

                     if (lnew) then
                        xvec(iu,ju) = rxnew(ju) - rxnew(iu)
                        yvec(iu,ju) = rynew(ju) - rynew(iu)
                        zvec(iu,ju) = rznew(ju) - rznew(iu)
                     else
                        xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
                        yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
                        zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
                     endif

                     xvec(ju,iu) = - xvec(iu,ju)
                     yvec(ju,iu) = - yvec(iu,ju)
                     zvec(ju,iu) = - zvec(iu,ju)
                     
                     distij(iu,ju) = dsqrt(xvec(iu,ju)**2
     &                    + yvec(iu,ju)**2 + zvec(iu,ju)**2)
                     distij(ju,iu) = distij(iu,ju)
                     
                       

                  enddo
               endif
               
            endif
         enddo

c ********************************
c * end of loop over trial units *
c ********************************
 200  continue

c      write(2,*) 'end ROSENBLUTH'

c ------------------------------------------------------------------

      return
      end


