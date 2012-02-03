      subroutine rigfix(lnew,i,ibox,imolty,lterm,wrig)

c rigfix
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
      include 'rosen.inc'
      include 'fix.inc'
      include 'connect.inc'
      include 'system.inc'

      logical lnew,ovrlap,lterm,lovra,lfind,lshit

      integer iw,i,ibox,imolty,iufrom,iuprev,ntogrow,count,iu,counta
     &     ,ilist,ja,max,num,inum,j,ju,iv,nlist,ichoi,ichtor,ip,itor
     &     ,it,jut2,jut3,jut4,jttor,iwalk,glist,ifrom,inuma

      parameter(max=10)


      double precision xub,yub,zub,lengtha,lengthb,dum,xfix,yfix,zfix
     &     ,phia,bendang,thetac,twopi,phidisp,phi,rlength,vdha,vtor
     &     ,vtorsion,phitors,bf_tor,random,ran_tor,bs,rxpa,rypa,rzpa
     &     ,bsuma,vtrya,vtrintraa,vtrexta,vtrelecta,vtrewalda
     &     ,vtrorienta,vtrintera,bsum,rbf,wrig


      dimension ilist(numax),inum(max),xfix(numax),yfix(numax)
     &     ,zfix(numax),lfind(numax),phia(numax),bendang(numax)
     &     ,rlength(numax),vtorsion(nchtor_max),phitors(nchtor_max)
     &     ,bf_tor(nchtor_max),rxpa(numax,nchmax),rypa(numax,nchmax)
     &     ,rzpa(numax,nchmax),vtrelecta(nchmax),vtrewalda(nchmax)
     &     ,bsuma(nchmax),vtrya(nchmax),vtrintraa(nchmax)
     &     ,vtrorienta(nchmax),vtrexta(nchmax),glist(max)
     &     ,lovra(nchmax),vtrintera(nchmax),ifrom(numax),inuma(max)

c     ----------------------------------------------------------

c      write(6,*) 'START RIGFIX'

      twopi = dacos(-1.0d0) * 2.0d0
      wrig = 1.0d0
      do j = 1, nunit(imolty)
         lfind(j) = .false.
      enddo
            
      ichoi = nchoi(imolty)
      ichtor = nchtor(imolty)

      do iw = 1, nrigi
         iufrom = rfrom(iw)
         iuprev = rprev(iw)
         ntogrow = rnum(iw)
         
         lfind(iufrom) = .true.

c     --- we must first set up cone for old configuration
         xfix(iufrom) = rxu(i,iufrom) - rxu(i,iuprev)
         yfix(iufrom) = ryu(i,iufrom) - ryu(i,iuprev)
         zfix(iufrom) = rzu(i,iufrom) - rzu(i,iuprev)

         lengthb = dsqrt(xfix(iufrom)**2 + yfix(iufrom)**2 
     &        + zfix(iufrom)**2)

         xub = xfix(iufrom) / lengthb
         yub = yfix(iufrom) / lengthb
         zub = zfix(iufrom) / lengthb
         
         call cone(1,xub,yub,zub,dum,dum,dum,dum,dum )
         
c     --- now we must cycle through all sites that we want to be rigid
         do count = 1, ntogrow
            iu = rlist(iw,count)
            
            xfix(iu) = rxu(i,iu) - rxu(i,iufrom)
            yfix(iu) = ryu(i,iu) - ryu(i,iufrom)
            zfix(iu) = rzu(i,iu) - rzu(i,iufrom)
            
            lengtha = dsqrt(xfix(iu)**2 + yfix(iu)**2
     &           + zfix(iu)**2)
                        
            rlength(iu) = lengtha
            
            thetac = -(xfix(iu)*xfix(iufrom) + yfix(iu)*yfix(iufrom)
     &           + zfix(iu)*zfix(iufrom)) / (lengtha*lengthb)
            
            if (abs(thetac).gt.1.0d0) stop 'screwup in rigfix'
            
            bendang(iu) = dacos(thetac)
            
            xub = xfix(iu) / lengtha
            yub = yfix(iu) / lengtha
            zub = zfix(iu) / lengtha
            
c     --- determine the phi value associated with this
            call cone(3,dum,dum,dum,bendang(iu),phia(iu),xub,yub,zub)
            
            lfind(iu) = .true.
            inum(count) = iu
         enddo

         counta = 0
         num = ntogrow
 5       continue
         ja = 0
         do count = 1, num
            iu = inum(count)
         
            do iv = 1, invib(imolty,iu)
               
               ju = ijvib(imolty,iu,iv)

               if (.not.lfind(ju)) then
                  ja = ja + 1
                  counta = counta + 1
                  ilist(counta) = ju

                  ifrom(counta) = iu
                  lfind(ju) = .true.
                  inuma(ja) = ju
                  
                  xfix(ju) = rxu(i,ju) - rxu(i,iufrom)
                  yfix(ju) = ryu(i,ju) - ryu(i,iufrom)
                  zfix(ju) = rzu(i,ju) - rzu(i,iufrom)

                  lengtha = dsqrt(xfix(ju)**2 + yfix(ju)**2
     &                 + zfix(ju)**2)

                  rlength(ju) = lengtha

                  thetac = -(xfix(ju)*xfix(iufrom) + yfix(ju)
     &                 *yfix(iufrom) + zfix(ju)*zfix(iufrom))
     &                 / (lengtha*lengthb)
                  
                  bendang(ju) = dacos(thetac)

                  xub = xfix(ju) / lengtha
                  yub = yfix(ju) / lengtha
                  zub = zfix(ju) / lengtha
c     --- determine the phi value associated with this
                  call cone(3,dum,dum,dum,bendang(ju)
     &                 ,phia(ju),xub,yub,zub)

               endif
            enddo
         enddo
         num = ja
                  
         if (ja.ne.0) then
            do j = 1, ja
               inum(j) = inuma(j)
            enddo
            goto 5
         endif
         nlist = counta
         
c     --- now that we determined the phi values for the old configuration,
c     --- let's set up cone for the new configuration

         if (lnew) then
            xub = xvec(iuprev,iufrom)
            yub = yvec(iuprev,iufrom)
            zub = zvec(iuprev,iufrom)
            
            lengthb = distij(iuprev,iufrom)

            xub = xub / lengthb
            yub = yub / lengthb
            zub = zub / lengthb
            
            call cone(1,xub,yub,zub,dum,dum,dum,dum,dum)
         endif
                 
         do ip = 1, ichoi
            bsum_tor(ip) = 0.0d0
            
            do itor = 1, ichtor
               vdha = 0.0d0
          
               lshit = .false.
               if (lnew.and.ip.eq.17.and.itor.eq.6) then
                  lshit = .true.
               endif

               if (.not.lnew.and.ip.eq.1.and.itor.eq.1) then
                  lshit = .true.
               endif


               do count = 1, ntogrow
                  iu = rlist(iw,count)
                  if (.not.lnew.and.itor.eq.1.and.ip.eq.1) then
                     xub = rxu(i,iu) - rxu(i,iufrom)
                     yub = ryu(i,iu) - ryu(i,iufrom)
                     zub = rzu(i,iu) - rzu(i,iufrom)
                  else
                     phidisp = random() * twopi
                    
                     phi = phia(iu) + phidisp

                     call cone(2,dum,dum,dum,bendang(iu),phi
     &                    ,xub,yub,zub)

                     lengtha = rlength(iu)

                     xub = xub * lengtha
                     yub = yub * lengtha
                     zub = zub * lengtha


                  endif
         
                  xvec(iufrom,iu) = xub
                  yvec(iufrom,iu) = yub
                  zvec(iufrom,iu) = zub

                  xvec(iu,iufrom) = -xub
                  yvec(iu,iufrom) = -yub
                  zvec(iu,iufrom) = -zub

                  do it = 1, intor(imolty,iu)
                     jut2 = ijtor2(imolty,iu,it)
                     jut3 = ijtor3(imolty,iu,it)
                     jut4 = ijtor4(imolty,iu,it)

                     if (jut2.eq.iufrom.and.jut3.eq.iuprev
     &                    .and..not.lplace(imolty,jut4)) then
c     --- check to see if jut4 exists
                        if (.not. lexist(jut4)) then
                           write(6,*) 'iu,jut2,jut3,jut4',iu
     &                          ,jut2,jut3,jut4
                           stop 'trouble, jut4 does not exist in rigfix'
                        endif
                        
                        jttor = ittor(imolty,iu,it)
                        
                        call calctor(iu,jut2,jut3,jut4,jttor,vtor)
                        vdha = vdha + vtor
                     endif
                  enddo
               enddo
            
               vtorsion(itor) = vdha
               phitors(itor) = phidisp
               bf_tor(itor) = dexp(-beta*vdha)
               bsum_tor(ip) = bsum_tor(ip) + bf_tor(itor)
            enddo
                
c     --- pick a torsion at random
            if (lnew .or. ip .ne. 1) then
               ran_tor = random()*bsum_tor(ip)
               bs = 0.0d0
               do itor = 1, ichtor
                  bs = bs + bf_tor(itor)
                  if (ran_tor .lt. bs) then
                     vtgtr(ip) = vtorsion(itor)
                     phidisp = phitors(itor)
                     goto 100
                  endif
               enddo
            else
               vtgtr(ip) = vtorsion(1)
               phidisp = phitors(1)
            endif

 100        continue

            bsum_tor(ip) = bsum_tor(ip) / dble(ichtor)
            
            do count = 1, ntogrow
               iu = rlist(iw,count)
               if (lnew .or. ip.ne.1) then
                  phi = phia(iu) + phidisp
                  
                  call cone(2,dum,dum,dum,bendang(iu),phi
     &                 ,xub,yub,zub)

                  lengtha = rlength(iu)
                  
                  xub = xub * lengtha
                  yub = yub * lengtha
                  zub = zub * lengtha
                  
                  if (lnew) then
                     rxpa(iu,ip) = xub + rxnew(iufrom)
                     rypa(iu,ip) = yub + rynew(iufrom)
                     rzpa(iu,ip) = zub + rznew(iufrom)
                  else
                     rxpa(iu,ip) = xub + rxu(i,iufrom)
                     rypa(iu,ip) = yub + ryu(i,iufrom)
                     rzpa(iu,ip) = zub + rzu(i,iufrom)
                  endif
               else
                  rxpa(iu,ip) = rxu(i,iu)
                  rypa(iu,ip) = ryu(i,iu)
                  rzpa(iu,ip) = rzu(i,iu)
               endif
            enddo
            
                 
c     --- we must determine the positions of the remaining sites
            do counta = 1, nlist
               iu = ilist(counta)
               
               if (lnew.or.ip.ne.1) then
                  phi = phia(iu) + phidisp
                  
                  call cone(2,dum,dum,dum,bendang(iu),phi
     &                 ,xub,yub,zub)
                  
                  lengtha = rlength(iu)

                  xub = xub * lengtha
                  yub = yub * lengtha
                  zub = zub * lengtha
               
                  if (lnew) then
                     rxpa(iu,ip) = xub + rxnew(iufrom)
                     rypa(iu,ip) = yub + rynew(iufrom)
                     rzpa(iu,ip) = zub + rznew(iufrom)
                  else
                     rxpa(iu,ip) = xub + rxu(i,iufrom)
                     rypa(iu,ip) = yub + ryu(i,iufrom)
                     rzpa(iu,ip) = zub + rzu(i,iufrom)
                  endif
               else
                  rxpa(iu,ip) = rxu(i,iu)
                  rypa(iu,ip) = ryu(i,iu)
                  rzpa(iu,ip) = rzu(i,iu)
               endif
               
            enddo

         enddo
    
c     --- now calculate the intramolecular energies of the beads
      
c     --- initialize rosenbluth weight
         do ip = 1, ichoi
            bsuma(ip) = 1.0d0
            vtrya(ip) = 0.0d0
            lovra(ip) = .false.
            vtrintraa(ip) = 0.0d0
            vtrexta(ip)   = 0.0d0
            vtrintera(ip) = 0.0d0
            vtrelecta(ip) =  0.0d0
            vtrewalda(ip) = 0.0d0
            vtrorienta(ip) = 0.0d0
         enddo

         do count = 1, ntogrow
            iu = rlist(iw,count)
            
            do ip = 1, ichoi
               rxp(count,ip) = rxpa(iu,ip)
               ryp(count,ip) = rypa(iu,ip)
               rzp(count,ip) = rzpa(iu,ip)
            enddo
            glist(count) = iu
         enddo


         call boltz(lnew,.false.,ovrlap,i,i,imolty,ibox,ichoi
     &        ,iufrom,ntogrow,glist)
         
         if (ovrlap) then
            lterm = .true.
            return
         endif

c     --- propagate rosenbluth weigth and energies
         do ip = 1, ichoi
            if (lovr(ip)) then
               lovra(ip) = .true.
            endif
            bsuma(ip) = bsuma(ip) * bfac(ip)
            vtrya(ip) = vtrya(ip) + vtry(ip)
            vtrintraa(ip) = vtrintraa(ip) + vtrintra(ip) 
	    vtrexta(ip)   = vtrexta(ip) + vtrext(ip)
	    vtrintera(ip) = vtrintera(ip) + vtrinter(ip)
            vtrelecta(ip) =  vtrelecta(ip) + vtrelect(ip)
            vtrewalda(ip) = vtrewalda(ip) + vtrewald(ip)
            vtrorienta(ip) = vtrorienta(ip) + vtrorient(ip)
         enddo
         
c     --- now run through all other sites
         do counta = 1, nlist
            count = 1
            iu = ilist(counta)

            do ip = 1, ichoi
               rxp(count,ip) = rxpa(iu,ip)
               ryp(count,ip) = rypa(iu,ip)
               rzp(count,ip) = rzpa(iu,ip)
            enddo

            glist(count) = iu

            call boltz(lnew,.false.,ovrlap,i,i,imolty,ibox,ichoi
     &           ,ifrom(counta),1,glist)

            if (ovrlap) then
               lterm = .true.
               return
            endif

c     --- propagate rosenbluth weigth and energies
            do ip = 1, ichoi
               if (lovr(ip)) then
                  lovra(ip) = .true.
               endif
               bsuma(ip) = bsuma(ip) * bfac(ip)
               vtrya(ip) = vtrya(ip) + vtry(ip)
               vtrintraa(ip) = vtrintraa(ip) + vtrintra(ip) 
               vtrexta(ip)   = vtrexta(ip) + vtrext(ip)
               vtrintera(ip) = vtrintera(ip) + vtrinter(ip)
               vtrelecta(ip) =  vtrelecta(ip) + vtrelect(ip)
               vtrewalda(ip) = vtrewalda(ip) + vtrewald(ip)
               vtrorienta(ip) = vtrorienta(ip) + vtrorient(ip)
            enddo         
         enddo

         bsum = 0.0d0
c     --- add up rosenbluth weight
         do ip = 1, ichoi
            if (.not. lovra(ip)) then
               bsum = bsum + bsuma(ip) * bsum_tor(ip)
            endif
         enddo

         if (lnew) then
            wrig = wrig * bsum
            
            if (wrig .lt. softlog) then
               lterm = .true.
               return
            endif

            rbf = bsum * random()
            bs = 0.0d0
            do ip = 1, ichoi
               if ( .not. lovra(ip) ) then
                  bs = bs + bsuma(ip) * bsum_tor(ip)
                  if (rbf .lt. bs ) then
                     iwalk = ip
                     goto 120
                  endif
               endif
            enddo
 120        continue
         else
            wrig = wrig * bsum
            if (wrig .lt. softlog) then
               lterm = .true.
               write(6,*) 'RIGFIX OLD REJECTED'
               return
            endif
         endif

c     --- now we must add up energies and record new positions
         if (lnew) then
            vnewt = vnewt + vtrya(iwalk) + vtgtr(iwalk)  
            vnewtg = vnewtg + vtgtr(iwalk)
            vnewext   = vnewext   + vtrext(iwalk)
            vnewintra = vnewintra + vtrintraa(iwalk)
            vnewinter = vnewinter + vtrintera(iwalk)
            vnewelect = vnewelect + vtrelecta(iwalk)
            vnewewald = vnewewald + vtrewalda(iwalk)
            vneworient = vneworient + vtrorienta(iwalk)
         else
            voldt = voldt + vtrya(1) + vtgtr(1)  
            voldtg = voldtg + vtgtr(1)
            voldext   = voldext   + vtrext(1)
            voldintra = voldintra + vtrintraa(1)
            voldinter = voldinter + vtrintera(1)
            voldelect = voldelect + vtrelecta(1)
            voldewald = voldewald + vtrewalda(1)
            voldorient = voldorient + vtrorienta(1)
         endif
         
         do count = 1, ntogrow
            iu = rlist(iw,count)
            
            lexist(iu) = .true.

            if (lnew) then
               rxnew(iu) = rxpa(iu,iwalk)
               rynew(iu) = rypa(iu,iwalk)
               rznew(iu) = rzpa(iu,iwalk)
            endif
            ju = iufrom

            if (lnew) then
               xvec(iu,ju) = rxnew(ju) - rxnew(iu)
               yvec(iu,ju) = rynew(ju) - rynew(iu)
               zvec(iu,ju) = rznew(ju) - rznew(iu)
            else
               xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
               yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
               zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
            endif
            
            distij(iu,ju) = dsqrt(xvec(iu,ju)**2
     &           + yvec(iu,ju)**2 + zvec(iu,ju)**2)
            
            distij(ju,iu) = distij(iu,ju)
            
            xvec(ju,iu) = -xvec(iu,ju)
            yvec(ju,iu) = -yvec(iu,ju)
            zvec(ju,iu) = -zvec(iu,ju)
            
         enddo

         do counta = 1, nlist
            iu = ilist(counta)

            lexist(iu) = .true.

            if (lnew) then
               rxnew(iu) = rxpa(iu,iwalk)
               rynew(iu) = rypa(iu,iwalk)
               rznew(iu) = rzpa(iu,iwalk)
            endif

            ju = ifrom(counta)


            if (lnew) then
               xvec(iu,ju) = rxnew(ju) - rxnew(iu)
               yvec(iu,ju) = rynew(ju) - rynew(iu)
               zvec(iu,ju) = rznew(ju) - rznew(iu)
            else
               xvec(iu,ju) = rxu(i,ju) - rxu(i,iu)
               yvec(iu,ju) = ryu(i,ju) - ryu(i,iu)
               zvec(iu,ju) = rzu(i,ju) - rzu(i,iu)
            endif
            
            distij(iu,ju) = dsqrt(xvec(iu,ju)**2
     &           + yvec(iu,ju)**2 + zvec(iu,ju)**2)
            
            distij(ju,iu) = distij(iu,ju)
            
            xvec(ju,iu) = -xvec(iu,ju)
            yvec(ju,iu) = -yvec(iu,ju)
            zvec(ju,iu) = -zvec(iu,ju)
            
         enddo
      enddo

c      write(6,*) 'END RIGFIX'
        
      return
      end





