      subroutine sumup( ovrlap, v, vinter,vtail,vintra,vvib,
     +                  vbend,vtg,vext,velect,vflucq,ibox, lvol)

c sumup     
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
c    **                                                               **
c    ** logical ovrlap            true for substantial atom overlap   **
c    *******************************************************************
 
      implicit none

c *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'poten.inc'
      include 'conver.inc' 
      include 'external.inc'
      include 'zeolite.inc'
      include 'connect.inc'
      include 'ewaldsum.inc'
      include 'fepsi.inc'
      include 'inputdata.inc'
      include 'qqlist.inc'
      include 'clusterbias.inc'
      include 'neigh.inc'
      include 'cell.inc'
      include 'nsix.inc'
      include 'peboco.inc'     
 
      logical ovrlap, lvol
      logical lexplt,lqimol,lqjmol,lcoulo,lij2,liji,lqchgi,lexclude
      integer i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij, iunit
     +     , ip1, ip2, ip3,ibox,nmcount,iii,jjj,ip
     +     ,iivib, jjvib, jjben, jjtor, it, ntj,itype,k
      double precision v, vinter, vintra, vtail, vvib, vbend, vtg, vext
     &     ,velect,vflucq,qqii
      double precision rcutsq,rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij
     &       ,rijsq,sr2, sr6, rho, thetac, theta 
     +     ,xaa1, yaa1, zaa1, xa1a2, ya1a2, za1a2, daa1, da1a2, dot
     +     ,vtorso, dzui, dz3, dz12, rhoz
     +     ,mmff,rij,vrecipsum,erfunc
     +     ,rvdwsq,rchgsq,rbcut,ninesix

      double precision sx,sy,sz

c      double precision vtemp
      double precision vsc

      double precision coru,coruz,xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,qave
      double precision ljsami,ljpsur,ljmuir,exsami,exmuir,exzeo,exsix
      double precision exgrph,vintera,velecta,vol
      double precision xvec(numax,numax),yvec(numax,numax)
     +                ,zvec(numax,numax),distij(numax,numax),epsilon2
     +                ,sigma2

      dimension lcoulo(numax,numax),lexclude(nmax)
c --------------------------------------------------------------------
      vintera = 0.0d0
      velecta = 0.0d0

c      write(6,*) 'start SUMUP'
      ovrlap = .false.

      if ( lpbc ) call setpbc (ibox)

      rvdwsq = rcut * rcut
      rchgsq = rcutchg(ibox)*rcutchg(ibox)

c     * rcutsq is not currently used
      if ( rvdwsq .gt. rchgsq .or. lchgall ) then
         rcutsq = rvdwsq
         rbcut = rcut
      else
         rcutsq = rchgsq
         rbcut = rcutchg(ibox) 
      endif

      rminsq = rmin * rmin

      v = 0.0d0
      vinter = 0.0d0
      vintra = 0.0d0
      vtail = 0.0d0
      vtg = 0.0d0
      vbend = 0.0d0
      vvib = 0.0d0
      vext = 0.0d0
      velect = 0.0d0
      vflucq = 0.0d0

c *** check the molecule count ***
      nmcount = 0
      do i = 1, nchain
         if ( nboxi(i) .eq. ibox ) then
            nmcount=nmcount+1
         endif
         neigh_cnt(i) = 0
      enddo
      if ( nmcount .ne. nchbox(ibox) ) then
         write(6,*) 'SUMUP: nmcount ne nchbox', nmcount, nchbox
         stop
      endif
 
c ###############################################################

C *******************************
C *** INTERCHAIN INTERACTIONS ***
C *******************************

c --- loop over all chains i 
c --- not if lgrand and ibox =2
      if (.not.(lgrand.and.ibox.eq.2)) then
       do 100 i = 1, nchain - 1
 
c ### check if i is in relevant box ###
          if ( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             lqimol = lelect(imolty)

             if ( lcutcm .and. lvol ) then
                xcmi = xcm(i)
                ycmi = ycm(i)
                zcmi = zcm(i)
                rcmi = rcmu(i)
             else
                lij2 = .true.
             endif

             if ( nugrow(imolty) .eq. nunit(imolty) ) then
                lexplt = .false.
             else
                lexplt = .true.
             endif
             
c --- loop over all chains j with j>i 
             do 99 j = i + 1, nchain
c ### check for simulation box ###
                if ( nboxi(j) .eq. ibox ) then


                   jmolty = moltyp(j)
                   lqjmol = lelect(jmolty)

                   if (lcutcm .and. lvol ) then
c                     --- check if ctrmas within rcmsq
                      rxuij = xcmi - xcm(j)
                      ryuij = ycmi - ycm(j)
                      rzuij = zcmi - zcm(j)
c                     --- minimum image the ctrmas pair separations
c                      if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
                      if ( lpbc ) then 
ccccccccccccccccccccccccccccccccccccc
                         if (lsolid(ibox) .and. .not. lrect(ibox)) then
c     ******************
c     *** Collin's Code:
c     *** non-rectangular box
                            sx = rxuij*hmati(ibox,1)+ryuij*hmati(ibox
     &                           ,4)  +rzuij*hmati(ibox,7)
                            sy = rxuij*hmati(ibox,2)+ryuij*hmati(ibox
     &                           ,5)  +rzuij*hmati(ibox,8)
                            sz = rxuij*hmati(ibox,3)+ryuij*hmati(ibox
     &                           ,6)  +rzuij*hmati(ibox,9)
                            
                            if ( sx .gt. 0.5d0 ) then
                               sx = sx-1d0
                            elseif ( sx .lt. -0.5d0 ) then
                               sx = sx+1d0
                            endif
                            if ( sy .gt. 0.5d0 ) then
                               sy = sy-1d0
                            elseif ( sy .lt. -0.5d0 ) then
                               sy = sy+1d0
                            endif
                            if ( sz .gt. 0.5d0 ) then
                               sz = sz-1d0
                            elseif ( sz .lt. -0.5d0 ) then
                               sz = sz+1d0
                            endif
                            rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)+
     &                           sz*hmat(ibox,7)
                            ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)+
     &                           sz*hmat(ibox,8)
                            rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)+
     &                           sz*hmat(ibox,9)
                            
                         else
                            
                               if ( lpbcx ) then
                                  if ( lfold ) then
                                     if ( rxuij .gt. hbx ) then
                                        rxuij=rxuij-bx
                                     else
                                        if (rxuij.lt.-hbx) rxuij=rxuij+
     &                                       bx
                                     endif
                                  else
                                     rxuij = rxuij - bx*dint(rxuij*bxi+
     &                                    dsign(0.5d0,rxuij))
                                  endif
                               endif
                               
                               if ( lpbcy ) then
                                  if ( lfold ) then
                                     if ( ryuij .gt. hby ) then
                                        ryuij=ryuij-by
                                     else
                                        if (ryuij.lt.-hby) ryuij=ryuij+
     &                                       by
                                     endif
                                  else
                                     ryuij = ryuij - by*dint(ryuij*byi+
     &                                    dsign(0.5d0,ryuij))
                                  endif
                               endif
                               
                               if ( lpbcz ) then
                                  if ( lfold ) then
                                     if (rzuij.gt.hbz) then
                                        rzuij=rzuij-bz
                                     else
                                        if (rzuij.lt.-hbz) rzuij=rzuij+
     &                                       bz
                                     endif
                                  else
                                     rzuij = rzuij - bz*dint(rzuij*bzi+
     &                                    dsign(0.5d0,rzuij))
                                  endif
                               endif
                            endif
                            
                            
                         endif
                         
cccccccccccccccccccccccccccccccccccc

                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rcm = rbcut + rcmi + rcmu(j)
                      rcmsq = rcm*rcm

c                      if ( lneighbor .and. rcmsq .lt. rbsmax**2 
c     &                     .and. rcmsq .gt. rbsmin**2 ) then
c                         neigh_cnt(i,jmolty)=neigh_cnt(i,jmolty)+1
c                         neighbor(neigh_cnt(i,jmolty),i,jmolty)=j
c                         neigh_cnt(j,imolty)=neigh_cnt(j,imolty)+1
c                         neighbor(neigh_cnt(j,imolty),j,imolty)=i
c                      endif

                      if ( rijsq .gt. rcmsq ) then
                         if ( lqimol .and. lqjmol .and. lchgall ) then
                            lij2 = .false.
                            goto 98
                         else
                            goto 99
                         endif
                      else
                         lij2 = .true.
                      endif
                   endif

 98                do ii = 1,nunit(imolty)
                      ntii = ntype(imolty,ii)
                      liji = lij(ntii)
                      lqchgi = lqchg(ntii)
                      rxui = rxu(i,ii)
                      ryui = ryu(i,ii)
                      rzui = rzu(i,ii)

c --- loop over all beads jj of chain j 
                      do 97 jj = 1, nunit(jmolty) 
c --- check exclusion table
                         if ( lexclu(imolty,ii,jmolty,jj) ) goto 97
                         
                         ntjj = ntype(jmolty,jj)
                         if ( lij2 ) then
                            if ( (.not. (liji .and. lij(ntjj))) 
     &                           .and. (.not. (lqchgi .and. 
     &                           lqchg(ntjj)))) goto 97
                         else
                            if (.not. (lqchgi .and. lqchg(ntjj)))
     &                           goto 97
                         endif
                         if (lexpsix .or. lmmff) then
                            ntij = (ntii+ntjj)/2
                         elseif (lninesix) then
                            ntij = (ntii-1)*nxatom + ntjj
                         else
                            ntij = (ntii-1)*nntype + ntjj
                         endif
                         
                         rxuij = rxui - rxu(j,jj)
                         ryuij = ryui - ryu(j,jj)
                         rzuij = rzui - rzu(j,jj)
c *** minimum image the pair separations ***
c                         if (lpbc) call mimage (rxuij,ryuij,rzuij,ibox)
                         if (lpbc) then 
ccccccccccccccccccccccccccccccccccccc
                            
                            if (lsolid(ibox) .and. .not. lrect(ibox))
     &                           then
c     ******************
c     *** Collin's Code:
c     *** non-rectangular box
                               sx = rxuij*hmati(ibox,1)+ryuij*hmati(ibox
     &                              ,4)  +rzuij*hmati(ibox,7)
                               sy = rxuij*hmati(ibox,2)+ryuij*hmati(ibox
     &                              ,5)  +rzuij*hmati(ibox,8)
                               sz = rxuij*hmati(ibox,3)+ryuij*hmati(ibox
     &                              ,6)  +rzuij*hmati(ibox,9)
                            
                               if ( sx .gt. 0.5d0 ) then
                                  sx = sx-1d0
                               elseif ( sx .lt. -0.5d0 ) then
                                  sx = sx+1d0
                               endif
                               if ( sy .gt. 0.5d0 ) then
                                  sy = sy-1d0
                               elseif ( sy .lt. -0.5d0 ) then
                                  sy = sy+1d0
                               endif
                               if ( sz .gt. 0.5d0 ) then
                                  sz = sz-1d0
                               elseif ( sz .lt. -0.5d0 ) then
                                  sz = sz+1d0
                               endif
                               rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)+
     &                              sz*hmat(ibox,7)
                               ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)+
     &                              sz*hmat(ibox,8)
                               rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)+
     &                              sz*hmat(ibox,9)
                               
                                              
                            else
                               
                               if ( lpbcx ) then
                                  if ( lfold ) then
                                     if ( rxuij .gt. hbx ) then
                                        rxuij=rxuij-bx
                                     else
                                        if (rxuij.lt.-hbx) rxuij=rxuij+
     &                                       bx
                                     endif
                                  else
                                     rxuij = rxuij - bx*dint(rxuij*bxi+
     &                                    dsign(0.5d0,rxuij))
                                  endif
                               endif
                            
                               if ( lpbcy ) then
                                  if ( lfold ) then
                                     if ( ryuij .gt. hby ) then
                                        ryuij=ryuij-by
                                     else
                                        if (ryuij.lt.-hby) ryuij=ryuij+
     &                                       by
                                     endif
                                  else
                                     ryuij = ryuij - by*dint(ryuij*byi+
     &                                    dsign(0.5d0,ryuij))
                                  endif
                               endif
                               
                               if ( lpbcz ) then
                                  if ( lfold ) then
                                     if (rzuij.gt.hbz) then
                                        rzuij=rzuij-bz
                                     else
                                        if (rzuij.lt.-hbz) rzuij=rzuij+
     &                                       bz
                                     endif
                                  else
                                     rzuij = rzuij - bz*dint(rzuij*bzi+
     &                                    dsign(0.5d0,rzuij))
                                  endif
                               endif
                            endif
                            
                            
                         endif


cccccccccccccccccccccccccccccccccccc

                         rijsq = (rxuij*rxuij)+(ryuij*ryuij)
     &                        + (rzuij*rzuij)
                         
c                           if ( i .eq. 12 .and. ii .eq. 6 .and.
c     &                        j .eq. 95 .and. jj .eq. 1 ) then
c                           write(6,*) 'CONTROL CONTROL CONTROL'
c                           write(6,*) 'box',ibox,nboxi(i),nboxi(j)
c                           write(6,*) 'i xyz',rxui,ryui,rzui
c                           write(6,*) 'j xyz',rxu(j,jj),ryu(j,jj),
c     &                          rzu(j,jj)
c                           write(6,*) 'r*uij',rxuij,ryuij,rzuij
c                           write(6,*) 'dist2',rijsq
c                           write(6,*) 'distance', dsqrt(rijsq)
c                        endif

                         if ( rijsq .lt. rminsq .and. .not.
     &                        (lexpand(imolty) .or. 
     &                        lexpand(jmolty))) then
                            if ( .not. lvol ) then
                               write(6,*) 'overlap inter'
                               write(6,*) 'rijsq rminsq', rijsq, rminsq
                               write(6,*) 'i ii', i, ii
                               write(6,*) 'i-pos', rxui,ryui,rzui
                               write(6,*) 'j jj', j, jj
                               write(6,*) 'j-pos', 
     &                              rxu(j,jj),ryu(j,jj),rzu(j,jj)
                            endif
                            ovrlap = .true.
                            return
                         elseif ( rijsq .lt. rvdwsq .or. lijall) then
                            if (llj.and.(.not.(lexpand(imolty).or.
     &                           lexpand(jmolty)))) then
                               if ( lij(ntii) .and. lij(ntjj) ) then
                                  sr2 = sig2ij(ntij) / rijsq
                                  epsilon2=epsij(ntij)        
                                  sr6 = sr2 * sr2 * sr2
                                  vinter = vinter 
     &                                 + sr6*(sr6-1.0d0)*epsilon2
                               endif
                            elseif ( lsami ) then
                               vinter = vinter + ljsami(rijsq,ntij)
                            elseif (lexpsix) then
                               vinter = vinter + exsix(rijsq,ntij)
                            elseif (lmmff) then
                               vinter = vinter + mmff(rijsq,ntij)
                            elseif (lninesix) then
                               vinter = vinter + ninesix(rijsq,ntij) 
                            elseif ( lmuir ) then
                               vinter = vinter + ljmuir(rijsq,ntij)
                            elseif ( lpsurf ) then
                               vinter = vinter + ljpsur(rijsq,ntij)
                            else if (lshift) then
                               sr2 = sig2ij(ntij) / rijsq
                               sr6 = sr2 * sr2 * sr2
                               vinter = vinter + 
     +                           sr6*(sr6-1.0d0)*epsij(ntij)-ecut(ntij) 

                            elseif ( lfepsi ) then
                               if ( lij(ntii) .and. lij(ntjj) ) then
                                  sr6 = rijsq*rijsq*rijsq
                                  
                                  if ( (.not. lqchg(ntii)) .and. 
     &                                 (.not. lqchg(ntjj)) ) then
                                     if ( nunit(imolty) .eq. 4 ) then
c *** TIP-4P structure (temperary use ???)
                                        qave=(qqu(i,4)+qqu(j,4))/2.0d0
                                     else
                                        qave=(qqu(i,4)+qqu(i,5)+qqu(j,4)
     &                                       +qqu(j,5))*0.85d0
                                     endif
                                  else
                                     qave = (qqu(i,ii)+qqu(j,jj))/2.0d0
                                  endif
                                  if ( lexpand(imolty) 
     &                                 .and. lexpand(jmolty)) then
                                     epsilon2=dsqrt(epsilon(imolty,ii)
     &                                    *epsilon(jmolty,jj))
                                  elseif (lexpand(imolty)) then
                                     epsilon2=dsqrt(epsilon(imolty,ii)
     &                                    *epsi(ntjj))
                                  elseif ( lexpand(jmolty) ) then
                                     
                                     epsilon2=dsqrt(epsi(ntii)*
     &                                    epsilon(jmolty,jj))
                                  else
                                     epsilon2=epsij(ntij)
                                  endif
                                  vinter = vinter + 
     &                                 ((aslope*(qave-a0)*(qave-a0)
     &                                 +ashift)/sr6 - (bslope*(qave-
     &                                 b0)*(qave-b0)+bshift))/
     &                                 sr6*epsilon2  
                               
                               endif      
                            else
                               if ( lij(ntii) .and. lij(ntjj) ) then
                                  if ( lexpand(imolty) 
     &                                 .and. lexpand(jmolty)) then
                                     sigma2=(sigma(imolty,ii)+
     &                                    sigma(jmolty,jj))/2.0d0
                                     sr2 = sigma2*sigma2/rijsq
                                     epsilon2=dsqrt(epsilon(imolty,ii)
     &                                    *epsilon(jmolty,jj))
                                  elseif ( lexpand(imolty) ) then
                                     sigma2=(sigma(imolty,ii)+
     &                                    sigi(ntjj))/2.0d0
                                     sr2 = sigma2*sigma2/rijsq
                                     epsilon2=dsqrt(epsilon(imolty,ii)
     &                                    *epsi(ntjj))
                                  elseif ( lexpand(jmolty) ) then
                                     sigma2=(sigma(jmolty,jj)+
     &                                    sigi(ntii))/2.0d0
                                     sr2 = sigma2*sigma2/rijsq
                                     epsilon2=dsqrt(epsi(ntii)*
     &                                    epsilon(jmolty,jj))
                                  else
                                     sr2 = sig2ij(ntij) / rijsq
                                     epsilon2=epsij(ntij)
                                  endif
                                  sr6 = sr2 * sr2 * sr2
                                  vinter = vinter 
     &                                 + sr6*(sr6-1.0d0)*epsilon2
                               endif
                               
                            endif

                         endif

                         if ( lchgall.and. lqchg(ntii) 
     &                        .and. lqchg(ntjj)) then
                            if ( lewald ) then
                               velect = velect + qqu(i,ii)*qqu(j,jj)*
     &                              erfunc(calp(ibox)*dsqrt(rijsq))/
     &                              dsqrt(rijsq)

                            else
                               velect = velect + qqu(i,ii)*qqu(j,jj)/
     &                              dsqrt(rijsq)
                            endif
                         
                         elseif ( lqimol .and. lqjmol .and. 
     &                           lqchgi .and. lqchg(ntjj) ) then

c --- All-Atom charges (charge-group look-up table)

                            iii = leaderq(imolty,ii)
                            jjj = leaderq(jmolty,jj)

                            if ( iii .eq. ii .and. jjj .eq. jj )then
c --- set up the table
                               if ( rijsq .lt. rchgsq ) then
                                  lcoulo(iii,jjj) = .true.
                               else
                                  lcoulo(iii,jjj) = .false.
                               endif
                            endif

                            if ( lcoulo(iii,jjj) ) then
                               if ( lewald ) then

                                  velect = velect + qqu(i,ii)*
     &                                 qqu(j,jj)*erfunc(calp(ibox)*
     &                                 dsqrt(rijsq))/dsqrt(rijsq)
                                  
                               else

                                  velect = velect + qqu(i,ii)
     &                                 *qqu(j,jj)/dsqrt(rijsq)
                               endif
                            endif
                            
                         endif
                         
                         if ( lneighbor .and. ii .eq. 1 .and. 
     &                        jj .eq. 1 .and. rijsq .lt. rbsmax**2 
     &                        .and. rijsq .gt. rbsmin**2 ) then
c                            neigh_cnt(i,jmolty)=neigh_cnt(i,jmolty)+1
c                            neighbor(neigh_cnt(i,jmolty),i,jmolty)=j
c                            neigh_cnt(j,imolty)=neigh_cnt(j,imolty)+1
c                            neighbor(neigh_cnt(j,imolty),j,imolty)=i
                            
                            neigh_cnt(i)=neigh_cnt(i)+1
                            neighbor(neigh_cnt(i),i)=j
                            neigh_cnt(j)=neigh_cnt(j)+1
                            neighbor(neigh_cnt(j),j)=i
                         endif
                         
 97                   continue

                   enddo
                endif
 99          continue
          endif

 100   continue
       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &      .and. .not. lninesix ) then 
              vinter = 4.0d0 * vinter
       endif
c       write(6,*) 'ltailc',ltailc
       if (ltailc) then
c--     add tail corrections for the Lennard-Jones energy
          if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
             vol = (hmat(ibox,1) * (hmat(ibox,5) * hmat(ibox,9) -
     &            hmat(ibox,8) * hmat(ibox,6)) + hmat(ibox,4)
     &            * (hmat(ibox,8) * hmat(ibox,3) - hmat(ibox,2)
     &            * hmat(ibox,9)) + hmat(ibox,7) * (hmat(ibox,2)
     &            * hmat(ibox,6) - hmat(ibox,5)*hmat(ibox,3)))
          else
             vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
          endif
          do imolty = 1, nmolty
             do jmolty = 1, nmolty
c                rho = ncmt(ibox,jmolty) / 
c     &               ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
                rho = ncmt(ibox,jmolty) / vol
                vtail = vtail + 
     &               ncmt(ibox,imolty) * coru(imolty,jmolty,rho)
c                write(6,*) 'vtail',vtail
             enddo
          enddo
c-----
          vinter = vinter + vtail
c----
       endif
      endif

c$$$      write(6,*)
c$$$      write(6,*) '+++++++'
c$$$      vtemp = velect
c$$$      write(6,*) 'direct space part:',velect*qqfact

      if ( lewald ) then
         call recipsum(ibox,vrecipsum)
c---- update self terms and correction terms
c         sself = 0.0d0
c         correct = 0.0d0
c * combine to reduce numerical error
         vsc = 0.0d0
         if ( ldielect ) then
            call dipole(ibox,0)
         endif
         do i = 1,nchain
            if (nboxi(i) .eq. ibox) then
               imolty = moltyp(i)
               do ii = 1,nunit(imolty)
c * correct but numerically inprecise
c                  sself = sself + qqu(i,ii)*qqu(i,ii)
c * 1.772.. is the square root of pi
                  vsc = vsc - qqu(i,ii)*qqu(i,ii)*
     &                 calp(ibox)/1.772453851d0

                  do jj = ii+1,nunit(imolty)

c * correct should only be calculated if ii and jj should NOT interact,
c * so only calculating it if lqinclu is false
c                    * this part is 1,2 and 1,3
                     if (.not. lqinclu(imolty,ii,jj)) then
                        rxuij = rxu(i,ii) - rxu(i,jj)
                        ryuij = ryu(i,ii) - ryu(i,jj)
                        rzuij = rzu(i,ii) - rzu(i,jj)
                        rij = dsqrt(rxuij*rxuij + ryuij*ryuij + 
     &                       rzuij*rzuij)
c                        correct = correct + qqu(i,ii)*qqu(i,jj)*
c     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
                        vsc = vsc + qqu(i,ii)*qqu(i,jj)*
     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij

c                    * this part is 1,4 which we want to scale by qscale
c                    * qscale = 0 means no 1,4 charge interaction
                     elseif (lqinclu(imolty,ii,jj) .and. 
     &                       (.not. linclu(imolty,ii,jj))) then
                        rxuij = rxu(i,ii) - rxu(i,jj)
                        ryuij = ryu(i,ii) - ryu(i,jj)
                        rzuij = rzu(i,ii) - rzu(i,jj)
                        rij = dsqrt(rxuij*rxuij + ryuij*ryuij + 
     &                       rzuij*rzuij)
c                        correct = correct + 
c     &                       (1.0d0 - qscale)*qqu(i,ii)*qqu(i,jj)*
c     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
                        vsc = vsc + 
     &                       (1.0d0 - qscale)*qqu(i,ii)*qqu(i,jj)*
     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
                     endif
                  enddo
               enddo
            endif
         enddo
c            vdipole = (dipolex*dipolex+dipoley*dipoley+
c     &           dipolez*dipolez)*(2.0d0*onepi)/(3.0d0*
c     &           boxlx(ibox)**3.0d0)
c            write(6,*) dipolex,dipoley,dipolez

c * old correct way but numerically inprecise
c         sself = -sself*calp(ibox)/dsqrt(onepi)
c         velect = velect + sself + correct + vrecipsum/qqfact
         velect = velect + vsc + vrecipsum/qqfact

      endif

c$$$c at this point velect contains all intermolecular charge interactions,
c$$$c plus the ewald self term and intramolecular corrections 
c$$$
c       write(6,*)
c       write(6,*) '== After Inter === velect is:',velect*qqfact
c$$$
c$$$       vtemp = velect

c ################################################################

c * have to recalculate ewald terms if volume changes
      if ( .not. lvol .or. (lvol .and. lewald) ) then

C *******************************
C *** INTRACHAIN INTERACTIONS ***
C *******************************
c         write(6,*) 'starting intrachain'
c --- loop over all chains i 
c         write(6,*) 'nchain',nchain
         do i = 1, nchain

c            lcoulo(1,5) = .false.


c ### check if i is in relevant box ###
c          write(6,*) 'nboxi(i),i,ibox',nboxi(i),i,ibox
            if ( nboxi(i) .eq. ibox ) then

               imolty = moltyp(i)
c             write(6,*) 'imolty',imolty

c             write(6,*) 'nunit(imolty)',nunit(imolty)
               do ii = 1, nunit(imolty)-1

c                write(6,*) 'ntype(imolty,ii),ii',ntype(imolty,ii),ii

                  ntii = ntype(imolty,ii)

                  rxui = rxu(i,ii)
                  ryui = ryu(i,ii)
                  rzui = rzu(i,ii)

                  do jj = ii+1, nunit(imolty)
                  
                   if ( linclu(imolty,ii,jj) .or. 
     &                  lqinclu(imolty,ii,jj)) then

                      ntjj = ntype(imolty,jj)
c                    write(6,*) 'sumup interaction',ii,jj
                      if (lexpsix .or. lmmff) then
                         ntij = (ntii+ntjj)/2
                      elseif (lninesix) then
                         ntij = (ntii-1)*nxatom + ntjj
                      else
                         ntij = (ntii-1)*nntype + ntjj
                      endif

                      rxuij = rxu(i,ii) - rxu(i,jj)
                      ryuij = ryu(i,ii) - ryu(i,jj)
                      rzuij = rzu(i,ii) - rzu(i,jj)

                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     
                      if (linclu(imolty,ii,jj)) then

                       if ( rijsq .lt. rminsq .and. 
     &                      .not. lexpand(imolty)) then
                          if ( .not. lvol ) then
                             write(6,*) 'overlap intra'
                             write(6,*) 'rijsq rminsq', rijsq, rminsq
                             write(6,*) 'i ii', i, ii
                             write(6,*) 'i-pos', rxui,ryui,rzui
                             write(6,*) 'jj', jj
                             write(6,*) 'j-pos'
     +                             ,rxu(i,jj),ryu(i,jj),rzu(i,jj)
                          endif
                          ovrlap = .true.
                          return
                       elseif ( rijsq .lt. rvdwsq .or. lijall) then
                        
                          if (llj.and.(.not.(lexpand(imolty)
     &                         ))) then
                             sr2 = sig2ij(ntij) / rijsq
                             epsilon2=epsij(ntij)
                             sr6 = sr2 * sr2 * sr2
                             vintra = vintra 
     &                            + sr6*(sr6-1.0d0)*epsilon2
                             
                             if (lainclu(imolty,ii,jj)) then
                                vintra = vintra + 0.25d0 * 
     &                               a15(a15type(imolty,ii,jj)) /
     &                               ((rijsq**2)*(rijsq**2)*
     &                               (rijsq**2))
                             endif
                          elseif ( lsami ) then
                             vintra = vintra + ljsami(rijsq,ntij)
                          elseif (lexpsix) then
                             vintra = vintra + exsix(rijsq,ntij)
                          elseif ( lmuir ) then
                             vintra = vintra + ljmuir(rijsq,ntij)
                          elseif (lmmff) then
                             vintra = vintra + mmff(rijsq,ntij)
                          elseif (lninesix) then
                             vintra = vintra + ninesix(rijsq,ntij)
                          elseif ( lpsurf ) then
                             vintra = vintra + ljpsur(rijsq,ntij)
                          else if (lshift) then
                             sr2 = sig2ij(ntij) / rijsq
                             sr6 = sr2 * sr2 * sr2
                             vintra = vintra + 
     +                           sr6*(sr6-1.0d0)*epsij(ntij)-ecut(ntij) 
                          else
                             if ( lexpand(imolty) ) then
                                sigma2=(sigma(imolty,ii)+
     &                               sigma(imolty,jj))/2.0d0
                                sr2 = sigma2*sigma2/rijsq
                                epsilon2 = dsqrt(epsilon(imolty,ii)
     &                               *epsilon(imolty,jj))
                             else
                                sr2 = sig2ij(ntij) / rijsq
                                epsilon2 = epsij(ntij)
                             endif
                             sr6 = sr2 * sr2 * sr2
                             vintra = vintra + sr6*(sr6-1.0d0)
     &                            *epsilon2

c * OH 1-5 interaction
                             if (lainclu(imolty,ii,jj)) then
                                vintra = vintra + 0.25d0 * 
     &                               a15(a15type(imolty,ii,jj)) /
     &                               ((rijsq**2)*(rijsq**2)*(rijsq**2))
                             endif

                          endif
                       endif

                      endif

c * calculate intramolecular charge interaction
                      if ( lchgall .and. lqchg(ntii)
     &                     .and. lqchg(ntjj)) then
                         if ( lewald ) then

                            if (lqinclu(imolty,ii,jj)) then
c                               write(6,*) 'sumup including ii,jj:',ii,jj
c     &                              ,qqfact*qqu(i,ii)*qqu(i,jj)*
c     &                              erfunc(calp(ibox)*dsqrt(rijsq))/
c     &                              dsqrt(rijsq)

c                              * 1,5 and beyond, full interaction
                               if (linclu(imolty,ii,jj)) then
                                  velect = velect + qqu(i,ii)*qqu(i,jj)*
     &                                 erfunc(calp(ibox)*dsqrt(rijsq))/
     &                                 dsqrt(rijsq)
c                              * 1,4 which we scale by qscale
                               else
                                  velect = velect + qscale*qqu(i,ii)*
     &                                 qqu(i,jj)*erfunc(calp(ibox)*
     &                                 dsqrt(rijsq))/dsqrt(rijsq)
                               endif
                            endif

                         else

                            if (lqinclu(imolty,ii,jj)) then
c                              * 1,5 and beyond, full interaction
                               if (linclu(imolty,ii,jj)) then
                                  velect = velect + qqu(i,ii)*qqu(i,jj)/
     &                                 dsqrt(rijsq)
                               else
c                              * 1,4 which we scale by qscale
                                  velect = velect + qscale*qqu(i,ii)*
     &                                 qqu(i,jj)/dsqrt(rijsq)
                               endif
                            endif

                         endif

                      elseif ( lelect(imolty) .and.
     &                        (lqchg(ntii) .and. lqchg(ntjj)) ) then

c --- All-Atom charges (charge-group look-up table)
                         iii = leaderq(imolty,ii)
                         jjj = leaderq(imolty,jj)

                         if ( iii .eq. ii .and. jjj .eq. jj ) then
c --- set up the table for neutral groups
                            if ( rijsq .lt. rchgsq ) then
                               lcoulo(iii,jjj) = .true.
                            else
                               lcoulo(iii,jjj) = .false.
                            endif
                         endif

c *** set up table for neighboring groups- make sure they interact when 
c *** leaderqs are only 2 bonds apart.
                         if (.not. lqinclu(imolty,iii,jjj)) then
                            lcoulo(iii,jjj)  = .true.
                         endif

c                         if (ibox .eq. 2) then
c                            write(6,*) 'ii,iii,jj,jjj',ii,iii,jj,jjj,
c     &                           lcoulo(iii,jjj),dsqrt(rijsq),
c     &                           dsqrt(rchgsq)
c                         endif

                         if ( lcoulo(iii,jjj) ) then

                            if (lqinclu(imolty,ii,jj)) then
                               if (linclu(imolty,ii,jj)) then
c                              * 1,5 and beyond, full interaction 
                                  if ( lewald ) then
                                     velect = velect + qqu(i,ii)*
     &                                    qqu(i,jj)*erfunc(calp(ibox)*
     &                                    dsqrt(rijsq))/dsqrt(rijsq)
                                  else
                                     velect = velect + qqu(i,ii)
     &                                    *qqu(i,jj)/dsqrt(rijsq)
                                  endif
                               else
c                              * 1,4 which we scale by qscale
                                  if ( lewald ) then
                                     velect = velect + qscale*qqu(i,ii)*
     &                                    qqu(i,jj)*erfunc(calp(ibox)*
     &                                    dsqrt(rijsq))/dsqrt(rijsq)
                                  else
                                     velect = velect + qscale*qqu(i,ii)
     &                                    *qqu(i,jj)/dsqrt(rijsq)
                                  endif
                               endif
                            endif

                         endif

                      endif
c end intramolecular charge interaction
c                        write(6,*) 'rijsq',rijsq,ii,jj,
c     &                     epsij(ntij),sig2ij(ntij)
c                        write(6,*) 'intra',sr6*(sr6-1.0d0)*epsij(ntij)
c     &                       ,'vintra',vintra
                   endif
                enddo
             enddo
          endif
       enddo

c$$$       vtemp = velect - vtemp
c$$$
c$$$c       write(6,*) '== Intra Velect ===',vtemp*qqfact
c$$$c       write(6,*) '== After Intra  === velect is:',velect*qqfact

       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &      .and. .not. lninesix ) then
               vintra = 4.0d0 * vintra
       endif
 
c ################################################################

C *************************************
C *** INTRACHAIN FLUCQ INTERACTIONS ***
C *************************************
       do i = 1,nchain
c --- calculate intramolecular flucq energy for chain i 
          if( nboxi(i) .eq. ibox ) then
             imolty = moltyp(i)
             if ( lelect(imolty) ) then

                if ( lflucq(imolty) ) then
                   iunit = nunit(imolty)
                   do ii = 1, iunit
                      
                      ntii = ntype(imolty,ii)
                      qqii = qqu(i,ii)
                      do jj = ii, iunit
                         
                         if ( ii .eq. jj) then
                            vflucq = vflucq + xiq(ntii)*qqii 
     &                           + jayself(ntii)*qqii*qqii
                         else
                            ntjj = ntype(imolty,jj)
                            ntij = (ntii-1)*nntype + ntjj
                            
                            vflucq = vflucq 
     &                           + jayq(ntij)*qqii*qqu(i,jj)
                         endif
                      enddo
                   enddo
                   vflucq = vflucq - fqegp(imolty)
                else
                   vflucq = 0.0d0
                endif
             endif
          endif
       enddo

C **************************************************
C *** CALCULATION OF VIB. + BEND. + TORS. ENERGY ***
C **************************************************

C NOTE here virtual coordinates can be used!!!
c
       do i = 1, nchain
          
          imolty = moltyp(i)

c ### check if i is in relevant box ###
          if ( nboxi(i) .eq. ibox ) then

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
             if ( brvibk(1) .gt. 0.01d0 .or. lninesix) then
                do j = 2, nunit(imolty)
                   do jjvib = 1, invib(imolty,j)
                      ip1 = ijvib(imolty,j,jjvib)
                      it  = itvib(imolty,j,jjvib)
                      if ( ip1 .lt. j ) vvib = vvib + 
     +             brvibk(it) * ( distij(ip1,j) - brvib(it) )**2
                   enddo
                enddo
             endif

c - bending -
c ### molecule with bond bending 
             do j = 3, nunit(imolty)
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

c                      write(6,*) 'ip2,ip1,j',ip2,ip1,j
c                      write(6,*) 'bend energy, theta '
c     &                     ,brbenk(it) * (theta-brben(it))**2,theta
                   endif
                enddo
             enddo

c - torsions -
c ### molecule with dihedral potenials ###
             do j = 4, nunit(imolty)
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

                      vtg = vtg + vtorso( thetac, it )
c                      write(6,*) 'vtg:',vtg,it
                   endif
                enddo
             enddo
          endif
       enddo
      endif
 
c ################################################################

C ***************************************************************
C *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
C ***************************************************************
 
c --- not if lgrand and ibox =2
c      if (.not.(lgrand.and.ibox.eq.2)) then
      if ( ljoe .or. lsami .or. lmuir .or. 
     +               lexzeo .or. lgraphite ) then
         do i = 1, nchain

c ### check if i is in relevant box ###
            if ( nboxi(i) .eq. ibox ) then

               do j = 1, nunit(imolty)

                  ntj = ntype(imolty,j)
 
                  if ( ljoe ) then
                     if ( extc12(ntj) .gt. 0.1d0 ) then
                        dzui = rzu(i,j) - extz0(ntj)
                        dz3  = dzui * dzui * dzui
                        dz12 = dz3**4
                        vext = vext + 
     +                    (extc12(ntj)/dz12) - (extc3(ntj)/dz3)  
                     endif
                  endif
		  
		  if( lgraphite ) then
			ntij = (ntj-1)*nntype + ntsubst
c			vext = vext + exgrph(rxu(i,j),ryu(i,j),
c     &				      rzu(i,j),ntij)
	       	  endif

                  if ( lsami ) 
     &                 vext = vext + exsami(rzu(i,j),ntj)
                  if ( lmuir ) 
     &                 vext = vext + exmuir(rzu(i,j),ntj)
                  if ( lexzeo ) vext = vext + 
     &                 exzeo(rxu(i,j),ryu(i,j),rzu(i,j),ntj)

               enddo
            endif
         enddo
c      endif
      if (ltailc .and. lexzeo) then
c     ---    add tailcorrections for the zeolite
         rhoz=nzeo/(zeorx*zeory*zeorz)
         vext=vext+nchbox(ibox)*coruz(iunit,rhoz)
      endif
      endif
C ----------------------------------------------------------------------------

c$$$      write(6,*) 'self,corr:',
c$$$     &     (velect-vtemp-vrecipsum/qqfact)*qqfact
c      write(6,*) 'vsc, new self cor:',vsc*qqfact
c      write(6,*) 'recip space part :',vrecipsum
c$$$      write(6,*) 'sc and recip:',vsc*qqfact + vrecipsum

      velect = qqfact*velect 
      v = vinter + vintra + vvib + vbend + vtg + vext + velect 
     &     + vflucq

      if ( .not. lvol ) then
         write(6,*)
         write(6,*) 'sumup control'
         write(6,*) 'number of chains', nmcount
         do i = 1, nmolty
            write(6,*) 'number of chains of type',i,ncmt(ibox,i)
         enddo
         write(6,*) 'inter lj energy ', vinter
         write(6,*) 'intra lj energy ', vintra
	 if (ltailc) write(6,*)
     +              'Tail correction ', vtail
     
	 if (ltailc.and.lexzeo) write(6,*)
     +              'Tail corr. zeol', nchbox(ibox)*coruz(iunit,rhoz) 
         write(6,*) 'bond vibration  ', vvib
         write(6,*) 'bond bending    ', vbend
         write(6,*) 'torsional       ', vtg
         write(6,*) 'external surf.  ', vext
         write(6,*) 'coulombic energy', velect
c         write(6,*) 'exact energy    ',
c     +        1.74756*1.67*831.441/3.292796
         write(6,*) 'fluc Q energy   ', vflucq
         write(6,*) 'total energy    ', v
      endif

c      write(6,*) 'end SUMUP'

      return
      end
