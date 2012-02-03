      subroutine sumup( ovrlap, v, vinter,vtail,vintra,vvib,
     +                  vbend,vtg,vext,velect,vflucq,ibox,lvol)
                       

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
      include 'ipswpar.inc'
      include 'eepar.inc'
c kea include for garofalini potential
      include 'garofalini.inc'
      include 'tabulated.inc'
 
      logical ovrlap, lvol 
      logical lexplt,lqimol,lqjmol,lcoulo,lij2,liji,lqchgi,lexclude
      integer i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij, iunit
     +     , ip1, ip2, ip3,ibox,nmcount,iii,jjj,ip
     +     ,iivib, jjvib, jjben, jjtor, it, ntj,itype,k, mmm
      double precision v, vinter, vintra, vtail, vvib, vbend, vtg, vext
     &     ,velect,vflucq,qqii
      double precision rcutsq,rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij
     &       ,rijsq,sr2, sr6, rho, thetac, theta 
     +     ,xaa1, yaa1, zaa1, xa1a2, ya1a2, za1a2, daa1, da1a2, dot
     +     ,vtorso, dzui, dz3, dz12, rhoz,xcc,ycc,zcc,tcc,spltor
     +     ,mmff,rij,vrecipsum,erfunc
     +     ,rbcut,ninesix,vwell, genlj
c tabulated potential variables
      double precision tabulated_vib, tabulated_bend, tabulated_vdW, 
     &     tabulated_elect, rbend, rbendsq

      double precision sx,sy,sz

c      double precision vtemp
      double precision vsc

      double precision coru,coruz,xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,qave
      double precision ljsami,ljpsur,ljmuir,exsami,exmuir,exzeo,exsix
      double precision exgrph,vintera,velecta,vol
      double precision xvec(numax,numax),yvec(numax,numax)
     +                ,zvec(numax,numax),distij(numax,numax),epsilon2
     +                ,sigma2, distij2(numax,numax)
      double precision slitpore, v_elect_field, field
      dimension lcoulo(numax,numax),lexclude(nmax)
c --------------------------------------------------------------------
      vintera = 0.0d0
      velecta = 0.0d0

c      write(iou,*) 'start SUMUP'
      ovrlap = .false.

      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut = rcut(ibox)
      field = Elect_field(ibox)

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
ckea - 3body garofalini term
      v3garo = 0.0d0

c *** check the molecule count ***
      nmcount = 0
      do i = 1, nchain
         if ( nboxi(i) .eq. ibox ) then
            nmcount=nmcount+1
         endif
         neigh_cnt(i) = 0
      enddo
      if ( nmcount .ne. nchbox(ibox) ) then
         write(iou,*) 'SUMUP: nmcount ne nchbox', nmcount, nchbox
         stop
      endif
 
c ###############################################################

C *******************************
C *** INTERCHAIN INTERACTIONS ***
C *******************************

c --- loop over all chains i 
c --- not if lgrand and ibox =2
c --- JLR 11-24-09 don't loop if box is ideal
c      if (.not.(lgrand.and.(ibox.eq.2))) then
      if (.not.(lgrand.and.ibox.eq.2) .and.
     &  .not.lideal(ibox) ) then
c --- END JLR 11-24-09
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
                      if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rij  = dsqrt(rijsq)
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
c Generalized Lennard Jones
                        elseif (lgenlj) then
                           ntij = (ntii-1)*nntype + ntjj
c KEA garofalini
                         elseif (lgaro) then
                            if (ntii.eq.ntjj) then
                               ntij = ntii
                            else
                               ntij = ntii+ntjj+1
                            endif
                         else
                            ntij = (ntii-1)*nntype + ntjj
                         endif
                         if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
                         
                         rxuij = rxui - rxu(j,jj)
                         ryuij = ryui - ryu(j,jj)
                         rzuij = rzui - rzu(j,jj)
c *** minimum image the pair separations ***
                         if (lpbc) call mimage (rxuij,ryuij,rzuij,ibox)

                         rijsq = (rxuij*rxuij)+(ryuij*ryuij)
     &                        + (rzuij*rzuij)
                         rij   = dsqrt(rijsq)
                         
c                           if ( i .eq. 12 .and. ii .eq. 6 .and.
c     &                        j .eq. 95 .and. jj .eq. 1 ) then
c                           write(iou,*) 'CONTROL CONTROL CONTROL'
c                           write(iou,*) 'box',ibox,nboxi(i),nboxi(j)
c                           write(iou,*) 'i xyz',rxui,ryui,rzui
c                           write(iou,*) 'j xyz',rxu(j,jj),ryu(j,jj),
c     &                          rzu(j,jj)
c                           write(iou,*) 'r*uij',rxuij,ryuij,rzuij
c                           write(iou,*) 'dist2',rijsq
c                           write(iou,*) 'distance', dsqrt(rijsq)
c                        endif

                         if ( rijsq .lt. rminsq .and. .not.
     &                        (lexpand(imolty) .or. 
     &                        lexpand(jmolty))) then
                            if ( .not. lvol ) then
                               write(iou,*) 'overlap inter'
                               write(iou,*) 'rijsq rminsq', rijsq,rminsq
                               write(iou,*) 'i ii', i, ii
                               write(iou,*) 'i-pos', rxui,ryui,rzui
                               write(iou,*) 'j jj', j, jj
                               write(iou,*) 'j-pos', 
     &                              rxu(j,jj),ryu(j,jj),rzu(j,jj)
                            endif
                            ovrlap = .true.
                            return
                         elseif ( rijsq .lt. rcutsq .or. lijall) then

                            if (L_vdW_table.and.(.not.(
     &                           lexpand(imolty).or.lexpand(jmolty)))
     &                           )then
                               call lininter_vdW(rij, 
     &                              tabulated_vdW, ntii, ntjj)

                               vinter = vinter + tabulated_vdW

                            elseif (llj.and.(.not.(lexpand(imolty).or.
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
c Generalized Lennard Jones potential
                            elseif (lgenlj) then
                               sr2 = sig2ij(ntij) / rijsq
                               epsilon2=epsij(ntij)
                         vinter = vinter + genlj(rijsq,sr2,epsilon2)
                            elseif ( lmuir ) then
                               vinter = vinter + ljmuir(rijsq,ntij)
                            elseif ( lpsurf ) then
                               vinter = vinter + ljpsur(rijsq,ntij)
c KEA garofalini potential
                            elseif ( lgaro) then
                               vinter = vinter + garofalini(rijsq,ntij
     &                              ,qqu(i,ii),qqu(j,jj),i,j)
                               if(lshift) then
                                  vinter = vinter-ecut(ntij)
                               endif
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

! charge interactions
ckea - skip for garofalini; included in vinter
                         if(lgaro) then
                         elseif ( lchgall.and. lqchg(ntii) 
     &                        .and. lqchg(ntjj)) then
                            if ( lewald ) then
                               velect = velect + qqu(i,ii)*qqu(j,jj)*
     &                              erfunc(calp(ibox)*rij)/
     &                              rij
                            else
                               velect = velect + qqu(i,ii)*qqu(j,jj)/
     &                              rij
                            endif
                         elseif ( lqimol .and. lqjmol .and. 
     &                           lqchgi .and. lqchg(ntjj) ) then

                            if (lewald) then               
                               if (rijsq.lt.rcutsq) then
                                  velect = velect + qqu(i,ii)*
     &                                  qqu(j,jj)*erfunc(calp(ibox)*
     &                                 rij)/rij
                               endif
                            else
c --- All-Atom charges (charge-group look-up table)
                              iii = leaderq(imolty,ii)
                              jjj = leaderq(jmolty,jj)
                              if ( iii .eq. ii .and. jjj .eq. jj )then
c --- set up the table
                                if ( rijsq .lt. rcutsq ) then
                                   lcoulo(iii,jjj) = .true.
                                else
                                   lcoulo(iii,jjj) = .false.
                                endif
                              endif
                              if ( lcoulo(iii,jjj) ) then
                                 if (L_elect_table) then 
                                    call lininter_elect(rij, 
     &                                   tabulated_elect, ntii, ntjj)
                                    velect = velect + qqu(i,ii)*
     &                                   qqu(j,jj)*tabulated_elect
                                 else
                                    velect = velect + qqu(i,ii)
     &                                   *qqu(j,jj)/rij
                                 endif
                              endif
                            endif
                         endif

! End inter-charge loop

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
                         elseif(lgaro) then
                            if((ntij.eq.4.and.rijsq.lt.grijsq(2,1)).or.
     &                           (ntij.eq.6.and.rijsq.lt.grijsq(3,1)))
     &                           then
                               write(64,*) 'neighbor',i,' (',
     &                              neigh_cnt(i)+1,')',j,' (',
     &                              neigh_cnt(j)+1,')'
                               neigh_cnt(i)=neigh_cnt(i)+1
                               neighbor(neigh_cnt(i),i)=j
                               neigh_cnt(j)=neigh_cnt(j)+1
                               neighbor(neigh_cnt(j),j)=i
                               ndij(neigh_cnt(i),i) = rij
                             ndij(neigh_cnt(j),j) = ndij(neigh_cnt(i),i)
                               nxij(neigh_cnt(i),i) = rxuij
                               nyij(neigh_cnt(i),i) = ryuij
                               nzij(neigh_cnt(i),i) = rzuij
                               nxij(neigh_cnt(j),j) = -rxuij
                               nyij(neigh_cnt(j),j) = -ryuij
                               nzij(neigh_cnt(j),j) = -rzuij
                            endif
                         endif
 97                   continue
                   enddo
                endif
 99          continue
          endif

 100   continue
       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     & .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro
     &      .and..not.L_vdW_table) then 
              vinter = 4.0d0 * vinter
       endif

c KEA garofalini 3 body potential
         if (lgaro) then
            call triad
            call vthreebody(v3garo)
         endif


c       write(iou,*) 'ltailc',ltailc
       if (ltailc) then
c--     add tail corrections for the Lennard-Jones energy
          if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
             vol = cell_vol(ibox)
          else
             vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
          endif
          do imolty = 1, nmolty
             do jmolty = 1, nmolty
c                rho = ncmt(ibox,jmolty) / 
c     &               ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
                rho = ncmt(ibox,jmolty) / vol
                vtail = vtail + 
     &               ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
c                write(iou,*) 'vtail',vtail
             enddo
          enddo
c-----
  
          vinter = vinter + vtail
c----
       endif
      endif

c$$$      write(iou,*)
c$$$      write(iou,*) '+++++++'
c$$$      vtemp = velect
c$$$      write(iou,*) 'direct space part:',velect*qqfact

      if ( ldielect ) then
           call dipole(ibox,0)
      endif

      if ( lewald ) then
         call recipsum(ibox,vrecipsum)
c---- update self terms and correction terms
         sself = 0.0d0
         correct = 0.0d0
c * combine to reduce numerical error
c         vsc = 0.0d0
         do i = 1,nchain
            if (nboxi(i) .eq. ibox) then
               imolty = moltyp(i)
               do ii = 1,nunit(imolty)
                  sself = sself + qqu(i,ii)*qqu(i,ii)
c * 1.772.. is the square root of pi
c                  vsc = vsc - qqu(i,ii)*qqu(i,ii)*
c     &                 calp(ibox)/1.772453851d0
                  do jj = ii+1,nunit(imolty)
c * correct should only be calculated if ii and jj should NOT interact,
c * so only calculating it if lqinclu is false
c                    * this part is 1,2 and 1,3
                     if (.not. lqinclu(imolty,ii,jj)) then
                        rxuij = rxu(i,ii) - rxu(i,jj)
                        ryuij = ryu(i,ii) - ryu(i,jj)
                        rzuij = rzu(i,ii) - rzu(i,jj)
c --- JLR 11-17-09  need call to mimage for intrachain
                        if (lpbc) call mimage ( rxuij,ryuij,rzuij,ibox )
c --- END JLR 11-17-09
                        rij = dsqrt(rxuij*rxuij + ryuij*ryuij + 
     &                       rzuij*rzuij)
                        correct = correct + qqu(i,ii)*qqu(i,jj)*
     &                            (erfunc(calp(ibox)*rij)-1.0d0)/rij
c                        vsc = vsc + qqu(i,ii)*qqu(i,jj)*
c     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
                     elseif (lqinclu(imolty,ii,jj)) then
                        rxuij = rxu(i,ii) - rxu(i,jj)
                        ryuij = ryu(i,ii) - ryu(i,jj)
                        rzuij = rzu(i,ii) - rzu(i,jj)
c --- JLR 11-17-09  need call to mimage for intrachain
                        if (lpbc) call mimage ( rxuij,ryuij,rzuij,ibox )
c --- END JLR 11-17-09
                        rij = dsqrt(rxuij*rxuij + ryuij*ryuij + 
     &                       rzuij*rzuij)
                        correct=correct+(1.0d0 - qscale2(imolty,ii,jj))
     &                           *qqu(i,ii)*qqu(i,jj)*
     &                           (erfunc(calp(ibox)*rij)-1.0d0)/rij
c                        vsc = vsc +
c     &                       (1.0d0 - qscale2(imolty,ii,jj))*qqu(i,ii)
c     *                             *qqu(i,jj)*
c     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
                     endif
                  enddo
               enddo
            endif
         enddo
c            vdipole = (dipolex*dipolex+dipoley*dipoley+
c     &           dipolez*dipolez)*(2.0d0*onepi)/(3.0d0*
c     &           boxlx(ibox)**3.0d0)
c            write(iou,*) dipolex,dipoley,dipolez

         sself = -sself*calp(ibox)/dsqrt(onepi)
         velect = velect + sself + correct + vrecipsum/qqfact

C         velect = velect + vsc + vrecipsum/qqfact

      endif

c$$$c at this point velect contains all intermolecular charge interactions,
c$$$c plus the ewald self term and intramolecular corrections 
c$$$
c       write(iou,*)
c       write(iou,*) '== After Inter === velect is:',velect*qqfact
c$$$
c$$$       vtemp = velect

c ################################################################

c * have to recalculate ewald terms if volume changes
      if ( .not. lvol .or. (lvol .and. lewald) ) then

C *******************************
C *** INTRACHAIN INTERACTIONS ***
C *******************************
c         write(iou,*) 'starting intrachain'
c --- loop over all chains i 
c         write(iou,*) 'nchain',nchain
         do i = 1, nchain

c            lcoulo(1,5) = .false.


c ### check if i is in relevant box ###
c          write(iou,*) 'nboxi(i),i,ibox',nboxi(i),i,ibox
            if ( nboxi(i) .eq. ibox ) then

               imolty = moltyp(i)
c             write(iou,*) 'imolty',imolty

c             write(iou,*) 'nunit(imolty)',nunit(imolty)
               do ii = 1, nunit(imolty)-1

c                write(iou,*) 'ntype(imolty,ii),ii',ntype(imolty,ii),ii

                  ntii = ntype(imolty,ii)

                  rxui = rxu(i,ii)
                  ryui = ryu(i,ii)
                  rzui = rzu(i,ii)

                  do jj = ii+1, nunit(imolty)
                  
                   if ( linclu(imolty,ii,jj) .or. 
     &                  lqinclu(imolty,ii,jj)) then

                      ntjj = ntype(imolty,jj)
c                    write(iou,*) 'sumup interaction',ii,jj
                      if (lexpsix .or. lmmff) then
                         ntij = (ntii+ntjj)/2
                      elseif (lninesix) then
                         ntij = (ntii-1)*nxatom + ntjj
                      elseif (lgenlj) then
                         ntij = (ntii-1)*nntype + ntjj
                      else
                         ntij = (ntii-1)*nntype + ntjj
                      endif
                      if (lexpee) rminsq = rminee(ntij)*rminee(ntij)

                      rxuij = rxu(i,ii) - rxu(i,jj)
                      ryuij = ryu(i,ii) - ryu(i,jj)
                      rzuij = rzu(i,ii) - rzu(i,jj)
c --- JLR 11-17-09  need call to mimage for intrachain
                      if (lpbc) call mimage ( rxuij,ryuij,rzuij,ibox )
c --- END JLR 11-17-09
                      rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                      rij  = dsqrt(rijsq) 
                      if (linclu(imolty,ii,jj)) then

                       if ( rijsq .lt. rminsq .and. 
     &                      .not. lexpand(imolty)) then
                          if ( .not. lvol ) then
                             write(iou,*) 'overlap intra'
                             write(iou,*) 'rijsq rminsq', rijsq, rminsq
                             write(iou,*) 'i ii', i, ii
                             write(iou,*) 'i-pos', rxui,ryui,rzui
                             write(iou,*) 'jj', jj
                             write(iou,*) 'j-pos'
     +                             ,rxu(i,jj),ryu(i,jj),rzu(i,jj)
                          endif
                          ovrlap = .true.
                          return
                       elseif ( rijsq .lt. rcutsq .or. lijall) then
                          
                          

                          if (L_vdW_table.or.L_bend_table.and.(.not.
     &                         (lexpand(imolty)))) then

                                do mmm=1,inben(imolty,ii)
                                   if (ijben3(imolty,ii,mmm).eq.jj) then

                                      call lininter_bend(rij,
     &                                     tabulated_bend, 
     &                                     itben(imolty,ii,mmm))
                                      vintra = vintra + tabulated_bend

                                      goto 94
                                   endif
                                enddo

                             call lininter_vdW(rij, tabulated_vdW, 
     &                            ntii, ntjj)
                             vintra = vintra + tabulated_vdW

                          elseif (llj.and.(.not.(lexpand(imolty)
     &                         ))) then
                             sr2 = sig2ij(ntij) / rijsq
                             epsilon2=epsij(ntij)
                             sr6 = sr2 * sr2 * sr2
                             vintra = vintra 
     &                            + sr6*(sr6-1.0d0)*epsilon2
     &				    *ljscale(imolty,ii,jj)
                             
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
                          elseif (lgenlj) then
                             sr2 = sig2ij(ntij) / rijsq
                             epsilon2=epsij(ntij)
                             vintra = vintra + genlj(rijsq,sr2,epsilon2)
                          elseif ( lpsurf ) then
                             vintra = vintra + ljpsur(rijsq,ntij)
                          else if (lshift) then
                             sr2 = sig2ij(ntij) / rijsq
                             sr6 = sr2 * sr2 * sr2
                             vintra = vintra + 
     +                       (sr6*(sr6-1.0d0)*epsij(ntij)-ecut(ntij))
     +				  *ljscale(imolty,ii,jj) 
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
     &                            *epsilon2*ljscale(imolty,ii,jj)

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

 94                   if ( lchgall .and. lqchg(ntii)
     &                     .and. lqchg(ntjj)) then
                         if (lqinclu(imolty,ii,jj)) then
                             if ( lewald ) then
                                 velect = velect + qscale2(imolty,ii,jj)
     &                                   * qqu(i,ii)*qqu(i,jj)*
     &                                 erfunc(calp(ibox)*rij)/
     &                                 rij
                             else
                                  velect = velect + 
     &                              qscale2(imolty,ii,jj)*qqu(i,ii)*
     &                                 qqu(i,jj)/rij
                             endif
                         endif
 
                      elseif ( lelect(imolty) .and.
     &                        (lqchg(ntii) .and. lqchg(ntjj)) ) then
                         if (lewald) then
                             if (lqinclu(imolty,ii,jj).and.
     &                            (rijsq.lt.rcutsq)) then
                                 velect = velect +
     &                              qscale2(imolty,ii,jj)*qqu(i,ii)*
     &                              qqu(i,jj)*erfunc(calp(ibox)*
     &                              rij)/rij
                             endif 
                         else
c --- All-Atom charges (charge-group look-up table)
                            iii = leaderq(imolty,ii)
                            jjj = leaderq(imolty,jj)
                            if ( iii .eq. ii .and. jjj .eq. jj ) then
c --- set up the table for neutral groups
                               if ( rijsq .lt. rcutsq ) then
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
                            if ( lcoulo(iii,jjj) ) then
                               if (lqinclu(imolty,ii,jj)) then
                                  if (L_elect_table) then
c       check this                   if (rij.gt.rmin) then 
                                     call lininter_elect(rij, 
     &                                    tabulated_elect, ntii,ntjj)
                                     velect = velect + 
     &                                    qscale2(imolty,ii,jj)*
     &                                    qqu(i,ii)*qqu(i,jj)*
     &                                    tabulated_elect
c                                     endif
                                  else
                                     velect = velect + 
     &                                    qscale2(imolty,ii,jj)
     &                                    *qqu(i,ii)*qqu(i,jj)/rij
                                  endif
                               endif
                            endif
                         endif
                      endif
                   endif
c end intramolecular charge
                enddo
             enddo
          endif
       enddo

c$$$       vtemp = velect - vtemp
c$$$
c$$$c       write(iou,*) '== Intra Velect ===',vtemp*qqfact
c$$$c       write(iou,*) '== After Intra  === velect is:',velect*qqfact

       if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     & .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro
     &      .and..not.L_vdW_table) then
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
                   distij2(ii,jj) = ( xvec(ii,jj)**2
     +                 + yvec(ii,jj)**2 + zvec(ii,jj)**2 )
                   distij(ii,jj) = dsqrt(distij2(ii,jj))

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
                      if ( ip1. lt. j .and. L_vib_table) then
                         call lininter_vib(distij(ip1,j), 
     &                        tabulated_vib, it)
                         vvib = vvib + tabulated_vib
c                         write(2,*) 'TABULATED VVIB: ', tabulated_vib, 
c     &                        distij(ip1,j), ip1, j
                      endif
                      if ( ip1 .lt. j .and..not.L_vib_table) vvib = vvib
     +                 + brvibk(it) * ( distij(ip1,j) - brvib(it) )**2
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

c$$$                      if (L_bend_table) then
c$$$                         rbendsq=distij2(ip1,j)+distij2(ip1,ip2)
c$$$     &                           -2.0d0*distij(ip1,j)*distij(ip1,ip2)
c$$$     &                           *thetac
c$$$                         rbend = dsqrt(rbendsq)
c$$$                         call lininter_bend(rbend, tabulated_bend, it)
c$$$                         vbend = vbend + tabulated_bend
c$$$                      else
                        vbend = vbend + 
     +                        brbenk(it) * (theta-brben(it))**2
c                      endif

c                        write(iou,*) 'j,ip1,ip2, it',j,ip1,ip2, it
c                        write(iou,*) 'bend energy, theta '
c     &                       ,brbenk(it) * (theta-brben(it))**2,theta

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

                      if (thetac.gt.1.0d0) thetac=1.0d0
                      if (thetac.lt.-1.0d0) thetac=-1.0d0
c     KEA -- added for extending range to +/- 180
c     and for asymmetric potentials
                      if (L_tor_table) then
c     *** calculate cross product of cross products ***
                         xcc = yaa1*za1a2 - zaa1*ya1a2
                         ycc = zaa1*xa1a2 - xaa1*za1a2
                         zcc = xaa1*ya1a2 - yaa1*xa1a2
c     *** calculate scalar triple product ***
                         tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2)
     &                        + zcc*zvec(ip1,ip2)
                         theta = dacos (thetac)

                         if (tcc .lt. 0.0d0) theta = -theta
                         if (L_spline) then
                            call splint(theta,spltor,it)
                         elseif(L_linear) then
                            call lininter(theta,spltor,it)
                         endif

                         vtg = vtg+spltor
                      else
                         vtg = vtg + vtorso( thetac, it )
                      endif
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
c --- for adsorption isotherms, don't calculate energy w/surface
c --- in box 2
      if ( .not. (lslit.and.ibox.eq.2)) then
      if ( ljoe .or. lsami .or. lmuir .or. 
     +               lexzeo .or. lgraphite .or. lslit ) then
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
		  
		  if (lslit) then
	             ntij = (ntj-1)*nntype + ntsubst
c -- calculate interaction with surface at the bottom of the box		  		  
		     vext = vext + slitpore(rzu(i,j),ntij)
c -- calculate interaction with the surface at the top of the box
		     dzui = boxlz(ibox)-rzu(i,j)
		     vext = vext +slitpore(dzui,ntij)
	          endif  
		  
		  if( lgraphite ) then
			ntij = (ntj-1)*nntype + ntsubst
			vext = vext + exgrph(rxu(i,j),ryu(i,j),
     &				      rzu(i,j),ntij)
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
      endif
      if (ltailc .and. lexzeo) then
c     ---    add tailcorrections for the zeolite
         rhoz=nzeo/(zeorx*zeory*zeorz)
         vext=vext+nchbox(ibox)*coruz(iunit,rhoz,ibox)
      endif
      endif

c **********************************************************************
c *** calculation of interaction energy with external electric field ***
c *** added 06/24/07 by KM
c **********************************************************************

      if (lelect_field) then
         do i = 1, nchain
           if (nboxi(i).eq.ibox) then
              if(lelect(moltyp(i))) then
                 do j = 1,nunit(moltyp(i)) 
                    vext = vext + v_elect_field(i,j,rzu(i,j),field)
                 enddo
              endif 
           endif 
         enddo
         vext = vext * eXV_to_K 
      endif  

c --------------------------------------------------------------------
c calculation of additional gaussian potential needed in thermodynamic
c integration in stages b and c
c --------------------------------------------------------------------
      vwell = 0.0d0
                                                                                
      if (lmipsw) then
         do i = 1, nchain
            imolty = moltyp(i)
            if (lwell(imolty)) then
            rxui = xcm(i)
            ryui = ycm(i)
            rzui = zcm(i)
            do j = 1, nwell(imolty)*nunit(imolty)
               k = j-int(j/nunit(imolty))*nunit(imolty)
               if (k.eq.0) k = nunit(imolty)
               rxuij = rxui-rxwell(j,imolty)
               ryuij = ryui-rywell(j,imolty)
               rzuij = rzui-rzwell(j,imolty)
               call mimage(rxuij,ryuij,rzuij,ibox)
               rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
               rcm = rcut(ibox)+rcmu(i)
               rcmsq = rcm*rcm
               if (rijsq.lt.rcmsq) then
               do ii = 1, nunit(imolty)
                  if (awell(ii,k,imolty).lt.1.0d-6) goto 666
                  rxui = rxu(i,ii)
                  ryui = ryu(i,ii)
                  rzui = rzu(i,ii)
                  rxuij = rxui-rxwell(j,imolty)
                  ryuij = ryui-rywell(j,imolty)
                  rzuij = rzui-rzwell(j,imolty)
                  call mimage(rxuij,ryuij,rzuij,ibox)
                  rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                  vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
 666           enddo
               endif
            enddo
            endif
         enddo
      endif

C ----------------------------------------------------------------------------

!      write(iou,*) 'self,corr:',
!     &     (velect-vrecipsum/qqfact)*qqfact
!      write(iou,*) 'vsc, new self cor:',vsc*qqfact
!      write(iou,*) 'recip space part :',vrecipsum
!      write(iou,*) 'sc and recip:',vsc*qqfact + vrecipsum
      
      if (.not.L_elect_table) then
         velect = velect*qqfact
      endif

 
c      velect = qqfact*velect 

      v = vinter + vintra + vext + velect + vflucq + v3garo
     
!      write(iou,*) 'v in sumup',v
                                                                                
      vipsw = v
      vwellipsw = vwell
                                                                                
      if (lstagea) then
         v = (1.0d0-lambdais*(1.0d0-etais))*v
      elseif (lstageb) then
         v = etais*v+lambdais*vwell
      elseif (lstagec) then
         v = (etais+(1.0d0-etais)*lambdais)*v+(1.0d0-lambdais)*vwell
      endif
                                                                                
      v = v + vvib + vbend + vtg

      if ( .not. lvol ) then
         write(iou,*)
         write(iou,*) 'sumup control'
         write(iou,*) 'number of chains', nmcount
         do i = 1, nmolty
            write(iou,*) 'number of chains of type',i,ncmt(ibox,i)
         enddo
         write(iou,*) 'inter lj energy ', vinter
         write(iou,*) 'intra lj energy ', vintra
         if (ltailc) write(iou,*)
     +              'Tail correction ', vtail
     
         if (ltailc.and.lexzeo) write(iou,*)
     +              'Tail corr. zeol', nchbox(ibox)*coruz(iunit,rhoz,
     +                                 ibox) 
         write(iou,*) 'bond vibration  ', vvib
         write(iou,*) 'bond bending    ', vbend
         write(iou,*) 'torsional       ', vtg
         write(iou,*) 'external        ', vext
         write(iou,*) 'coulombic energy', velect
c         write(iou,*) 'exact energy    ',
c     +        1.74756*1.67*831.441/3.292796
         write(iou,*) 'fluc Q energy   ', vflucq
         write(iou,*) 'well energy     ', vwellipsw
         if(lgaro) write(iou,*) '3-body garo     ', v3garo
         write(iou,*) 'total energy    ', v
      endif

c      write(iou,*) 'end SUMUP'

      return
      end
