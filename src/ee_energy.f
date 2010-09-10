      subroutine ee_energy ( i,imolty, v, vintra, vinter,vext
     &     ,velect,vewald,vtail,flagon,ibox, istart, iend,lljii,ovrlap
     &     ,ltors,vtors,lcharge_table,lfavor)

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
      include 'cell.inc'
      include 'external.inc'
      include 'connect.inc'
      include 'ewaldsum.inc'
      include 'fepsi.inc'
      include 'qqlist.inc'
      include 'clusterbias.inc'
      include 'nsix.inc'
      include 'eepar.inc'

      logical::lqimol,lqjmol,lexplt,lcoulo,lfavor,lij2,liji,lqchgi
      logical::lljii,ovrlap,ltors,lcharge_table,lt,lfound

      integer::growii,growjj,k,cellinc,jcell,ic,nmole
      integer::i,ibox, istart, iend,ii,ntii,flagon,jjj,iii
     +       ,j,jj,ntjj,ntij,ntj,imolty,jmolty,ncell
      integer::iivib,jjtor,ip1,ip2,ip3,it,nchp2,acellinc,kmolty

      real(8)::ljsami,ljpsur,ljmuir,v,vintra, vinter,vext 
     +                ,rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq
     +                ,ryuij,rzuij,sr2,sr6,rij,rijsq,dzui,dz3,dz12
     +                ,exgrph,exsami,exmuir,exzeo,vtors,exsix,velect
     +                ,vewald,mmff,rbcut,ninesix,vharo,genlj
      real(8)::erfunc,qave,rho,vol,vtail
      real(8)::xvec,yvec,zvec,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2
     &     ,daa1,da1a2,dot,thetac,vtorso,coru
      real(8)::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,epsilon2,sigma2

      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)
      dimension lcoulo(numax,numax),cellinc(cmax),jcell(nmax)
      dimension acellinc(numax,27)

C --------------------------------------------------------------------

c      write(iou,*) 'start ENERGY'
      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut  = rcut(ibox)
 
      if (ldual) rcinsq = rcutin*rcutin

c      rminsq = rmin * rmin
      v = 0.0d0
      vinter = 0.0d0
      vintra = 0.0d0
      vext = 0.0d0
      velect = 0.0d0
      vewald = 0.0d0
      vtail = 0.0d0
      vtors = 0.0d0
      vharo = 0.0d0
      ovrlap = .false.
      if ( istart .eq. 1 .and. flagon .eq. 2) then
         neigh_icnt = 0
      endif

C *******************************
C *** INTERCHAIN INTERACTIONS ***
C *******************************
      lqimol = lelect(imolty)

      if (nugrow(imolty) .eq. nunit(imolty)) then
         lexplt = .false.
      else
         lexplt = .true.
         growii = nugrow(imolty)
      endif

      if ( lcutcm .or. lfavor ) then
c --- calculate the center of mass of chain i and give it a dummy #
         nchp2 = nchain + 2
         do ii = 1,nunit(imolty)
            rxu(nchp2,ii) = rxuion(ii,flagon) 
            ryu(nchp2,ii) = ryuion(ii,flagon) 
            rzu(nchp2,ii) = rzuion(ii,flagon) 
c	write(iou,*) ii,flagon,rxu(nchp2,ii),ryu(nchp2,ii),rzu(nchp2,ii)
         enddo
         nboxi(nchp2) = ibox
         moltyp(nchp2) = imolty
         call ctrmas(.false.,ibox,nchp2,9)
         xcmi = xcm(nchp2)
         ycmi = ycm(nchp2)
         zcmi = zcm(nchp2)
         rcmi = rcmu(nchp2)
c         write(iou,*) 'rcmi:',rcmi
      else
         lij2 = .true.
      endif

      if (licell.and.(ibox.eq.boxlink)) then
         do ii = istart, iend

            rxui = rxuion(ii,flagon)
            ryui = ryuion(ii,flagon)
            rzui = rzuion(ii,flagon)

c     --- check perodic boundaries
            if (rxui.gt.boxlx(ibox)) then
               rxui = rxui - boxlx(ibox)
            elseif (rxui.lt.0) then
               rxui = rxui + boxlx(ibox)
            endif

            if (ryui.gt.boxly(ibox)) then
               ryui = ryui - boxly(ibox)
            elseif (ryui.lt.0) then
               ryui = ryui + boxly(ibox)
            endif

            if (rzui.gt.boxlz(ibox)) then
               rzui = rzui - boxlz(ibox)
            elseif (rzui.lt.0) then
               rzui = rzui + boxlz(ibox)
            endif
                        
            call linkcell(3,i,rxui,ryui,rzui,cellinc)
            
            do j = 1, 27
               acellinc(ii,j) = cellinc(j)
            enddo
         enddo
                  
         ncell = 0
         
         do j = 1, 27
            ncell = ncell + 1
            cellinc(j) = acellinc(1,j)
         enddo

         if (abs(iend-istart).gt.0) then
            do ii = istart+1, iend
               
               do j = 1, 27
                  
                  ic = acellinc(ii,j)

                  lfound = .false.
                  do jj = 1, ncell

                     if (ic.eq.cellinc(jj)) then
                        lfound = .true.
                        goto 92
                     endif
                                          
                  enddo
 92               continue
                  
                  if (.not.lfound) then
                     ncell = ncell + 1
                     cellinc(ncell) = ic
                  endif
               enddo
            enddo
         endif
         
c         lt = .true.
         nmole = 0
         do j = 1, ncell
            ic = cellinc(j)
            
            do k = 1, nicell(ic)

               nmole = nmole + 1
               jcell(nmole) = iucell(ic,k)

c *** what is this??? always compute interactions with chain number 1?
c solute? lt = solute? removing...
c               if (jcell(nmole).eq.1) then
c                  lt = .false.
c               endif
                  
            enddo
         enddo

c         if (lt) then
c            nmole = nmole + 1
c            jcell(nmole) = 1
c         endif
      else
         nmole = nchain
      endif

c     --- loop over all chains except i - not for grand can. with ibox=2 !
      if (.not.(lgrand.and.(ibox.eq.2))) then    
         do 98 k = 1, nmole

            if (licell.and.(ibox.eq.boxlink)) then
               j = jcell(k)
            else
               j = k
            endif
            
            jmolty = moltyp(j)
            lqjmol = lelect(jmolty)
            growjj = nugrow(jmolty)

c ### check for simulation box ###
            if ( ( ibox .eq. nboxi(j) ) .and. (i .ne. j )) then

               if ( lneigh ) then
                  if ( .not. lnn(j,i) ) goto 98
               endif

               if (lcutcm .or. lfavor) then
c              --- check if ctrmas within rcmsq
                  rxuij = xcmi - xcm(j)
                  ryuij = ycmi - ycm(j)
                  rzuij = zcmi - zcm(j)
c              --- minimum image the ctrmas pair separations ***
                  if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  rcm = rbcut + rcmi + rcmu(j)
                  rcmsq = rcm*rcm
                  if ( lfavor ) then
                     favor(j) = (rminsq/rijsq)**2*5.0d0
                     favor2(j) = rminsq/rijsq
                  endif
                  if ( rijsq .gt. rcmsq .and. lcutcm) then 
                     if ( lqimol .and. lqjmol .and. lchgall ) then
                        lij2 = .false.
                        goto 108
                     else
                        goto 98
                     endif
                  else
                     lij2 = .true.
                  endif
               endif

               if ( lcharge_table .and. (.not. lchgall) ) then
c --- called from CBMC and must set up charge-interaction table ---
                  do ii = 1,nugrow(imolty)
                     do jj = 1,nugrow(jmolty)
                        iii = leaderq(imolty,ii)
                        jjj = leaderq(jmolty,jj)
                        if ( iii .eq. ii .and. jjj .eq. jj ) then
                           rxuij = rxuion(ii,flagon) - rxu(j,jj)
                           ryuij = ryuion(ii,flagon) - ryu(j,jj)
                           rzuij = rzuion(ii,flagon) - rzu(j,jj)
                           if ( lpbc ) 
     &                          call mimage(rxuij,ryuij,rzuij,ibox)
                           rijsq = rxuij*rxuij + ryuij*ryuij 
     +                          + rzuij*rzuij
                           if ((rijsq .lt. rcutsq) .or. lijall) then
                              lcoulo(ii,jj) = .true.
                           else
                              lcoulo(ii,jj) = .false.
                           endif
                        endif
                     enddo
                  enddo
               endif

c              --- loop over all beads ii of chain i 
 108           do ii = istart, iend
                  ntii = ntype(imolty,ii)
                  liji = lij(ntii)
                  lqchgi = lqchg(ntii)
                  rxui = rxuion(ii,flagon)
                  ryui = ryuion(ii,flagon)
                  rzui = rzuion(ii,flagon)
                  
c                 --- loop over all beads jj of chain j 
                  do 97 jj = 1, nunit(jmolty) 
c                    --- check exclusion table
                     if ( lexclu(imolty,ii,jmolty,jj) ) goto 97
                     
                     ntjj = ntype(jmolty,jj)
                     if ( lij2 ) then
                        if ( (.not. (liji .and. lij(ntjj))) 
     &                       .and. 
     &                       (.not. (lqchgi .and. lqchg(ntjj)))) 
     &                       goto 97
                     else
                        if (.not. (lqchgi .and. lqchg(ntjj)))
     &                       goto 97
                     endif
                     if ( lexpsix .or. lmmff ) then
                        ntij = (ntii+ntjj)/2
                     elseif (lninesix) then
                        ntij = (ntii-1)*nxatom + ntjj
                     elseif (lgenlj) then
                        ntij = (ntii-1)*nntype + ntjj 
                     else
                        ntij = (ntii-1)*nntype + ntjj
                     endif
                     rminsq = rminee(ntij)*rminee(ntij)
                     
                     rxuij = rxui - rxu(j,jj)
                     ryuij = ryui - ryu(j,jj)
                     rzuij = rzui - rzu(j,jj)
                   
c *** minimum image the pair separations ***
                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox)
                   
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   
                     if ( rijsq .lt. rminsq .and. .not. 
     &                    (lexpand(imolty) .or. lexpand(jmolty))) then
                        ovrlap = .true.
c	if (leemove) then
c                        write(iou,*) 'inter ovrlap:',i,j
c                        write(iou,*) 'i xyz',rxui,ryui,rzui
c                        write(iou,*) 'j xyz',rxu(j,jj),ryu(j,jj),rzu(j,jj) 
c                        write(iou,*) 'ii:',ii,'jj:',jj
c                        write(iou,*) 'distance', dsqrt(rijsq)
c	endif
                        return
                     endif
                     if ( (rijsq .lt. rcutsq) .or. lijall) then
                        if ( lsami ) then
                           vinter = vinter + ljsami(rijsq,ntij)
                        elseif (lexpsix) then
                           vinter = vinter + exsix(rijsq,ntij)
                        elseif (lmmff) then
                           vinter = vinter + mmff(rijsq,ntij)
                        elseif (lninesix) then
                           vinter = vinter + ninesix(rijsq,ntij)
                        elseif (lgenlj) then
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2=epsij(ntij)
                           vinter = vinter + genlj(rijsq,sr2,epsilon2)
                        elseif ( lmuir ) then
                           vinter = vinter + ljmuir(rijsq,ntij)
                        elseif ( lpsurf ) then
                           vinter = vinter + ljpsur(rijsq,ntij)
                        else if (lshift) then
                           sr2 = sig2ij(ntij) / rijsq
                           sr6 = sr2 * sr2 * sr2
                           vinter = vinter + sr6*(sr6-1.0d0)
     &                          *epsij(ntij)-ecut(ntij) 
                        elseif ( lij(ntii) .and. lij(ntjj) ) then
                           if ( lfepsi ) then
                              sr6 = rijsq*rijsq*rijsq
                              if ( (.not. lqchg(ntii)) .and. 
     &                             (.not. lqchg(ntjj)) ) then
                                 if ( nunit(imolty) .eq. 4 ) then
c *** TIP-4P structure (temperary use ???)
                                    qave = (qquion(4,flagon)
     &                                   +qqu(j,4))/2.0d0
                                 else
                                    qave=(qquion(4,flagon)
     &                                   +qquion(5,flagon)
     &                                   +qqu(j,4)+qqu(j,5))*0.85d0
                                 endif
                              else
                                 qave = (qquion(ii,flagon)
     &                                +qqu(j,jj))/2.0d0
                              endif
                              if ( lexpand(imolty) 
     &                             .and. lexpand(jmolty)) then
                                 epsilon2=dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif (lexpand(imolty)) then
                                 epsilon2=dsqrt(epsilon(imolty,ii)
     &                                *epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 epsilon2=dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 epsilon2 = epsij(ntij)
                              endif
                              vinter = vinter + 
     &                             ((aslope*(qave-a0)*(qave-a0)
     &                             +ashift)/sr6 - (bslope*(qave-
     &                             b0)*(qave-b0)+bshift))/
     &                             sr6*epsilon2
                           else
                              if ( lexpand(imolty) 
     &                             .and. lexpand(jmolty)) then
                                 sigma2=(sigma(imolty,ii)+
     &                                sigma(jmolty,jj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif ( lexpand(imolty) ) then
                                 sigma2=(sigma(imolty,ii)+
     &                                sigi(ntjj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsilon(imolty,ii)*
     &                                epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 sigma2=(sigma(jmolty,jj)+
     &                                sigi(ntii))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2=dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2 = epsij(ntij)
                              endif
                              sr6 = sr2 * sr2 * sr2
                              vinter = vinter 
     &                             + sr6*(sr6-1.0d0)*epsilon2
                           endif
                           
                        endif
                     endif

c                    --- electrostatics
                     if ( lchgall .and. lqchg(ntii) 
     &                    .and. lqchg(ntjj) ) then
                        if ( lewald ) then
                           velect = velect + qquion(ii,flagon)*qqu(j,jj)
     &                          *erfunc(calp(ibox)*dsqrt(rijsq))
     &                          /dsqrt(rijsq)
                        else
                           velect = velect + qquion(ii,flagon)
     &                          *qqu(j,jj)/dsqrt(rijsq)
                        endif
                     elseif ( lqimol .and. lqjmol .and. lqchg(ntii) 
     &                       .and. lqchg(ntjj) ) then

                        iii = leaderq(imolty,ii)
                        jjj = leaderq(jmolty,jj)

                        if ( ii .eq. iii .and. jj .eq. jjj ) then
c --- set up the charge-interaction table
                           if ( rijsq .lt. rcutsq ) then
                              lcoulo(ii,jj) = .true.
                           else
                              lcoulo(ii,jj) = .false.
                           endif
                        endif

                        if ( lcoulo(iii,jjj) ) then
                           if ( lewald ) then
                              velect = velect + qquion(ii,flagon)*
     &                             qqu(j,jj)*erfunc(calp(ibox)*
     &                             dsqrt(rijsq))/dsqrt(rijsq)
                           else
                              velect = velect + qquion(ii,flagon)
     &                             *qqu(j,jj)/dsqrt(rijsq)
                           endif
                        endif
                     endif

                     if ( lneighbor .and. ii .eq. 1 .and. 
     &                    jj .eq. 1 .and. flagon .eq. 2
     &                    .and. rijsq .lt. rbsmax**2 
     &                    .and. rijsq .gt. rbsmin**2) then
c                           neigh_icnt1(jmolty)=neigh_icnt1(jmolty)+1
c                           neighi1(neigh_icnt1(jmolty),jmolty)=j
                        neigh_icnt=neigh_icnt+1
                        neighi(neigh_icnt)=j
                     endif

c *** additional repulsion between negative charge on the aromatic
c *** ring and the positive H in e.g. phenol

 97               continue
               enddo
           endif
 98      continue
      endif

      if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff
     & .and. .not. lgenlj .and. .not. lninesix ) then
         vinter = 4.0d0 * vinter
      endif
      vinter = vinter+vharo

c ################################################################

c * the intramolecular van der waals and ewald terms have to be calculated 
c for the explicit atom placement models
C *******************************
C *** INTRACHAIN INTERACTIONS ***
C *******************************
      
c *** for expanded ensemble

   
c --- calculate intramolecular energy correction for chain i 
      do ii = istart, iend

         ntii = ntype(imolty,ii)
         rxui = rxuion(ii,flagon)
         ryui = ryuion(ii,flagon)
         rzui = rzuion(ii,flagon)
            
         do jj = 1,ii-1
               
            ntjj = ntype(imolty,jj)
            if ( lexpsix .or. lmmff ) then
               ntij = (ntii+ntjj)/2
            elseif (lninesix) then
               ntij = (ntii-1)*nxatom + ntjj
            elseif (lgenlj) then
               ntij = (ntii-1)*nntype + ntjj
            else
               ntij = (ntii-1)*nntype + ntjj
            endif
            rminsq = rminee(ntij)*rminee(ntij)
            
            rxuij = rxuion(ii,flagon) - rxuion(jj,flagon)
            ryuij = ryuion(ii,flagon) - ryuion(jj,flagon)
            rzuij = rzuion(ii,flagon) - rzuion(jj,flagon)
               
            rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  
c * calculation of intramolecular electrostatics
            if ( lqinclu(imolty,ii,jj) ) then
                  
               if ( lchgall .and. lqchg(ntii)
     &              .and. lqchg(ntjj) ) then
                  if ( lewald ) then

c                        write(iou,*) 'energy including ii,jj',ii,jj,
c     &                       qqfact*qquion(ii,flagon)*qquion(jj,flagon)
c     &                       *erfunc(calp(ibox)*dsqrt(rijsq))
c     &                       /dsqrt(rijsq)

c                    * 1,5 and beyond which gets full interaction
                     if (linclu(imolty,ii,jj)) then
                        velect = velect + 
     &                       qquion(ii,flagon)*qquion(jj,flagon)
     &                       *erfunc(calp(ibox)*dsqrt(rijsq))
     &                       /dsqrt(rijsq)
c                    * scale 1,4 interactions by qscale2(imolty,ii,jj)
                     else
                        velect = velect + qscale2(imolty,ii,jj)*
     &                       qquion(ii,flagon)
     &                       *qquion(jj,flagon)*erfunc(calp(ibox)
     &                       *dsqrt(rijsq))/dsqrt(rijsq)
                     endif
                  else
c                    * 1,5 and beyond which gets full interaction
                     if (linclu(imolty,ii,jj)) then
                        velect = velect + qquion(ii,flagon)
     &                       *qquion(jj,flagon)/dsqrt(rijsq)
c                    * scale 1,4 interactions by qscale2(imolty,ii,jj)
                     else
                        velect = velect + qscale2(imolty,ii,jj)*
     &                       qquion(ii,flagon)
     &                       *qquion(jj,flagon)/dsqrt(rijsq)
                     endif
                  endif

               elseif ( lqimol .and. lqchg(ntii)
     &                 .and. lqchg(ntjj) ) then

                  iii = leaderq(imolty,ii)
                  jjj = leaderq(imolty,jj)

                  if ( ii .eq. iii .and. jj .eq. jjj ) then
c --- set up the charge-interaction table
                     if ( rijsq .lt. rcutsq ) then
                        lcoulo(ii,jj) = .true.
                     else
                        lcoulo(ii,jj) = .false.
                     endif
                  endif

c *** set up table for neighboring groups- make sure they interact when 
c *** leaderqs are only 2 bonds apart.
                  if (.not. lqinclu(imolty,iii,jjj)) then
                     lcoulo(iii,jjj)  = .true.
                  endif

                  if ( lcoulo(iii,jjj) ) then
c                    * 1,5 and beyond which gets full interaction
                     if (linclu(imolty,ii,jj)) then

                        if ( lewald ) then
                           velect = velect + qquion(ii,flagon)*
     &                          qquion(jj,flagon)*erfunc(calp(ibox)*
     &                          dsqrt(rijsq))/dsqrt(rijsq)
                        else
                           velect = velect + qquion(ii,flagon)
     &                          *qquion(jj,flagon)/dsqrt(rijsq)
                        endif

c                    * scale 1,4 interactions by qscale2(imolty,ii,jj)
                     else

                        if ( lewald ) then
                           velect = velect + qscale2(imolty,ii,jj)*
     &                          qquion(ii,flagon)*
     &                          qquion(jj,flagon)*erfunc(calp(ibox)*
     &                          dsqrt(rijsq))/dsqrt(rijsq)
                        else
                           velect = velect + qscale2(imolty,ii,jj)*
     &                          qquion(ii,flagon)
     &                          *qquion(jj,flagon)/dsqrt(rijsq)
                        endif

                     endif
                  endif
               endif
               
            endif

c * calculation of other non-bonded interactions
            if ( linclu(imolty,ii,jj) ) then

               if (lljii) then

                  if ( rijsq .lt. rminsq .and. .not. 
     &                 lexpand(imolty)) then
                     ovrlap = .true.
                
c                     write(iou,*) 'intra ovrlap:',ii,jj
                     return
                  elseif ( rijsq .lt. rcutsq .or. lijall) then
                     if ( lsami ) then
                        vintra = vintra + ljsami(rijsq,ntij)
                     elseif (lexpsix) then
                        vintra = vintra + exsix(rijsq,ntij)
                     elseif (lmmff) then
                        vintra = vintra + mmff(rijsq,ntij)
                     elseif (lninesix) then
                        vintra = vintra + ninesix(rijsq,ntij)
                     elseif (lgenlj) then
                        sr2 = sig2ij(ntij) / rijsq
                        epsilon2=epsij(ntij)
                        vintra = vintra + genlj(rijsq,sr2,epsilon2)
                     elseif ( lmuir ) then
                        vintra = vintra + ljmuir(rijsq,ntij)
                     elseif ( lpsurf ) then
                        vintra = vintra + ljpsur(rijsq,ntij)
                     else if (lshift) then
                        sr2 = sig2ij(ntij) / rijsq
                        sr6 = sr2 * sr2 * sr2
                        vintra = vintra + 
     +                       sr6*(sr6-1.0d0)*epsij(ntij)-ecut(ntij) 
                     else
                        if ( lexpand(imolty) ) then
                           sigma2=(sigma(imolty,ii)+
     &                          sigma(imolty,jj))/2.0d0
                           sr2 = sigma2*sigma2/rijsq
                           epsilon2 = dsqrt(epsilon(imolty,ii)
     &                          *epsilon(imolty,jj))
                        else
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2 = epsij(ntij)
                        endif
                        sr6 = sr2 * sr2 * sr2
                        vintra = vintra + sr6*(sr6-1.0d0)
     &                       *epsilon2

c * OH 1-5 interaction
                             if (lainclu(imolty,ii,jj)) then
                                vintra = vintra + 0.25d0 * 
     &                               a15(a15type(imolty,ii,jj)) /
     &                               ((rijsq**2)*(rijsq**2)*(rijsq**2))
                             endif

                     endif
                  endif
               endif

            endif

            if (lewald ) then
c     compute the ewald intramolecular (self and correction) terms for 
c     the interactions of the placed atoms with themselves, and with the
c     rest of their own molecule, if there's no interaction
               rij = dsqrt(rijsq)
c              * these are 1,2 and 1,3
               if (.not. lqinclu(imolty,ii,jj)) then
                  vewald = vewald + qquion(ii,flagon)*qquion(jj,flagon)*
     &                 (erfunc(calp(ibox)*rij)-1.0d0)/rij
c              * 1,4 interaction which we scale by qscale2(imolty,ii,jj)
               elseif (lqinclu(imolty,ii,jj) .and. 
     &                 (.not. linclu(imolty,ii,jj))) then
                  vewald = vewald + (1.0d0-qscale2(imolty,ii,jj))*
     &                 qquion(ii,flagon)*
     &                 qquion(jj,flagon)*
     &                 (erfunc(calp(ibox)*rij)-1.0d0)/rij
               endif
            endif
         enddo
         if ( lewald ) then
c           -- self interaction term, 1.772 is sqrt of pi
            vewald = vewald - qquion(ii,flagon)*qquion(ii,flagon)*
     &           calp(ibox)/1.772453851d0
         endif
      enddo
      if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     & .and. .not. lgenlj  .and. .not. lninesix ) then 
           vintra = 4.0d0 * vintra 
      endif

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
          do kmolty = 1, nmolty
             do jmolty = 1, nmolty
c                rho = ncmt(ibox,jmolty) /
c     &               ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
                if (flagon.eq.1) then
                   rho = ncmt(ibox,jmolty) / vol
                   vtail = vtail +
     &               ncmt(ibox,kmolty) * coru(kmolty,jmolty,rho,ibox)
c                write(iou,*) 'vtail',vtail
                else
                   if (jmolty.eq.ee_moltyp(mstate)) then
                      rho = (ncmt(ibox,jmolty)-1) / vol
                   elseif (jmolty.eq.imolty) then
                      rho = (ncmt(ibox,jmolty)+1) / vol
                   else
                      rho = ncmt(ibox,jmolty) / vol
                   endif
                   if (kmolty.eq.imolty) then
                      vtail = vtail +
     &                  (ncmt(ibox,kmolty)+1) * coru(kmolty,jmolty,rho
     *                   ,ibox)
                   elseif (kmolty.eq.ee_moltyp(mstate)) then
                      vtail = vtail +
     &                  (ncmt(ibox,kmolty)-1) * coru(kmolty,jmolty,rho
     &                    ,ibox)
                   else
                      vtail = vtail +
     &                  ncmt(ibox,kmolty) * coru(kmolty,jmolty,rho,ibox)
                   endif
                endif
             enddo
          enddo
c-----
          vinter = vinter + vtail
c----
       endif

c ################################################################

C ***************************************************************
C *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
C ***************************************************************
      
      if ( ljoe .or. lsami .or. lmuir .or. lexzeo
     +				      .or. lgraphite) then
c ---  not for grand can. with ibox=2 !
         if (.not.(lgrand.and.(ibox.eq.2))) then    
            do 399 j = istart,iend

               ntj = ntype(imolty,j)
 
               if ( ljoe ) then
                  if ( extc12(ntj) .gt. 0.1d0 ) then
                     dzui = rzuion(j,flagon) - extz0(ntj)
                     dz3  = dzui * dzui * dzui
                     dz12 = dz3**4
                     vext = vext + 
     +                    (extc12(ntj)/dz12) - (extc3(ntj)/dz3)
                  endif
               endif
	       
 	       if( lgraphite ) then
	       		ntij = (ntj-1)*nntype + ntsubst
			vext = vext + exgrph(rxuion(j,flagon),
     +		    	       ryuion(j,flagon),rzuion(j,flagon),ntij)
	       endif
	       
               if ( lsami ) vext = vext + exsami(rzuion(j,flagon),ntj)
               if ( lmuir ) vext = vext + exmuir(rzuion(j,flagon),ntj)

               if ( lexzeo ) vext = vext + exzeo(rxuion(j,flagon)
     &              ,ryuion(j,flagon),rzuion(j,flagon),ntj)
               
 399        continue

         endif
      endif

c *********************************************************************
c *** calculation of torsion energy for explicit atom methy groups ****
c *********************************************************************

      if ( ltors ) then
         do ii = 1, nunit(imolty)
            rxui=rxuion(ii,flagon)
            ryui=ryuion(ii,flagon)
            rzui=rzuion(ii,flagon)
            do iivib = 1, invib(imolty,ii)
               jj = ijvib(imolty,ii,iivib)
               xvec(ii,jj) = rxuion(jj,flagon) - rxui
               yvec(ii,jj) = ryuion(jj,flagon) - ryui
               zvec(ii,jj) = rzuion(jj,flagon) - rzui
c              --- account for explicit atoms in opposite direction
               xvec(jj,ii)   = -xvec(ii,jj)
               yvec(jj,ii)   = -yvec(ii,jj)
               zvec(jj,ii)   = -zvec(ii,jj)
            enddo
         enddo

         do j = nugrow(imolty)+1, nunit(imolty)
            do jjtor = 1, intor(imolty,j)
               ip3 = ijtor4(imolty,j,jjtor)
               if ( ip3 .lt. j ) then
                  ip1 = ijtor2(imolty,j,jjtor)
                  ip2 = ijtor3(imolty,j,jjtor)
                  it  = ittor(imolty,j,jjtor)
c*** calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                  xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     +                 zvec(ip1,j) * yvec(ip1,ip2)
                  yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     +                 xvec(ip1,j) * zvec(ip1,ip2)
                  zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     +                 yvec(ip1,j) * xvec(ip1,ip2)
                  xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     +                 zvec(ip1,ip2) * yvec(ip3,ip2)
                  ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     +                 xvec(ip1,ip2) * zvec(ip3,ip2)
                  za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     +                 yvec(ip1,ip2) * xvec(ip3,ip2)
c *** calculate lengths of cross products ***
                  daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                  da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
c *** calculate dot product of cross products ***
                  dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                  thetac = - dot / ( daa1 * da1a2 )
                  vtors = vtors + vtorso( thetac, it )
               endif
            enddo
         enddo
      endif

C ----------------------------------------------------------------------------
      
      velect = velect * qqfact
      vewald = vewald * qqfact

c     note that vintra is only computed when the flag lljii is true
      v = vinter + vext + vintra + velect + vewald 

c      write(iou,*) 'end ENERGY'

      return
      end





