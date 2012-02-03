      subroutine boltz( lnew,lfirst,ovrlap,i,icharge,imolty,ibox
     &     ,ichoi,iufrom,ntogrow,glist,maxlen)

c boltz
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
c    ** calculates the potential energy and the boltzmann factor      **
c    ** for ichoi trial positions.                                    **
c    ** logical ovrlap            true for walk termination           **
c    *******************************************************************

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'poten.inc'
      include 'neigh.inc'
      include 'cbmc.inc' 
      include 'external.inc'
      include 'connect.inc'
      include 'ewaldsum.inc'
      include 'fepsi.inc'
      include 'qqlist.inc'
      include 'nsix.inc'
      include 'cell.inc'
      include 'peboco.inc'

      logical lnew,ovrlap,lcmno,lfirst,lcompute,lcoulo
      logical lqimol,lqjmol,liji,lqchgi
      integer ichoi,growjj,igrow,count,glist,icharge,cnt,jcell,ic
      integer i,imolty,ibox,ntogrow,itrial,ntii,j,jj,ntjj,ntij
     +       ,iu,jmolty,jjj,iufrom,ii,zz,bdmol_b,cellinc,k,nmole

      double precision ljsami,rminsq,rxui,sr6,ryui,rzui
     +     ,rxuij,ryuij,rzuij,rijsq,sr2,dzui,dz3,dz12
     +     ,exzeo,exsami,exmuir,exgrph,ljpsur,ljmuir,exsix
     +     ,mmff,maxlen,rcm,rcmsq,rcinvdsq,rcinchsq
     +     ,corr,erfunc,rcutmax,rcvdwsq,rcchgsq,ninesix
      double precision vinter,vintra,vext,velect,vewald,qave,
     &     epsilon2,sigma2

      double precision sx,sy,sz

      dimension lcmno(nmax),lcoulo(numax,numax)
      dimension glist(numax),cellinc(27),jcell(nmax)

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
c      write(6,*) 'start BOLTZ'

c      lcompute = .false.

      if ( lpbc ) call setpbc (ibox)

c     --- determine the potential cutoffs
      rcvdwsq = rcut*rcut
      rcchgsq = rcutchg(ibox)*rcutchg(ibox)
      
      if ( ldual ) then
c        --- use rcutin for both types of interactions (except intra)
         rcinvdsq = rcutin*rcutin
         rcinchsq = rcinvdsq
      else
c        --- compute the cutoffs squared for each interaction
         rcinvdsq = rcvdwsq
         rcinchsq = rcchgsq

         if ( lcutcm ) then
c           --- decide which cut-off is longer and use that one
c           --- not needed when ldual is true since will use rcutin then
            if ( rcut .gt. rcutchg(ibox) ) then
               rcutmax = rcut
            else
               rcutmax = rcutchg(ibox)
            endif
         endif
      endif

c     --- compute minimum cutoff squared
      rminsq = rmin * rmin

      lqimol = lelect(imolty)
      igrow = nugrow(imolty)

      if ( lcutcm .and. (.not. lfirst) ) then
c     --- check previous bead (iufrom) COM for each molecule in the box
         if ( lnew ) then
c           ### for trial chain ###
            rxui  = rxnew(iufrom)
            ryui  = rynew(iufrom)
            rzui  = rznew(iufrom)
         else
c           ### for old chain ###
            rxui  = rxu(i,iufrom)
            ryui  = ryu(i,iufrom)
            rzui  = rzu(i,iufrom)
         endif
         
         do j = 1,nchain
            lcmno(j) = .false.
            if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
               rxuij = rxui-xcm(j)
               ryuij = ryui-ycm(j)
               rzuij = rzui-zcm(j)
c                --- minimum image the pseudo-ctrmas pair separation
c               if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
               
ccccccccccccccccccccccccccccccccccccccccccc


                           if (lpbc) then 
                              
                              if (lsolid(ibox) .and. .not. lrect(ibox))
     &                             then
c     ******************
c     *** Collin's Code:
c     *** non-rectangular box
                                 sx = rxuij*hmati(ibox,1)+ryuij*hmati
     &                                (ibox,4)  +rzuij*hmati(ibox,7)
                                 sy = rxuij*hmati(ibox,2)+ryuij*hmati
     &                                (ibox,5)  +rzuij*hmati(ibox,8)
                                 sz = rxuij*hmati(ibox,3)+ryuij*hmati
     &                                (ibox,6)  +rzuij*hmati(ibox,9)
                                 
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
                                 rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)
     &                                + sz*hmat(ibox,7)
                                 ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)
     &                                + sz*hmat(ibox,8)
                                 rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)
     &                                +sz*hmat(ibox,9)
                                 
                              else
                                 
                                 if ( lpbcx ) then
                                    if ( lfold ) then
                                       if ( rxuij .gt. hbx ) then
                                          rxuij=rxuij-bx
                                       else
                                          if (rxuij.lt.-hbx) rxuij=rxuij
     &                                         +  bx
                                       endif
                                    else
                                       rxuij = rxuij - bx*dint(rxuij*bxi
     &                                      + dsign(0.5d0,rxuij))
                                    endif
                                 endif
                                 
                                 if ( lpbcy ) then
                                    if ( lfold ) then
                                       if ( ryuij .gt. hby ) then
                                          ryuij=ryuij-by
                                       else
                                          if (ryuij.lt.-hby) ryuij=ryuij
     &                                         + by
                                       endif
                                    else
                                       ryuij = ryuij - by*dint(ryuij*byi
     &                                      +dsign(0.5d0,ryuij))
                                    endif
                                 endif
                                 
                                 if ( lpbcz ) then
                                    if ( lfold ) then
                                       if (rzuij.gt.hbz) then
                                          rzuij=rzuij-bz
                                       else
                                          if (rzuij.lt.-hbz) rzuij=rzuij
     &                                         + bz
                                       endif
                                    else
                                       rzuij = rzuij - bz*dint(rzuij*bzi
     &                                      +dsign(0.5d0,rzuij))
                                    endif
                                 endif
                              endif
                              
                              
                           endif

ccccccccccccccccccccccccccccccccccccccccccccccc

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
               

               if ( ldual ) then
                  rcm = rcutin + rcmu(j) + maxlen
                  rcmsq = rcm*rcm
               else
                  rcm = rcutmax + rcmu(j) + maxlen
                  rcmsq = rcm*rcm
               endif
               if (rijsq .gt. rcmsq ) lcmno(j) = .true.
            endif
         enddo
      endif

      do 20 itrial = 1, ichoi
         lovr(itrial) = .false.
         vinter = 0.0d0
         vintra = 0.0d0
         vext = 0.0d0
         velect = 0.0d0
         vewald = 0.0d0 

         do count = 1,ntogrow
            ii = glist(count)

            if (lewald) then
c              -- This part does not change for fixed charge moves, but is
c              -- used in the swap rosenbluth weight. - ewald self term
c              -- 1.772 is sqrt of pi
               vewald = vewald - qqu(icharge,ii)*qqu(icharge,ii)
     &              *calp(ibox)/1.772453851d0
            endif

         enddo

c        --- no intramolecular interactions if this is the first bead
         if ( .not. lfirst ) then

C *****************************************
C *** INTRACHAIN BEAD-BEAD INTERACTIONS ***
C *****************************************

c        --- cycle through molecule and check bead by bead
         do 17 iu = 1, igrow

c           --- see if iu exists in the new chain yet
            if (.not. lexist(iu)) goto 17
            
c           --- loop over all the grown beads
            do count = 1,ntogrow
               ii = glist(count)

c           --- see if iu has nonbonded intramolecular interaction with ii
               if ( linclu(imolty,ii,iu) .or. lewald ) then

c                 --- assign bead type for ii,iu, and the cross term
                  ntii = ntype(imolty,ii)
                  ntjj = ntype(imolty,iu)

                  if (lexpsix .or. lmmff ) then
                     ntij = (ntii+ntjj)/2
                  elseif (lninesix) then
                     ntij = (ntii-1)*nxatom + ntjj
                  else
                     ntij = (ntii-1)*nntype + ntjj
                  endif
            
c                 --- determine distances
                  if ( lnew ) then
c                    --- use new trial chain coordinates
                     rxuij  = rxnew(iu) - rxp(count,itrial)
                     ryuij  = rynew(iu) - ryp(count,itrial)
                     rzuij  = rznew(iu) - rzp(count,itrial)
                  else
c                    --- use old chain coordinates
                     rxuij  = rxu(i,iu) - rxp(count,itrial)
                     ryuij  = ryu(i,iu) - ryp(count,itrial)
                     rzuij  = rzu(i,iu) - rzp(count,itrial)
                  endif

                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  
               endif

               if ( linclu(imolty,ii,iu) .or.
     &              lqinclu(imolty,ii,iu)) then

                  if ( linclu(imolty,ii,iu) ) then

                     if ( rijsq .lt. rminsq .and. 
     &                    .not. lexpand(imolty) ) then
                        lovr(itrial) = .true.
c                     write(6,*) 'intra overlap'
                        goto 19
                     elseif ( rijsq .lt. rcvdwsq .or. lijall) then
                        if (llj.and.(.not.(lexpand(imolty)
     &                       ))) then
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2=epsij(ntij)
                           sr6 = sr2 * sr2 * sr2
                           vintra = vintra 
     &                          + sr6*(sr6-1.0d0)*epsilon2

c * OH 1-5 interaction
                           if (lainclu(imolty,ii,iu)) then
                              vintra = vintra + 0.25d0 * 
     &                             a15(a15type(imolty,ii,iu)) /
     &                             ((rijsq**2)*(rijsq**2)*(rijsq**2))
                           endif                   
                        elseif ( lsami ) then
                           vintra = vintra + ljsami(rijsq,ntij)
                        elseif (lexpsix) then
                           vintra = vintra + exsix(rijsq,ntij)
                        elseif (lmmff) then
                           vintra = vintra + mmff(rijsq,ntij)
                        elseif (lninesix) then
                           vintra = vintra + ninesix(rijsq,ntij)
                        elseif ( lmuir ) then
                           vintra = vintra + ljmuir(rijsq,ntij)
                        elseif ( lpsurf ) then
                           vintra = vintra + ljpsur(rijsq,ntij)
                        else if (lshift) then
                           sr2 = sig2ij(ntij) / rijsq
                           sr6 = sr2 * sr2 * sr2
                           vintra = vintra + sr6 
     &                          *(sr6-1.0d0)*epsij(ntij) - ecut(ntij) 
                        else
                           if ( lexpand(imolty) ) then
                              sigma2 = (sigma(imolty,ii)
     &                             +sigma(imolty,iu))/2.0d0
                              sr2 = sigma2*sigma2 / rijsq
                              epsilon2 = dsqrt(epsilon(imolty,ii)
     &                             *epsilon(imolty,iu))
                           else
                              sr2 = sig2ij(ntij) / rijsq
                              epsilon2 = epsij(ntij)
                           endif
                           sr6 = sr2 * sr2 * sr2
                           vintra = vintra + sr6*(sr6-1.0d0)*epsilon2
c * OH 1-5 interaction
                           if (lainclu(imolty,ii,iu)) then
                              vintra = vintra + 0.25d0 * 
     &                             a15(a15type(imolty,ii,iu)) /
     &                             ((rijsq**2)*(rijsq**2)*(rijsq**2))
                           endif

                        endif
                     endif
                  endif

c intramolecular charge interaction 
c                    --- compute velect (coulomb and ewald)
                  if ( lelect(imolty) .and.
     &                 lqinclu(imolty,ii,iu)) then

c *** boltz.f has problem to compute the electrostatic interactions
c *** in a group-based way because the leader q might not be grown at
c *** present, so it calculates electrostatic interaction not based on 
c *** group but on its own distance in SC, but should be corrected 
c *** later by calling energy subroutine.

                     lcompute = .false.

                     if ( rijsq .lt. rcinchsq )
     &                    lcompute = .true.

c     --- if lcompute then compute the electrostatic energy 
                     if ( lcompute ) then

                        if (lewald) then
c                          --- compute real space term of vewald

                           if (linclu(imolty,ii,iu)) then
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(icharge,iu)*
     &                             erfunc(calp(ibox)*dsqrt(rijsq))/
     &                             dsqrt(rijsq)
c                          * scale 1,4 by qscale
                           else
                              velect = velect + qscale*qqu(icharge,ii)
     &                             *qqu(icharge,iu)*
     &                             erfunc(calp(ibox)*dsqrt(rijsq))/
     &                             dsqrt(rijsq)

c                 --- ewald sum correction term
                              corr = (1.0d0 - qscale)*qqu(icharge,ii)
     &                             *qqu(icharge,iu)*(erfunc(calp(ibox)
     &                             * dsqrt(rijsq))-1.0d0) /dsqrt(rijsq)
                              vewald = vewald + corr

                           endif

c                                 write(6,*) '      ',iu,ii,qqfact*
c     &                                qqu(icharge,ii)*qqu(icharge,iu)*
c     &                                erfunc(calp(ibox)*dsqrt(rijsq))/
c     &                                dsqrt(rijsq)

                        else
c                             --- compute all electrostatic interactions
                           if (linclu(imolty,ii,iu)) then
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(i,iu)/dsqrt(rijsq)
c                          * scale 1,4 by qscale
                           else
                              velect = velect + qscale*qqu(icharge,ii)
     &                             *qqu(i,iu)/dsqrt(rijsq)
                           endif
                        endif
                     endif
                  endif
c end charge calculation 

c will only add correction if lqinclu is false.
 16            elseif ( lewald ) then

c                 --- ewald sum correction term
                  corr = qqu(icharge,ii)*qqu(icharge,iu)
     &                 *(erfunc(calp(ibox)
     &                 * dsqrt(rijsq))-1.0d0) /dsqrt(rijsq)
                  vewald = vewald + corr
               endif

            enddo

 17      continue

         if ( lewald .and. ntogrow .gt. 1) then
c          --- ewald sum correction term for interactions of the 
c          --- growing beads with each other
c * this is 1-3, so don't need to consult lqinclu
c * should change this since it corrects for all currently grown beads in this
c * step, which could at somepoint be further than 1-3 apart!!! (say for rigrot...)
            do cnt = 1,ntogrow-1
               iu = glist(cnt)
               do count = cnt+1,ntogrow
                  ii = glist(count)
                  ntii = ntype(imolty,ii)
                  
c     --- assign bead type for ii,iu, and the cross term
                  ntjj = ntype(imolty,iu)
                  ntij = (ntii-1)*nntype + ntjj
c     --- determine distances - use trial chain coordinates
                  rxuij  = rxp(cnt,itrial) - rxp(count,itrial)
                  ryuij  = ryp(cnt,itrial) - ryp(count,itrial)
                  rzuij  = rzp(cnt,itrial) - rzp(count,itrial)
                  
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  
c     --- ewald sum correction term
                  corr = qqu(icharge,ii)*qqu(icharge,iu)
     &                 *(erfunc(calp(ibox)
     &                 * dsqrt(rijsq))-1.0d0) /dsqrt(rijsq)
                  vewald = vewald + corr
               enddo
            enddo
         endif

         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &         .and. .not. lninesix ) then 
                 vintra = 4.0d0 * vintra
         endif
 
c 
        end if

c     grand-canonical: if ibox = 2 (ideal gas box) only intra-chain
         if ( .not. (lgrand .and. ibox .eq. 2) ) then

            if (licell.and.(ibox.eq.boxlink)) then
c     --- we only use count = 1, the rest should be taken care of
c     --- in rintramax
               count = 1
               rxui = rxp(count,itrial)
               ryui = ryp(count,itrial)
               rzui = rzp(count,itrial)

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
                  
               nmole = 0
               do j = 1, 27
                  ic = cellinc(j)
                     
                  do k = 1, nicell(ic)

                     nmole = nmole + 1
                     jcell(nmole) = iucell(ic,k)
                             
                  enddo
               enddo
            else
               nmole = nchain                  
            endif

C *******************************
C *** INTERCHAIN INTERACTIONS ***
C *******************************
c        ---    loop over all chains except i
         do 98 k = 1, nmole

            if (licell.and.(ibox.eq.boxlink)) then
               j = jcell(k)
            else
               j = k
            endif

c     --- check for simulation box 
            if ( ( nboxi(j) .eq. ibox ) .and. ( i .ne. j ) ) then
               if ( lneigh ) then
c                 --- check neighbor list
                  if ( .not. lnew ) then
                     if ( .not. lnn(j,i) ) goto 98
                  endif
               endif

               jmolty = moltyp(j)
               lqjmol = lelect(jmolty)
               growjj = nugrow(jmolty)

               if ( .not. lfirst ) then
c                 --- check COM table calculated above
                  if (lcmno(j) .and. lcutcm) 
     &                 goto 98
               endif
               
c              --- loop over all beads of molecule i grown this step
 108           do count = 1,ntogrow
c                 --- assign bead type for ii
                  ii = glist(count)
                  ntii = ntype(imolty,ii)
                  liji = lij(ntii)
                  lqchgi = lqchg(ntii)
                  
c                 --- assign positions to r*ui
                  rxui = rxp(count,itrial)
                  ryui = ryp(count,itrial)
                  rzui = rzp(count,itrial)
            
                  if ( lfirst .and. lcutcm ) then
c                    --- check if ctrmas within rcmsq
    
                     rxuij = rxui-xcm(j)
                     ryuij = ryui-ycm(j)
                     rzuij = rzui-zcm(j)

c                    --- minimum image the ctrmas pair separations
c                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
cccccccccccccccccccccccccccccccccccccccccccccc

                           if (lpbc) then 
                              
                              if (lsolid(ibox) .and. .not. lrect(ibox))
     &                             then
c     ******************
c     *** Collin's Code:
c     *** non-rectangular box
                                 sx = rxuij*hmati(ibox,1)+ryuij*hmati
     &                                (ibox,4)  +rzuij*hmati(ibox,7)
                                 sy = rxuij*hmati(ibox,2)+ryuij*hmati
     &                                (ibox,5)  +rzuij*hmati(ibox,8)
                                 sz = rxuij*hmati(ibox,3)+ryuij*hmati
     &                                (ibox,6)  +rzuij*hmati(ibox,9)
                                 
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
                                 rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)
     &                                + sz*hmat(ibox,7)
                                 ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)
     &                                + sz*hmat(ibox,8)
                                 rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)
     &                                +sz*hmat(ibox,9)
                                 
                              else
                                 
                                 if ( lpbcx ) then
                                    if ( lfold ) then
                                       if ( rxuij .gt. hbx ) then
                                          rxuij=rxuij-bx
                                       else
                                          if (rxuij.lt.-hbx) rxuij=rxuij
     &                                         +  bx
                                       endif
                                    else
                                       rxuij = rxuij - bx*dint(rxuij*bxi
     &                                      + dsign(0.5d0,rxuij))
                                    endif
                                 endif
                                 
                                 if ( lpbcy ) then
                                    if ( lfold ) then
                                       if ( ryuij .gt. hby ) then
                                          ryuij=ryuij-by
                                       else
                                          if (ryuij.lt.-hby) ryuij=ryuij
     &                                         + by
                                       endif
                                    else
                                       ryuij = ryuij - by*dint(ryuij*byi
     &                                      +dsign(0.5d0,ryuij))
                                    endif
                                 endif
                                 
                                 if ( lpbcz ) then
                                    if ( lfold ) then
                                       if (rzuij.gt.hbz) then
                                          rzuij=rzuij-bz
                                       else
                                          if (rzuij.lt.-hbz) rzuij=rzuij
     &                                         + bz
                                       endif
                                    else
                                       rzuij = rzuij - bz*dint(rzuij*bzi
     &                                      +dsign(0.5d0,rzuij))
                                    endif
                                 endif
                              endif
                              
                              
                           endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc


                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     
c                    --- determine cutoff
                     if ( ldual ) then
c                       --- must be lfirst so no previous bead
                        rcm = rcutin + rcmu(j)
                     else
c                       --- standard lcutcm cutoff
                        rcm = rcutmax + rcmu(j)
                     endif

c                    --- check if interaction distance is greater than cutoff
                     rcmsq = rcm*rcm
                     if ( rijsq .gt. rcmsq )
     &                    goto 98
                  endif

c                 --- loop over all beads jj of chain j
                  do 97 jj = 1,nunit(jmolty)
               
c                    --- check exclusion table
                     if ( lexclu(imolty,ii,jmolty,jj) ) goto 97
                     
c *** start iswatch add-on ***
c *** is there a way to pull this out of the loops? ***
                     if ( liswatch ) then
                        if (j .eq. other) then
                              if ( .not. liswinc(jj,jmolty) ) then
c                                 write(6,*) 'iSwatch-skipping:',jj
                                 goto 97
                              endif
                        endif
                     endif
c *** end iswatch add-on ***

                     ntjj = ntype(jmolty,jj)
                     if ( (.not. (liji .and. lij(ntjj))) 
     &                    .and. 
     &                    (.not. (lqchgi .and. lqchg(ntjj)))) 
     &                    goto 97
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

c                    --- minimum image the pair separations ***
c                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
                     
cccccccccccccccccccccccccccccccccccccccccc
                     
                     if (lpbc) then 
                              
                        if (lsolid(ibox) .and. .not. lrect(ibox))
     &                       then
c     ******************
c     *** Collin's Code:
c     *** non-rectangular box
                           sx = rxuij*hmati(ibox,1)+ryuij*hmati
     &                          (ibox,4)  +rzuij*hmati(ibox,7)
                           sy = rxuij*hmati(ibox,2)+ryuij*hmati
     &                          (ibox,5)  +rzuij*hmati(ibox,8)
                           sz = rxuij*hmati(ibox,3)+ryuij*hmati
     &                          (ibox,6)  +rzuij*hmati(ibox,9)
                           
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
                           rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)
     &                          + sz*hmat(ibox,7)
                           ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)
     &                          + sz*hmat(ibox,8)
                           rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)
     &                          +sz*hmat(ibox,9)
                           
                        else
                           
                           if ( lpbcx ) then
                              if ( lfold ) then
                                 if ( rxuij .gt. hbx ) then
                                    rxuij=rxuij-bx
                                 else
                                    if (rxuij.lt.-hbx) rxuij=rxuij
     &                                   +  bx
                                 endif
                              else
                                 rxuij = rxuij - bx*dint(rxuij*bxi
     &                                      + dsign(0.5d0,rxuij))
                              endif
                           endif
                           
                           if ( lpbcy ) then
                              if ( lfold ) then
                                 if ( ryuij .gt. hby ) then
                                    ryuij=ryuij-by
                                 else
                                    if (ryuij.lt.-hby) ryuij=ryuij
     &                                   + by
                                 endif
                              else
                                 ryuij = ryuij - by*dint(ryuij*byi
     &                                +dsign(0.5d0,ryuij))
                              endif
                           endif
                           
                           if ( lpbcz ) then
                              if ( lfold ) then
                                 if (rzuij.gt.hbz) then
                                    rzuij=rzuij-bz
                                 else
                                    if (rzuij.lt.-hbz) rzuij=rzuij
     &                                   + bz
                                 endif
                              else
                                 rzuij = rzuij - bz*dint(rzuij*bzi
     &                                +dsign(0.5d0,rzuij))
                              endif
                           endif
                        endif
                     endif
                     
cccccccccccccccccccccccccccccccccccccccccc
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

c                    --- compute vinter (ex. lennard-jones)
                     if ( rijsq .lt. rminsq .and. .not. 
     &                    (lexpand(imolty) .or. lexpand(jmolty))) then
                        lovr(itrial) = .true.
c                        write(6,*) 'j:',j,jj
c                        write(6,*) 'rjsq:',rijsq,rminsq
                        goto 19
                     elseif (rijsq .lt. rcinvdsq .or. lijall) then
                        if (llj.and.(.not.(lexpand(imolty).or.
     &                       lexpand(jmolty)))) then
                           if ( lij(ntii) .and. lij(ntjj) ) then
                              sr2 = sig2ij(ntij) / rijsq
                              epsilon2=epsij(ntij)
                              sr6 = sr2 * sr2 * sr2
                              vinter = vinter 
     &                             + sr6*(sr6-1.0d0)*epsilon2
                           endif
                        elseif ( lsami) then
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
                           vinter = vinter + sr6*(sr6-1.0d0)*epsij(ntij) 
     +                          - ecut(ntij) 
                        else if ( lij(ntii) .and. lij(ntjj) ) then
                           if (lfepsi) then
                              sr6 = rijsq*rijsq*rijsq
                              if ( (.not. lqchg(ntii)) .and. 
     &                             (.not. lqchg(ntjj)) ) then
                                 if ( nunit(imolty) .eq. 4 ) then
c *** TIP-4P structure (temperary use ???)
                                    qave = (qqu(i,4)+qqu(j,4))/2.0d0
                                 else
                                    qave = (qqu(i,4)+qqu(i,5)+
     &                                   qqu(j,4)+qqu(j,5))*0.85d0
                                 endif
                              else
                                 qave = (qqu(i,ii)+qqu(j,jj))/2.0d0
                              endif
                              if ( lexpand(imolty) 
     &                             .and. lexpand(jmolty)) then
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif (lexpand(imolty)) then
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 epsilon2 = dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 epsilon2 = epsij(ntij)
                              endif
                              vinter = vinter + ((aslope*(qave-a0)*
     &                             (qave-a0)+ashift)/sr6 - (bslope*
     &                             (qave-b0)*(qave-b0)+bshift))/
     &                             sr6*epsilon2
                           else
                              if ( lexpand(imolty) .and. 
     &                             lexpand(jmolty)) then
                                 sigma2 = (sigma(imolty,ii)+
     &                                sigma(jmolty,jj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsilon(jmolty,jj))
                              elseif ( lexpand(imolty) ) then
                                 sigma2 = (sigma(imolty,ii)+
     &                                sigi(ntjj))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2 = dsqrt(epsilon(imolty,ii)*
     &                                epsi(ntjj))
                              elseif ( lexpand(jmolty) ) then
                                 sigma2 = (sigma(jmolty,jj)+
     &                                sigi(ntii))/2.0d0
                                 sr2 = sigma2*sigma2/rijsq
                                 epsilon2 = dsqrt(epsi(ntii)*
     &                                epsilon(jmolty,jj))
                              else
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2 = epsij(ntij)
                              endif
                              sr6 = sr2 * sr2 * sr2
                              vinter = vinter + 
     &                             sr6*(sr6-1.0d0)*epsilon2

                           endif
                        endif
                     endif

c                    --- compute velect (coulomb and ewald)
c * lcompute has not been set yet, this was wrong
c                     if ( lcompute.and.lqimol .and. lqjmol .and. 
c     &                    lqchg(ntii) .and. lqchg(ntjj) ) then
                     if ( lqimol .and. lqjmol .and. 
     &                    lqchg(ntii) .and. lqchg(ntjj) ) then

c *** boltz.f has problem to compute the electrostatic interactions
c *** in a group-based way because the leader q might not be grown at
c *** present, so it calculates electrostatic interaction not based on 
c *** group but on its own distance in SC, but should be corrected 
c *** later by calling energy subroutine.

                        lcompute = .false.

                        if ( rijsq .lt. rcinchsq ) 
     &                       lcompute = .true.

c     --- if lcompute then compute the electrostatic energy 
                        if ( lcompute ) then

                           if (lewald) then
c                             --- compute real space term of velect
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(j,jj)*
     &                             erfunc(calp(ibox)*dsqrt(rijsq))/
     &                             dsqrt(rijsq)
                           else
c                             --- compute all electrostatic interactions
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(j,jj)/
     &                             dsqrt(rijsq)
                           endif
                        endif
                     endif
                     
 97               continue 
               enddo
            endif
 98      continue
         
         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &         .and. .not. lninesix ) then
                 vinter = 4.0d0 * vinter
         endif
      end if         
c ################################################################

C **************************************************************
C *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
C ***************************************************************
 
c ---  not for grand can. with ibox=2 !
         if ( ljoe .or. lsami .or. lmuir .or. lexzeo
     +         .or. lgraphite ) then
            do count = 1,ntogrow
c              --- assign bead type for ii
               ii = glist(count)
               ntii = ntype(imolty,ii)
               rxui = rxp(count,itrial)
               ryui = ryp(count,itrial)
               rzui = rzp(count,itrial)
               if ( ljoe ) then
                  if ( extc12(ntii) .gt. 0.1d0 ) then
                     dzui = rzui - extz0(ntii)
                     dz3  = dzui * dzui * dzui
                     dz12 = dz3**4
                     vext = vext +
     +                    (extc12(ntii)/dz12) - (extc3(ntii)/dz3)    
                  endif
               endif 
               if( lgraphite ) then
                  ntij = (ntii-1)*nntype + ntsubst
                  vext = vext + exgrph(rxui,ryui,rzui,ntij)
               endif

               if ( lsami ) vext = vext + exsami(rzui,ntii)
               if ( lmuir ) vext = vext + exmuir(rzui,ntii)
               if ( lexzeo ) vext = vext + exzeo(rxui,ryui,rzui,ntii)
            enddo
         endif
C ----------------------------------------------------------------------------
 
C *********************************************
C *** CALCULATION OF TOTAL POTENTIAL ENERGY ***
C *********************************************
         
 19      if ( lovr(itrial) ) then
            bfac(itrial) = 0.0d0
         else
            velect = velect*qqfact
            vewald = vewald*qqfact
            vtry(itrial) = vinter + vintra + vext + velect + vewald
            vtrintra(itrial) = vintra
            vtrext(itrial)   = vext
            vtrinter(itrial) = vinter
            vtrelect(itrial) = velect
            vtrewald(itrial) = vewald
            if ( ( vtry(itrial) * beta ) .gt. (2.3d0*softcut) ) then
c     write(6,*) 'caught by softcut',vtry(itrial)*beta
               lovr(itrial) = .true.
               bfac(itrial) = 0.0d0
            elseif ( ( vtry(itrial) * beta ) .lt. -2.303d0*308 ) then
               write(6,*) '### warning: weight too big out of range'
               lovr(itrial) = .true.
               bfac(itrial) = 0.0d0
            else
               bfac(itrial) = dexp ( -(vtry(itrial)*beta) )
            endif
         endif

 20   continue



c ----------------------------------------------------------------------------
      ovrlap = .true.
      do itrial = 1, ichoi
         if ( .not. lovr(itrial)) ovrlap = .false.
      enddo

c      write(6,*) 'end BOLTZ'

      return

      end


