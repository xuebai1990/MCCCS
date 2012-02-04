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
      include 'ipswpar.inc'
      include 'eepar.inc'
      include 'conver.inc'
      include 'tabulated.inc'
c RP added for MPI
      include 'mpif.h'
      include 'mpi.inc'

      logical lnew,ovrlap,lcmno,lfirst,lcompute,lcoulo
      logical lqimol,lqjmol,liji,lqchgi
      integer ichoi,growjj,igrow,count,glist,icharge,cnt,jcell,ic
      integer i,imolty,ibox,ntogrow,itrial,ntii,j,jj,ntjj,ntij
     +       ,iu,jmolty,jjj,iufrom,ii,zz,bdmol_b,cellinc,k,nmole


!      integer NRtype 

      double precision ljsami,rminsq,rxui,sr6,ryui,rzui
     +     ,rxuij,ryuij,rzuij,rij,rijsq,sr2,dzui,dz3,dz12
     +     ,exzeo,exsami,exmuir,exgrph,ljpsur,ljmuir,exsix
     +     ,mmff,maxlen,rcm,rcmsq
     +     ,corr,erfunc,rcutmax,ninesix, genlj
      double precision vinter,vintra,vext,velect,vewald,qave,
     &     epsilon2,sigma2,vwell,v,rcutsq,rcinsq

      double precision sx,sy,sz,v_elect_field, field
      double precision slitpore
      	
      dimension lcmno(nmax),lcoulo(numax,numax)
      dimension glist(numax),cellinc(27),jcell(nmax)

      double precision tabulated_bend, tabulated_vdW, tabulated_elect
      integer mmm

c------------- RP added for MPI
      integer my_start,my_end,loops_per_proc,scount,my_itrial
      double precision my_vtry(nchmax),my_vtrintra(nchmax),
     &   my_vtrext(nchmax),my_vtrinter(nchmax),my_vtrelect(nchmax)
     &   ,my_vtrewald(nchmax),my_bfac(nchmax),my_vipswot(nchmax)
     & ,my_vwellipswot(nchmax),my_vipswnt(nchmax),my_vwellipswnt(nchmax)
      logical my_lovr(nchmax)
      integer ncount_arr(numprocmax+1),ncount_displs(numprocmax+1)
c ------------------------------------------
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
c      write(iou,*) 'start BOLTZ'

c      lcompute = .false.

      if ( lpbc ) call setpbc (ibox)

c     --- determine the potential cutoffs

      rcutsq = rcut(ibox)*rcut(ibox)
      field = Elect_field(ibox)

c KM
c initialize variables
      do j=1,ichoi
         my_lovr(j) = .false.
         lovr(j) = .false.
         my_vtry(j) = 0.0d0
         my_vtrintra(j) = 0.0d0
         my_vtrelect(j) = 0.0d0
         my_vtrext(j) = 0.0d0
         my_vtrinter(j) = 0.0d0
         my_vtrewald(j) = 0.0d0
         my_bfac(j) = 0.0d0
         my_vipswot(j) = 0.0d0
         my_vwellipswot(j) = 0.0d0
         my_vipswnt(j) = 0.0d0
         my_vwellipswnt(j) = 0.0d0
      enddo
      do j=1,numprocs
         ncount_arr(j) = 0
         ncount_displs(j) = 0
      enddo
      

      if ( ldual ) then
c        --- use rcutin for both types of interactions (except intra)
         rcinsq = rcutin*rcutin
      else
c        --- compute the cutoffs squared for each interaction
         rcinsq =  rcutsq
         if ( lcutcm ) then
c           --- not needed when ldual is true since will use rcutin then
               rcutmax = rcut(ibox)
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
               if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
               rij  = dsqrt(rijsq)                

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

c RP added for MPI

      loops_per_proc = ichoi/numprocs

      my_start = (myid*loops_per_proc)+1
      if(myid .eq. (numprocs-1))then
        my_end = ichoi
      else
         my_end   = (myid + 1)*loops_per_proc
      endif

c RP added for MPI
c      write(6,*)'170: boltz my_start=',my_start,'my_end=',my_end
c     &    ,'ichoi=',ichoi,'loops_per_proc=',loops_per_proc,'myid=',myid
c
      my_itrial  = 0
c      do 20 itrial = 1, ichoi
      do 20 itrial = my_start,my_end
         my_itrial  = my_itrial + 1
         my_lovr(my_itrial) = .false.
       
c         lovr(itrial) = .false.
         vinter = 0.0d0
         vintra = 0.0d0
         vext = 0.0d0
         velect = 0.0d0
         vewald = 0.0d0 

c -- if L_Coul_CBMC is true  only then compute electrostatic interactions/corrections
         if(L_Coul_CBMC) then
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
         endif

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
                  elseif (lgenlj) then
                     ntij = (ntii-1)*nntype + ntjj
                  else
                     ntij = (ntii-1)*nntype + ntjj
                  endif
                  if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
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
                  if (lpbc) call mimage ( rxuij,ryuij,rzuij,ibox )
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  rij = dsqrt(rijsq)
               endif
               if ( linclu(imolty,ii,iu) .or.
     &              lqinclu(imolty,ii,iu)) then
                  if ( linclu(imolty,ii,iu) ) then
                     if ( rijsq .lt. rminsq .and. 
     &                    .not. lexpand(imolty) ) then
c RP added for MPI
                        my_lovr(my_itrial) = .true.

c                     write(iou,*) 'intra overlap'
                        goto 19
                     elseif ( rijsq .lt. rcutsq .or. lijall) then
                        if (L_vdW_table.or. L_bend_table.and.
     &                       (.not.lexpand(imolty))) then

                           do mmm=1,inben(imolty,ii)
                              if (ijben3(imolty,ii,mmm).eq.jj) then
                                 
                                 call lininter_bend(rij,
     &                                tabulated_bend, 
     &                                itben(imolty,ii,mmm))
                                 vintra = vintra + tabulated_bend
                                 
                                 goto 96
                              endif
                           enddo

                           call lininter_vdW(rij, 
     &                          tabulated_vdW, ntii, ntjj)
                           vintra = vintra + tabulated_vdW

                        elseif (llj.and.(.not.(lexpand(imolty)
     &                       ))) then
                           sr2 = sig2ij(ntij) / rijsq
                           epsilon2=epsij(ntij)
                           sr6 = sr2 * sr2 * sr2
                           vintra = vintra 
     &                          + sr6*(sr6-1.0d0)*epsilon2
     &				* ljscale(imolty,ii,iu)
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
                           vintra = vintra + sr6 
     &                          *((sr6-1.0d0)*epsij(ntij)-ecut(ntij))
     &				*ljscale(imolty,ii,iu)
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
     &					*ljscale(imolty,ii,iu)		   
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
 96               if (L_Coul_CBMC) then
                  if ( lelect(imolty) .and.
     &                 lqinclu(imolty,ii,iu)) then

c *** boltz.f has problem to compute the electrostatic interactions
c *** in a group-based way because the leader q might not be grown at
c *** present, so it calculates electrostatic interaction not based on 
c *** group but on its own distance in SC, but should be corrected 
c *** later by calling energy subroutine.

                     lcompute = .false.
                     if ( rijsq .lt. rcinsq )
     &                    lcompute = .true.
c     --- if lcompute then compute the electrostatic energy 
                        if ( lcompute ) then
                           if (L_elect_table) then
                              call lininter_elect(rij, 
     &                             tabulated_elect, ntii, ntjj)
                              velect = velect + qscale2(imolty,ii,iu)*
     &                             qqu(icharge,ii)*
     &                             qqu(icharge,iu)*tabulated_elect
                           elseif (lewald) then
c                          --- compute real space term of vewald
                              velect = velect + 
     &                          qscale2(imolty,ii,iu)*qqu(icharge,ii)
     &                             *qqu(icharge,iu)*
     &                             erfunc(calp(ibox)*rij)/
     &                             rij
c                 --- ewald sum correction term
                              corr = (1.0d0 - qscale2(imolty,ii,iu))*
     &                           qqu(icharge,ii)
     &                             *qqu(icharge,iu)*(erfunc(calp(ibox)
     &                             * rij)-1.0d0) /rij
                              vewald = vewald + corr
                           else
                              velect = velect + qscale2(imolty,ii,iu)
     &                             *qqu(icharge,ii)
     &                             *qqu(i,iu)/rij
                            endif
                        endif
                     endif
                     endif
c end charge calculation 

c will only add correction if lqinclu is false.
                  elseif ( lewald .and.L_Coul_CBMC) then
c                 --- ewald sum correction term
                     corr = qqu(icharge,ii)*qqu(icharge,iu)
     &                 *(erfunc(calp(ibox)
     &                 * rij)-1.0d0) /rij
                     vewald = vewald + corr
                  endif
            enddo

 17      continue

         if (L_Coul_CBMC) then 
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
                  rij   = dsqrt(rijsq)                  
c     --- ewald sum correction term
                  corr = qqu(icharge,ii)*qqu(icharge,iu)
     &                 *(erfunc(calp(ibox)
     &                 * rij)-1.0d0) /rij
                  vewald = vewald + corr
               enddo
            enddo
         endif
         endif
 
         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     & .and. .not. lgenlj .and. .not. lninesix
     & .and..not.L_vdW_table.and..not.L_bend_table) then 
                 vintra = 4.0d0 * vintra
         endif
 
c 
        end if

c     grand-canonical: if ibox = 2 (ideal gas box) only intra-chain
c --- JLR 11-24-09 don't compute if lideal          
c         if ( .not. ((lgrand .and. ibox .eq. 2))) then 
         if ( .not. (lgrand .and. ibox .eq. 2) 
     &       .and. .not.lideal(ibox) ) then
c --- END JLR 11-24-09
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
                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     rij   = dsqrt(rijsq)
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
c                                 write(iou,*) 'iSwatch-skipping:',jj
                                 goto 97
                              endif
                        endif
                     endif
c *** end iswatch add-on ***

                     ntjj = ntype(jmolty,jj)
                     if ( (.not. (liji .and. lij(ntjj))) 
     &                    .and. 
     &                    (.not. (lqchgi .and. lqchg(ntjj)))
     &                    .and..not.L_vdW_table) 
     &                    goto 97
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

                     rxuij = rxui - rxu(j,jj)
                     ryuij = ryui - ryu(j,jj)
                     rzuij = rzui - rzu(j,jj)

c                    --- minimum image the pair separations ***
                     if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     rij   = dsqrt(rijsq)
c                    --- compute vinter (ex. lennard-jones)
                     if ( rijsq .lt. rminsq .and. .not. 
     &                    (lexpand(imolty) .or. lexpand(jmolty))) then
                        my_lovr(my_itrial) = .true.
c                        write(iou,*) 'j:',j,jj
c                        write(iou,*) 'rjsq:',rijsq,rminsq
                        goto 19
                     elseif (rijsq .lt. rcinsq .or. lijall) then
                        if (L_vdW_table.and.(.not.
     &                       (lexpand(imolty).or.lexpand(jmolty))))
     &                       then
                           call lininter_vdW(rij, 
     &                          tabulated_vdW, ntii, ntjj)
                           vinter = vinter + tabulated_vdW
                        elseif (llj.and.(.not.(lexpand(imolty).or.
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
                     if (L_Coul_CBMC) then
                     if ( lqimol .and. lqjmol .and. 
     &                    lqchg(ntii) .and. lqchg(ntjj) ) then
c *** boltz.f has problem to compute the electrostatic interactions
c *** in a group-based way because the leader q might not be grown at
c *** present, so it calculates electrostatic interaction not based on 
c *** group but on its own distance in SC, but should be corrected 
c *** later by calling energy subroutine.
                        lcompute = .false.
                        if ( rijsq .lt. rcinsq ) 
     &                       lcompute = .true.
c     --- if lcompute then compute the electrostatic energy 
                        if ( lcompute ) then
                           if (L_elect_table) then
                              call lininter_elect(rij, 
     &                             tabulated_elect, ntii, ntjj)
                              velect = velect + qqu(icharge,ii)*
     &                             qqu(j,jj)*tabulated_elect
                           elseif (lewald) then
c                             --- compute real space term of velect
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(j,jj)*
     &                             erfunc(calp(ibox)*rij)/
     &                             rij
                           else
c                             --- compute all electrostatic interactions
                              velect = velect + qqu(icharge,ii)
     &                             *qqu(j,jj)/
     &                             rij
                           endif
                        endif
                     endif
                     endif 
 97               continue 
               enddo
            endif
            
 98      continue
         
         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &   .and. .not. lgenlj   .and. .not. lninesix
     &   .and..not.L_vdW_table ) then
                 vinter = 4.0d0 * vinter
         endif
      end if         
c ################################################################

C **************************************************************
C *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
C ***************************************************************
 
c ---  not for grand can. with ibox=2 !
c -- required for histogram reweighting to work for monolayer 
c -- phase diagrams.
c -- not used for adsorption isotherms
         if (.not. (lslit .and. ibox.eq.2)) then
         if ( ljoe .or. lsami .or. lmuir .or. lexzeo
     +         .or. lgraphite .or. lslit ) then
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
	        
	       if (lslit) then
	          ntij = (ntii-1)*nntype + ntsubst		  
c -- calculate interaction with surface at the bottom of the box		  
		  vext = vext+slitpore(rzui,ntij)
c -- calculate interaction with the surface at the top of the box
		  dzui = boxlz(ibox)-rzui
		  vext = vext +slitpore(dzui,ntij)
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
	 endif


         if (.not.ldual) then
            if (lelect_field) then
               if(lelect(moltyp(i))) then
                  if (nboxi(i).eq.ibox) then
                     do count = 1,ntogrow
                        rzui = rzp(count,itrial)
                        vext = vext + v_elect_field(i,count,rzui,field)
                     enddo
                  endif
               endif
               vext = vext * eXV_to_K
            endif 
         endif
         
c --------------------------------------------------------------------------
c well potential for thermodynamic integration stages b and c
c --------------------------------------------------------------------------
                                                                                
      vwell = 0.0d0
      if (lwell(imolty).and.lmipsw) then
                                                                                
         rxui = xcm(i)
         ryui = ycm(i)
         rzui = zcm(i)
         do j = 1, nwell(imolty)*nunit(imolty)
            k = j - int(j/nunit(imolty))*nunit(imolty)
            if (k.eq.0) k = nunit(imolty)
            rxuij = rxui-rxwell(j,imolty)
            ryuij = ryui-rywell(j,imolty)
            rzuij = rzui-rzwell(j,imolty)
            call mimage(rxuij,ryuij,rzuij,ibox)
            rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
            rcm = rcut(ibox)+rcmu(i)+maxlen
            rcmsq = rcm*rcm
            if (rijsq.lt.rcmsq) then
            do count = 1, ntogrow
               ii = glist(count)
               if (awell(ii,k,imolty).lt.1.0d-6) goto 666
               rxui = rxp(count,itrial)
               ryui = ryp(count,itrial)
               rzui = rzp(count,itrial)
               rxuij = rxui-rxwell(j,imolty)
               ryuij = ryui-rywell(j,imolty)
               rzuij = rzui-rzwell(j,imolty)
               call mimage(rxuij,ryuij,rzuij,ibox)
               rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
               vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
 666        enddo
            endif
         enddo
                                                                                
      endif
C ----------------------------------------------------------------------------
 
C *********************************************
C *** CALCULATION OF TOTAL POTENTIAL ENERGY ***
C *********************************************
!         write(23,*) 
!         write(23,*) 'wirting out total energy'
!         write(23,*) 'Lnew', lnew 
!         if (NRtype.eq.1) then
!             write(23,*) 'called from swap'
!         else
!             write(23,*) 'called from rigrot'
!         endif 


 19      if ( my_lovr(my_itrial) ) then
            my_bfac(my_itrial) = 0.0d0
         else
            if (.not.L_elect_table) then
               velect = velect*qqfact
               vewald = vewald*qqfact
            endif
            v = vinter+vintra+vext+velect+vewald
            if (.not.lnew) then
               my_vipswot(my_itrial) = v
               my_vwellipswot(my_itrial) = vwell
            else
               my_vipswnt(my_itrial) = v
               my_vwellipswnt(my_itrial) = vwell
            endif
                                                                                
            if (lstagea) then
               v = (1.0d0-lambdais*(1.0d0-etais))*v
            elseif (lstageb) then
               v = etais*v+lambdais*vwell
            elseif (lstagec) then
               v = (etais+(1.0d0-etais)*lambdais)*v+
     &                         (1.0d0-lambdais)*vwell
            endif

            my_vtry(my_itrial) = v
            my_vtrintra(my_itrial) = vintra
            my_vtrext(my_itrial)   = vext
            my_vtrinter(my_itrial) = vinter
            my_vtrelect(my_itrial) = velect
            my_vtrewald(my_itrial) = vewald
!            write(23,*) 'itrial' ,itrial
!            write(23,*) vtry(itrial), vtrintra(itrial), vtrext(itrial),
!     &       vtrinter(itrial),vtrelect(itrial), vtrewald(itrial)


            if ((my_vtry(my_itrial)*beta).gt.(2.3d0*softcut))then
c               write(iou,*) 'caught by softcut',vtry(itrial)*beta
               my_lovr(my_itrial) = .true.
               my_bfac(my_itrial) = 0.0d0
            elseif((my_vtry(my_itrial)*beta).lt.-2.303d0*308)then
c               write(iou,*) '### warning: weight too big out of range'
               my_lovr(my_itrial) = .true.
               my_bfac(my_itrial) = 0.0d0
            else
               my_bfac(my_itrial) = dexp ( -(my_vtry(my_itrial)*beta) )
            endif
         endif
 20   continue

c      scount = loops_per_proc
       loops_per_proc = (my_end-my_start) + 1

       CALL MPI_ALLGATHER(loops_per_proc,1,MPI_INTEGER,ncount_arr,
     &       1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
       ncount_displs(1) = 0

c KM for MPI
c cannot loop over i, must use j
c i is an argument to the boltz subroutine!
       do j = 2,numprocs
           ncount_displs(j) = ncount_displs(j-1) + ncount_arr(j-1)
       enddo

      call MPI_ALLGATHERV(my_lovr,loops_per_proc,MPI_LOGICAL
     &     ,lovr,ncount_arr,ncount_displs,
     &     MPI_LOGICAL,MPI_COMM_WORLD,ierr)

      ovrlap = .true.
      do itrial = 1, ichoi
         if ( .not. lovr(itrial)) ovrlap = .false.
      enddo

      call MPI_ALLGATHERV(my_vtry, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtry, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrintra, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrintra, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrext,loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrext, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrinter, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrinter, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrelect,loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrelect, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vtrewald, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vtrewald, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_bfac, loops_per_proc, 
     &     MPI_DOUBLE_PRECISION, bfac, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

c       do my_itrial=my_start,my_end
c        write(6,*)"938:boltz my_bfac(",my_itrial,")=",my_bfac(my_itrial)
c     &      ,'myid=',myid
c       enddo

      call MPI_ALLGATHERV(my_vipswot, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vipswot, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vwellipswot, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vwellipswot, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vipswnt, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vipswnt, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

      call MPI_ALLGATHERV(my_vwellipswnt, loops_per_proc,
     &     MPI_DOUBLE_PRECISION, vwellipswnt, ncount_arr, ncount_displs,
     &     MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

c ----------------------------------------------------------------------------


c      write(iou,*) 'end BOLTZ'

      return

      end


