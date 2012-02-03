      subroutine pressure ( press, surf, ibox )

c pressure
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
 
c    *****************************************************
c    ** calculates the pressure for a configuration.    **
c    *****************************************************
 
      implicit none
 
c *** common blocks ***

      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'poten.inc' 
      include 'expsix.inc'
      include 'merck.inc' 
      include 'connect.inc'
      include 'ewaldsum.inc'
      include 'conver.inc'
      include 'inputdata.inc'
      include 'fepsi.inc'
      include 'qqlist.inc'
      include 'nsix.inc'
      include 'ipswpar.inc'
      include 'cell.inc'

      logical lcoulo,lexplt,lqimol,lqjmol,lij2
      integer ibox,itype
      integer i,ii,j,jj,ntii,ntjj,ntij,imolty,jmolty,iii,jjj,k
 
      double precision press,repress,erfunc

      double precision rxui,ryui,rzui,rxuij,ryuij,rzuij,rijsq,
     +                 rcutsq,sr2,sr6,rhosq,corp,rij,rs1,
     +                 sr1,rs2,sr7,rs7,rs6

      double precision surf,pxx,pyy,pzz,rpxx,rpyy,rpzz,pxy,pyx
     &                ,pxz,pzx,pyz,pzy,rpxy,rpyx,rpxz,rpzx,rpyz,rpzy

      double precision fxcmi,fycmi,fzcmi,fij,xcmi,ycmi,zcmi,flj
      double precision rcm,rcmsq,rcmi
      double precision rvdwsq,rchgsq,rbcut,qave,
     &     volsq,epsilon2,sigma2,pwell,vol 

      dimension lcoulo(numax,numax)
C --------------------------------------------------------------------
      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut  = rcut(ibox)

      press = 0.0d0
      pxx = 0.0d0
      pyy = 0.0d0
      pzz = 0.0d0

      do i = 1, 3
         do j = 1, 3
            pips(i,j) = 0.0d0
         enddo
      enddo

C *******************************
C *** INTERCHAIN INTERACTIONS ***
C *******************************

      

c      if(LSOLPAR.and.(ibox.eq.2)) then
c         press= 1.380662d4 * ( ( nchbox(ibox) / beta) -
c     +     ( press/3.0d0 ) ) /
c     +     ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
c         surf = 0.0d0
c         return
c      endif


c --- loop over all chains i 
      do 100 i = 1, nchain - 1
 
c ### check if i is in relevant box ###
         if ( nboxi(i) .eq. ibox ) then

            imolty = moltyp(i)
            lqimol = lelect(imolty)
            if ( nugrow(imolty) .eq. nunit(imolty) ) then
               lexplt = .false.
            else
               lexplt = .true.
            endif
            xcmi = xcm(i)
            ycmi = ycm(i)
            zcmi = zcm(i)
            if (lcutcm) then
               rcmi = rcmu(i)
            else
               lij2 = .true.
            endif

c --- loop over all chains j with j>i 
            do 99 j = i + 1, nchain
               
c ### check for simulation box ###
               if ( nboxi(j) .eq. ibox ) then
                  
                  jmolty = moltyp(j)
                  lqjmol = lelect(jmolty)
                  fxcmi = 0.0d0
                  fycmi = 0.0d0
                  fzcmi = 0.0d0
                  if ( lcutcm ) then
c                    --- check if ctrmas within rcmsq
                     rxuij = xcmi - xcm(j)
                     ryuij = ycmi - ycm(j)
                     rzuij = zcmi - zcm(j)
c                    --- minimum image the ctrmas pair separations
                     if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)
                     rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                     rcm = rbcut + rcmi + rcmu(j)
                     rcmsq = rcm*rcm
                     if ( rijsq .gt. rcmsq ) then
                        if (lqimol .and. lqjmol .and. lchgall) then
                           lij2 = .false.
                           goto 108
                        else
                           goto 99
                        endif
                     else
                        lij2 = .true.
                     endif
                  endif


c --- loop over all beads ii of chain i 
 108              do 98 ii = 1, nunit(imolty)

                     ntii = ntype(imolty,ii)

                     rxui = rxu(i,ii)
                     ryui = ryu(i,ii)
                     rzui = rzu(i,ii)

c --- loop over all beads jj of chain j 
                     do 97 jj = 1, nunit(jmolty)

                        ntjj = ntype(jmolty,jj)
                        if ( lij2 ) then
                           if ( (.not. (lij(ntii) .and. lij(ntjj))) 
     &                          .and. 
     &                          (.not. (lqchg(ntii) .and. lqchg(ntjj)))) 
     &                          goto 97
                        else
                           if (.not. (lqchg(ntii) .and. lqchg(ntjj)))
     &                          goto 97
                        endif
                        if ( lexpsix .or. lmmff ) then
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
                        if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)

c     write(2,*) 'bead ruij',rxuij,ryuij,rzuij
                        
                        rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                        
c --- compute whether the charged groups will interact & fij
                        fij = 0.0d0
                        if ( lchgall ) then
                           if ( lewald ) then
                              fij = fij - (2.0d0*calp(ibox)*
     &                             dexp(-calp(ibox)*calp(ibox)
     &                             *rijsq)/dsqrt(onepi)
     &                             +erfunc(calp(ibox)*dsqrt(rijsq))/
     &                             dsqrt(rijsq))*qqu(i,ii)*qqu(j,jj)*
     &                             qqfact/rijsq
                           endif
                        elseif ( lqimol .and. lqjmol .and. 
     &                          lqchg(ntii) .and. lqchg(ntjj) ) then
                           
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
                                 fij = fij - (2.0d0*calp(ibox)*
     &                                dexp(-calp(ibox)*calp(ibox)
     &                                *rijsq)/dsqrt(onepi)
     &                                +erfunc(calp(ibox)*dsqrt(rijsq))/
     &                                dsqrt(rijsq))*qqu(i,ii)*qqu(j,jj)*
     &                                qqfact/rijsq
                              else
                                 fij = -(qqfact*qqu(i,ii)*qqu(j,jj)
     &                                /dsqrt(rijsq))/rijsq
                              endif
                           endif
                        endif

                        if ( rijsq .lt. rcutsq .or. lijall) then
                           if ( lexpsix ) then
                              rij = dsqrt(rijsq)
                              fij = fij + (-6.0d0*aexsix(ntij)/(rijsq
     &                             *rijsq*rijsq) + bexsix(ntij)
     &                             *cexsix(ntij)*rij*dexp( cexsix(ntij)
     &				   *rij ))/rijsq
c                              write(2,*) 'rij,fij,ntij',rij,fij,ntij
                           elseif ( lmmff ) then
                              rs2 = rijsq/(sigisq(ntij))
                              rs1 = dsqrt(rs2)
                              rs6 = rs2*rs2*rs2
                              rs7 = rs1*rs6
                              sr1 = 1.07d0/(rs1+0.07d0)
                              sr7 = sr1**7.0d0
                              fij = fij - 7.0d0*epsimmff(ntij)*rs1*sr7*(
     +                           sr1*(1.12d0/(rs7+0.12d0)-2.0d0)/1.07d0+
     +                           1.12d0*rs6/((rs7+0.12d0)*(rs7+0.12d0)))
     +                             /rijsq
                           elseif (lninesix) then
                              rij = dsqrt(rijsq)
                              fij = fij + ( 72.0d0*epsnx(ntij)/
     &                            (rij*rzero(ntij)) ) * 
     &                            (rzero(ntij)/rij)**7 *
     &                            (1.0d0-(rzero(ntij)/rij)**3)
                           else
                              if ( lfepsi ) then
                                 sr2 = 1.0d0/rijsq
                                 sr6 = sr2 * sr2 * sr2
                                 if ( (.not. lqchg(ntii)) .and. 
     &                                (.not. lqchg(ntjj)) ) then
                                    if ( nunit(imolty) .eq. 4 ) then
c *** TIP-4P structure (temperary use ???)
                                       qave = (qqu(i,4)+qqu(j,4))/2.0d0
                                    else
                                       qave = (qqu(i,4)+qqu(i,5)+
     &                                      qqu(j,4)+qqu(j,5))*0.85d0
                                    endif
                                 else
                                    qave = (qqu(i,ii)+qqu(j,jj))/2.0d0
                                 endif

                                 if ( lexpand(imolty) 
     &                                .and. lexpand(jmolty)) then
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsilon(jmolty,jj))
                                 elseif ( lexpand(imolty)) then
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsi(ntjj))
                                 elseif ( lexpand(jmolty) ) then
                                    epsilon2=dsqrt(epsilon(jmolty,jj)
     &                                   *epsi(ntii))
                                 else
                                    epsilon2=epsij(ntij)
                                 endif
                                 flj = 48.0d0*epsilon2*
     &                                sr6*(-sr6*(aslope*(qave-a0)*
     &                                (qave-a0)+ashift)+0.5d0*
     &                                (bslope*(qave-b0)*(qave-b0)+
     &                                bshift)) / rijsq
                              else
                                 if ( lexpand(imolty) 
     &                                .and. lexpand(jmolty) ) then
                                    sigma2=(sigma(imolty,ii)+
     &                                   sigma(jmolty,jj))/2.0d0
                                    sr2 = sigma2*sigma2/rijsq
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsilon(jmolty,jj))
                                 elseif ( lexpand(imolty) ) then
                                    sigma2=(sigma(imolty,ii)+
     &                                   sigi(ntjj))/2.0d0
                                    sr2 = sigma2*sigma2/rijsq
                                    epsilon2=dsqrt(epsilon(imolty,ii)
     &                                   *epsi(ntjj))
                                 elseif ( lexpand(jmolty) ) then
                                    sigma2=(sigma(jmolty,jj)+
     &                                   sigi(ntii))/2.0d0
                                    sr2 = sigma2*sigma2/rijsq
                                    epsilon2=dsqrt(epsilon(jmolty,jj)
     &                                   *epsi(ntii))
                                 else
                                    sr2 = sig2ij(ntij) / rijsq
                                    epsilon2 = epsij(ntij)
                                 endif
                                 sr6 = sr2 * sr2 * sr2
                                 flj = 48.0d0*sr6*(-sr6+0.5d0)
     &                                *epsilon2 / rijsq
                              endif
                              fij = fij + flj
                           endif
                        endif
                        fxcmi = fxcmi + fij * rxuij 
                        fycmi = fycmi + fij * ryuij 
                        fzcmi = fzcmi + fij * rzuij 
 97                  continue
                      
 98               continue


c --- calculate distance between c-o-m ---
                  rxuij = xcmi - xcm(j)
                  ryuij = ycmi - ycm(j)
                  rzuij = zcmi - zcm(j)

c *** minimum image the pair separations ***
                  if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)

c                  write(2,*) 'COM  ruij',rxuij,ryuij,rzuij

                  press = press + 
     +                 fxcmi*rxuij + fycmi*ryuij + fzcmi*rzuij

c * for surface tension                  
c * this is correct for the coulombic part and for LJ.  Note sign difference!
                  pxx = pxx - fxcmi*rxuij
                  pyy = pyy - fycmi*ryuij
                  pzz = pzz - fzcmi*rzuij
                  pips(1,2) = pips(1,2) - rxuij*fycmi
                  pips(1,3) = pips(1,3) - rxuij*fzcmi
                  pips(2,1) = pips(2,1) - ryuij*fxcmi
                  pips(2,3) = pips(2,3) - ryuij*fzcmi
                  pips(3,1) = pips(3,1) - rzuij*fxcmi
                  pips(3,2) = pips(3,2) - rzuij*fycmi
                  
               endif

 99         continue


         endif

 100  continue

c ################################################################

      if ( lewald ) then

c *** Compute the reciprocal space contribution
c *** by using the thermodynamic definition
c *** 
         call recippress(ibox,repress,rpxx,rpyy,rpzz,rpxy,rpyx,rpxz,
     &                  rpzx,rpyz,rpzy)
         press = press - repress

         pxx = pxx + rpxx
         pyy = pyy + rpyy
         pzz = pzz + rpzz
         pips(1,2) = pips(1,2) + qqfact*rpxy
         pips(1,3) = pips(1,3) + qqfact*rpxz
         pips(2,1) = pips(2,1) + qqfact*rpyx
         pips(2,3) = pips(2,3) + qqfact*rpyz
         pips(3,1) = pips(3,1) + qqfact*rpzx
         pips(3,2) = pips(3,2) + qqfact*rpzy

      endif
      pips(1,1) = pxx
      pips(2,2) = pyy
      pips(3,3) = pzz

      press = 1.380662d4 * ( ( (nchbox(ibox)+ ghost_particles(ibox))
     +      / beta) -
     +     ( press/3.0d0 ) ) / 
     +     ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
c      write (6,*) ' press' , press

      surf = pzz - 0.5d0*(pxx + pyy)
c * divide by surface area and convert from K to put surf in mN/m 
      surf = 1.380658d0*surf / (2.0d0*boxlx(ibox)*boxly(ibox))

c----check pressure tail correction

      if (ltailc) then
c --     add tail corrections for the Lennard-Jones energy

c --     Not adding tail correction for the ghost particles
c --     as they are ideal (no interaction) Neeraj.

         volsq = ( boxlx(ibox) * boxly(ibox) * boxlz(ibox) )**2
         do imolty=1, nmolty        
            do jmolty=1, nmolty  
               rhosq = ncmt(ibox,imolty)*ncmt(ibox,jmolty)
     +              / volsq
               press=press + 1.380662d4 * corp(imolty,jmolty,rhosq,ibox)
            enddo
        enddo
      endif

c      write (6,*) ' press tail' ,  press

      return
      end











