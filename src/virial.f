      subroutine virial(binvir,binvir2,nvirial,starvir,stepvir)

c virial
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
c    ** computes the 2nd virial coefficient                           **
c    ** B2(T) = -2Pi Int 0toInf [ Exp[-beta*u(r)] -1] r^2 dr          **
c    ** Using the trapazoid method of numerical integration give      **
c    ** B2(T) = -2*Pi*stepvir* Sum(i=2,n-1)[ <Exp[-beta*u(r)]-1> ri^2 **
c    ** + 1/2( <Exp[-beta*u(r1)]-1> r1 + <Exp[-beta*u(rn)]-1> rn      **
c    ** Marcus Martin 1-15-97                                         **
c    *******************************************************************
 
      implicit none

c *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'poten.inc'
      include 'conver.inc' 
      include 'external.inc'
      include 'connect.inc'
      include 'fepsi.inc'
      include 'inputdata.inc'
      include 'eepar.inc'

      integer i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij
     +     ,nnn,nvirial,ip,itemp,iii
      double precision vinter,rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij
     &       ,rijsq,sr2, sr6 ,velect,mayer

      double precision ljsami,ljpsur,ljmuir,xdiff,ydiff,zdiff
     &     ,dvircm,stepvir,exsix,starvir
      double precision binvir
      dimension binvir(maxvir,maxntemp),mayer(maxntemp)

      logical ovrlap,lqimol,lqjmol,lmatrix
  
c *** only use for polarizable models
      integer chgmax
      parameter (chgmax=10)
      integer ip1,ip2,iunit,numchg,info,ipiv(chgmax),mainsite(2,2),
     &     lam1,lam2
      double precision a(chgmax,chgmax),b2(chgmax,1),mainxiq(2,2)

      double precision consa1,consa2,consb1,consb2,selfadd1,selfadd2,
     &     vtotal,epsilon2,vmin

      double precision mass_t,binvir2(maxvir,maxntemp),
     &     factor,corr,vold,deri_u

C --------------------------------------------------------------------

c      write(2,*) 'start VIRIAL'
c      write(2,*) 'binvir',binvir
      rminsq = rmin * rmin
      vmin = -2000.0d0

      if ( lfepsi ) then
         consa1 = 4.0d0*aslope*a0*(2.0d0/1.6d0)
         consb1 = 4.0d0*bslope*b0*(2.0d0/1.6d0)
         consa2 = 2.0d0*aslope*(4.0d0/2.56d0)
         consb2 = 2.0d0*bslope*(4.0d0/2.56d0)
      endif
      

      if ( nboxi(1) .eq. nboxi(2) ) stop 'particles found in same box'

 
c ################################################################

C *******************************
C *** INTERCHAIN INTERACTIONS ***
C *******************************

      xdiff = xcm(2) - xcm(1)
      ydiff = ycm(2) - ycm(1)
      zdiff = zcm(2) - zcm(1)
      dvircm = starvir
      i = 1
      imolty = moltyp(i)
      iunit = nunit(imolty)
      mass_t = 0.0d0
      do ii = 1, iunit
         mass_t = mass_t + mass(ntype(imolty,ii))
      enddo
      mass_t = mass_t/1000d0
      factor = -(6.6260755d-34)**2*6.0221367d23*1d20 / 
     &     (24.0d0*onepi*mass_t*1.380658d-23*twopi)

      lqimol = lelect(imolty)

      if ( lflucq(imolty) .and. lflucq(moltyp(2))) then
         lmatrix = .true.

c *** count charge site for type 1 molecule
         numchg = 0
         do ii = 1, iunit
            ntii = ntype(imolty,ii)
            if ( lqchg(ntii) ) numchg = numchg + 1
         enddo
      else
         lmatrix = .false.
      endif

      do nnn = 1,nvirial

         if ( lmatrix ) then
c *** initialize matrix
            do ii = 1,chgmax
               do jj = 1,chgmax
                  a(ii,jj) = 0
               enddo
               b2(ii,1) = 0
            enddo
         endif

         ovrlap = .false.
         vinter = 0.0d0
         velect = 0.0d0
c --- loop over all beads ii of chain i
         ip1 = 0
         do 99 ii = 1, nunit(imolty) 

            ntii = ntype(imolty,ii)
            if ( lqchg(ntii) ) ip1 = ip1 + 1
            rxui = rxu(i,ii)  + xdiff - dvircm
            ryui = ryu(i,ii)  + ydiff
            rzui = rzu(i,ii)  + zdiff

c --- loop over chain 2 
            j = 2
            jmolty = moltyp(j)
            lqjmol = lelect(jmolty)
c     --- loop over all beads jj of chain j 
            ip2 = numchg
            do 97 jj = 1, nunit(jmolty) 
c     --- check exclusion table
               if ( lexclu(imolty,ii,jmolty,jj) ) goto 97
               
               ntjj = ntype(jmolty,jj)
               if ( lqchg(ntjj) ) ip2 = ip2 + 1
               
               if (lexpsix) then
                  ntij = (ntii+ntjj)/2
               else
                  ntij = (ntii-1)*nntype + ntjj
               endif
               if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
               
               rxuij = rxui - rxu(j,jj)
               ryuij = ryui - ryu(j,jj)
               rzuij = rzui - rzu(j,jj)

               rijsq = (rxuij*rxuij)+(ryuij*ryuij)
     &              + (rzuij*rzuij)
               if ( lmatrix .and. lqchg(ntii) .and. lqchg(ntjj) ) then
                  a(ip1,ip2) = qqfact/dsqrt(rijsq)
                  a(ip2,ip1) = a(ip1,ip2)
               endif

               if ( rijsq .lt. rminsq ) then
                  ovrlap = .true.
                  goto 100
               else
                  if ( lsami ) then
                     vinter = vinter + ljsami(rijsq,ntij)
                  elseif (lexpsix) then
                     vinter = vinter + exsix(rijsq,ntij)
                  elseif ( lmuir ) then
                     vinter = vinter + ljmuir(rijsq,ntij)
                  elseif ( lpsurf ) then
                     vinter = vinter + ljpsur(rijsq,ntij)
                  else if (lshift) then
                     sr2 = sig2ij(ntij) / rijsq
                     sr6 = sr2 * sr2 * sr2
                     vinter = vinter + 
     +                    sr6*(sr6-1.0d0)*epsij(ntij)-ecut(ntij) 
                  elseif ( lfepsi ) then
                     if ( lij(ntii) .and. lij(ntjj) ) then
                        sr6 = rijsq*rijsq*rijsq
                        epsilon2 = epsij(ntij)
                        selfadd1 = epsilon2*(consa1/sr6-consb1)/sr6
                        selfadd2 = epsilon2*(consa2/sr6-consb2)/sr6
                     endif
                  else
                     sr2 = sig2ij(ntij) / rijsq
                     sr6 = sr2 * sr2 * sr2
                     vinter = vinter + 
     +                    sr6*(sr6-1.0d0)*epsij(ntij)
                     
                  endif
                  
                  if ( lqimol .and. lqjmol .and. .not. lmatrix ) then
                     velect = velect + qqfact*qqu(i,ii)*qqu(j,jj)
     &                    /dsqrt(rijsq)
                  endif
                  
               endif
               
 97         continue
 99      continue

c *** use Matrix Minimization to solve the equilibrium charge

         if ( lmatrix ) then

            iii = 2
            do ip = 1,2
               imolty = moltyp(ip)
               iunit = nunit(imolty)
               ip1 = ( ip-1)*numchg
               do ii = 1, iunit
                  ntii = ntype(imolty,ii)
                  if ( lqchg(ntii) ) then
                     ip1 = ip1 + 1
                     a(ip1,ip1) = 2.0d0*jayself(ntii)
                     b2(ip1,1) = -xiq(ntii)
                     if ( dabs(xiq(ntii)) .gt. 1.0d-10 ) then
                        if ( iii .eq. 1 ) then
                           iii = 2
                        elseif ( iii .eq. 2 ) then
                           iii = 1
                        endif
                        mainsite(ip,iii) = ip1
                        mainxiq(ip,iii) = -xiq(ntii)
                     endif
                     ip2 = ip1
                     do jj = ii+1,iunit
                        ntjj = ntype(imolty,jj)
                        if ( lqchg(ntjj) ) then
                           ip2 = ip2 + 1
                           ntij = (ntii-1)*nntype + ntjj
                           a(ip1,ip2) = jayq(ntij)
                           a(ip2,ip1) = jayq(ntij)
                        endif
                     enddo
                  endif
               enddo
            enddo

c *** undetermined multiplier
            if ( lqtrans(imolty) ) then
               lam1 = ip2 + 1
               do ii = 1, numchg*2
                  a(lam1,ii) = 1
                  a(ii,lam2) = 1
               enddo
            else
               lam1 = ip2 + 1
               lam2 = ip2 + 2
               do ii = 1, numchg
                  a(lam1,ii) = 1
                  a(ii,lam1) = 1
               enddo
               do ii = numchg+1,ip2
                  a(lam2,ii) = 1
                  a(ii,lam2) = 1
               enddo
            endif
            
            if ( lfepsi ) then
               ip1 = mainsite(1,1)
               ip2 = mainsite(2,1)
               b2(ip1,1) = b2(ip1,1) + selfadd1
               b2(ip2,1) = b2(ip2,1) + selfadd1
               a(ip1,ip1) = a(ip1,ip1) + selfadd2
               a(ip1,ip2) = a(ip1,ip2) + selfadd2
               a(ip2,ip1) = a(ip1,ip2)
               a(ip2,ip2) = a(ip2,ip2) + selfadd2
               if ( nunit(imolty) .eq. 5 ) then
                  ip1 = mainsite(1,2)
                  ip2 = mainsite(2,2)
                  b2(ip1,1) = b2(ip1,1) + selfadd1
                  b2(ip2,1) = b2(ip2,1) + selfadd1
                  a(ip1,ip1) = a(ip1,ip1) + selfadd2
                  a(ip1,ip2) = a(ip1,ip2) + selfadd2
                  a(ip2,ip1) = a(ip1,ip2)
                  a(ip2,ip2) = a(ip2,ip2) + selfadd2
                  ip1 = mainsite(1,1)
                  ip2 = mainsite(1,2)
                  a(ip1,ip2) = a(ip1,ip2) + selfadd2
                  a(ip2,ip1) = a(ip1,ip2)
                  ip2 = mainsite(2,2)
                  a(ip1,ip2) = a(ip1,ip2) + selfadd2
                  a(ip2,ip1) = a(ip1,ip2)
                  ip1 = mainsite(2,1)
                  a(ip1,ip2) = a(ip1,ip2) + selfadd2
                  a(ip2,ip1) = a(ip1,ip2)
                  ip2 = mainsite(1,2)
                  a(ip1,ip2) = a(ip1,ip2) + selfadd2
                  a(ip2,ip1) = a(ip1,ip2)
               endif
            endif

c *** commenting this out JMS 6-20-00 ***
c            call dgesv(chgmax,1,a,chgmax,ipiv,b2,chgmax,info)
c *** 


c            write(21,*) b2
c            do ii = 1, nunit(imolty)
c               write(21,*) rxu(i,ii)+xdiff-dvircm,ryu(i,ii)+ydiff
c     &              ,rzu(i,ii)+zdiff
c            enddo
c            do ii = 1, nunit(imolty)
c               write(21,*) rxu(2,ii),ryu(2,ii),rzu(2,ii)
c            enddo
c            write(11,*) sr6**(1.0/6.0d0),b2(1,1),b2(2,1)

            velect = 0.0d0
            do ip = 1, 2
               ii = mainsite(ip,1)
               if ( lfepsi ) then
                  velect = velect-0.5d0*b2(ii,1)*(mainxiq(ip,1)
     &                 +selfadd1)
c                  write(21,*) b2(ii,1),mainxiq(ip,1)
                  vinter = epsilon2*((aslope*a0*a0+ashift)/sr6
     &                 -(bslope*b0*b0+bshift))/sr6
               else
                  velect = velect-0.5d0*b2(ii,1)*mainxiq(ip,1)
               endif
               if ( nunit(imolty) .eq. 5 ) then
                  ii = mainsite(ip,2)
                  if ( lfepsi ) then
                     velect = velect-0.5d0*b2(ii,1)*(mainxiq(ip,2)
     &                    +selfadd1)
c                     write(21,*) b2(ii,1),mainxiq(ip,2)
                     vinter = epsilon2*((aslope*a0*a0+ashift)/sr6
     &                    -(bslope*b0*b0+bshift))/sr6
                  else
                     velect = velect-0.5d0*b2(ii,1)*mainxiq(ip,2)
                  endif
               endif
            enddo
            velect = velect - fqegp(moltyp(1)) - fqegp(moltyp(2))
            if ( velect + vinter*4.0d0 .lt. vmin ) then
               vmin = velect + vinter*4.0d0
               write(11,*) 'vmin:',vmin
               do ii = 1, nunit(imolty)
                  write(11,*) rxu(i,ii)+xdiff-dvircm,ryu(i,ii)+ydiff
     &                 ,rzu(i,ii)+zdiff
               enddo
               do ii = 1, nunit(imolty)
                  write(11,*) rxu(2,ii),ryu(2,ii),rzu(2,ii)
               enddo
               rewind(11)
            endif
c            write(11,*) 'energy:',velect+vinter*4.0d0

         endif
                     

 100     if ( ovrlap ) then
            do itemp = 1, ntemp
               mayer(itemp) = -1.0d0
            enddo
         else
            if (.not.lsami .and. .not.lexpsix) vinter = 4.0d0*vinter
            do itemp = 1,ntemp
               mayer(itemp)
     &              = dexp(-(vinter+velect)/virtemp(itemp))-1.0d0
            enddo
         endif
c         write(2,*) 'mayer',mayer
c         write(2,*) 'nnn',nnn
         do itemp = 1,ntemp
            binvir(nnn,itemp) = binvir(nnn,itemp) + mayer(itemp)
         enddo
         if ( nnn .eq. 1 ) then
            corr = 0.0d0
            vold = vinter + velect
         else
            deri_u = (vinter+velect-vold)/stepvir
            corr = deri_u * deri_u * factor
            vold = vinter + velect
         endif
         do itemp = 1, ntemp
            binvir2(nnn,itemp) = binvir2(nnn,itemp) +
     &           (mayer(itemp)+1.0d0)*corr/(virtemp(itemp)**3)
         enddo
         dvircm = dvircm + stepvir
      enddo


c ################################################################

c      write(2,*) 'binvir',binvir

c      write(2,*) 'end VIRIAL'

      return
      end









