      subroutine nrtab 

c nrtab
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

C    WILL ONLY WORK FOR PURE SYSTEMS OF HOMOPOLYMERS !!!!! 
c    *******************************************************************
c    **     calculates a look-up table for the distribution of        **
c    **          internal bead-bead distances as needed by            **
c    **           fixed-end-point configurational-bias MC             **
c    *******************************************************************
c    ** not implemented for new version of CBMC

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'poten.inc'
      include 'connect.inc'
      include 'nrtab.inc'

      logical lpoly,lcone
c     lpoly used to have a meaning in iclude files but does not any longer
      integer i,j,iunit,iutry,iincre,istart,iend,idiff,nrtc,iu,
     +        iulast,iuprev,iuppre,ivib,iben,itor,ibin,it

      double precision xprev,yprev,zprev,dprev
     +                ,xpprev,ypprev,zpprev,xa1a2,ya1a2,za1a2
     +                ,da1a2,xub,yub,zub,va,vbba,vdha,vtorso
     +                ,x,y,z,dlast,xi1,xi2, xisq,thetac,theta
     +                ,random,xaa1,yaa1,zaa1,daa1,dot,bf,rij,dum

      double precision rxnew(numax),rynew(numax),rznew(numax)
      double precision nnrtab(nrtmax),hnrtab(nrtmax,nrtbin)

c     lcone has no meaning any longer 2-16-98

c ------------------------------------------------------------------

c we will generate NRTCON configurations 
c    of maximum length NRTMAX
c    which will be binned into NRTBIN bins

c ------------------------------------------------------------------

c zero the bins
      do i = 1, nrtmax
         nnrtab(i) = 0
         do j = 1, nrtbin
            hnrtab(i,j) = 0
         enddo
      enddo

c ------------------------------------------------------------------

c calculate the bin size
      do i = 2, nrtmax
         if ( lpoly ) then
            dmrtab(i) = dble(i) * brvib(1)
         else
            dmrtab(i) = dble(i) * brvib(1) * dsin( 0.5d0 * brben(1) )
            dmrtab(i) = dmrtab(i) * 1.1d0
         endif 
         dbrtab(i) = dmrtab(i) / dble(nrtbin)
      enddo

c ------------------------------------------------------------------

c *** store number of units in iunit ***
      iunit = nunit(1)
c *** select starting unit ***
      iutry = 1
      rxnew(1) = 0.0d0
      rynew(1) = 0.0d0
      rznew(1) = 0.0d0
      iincre = 1
      istart = iutry + 1
      iend = iunit

c ------------------------------------------------------------------

      do 999 nrtc = 1, nrtcon

c    *******************************************************************
c    **     performs a configurational-bias move for ideal chain      **
c    **         with all bonded-intramolecular interactions           **
c    *******************************************************************
 
      do 200 iu = istart, iend, iincre

         iulast = iu - iincre
         iuprev = iulast - iincre
         iuppre = iuprev - iincre

         ivib = iulast
         if ( iuprev .ge. 1 ) then
            iben = iuprev
            if ( iuppre .ge. 1 ) then
               itor = iuppre
            else
               itor = 0
            endif
         else
            iben = 0
            itor = 0
         endif

         if ( iben .gt. 0 .and. .not. lpoly ) then
c ---       vector from last to previous unit ---
            xprev = rxnew(iuprev) - rxnew(iulast)
            yprev = rynew(iuprev) - rynew(iulast)
            zprev = rznew(iuprev) - rznew(iulast)
            dprev = brvib(1)
            if ( itor .gt. 0 ) then
c ---          vector from previous to preprevious unit ---
               xpprev = rxnew(iuppre) - rxnew(iuprev)
               ypprev = rynew(iuppre) - rynew(iuprev)
               zpprev = rznew(iuppre) - rznew(iuprev)
c ***          calculate cross products d_a-1 x d_a-2 ***
               xa1a2 = yprev*zpprev - zprev*ypprev
               ya1a2 = zprev*xpprev - xprev*zpprev
               za1a2 = xprev*ypprev - yprev*xpprev
c ***          calculate lengths of cross products ***
               da1a2 = dsqrt ( xa1a2**2 + ya1a2**2 + za1a2**2 )
            endif
         endif
	 
c ---set up the cone
         if (lcone) then
            xub = -(xprev / dprev)
            yub = -(yprev / dprev)
            zub = -(zprev / dprev)
            call cone (1,xub,yub,zub,brben(1),dum,dum,dum )
	 endif
	    
c *** select ONE trial position ***
c --- set energies of trial position to zero ---
         vbba = 0.0d0
         vdha = 0.0d0

 108     if ( iben .gt. 0 .and. lcone ) then
            call cone (2,dum,dum,dum,dum,x,y,z )
            x = brvib(1) * x
            y = brvib(1) * y
            z = brvib(1) * z
            dlast = brvib(1)
         else
c --- calculate random vector on the unit sphere ---
 109        xi1 = ( 2.0d0 * random() ) - 1.0d0
            xi2 = ( 2.0d0 * random() ) - 1.0d0
            xisq = xi1**2 + xi2**2
            if ( xisq .lt. 1.0d0 ) then
               if ( brvibk(1) .gt. 0.1d0 ) 
     +              stop 'bond vibrations not implemented'
               x = brvib(1) * 2.0d0 * xi1 * dsqrt( 1.0d0 - xisq )
               y = brvib(1) * 2.0d0 * xi2 * dsqrt( 1.0d0 - xisq )
               z = brvib(1) * ( 1.0d0 - 2.0d0 * xisq )
               dlast = brvib(1)
            else
               goto 109
            endif
         endif

         if ( .not. lpoly ) then
            if ( iben .gt. 0 ) then
               if ( lcone ) then
                  vbba = 0.0d0
               else
c *** calculate the bond bending potential energy ***
                  thetac = ( x*xprev+y*yprev+z*zprev ) / (dlast*dprev)
                  theta = dacos(thetac)
                  vbba = brbenk(1) * ( theta - brben(1) )**2
               endif

               if ( itor .gt. 0 ) then
c *** calculate the dihedral energy ***
c *** calculate cross products d_a x d_a-1 ***
                  xaa1 = y*(-zprev) + z*yprev
                  yaa1 = z*(-xprev) + x*zprev
                  zaa1 = x*(-yprev) + y*xprev
c *** calculate lengths of cross products ***
                  daa1 = dsqrt ( xaa1**2 + yaa1**2 + zaa1**2 )
c *** calculate dot product of cross products ***
                  dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                  thetac = -(dot / ( daa1 * da1a2 ))
                  it = 0
                  vdha = vtorso( thetac, it )
               endif
            endif

            va = vbba + vdha
            bf = dexp ( -(va * beta) )
            if ( random() .ge. bf ) go to 108 
c            write(6,*) 'new ip', ip, '   va', va

         else
c *** poly-bead molecule ***
            va = 0.0d0
         endif

         rxnew (iu) = rxnew(iulast) + x
         rynew (iu) = rynew(iulast) + y
         rznew (iu) = rznew(iulast) + z

200      continue

c ------------------------------------------------------------------

c ***************************
c * analyse the trial chain *
c ***************************

         do i = 1, iunit - 2
            do j = i + 2, iunit

               idiff = j - i

               if ( idiff .le. nrtmax ) then
                  nnrtab(idiff) = nnrtab(idiff) + 1
                  x = rxnew(i) - rxnew(j)
                  y = rynew(i) - rynew(j)
                  z = rznew(i) - rznew(j)
                  rij = dsqrt( x*x + y*y + z*z )
                  if ( rij .gt. dmrtab(idiff) ) then
                     write(6,*) 'WARNING NRTAB: rij .gt. dmrtab'
                     write(6,*) 'idiff',idiff,
     &                    'rij',rij,'dmrtab',dmrtab(idiff)
                  endif
                  ibin = idint( rij / dbrtab(idiff) ) + 1
                  if ( ibin .gt. nrtbin ) ibin = nrtbin
                  hnrtab(idiff,ibin) = hnrtab(idiff,ibin) + 1

                  if ( idiff .le. 3 .and. ibin .eq. 1 ) then
                     write(6,*) 'nrtc',nrtc,'i',i,'j',j
                     write(6,*) 'ri',rxnew(i),rynew(i),rznew(i)
                     write(6,*) 'rj',rxnew(j),rynew(j),rznew(j)
                     write(6,*) 'rij',rij
                  endif

               endif

            enddo
         enddo

c ------------------------------------------------------------------

 999  continue

c ------------------------------------------------------------------

c *************************
c * calculate the weights *
c *************************

      do i = 2, nrtmax
         iu = 20 + i
         do j = 1, nrtbin
            wnrtab(i,j) = dble(hnrtab(i,j)) / dble(nnrtab(i))
            rij = (dble(j)-0.5d0)*dbrtab(i)
            write(iu,'(2x,f10.4,f12.8)') rij, wnrtab(i,j)
         enddo
      enddo

      write(21,*) nrtmax, nrtbin
      write(21,*) dbrtab
      write(21,*) dmrtab
      write(21,*) wnrtab

c ------------------------------------------------------------------

      return
      end
