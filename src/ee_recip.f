      subroutine ee_recip(ibox,vrecipnew,vrecipold,type)

c recip
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

c    *********************************************************************
c    ** calculates the reciprocal ewald-sum term for trans, rot, flucq, **
c    ** swatch and swap moves, and update the reciprocal ewald-sum.     **
c    ** rewritten on June 25/99 by Bin Chen.                            **
c    *********************************************************************

      implicit none
      integer ic,zz,ii,imolty,ibox,ncount,type
      double precision vrecipnew,vrecipold,sumr(2),sumi(2),arg
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ewaldsum.inc'
      include 'poten.inc'

      ncount = numvect(ibox)

      if ( type .eq. 1 ) then
         
c *** recalculate the reciprocal space part for one-particle move, translation,
c *** rotation, swap, flucq, and swatch.
c *** old conformation zz = 1 (which is 0 for swap inserted molecule)
c *** new conformation zz = 2 (which is 0 for swap removed molecule)

c         write(2,*) 'in recip:',moltion(1),moltion(2)
c         do zz = 1,2
c            imolty = moltion(zz)
c            do ii = 1, nunit(imolty)
c               write(2,*) rxuion(ii,zz),ryuion(ii,zz),rzuion(ii,zz),
c     &              qquion(ii,zz)
c            enddo
c         enddo

         do 30 ic = 1, ncount
            do 20 zz = 1,2
c --- zz = 1: old configuration 
c --- zz = 2: new configuration

               sumr(zz) = 0.0d0
               sumi(zz) = 0.0d0
               imolty = moltion(zz)
               do ii = 1, nunit(imolty)
c                  if ( lqchg(ntype(imolty,ii)) ) then
                     arg = kx(ic,ibox)*rxuion(ii,zz) +
     &                    ky(ic,ibox)*ryuion(ii,zz) +
     &                    kz(ic,ibox)*rzuion(ii,zz)
                     sumr(zz) = sumr(zz) + 
     &                    qquion(ii,zz)*dcos(arg)
                     sumi(zz) = sumi(zz) + 
     &                    qquion(ii,zz)*dsin(arg)
c                  endif
               enddo
 20         continue
            ssumrn(ic,ibox) = ssumr(ic,ibox) - sumr(1)
     &           + sumr(2)
            ssumin(ic,ibox) = ssumi(ic,ibox) - sumi(1)
     &           + sumi(2)
 30      continue
         vrecipnew = 0.0d0
         vrecipold = 0.0d0
         do ic = 1,ncount
            vrecipnew = vrecipnew + (ssumrn(ic,ibox)*
     &           ssumrn(ic,ibox) + ssumin(ic,ibox)*
     &           ssumin(ic,ibox))*prefact(ic,ibox)
            vrecipold = vrecipold + (ssumr(ic,ibox)*
     &           ssumr(ic,ibox) + ssumi(ic,ibox)*
     &           ssumi(ic,ibox))*prefact(ic,ibox)
         enddo

         vrecipnew = vrecipnew*qqfact
         vrecipold = vrecipold*qqfact

      elseif (type .eq. 2) then

c *** update the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumrn(ic,ibox)
            ssumi(ic,ibox) = ssumin(ic,ibox)
         enddo

      elseif (type .eq. 3) then

c *** store the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumro(ic,ibox) = ssumr(ic,ibox)
            ssumio(ic,ibox) = ssumi(ic,ibox)
         enddo

      elseif (type .eq. 4) then

c *** restore the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumr(ic,ibox) = ssumro(ic,ibox)
            ssumi(ic,ibox) = ssumio(ic,ibox)
         enddo

      endif

c      write(2,*) 'in recip:',ssumr(100,ibox),ibox,ssumrn(100,ibox)

      return 
      end

