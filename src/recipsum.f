      subroutine recipsum(ibox,vrecip)

c recipsum
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

c    *********************************************************************
c    ** calculates the total reciprocal space ewald-sum term for volume **
c    ** moves, written in 1998 by Bin Chen.                             **
c    ** rewritten in 2001 by Bin Chen.                                  **
c    ** rewritten again, probably by Bin.                               **
c    *********************************************************************

      integer ibox,nkx,nky,nkz
     &     ,nkx_min,nkx_max,nky_min,nky_max,nkz_min,nkz_max
     &     ,i,ii,imolty,kmax1,ncount

c * from h-matrix formulation
      integer l,m,n,m_min,m_max,n_min,n_max,kmaxl,kmaxm,kmaxn

      double precision alpsqr4,vol,ksqr,sumr,sumi,arg,boxlen,vrecip,
     &     bx1,by1,bz1,bmin,xratio,yratio,zratio,
     &     constx,consty,constz,hmatik(9),kx1,ky1,kz1,hmaxsq,calpi

      include 'control.inc'
      include 'system.inc'
      include 'coord.inc'
      include 'ewaldsum.inc'
      include 'poten.inc'
      include 'conver.inc'
      include 'cell.inc'

c *** Set up the reciprocal space vectors ***

      ncount = 0
      vrecip = 0.0d0

      calpi = calp(ibox)

      if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
         bx1 = boxlx(ibox)
         by1 = boxly(ibox)
         bz1 = boxlz(ibox)
         hmat(ibox,1) = bx1
         hmat(ibox,5) = by1
         hmat(ibox,9) = bz1
         do i = 1,9
            hmatik(i) = 0.0d0
         enddo
         hmatik(1) = twopi/bx1
         hmatik(5) = twopi/by1
         hmatik(9) = twopi/bz1
         kmaxl = dint(bx1*calpi)+1
         kmaxm = dint(by1*calpi)+1
         kmaxn = dint(bz1*calpi)+1
      else
         do i = 1,9
            hmatik(i) = twopi*hmati(ibox,i)
         enddo
         kmaxl = dint(hmat(ibox,1)*calpi)+2
         kmaxm = dint(hmat(ibox,5)*calpi)+2
         kmaxn = dint(hmat(ibox,9)*calpi)+2
      endif
    
      alpsqr4 = 4.0d0*calpi*calpi
         
      vol = hmat(ibox,1)*
     &     (hmat(ibox,5)*hmat(ibox,9) - hmat(ibox,8)*hmat(ibox,6))
     &     + hmat(ibox,4)*
     &     (hmat(ibox,8)*hmat(ibox,3) - hmat(ibox,2)*hmat(ibox,9))
     &     + hmat(ibox,7)*
     &     (hmat(ibox,2)*hmat(ibox,6) - hmat(ibox,5)*hmat(ibox,3))

      vol = vol/(4.0d0*onepi)

      hmaxsq = alpsqr4*onepi*onepi
         
c *** generate the reciprocal-space 

      do l = 0,kmaxl
         if ( l .eq. 0 ) then
            m_min = 0
         else
            m_min = -kmaxm
         endif
         do m = m_min, kmaxm
            if (l .eq. 0 .and. m .eq. 0) then
               n_min = 1
            else
               n_min = -kmaxn
            endif
            do n = n_min, kmaxn
               kx1 = dble(l)*hmatik(1)+dble(m)*hmatik(2)+
     &              dble(n)*hmatik(3)
               ky1 = dble(l)*hmatik(4)+dble(m)*hmatik(5)+
     &              dble(n)*hmatik(6)
               kz1 = dble(l)*hmatik(7)+dble(m)*hmatik(8)+
     &              dble(n)*hmatik(9)
               ksqr = kx1*kx1+ky1*ky1+kz1*kz1
c               if ( ksqr .lt. hmaxsq ) then
c --- sometimes these are about equal, which can cause different
c --- behavior on 32 and 64 bit machines without this .and. statement
               if ( ksqr .lt. hmaxsq .and.
     &              abs(ksqr-hmaxsq) .gt. 1d-9 ) then
                  ncount = ncount + 1
                  kx(ncount,ibox) = kx1
                  ky(ncount,ibox) = ky1
                  kz(ncount,ibox) = kz1
                  prefact(ncount,ibox) =
     &                 dexp(-ksqr/alpsqr4)/(ksqr*vol)
c     *** sum up q*cos and q*sin ***
                  sumr = 0.0d0
                  sumi = 0.0d0
                  do 10 i = 1,nchain
                     imolty = moltyp(i)
                     if ( .not. lelect(imolty) ) goto 10
                     if ( nboxi(i) .eq. ibox) then
                        do ii = 1,nunit(imolty)
                           if ( lqchg(ntype(imolty,ii)) ) then
                              arg = kx1*rxu(i,ii) +
     &                             ky1*ryu(i,ii) +
     &                             kz1*rzu(i,ii)
                              sumr = sumr + dcos(arg)*qqu(i,ii)
                              sumi = sumi + dsin(arg)*qqu(i,ii)
                           endif
                        enddo
                     endif
 10               continue
                  ssumr(ncount,ibox) = sumr
                  ssumi(ncount,ibox) = sumi
c *** Potential energy ***
                  vrecip = vrecip + (sumr*sumr + sumi*sumi)
     &                 * prefact(ncount,ibox)

               endif
            enddo
         enddo
      enddo
      vrecip = vrecip*qqfact
c      write(iou,*) 'in recipsum:',ssumr(100,ibox),ibox
c *** safety check ***
c      write(iou,*) 'A total of ',ncount,' vectors are used'
      if ( ncount .gt. vectormax ) then
         write(iou,*)  'choose a larger vectormax'
         stop
      endif

      numvect(ibox) = ncount
         
      return

      end


