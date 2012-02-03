      subroutine inclus( inclnum,inclmol,inclbead,inclsign,ncarbon,
     &     ainclnum,ainclmol,ainclbead,a15t)

c inclus
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

      include 'control.inc'
      include 'coord.inc'
      include 'connect.inc'
      include 'poten.inc'

      integer m,n,nb,mb,imolty,ioffset
      integer inclnum,inclmol,inclbead,inclsign,ncarbon
      integer ainclnum,ainclmol,ainclbead,a15t
      dimension inclmol(ntmax*numax*numax),inclsign(ntmax*numax*numax)
      dimension inclbead(ntmax*numax*numax,2),ncarbon(ntmax)
      dimension ainclmol(ntmax*numax*numax)
      dimension ainclbead(ntmax*numax*numax,2)
      dimension a15t(ntmax*numax*numax)

c ----------------------------------------------------------------

c - triple loop over all types of molecules -
      do imolty = 1, nmolty

         do m = 1, nunit(imolty)
            do n = 1, nunit(imolty)
               linclu(imolty,m,n) = .true.
               lqinclu(imolty,m,n) = .true.
c * by default, don't want any 1-5 r^12 interactions
               lainclu(imolty,m,n) = .false.
            enddo
         enddo
         
c - double loop over all units -
         do m = 1, nunit(imolty)

c - exclude all self interactions -
            linclu(imolty,m,m) = .false.
            lqinclu(imolty,m,m) = .false.
c - exclude all directly bonded beads (vibrations) -
            do n = 1, invib(imolty,m)
               nb = ijvib(imolty,m,n)
               linclu(imolty,m,nb) = .false.
               lqinclu(imolty,m,nb) = .false.
            enddo

c - exclude carbons around a quaternary center for explct
            if (invib(imolty,m) .eq. 4) then
               do n = 1,4
                  do nb = 1,4
                     linclu(imolty,ijvib(imolty,m,n)
     &                    ,ijvib(imolty,m,nb))=.false.
                     lqinclu(imolty,ijvib(imolty,m,n)
     &                    ,ijvib(imolty,m,nb))=.false.
                  enddo
               enddo
            endif

c - exclude all next-nearest neighbor bonded beads (bending) -
            do n = 1, inben(imolty,m)
               nb = ijben3(imolty,m,n)
               linclu(imolty,m,nb) = .false.
               lqinclu(imolty,m,nb) = .false.
            enddo
c - exclude all third-nearest neighbor bonded beads (torsions) -
            do n = 1, intor(imolty,m)
               nb = ijtor4(imolty,m,n)
               linclu(imolty,m,nb) = .false.
c * don't set lqinclu since we want 1-4 interactions, unless q14scale is zero
               if (qscale .lt. 0.001d0) then
                  lqinclu(imolty,m,nb) = .false.
               endif
            enddo

         enddo

c     - include or exclude additional beads accoring to incl
         do n = 1,inclnum
            if ( inclmol(n) .eq. imolty ) then
               m = inclbead(n,1)
               nb = inclbead(n,2)
               if ( inclsign(n) .eq. 1 ) then
                  linclu(imolty,m,nb) = .true.
                  linclu(imolty,nb,m) = .true.
                  lqinclu(imolty,m,nb) = .true.
                  lqinclu(imolty,nb,m) = .true.
               elseif (inclsign(n) .eq. -1 ) then
                  linclu(imolty,m,nb) = .false.
                  linclu(imolty,nb,m) = .false.
                  lqinclu(imolty,m,nb) = .false.
                  lqinclu(imolty,nb,m) = .false.
               else
                  write(6,*) 'INCLUS: n,inclsign(n)',n,inclsign(n)
                  stop 'inclusign must be 1 or -1'
               endif
            endif

         enddo

c * add in 1-5 interactions according to aincl
         do m = 1,ainclnum
            if ( ainclmol(m) .eq. imolty ) then
               mb = ainclbead(m,1)
               nb = ainclbead(m,2)
                  lainclu(imolty,mb,nb) = .true.
                  lainclu(imolty,nb,mb) = .true.
                  a15type(imolty,mb,nb) = a15t(m)
                  a15type(imolty,nb,mb) = a15t(m)
            endif
         enddo

c - exclude all hydrogens that have their carbons excluded 
         if ( ncarbon(imolty) .lt. nunit(imolty) ) then
            if (ncarbon(imolty) .eq. 3 .and. nunit(imolty) .eq. 8) then
c - ethane with bead 3 being hydrogen
               ioffset = 0
            else
               ioffset = 1
            endif
            do m = ncarbon(imolty)+ioffset,nunit(imolty)
c - hydrogens only have one vibration and that is to the C atom               
               mb = ijvib(imolty,m,1)
               do nb = 1, ncarbon(imolty)
                  if ( .not. linclu(imolty,mb,nb) ) then
                     linclu(imolty,m,nb) = .false.
                     linclu(imolty,nb,m) = .false.
                     lqinclu(imolty,m,nb) = .false.
                     lqinclu(imolty,nb,m) = .false.
                  endif
               enddo
               do n=m+1,nunit(imolty)
                  nb = ijvib(imolty,n,1)
                  linclu(imolty,n,nb) = .false.
                  linclu(imolty,nb,n) = .false.
                  lqinclu(imolty,n,nb) = .false.
                  lqinclu(imolty,nb,n) = .false.
                  if ( .not. linclu(imolty,mb,nb) ) then
                     linclu(imolty,m,n) = .false.
                     linclu(imolty,n,m) = .false.
                     linclu(imolty,m,nb) = .false.
                     linclu(imolty,nb,m) = .false.
                     linclu(imolty,n,mb) = .false.
                     linclu(imolty,mb,n) = .false.
                     lqinclu(imolty,m,n) = .false.
                     lqinclu(imolty,n,m) = .false.
                     lqinclu(imolty,m,nb) = .false.
                     lqinclu(imolty,nb,m) = .false.
                     lqinclu(imolty,n,mb) = .false.
                     lqinclu(imolty,mb,n) = .false.
                  endif
               enddo
            enddo
         endif

c * exclude charge interactions if lqchg is false
         do m = 1,nunit(imolty)
            do n = m+1,nunit(imolty)
               if (.not. lqchg(ntype(imolty,m)) .or.
     &              .not. lqchg(ntype(imolty,n))) then
                  lqinclu(imolty,m,n) = .false.
                  lqinclu(imolty,n,m) = .false.
               endif
            enddo
         enddo

c * self consistency check *

         do m = 1,nunit(imolty)
            do n = m+1,nunit(imolty)

               if (linclu(imolty,m,n) .neqv. linclu(imolty,n,m)) then
                  linclu(imolty,m,n) = .false.
                  linclu(imolty,n,m) = .false.
               endif

               if (lqinclu(imolty,m,n) .neqv. lqinclu(imolty,n,m)) then
                  lqinclu(imolty,m,n) = .false.
                  lqinclu(imolty,n,m) = .false.
               endif

            enddo
         enddo

         if (lrigid(imolty)) then
c - don't include rigid beads

c - there will be no intramolecular forces between rigid beads
c - or beads connected one away from a rigid bead          
            do m = 1, invib(imolty,riutry(imolty,1))
               mb = ijvib(imolty,riutry(imolty,1),m)
               do n = riutry(imolty,1), nunit(imolty)
                  linclu(imolty,n,mb) = .false.
                  linclu(imolty,mb,n) = .false.
                  lqinclu(imolty,n,mb) = .false.
                  lqinclu(imolty,mb,n) = .false.
               enddo
            enddo

 
            do m = riutry(imolty,1), nunit(imolty)
               do n = riutry(imolty,1), nunit(imolty)
                  linclu(imolty,m,n) = .false.
                  lqinclu(imolty,m,n) = .false.
               enddo
            enddo
         endif

         write(6,*) 
         write(6,*) 'INCLUSION TABLE'

         do m = 1, nunit(imolty)
            write(6,*) m, (linclu(imolty,m,n),n=1,nunit(imolty))
         enddo
         write(6,*) 

         write(6,*) 
         write(6,*) 'CHARGE INCLUSION TABLE'

         do m = 1, nunit(imolty)
            write(6,*) m, (lqinclu(imolty,m,n),n=1,nunit(imolty))
         enddo
         write(6,*) 

c * not really that important to write out
c         write(6,*) 
c         write(6,*) '1-5 OH INTERACTION TABLE'
c         do m = 1, nunit(imolty)
c            write(6,*) m, (lainclu(imolty,m,n),n=1,nunit(imolty))
c         enddo
c         write(6,*) 

      enddo
c      write(6,*) 'finished inclus'
      return
      end






