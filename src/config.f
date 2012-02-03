      subroutine config

c config
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
c    ** performs a lengthconserving configurational bias move         **
c    ** for linear, branched, anisotropic, and explicit atom          **
c    ** molecules                                                     **
c    ** rewritten from old config and branch subroutines by           **
c    ** M.G. Martin 9-19-97                                           **
c    ** number of trial attempts starting at unit inb is stored in    **
c    **    bncb ( inb ).                                              **
c    ** number of successful generations of trial configuration is in **
c    **    bscb ( 1,inb ).                                            **
c    ** number of accepted trial configurations is in                 **
c    **    bscb ( 2,inb ).                                            **
c    *******************************************************************
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'cbmc.inc'
      include 'rosen.inc' 
      include 'inputdata.inc'
      include 'ewaldsum.inc'
      include 'poten.inc'
      include 'neigh.inc'
      
      logical lterm, ovrlap, ltors, lneighij,lfixnow

      integer i,j,k,iii,ibox,iunit,igrow,icbu,islen,imolty,iutry

      integer istt,iett,nchp1,ic,ncount,total,bin,count,findex,iw
      integer ddum,idum,ip

      dimension ddum(27)

      double precision v,vintra,vinter,vext,velect,vewald,vtorold
     & ,vtornew,delen,deleo,vdum,tofo,wplace,wrig,vorient

      double precision dchain,random,rchain,wnlog,wolog,wdlog,wratio

      double precision vrecipn,vrecipo,cwtorfo,cwtorfn,x,y,z

c ------------------------------------------------------------------

c      write(6,*) 'start CONFIG'
c ***    select a chain at random ***
      rchain  = random()
      do icbu = 1,nmolty
         if ( rchain .lt. pmcbmt(icbu) ) then
            imolty = icbu
            rchain = 2.0d0
         endif
      enddo

c     *** determine whether to use fecbmc or not ***
      if (random().lt.pmfix(imolty)) then
         lfixnow = .true.
      else
         lfixnow = .false.
      endif

      if (lgrand) then
c ---    select a chain in box 1
c         write(6,*) 'counters not implemented properly for grand'
         if (ncmt(1,imolty).eq.0) return
         i = idint( dble(ncmt(1,imolty))*random() ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(6,*) 'screwup config'
         ibox=1
      else 
         dchain = dble(temtyp(imolty))
         i = int( dchain*random() + 1 )
         i = parall(imolty,i)
         ibox = nboxi(i)
         if ( moltyp(i) .ne. imolty ) stop 'screwup config'
      endif

c *** store number of units in iunit and # to be grown in igrow ***
      iunit = nunit(imolty)
      igrow = nugrow(imolty)
 
c *** store position of trial chain in r x/y/z cbu ***
      do icbu = 1, igrow
         rxnew(icbu) = rxu(i,icbu)
         rynew(icbu) = ryu(i,icbu)
         rznew(icbu) = rzu(i,icbu)
      enddo

      if (lfixnow) then
         call safeschedule(igrow,imolty,islen,iutry,findex,1)
      else
         call schedule(igrow,imolty,islen,iutry,0,1)
      endif

c     --- determine how many beads are being regrown
      total = 0
      do icbu = 1,igrow
         if ( .not. lexshed(icbu) ) total = total + 1
      enddo

      if (lfixnow) then
         fbncb(imolty,findex-1) = fbncb(imolty,findex-1) + 1.0d0
      else
         bncb(imolty,total) = bncb(imolty,total) + 1.0d0
      endif
 
c      if ( lelect(imolty) ) then
c        ---  Call qqcheck to setup the group based qq cutoff
c         call qqcheck(i,ibox,rxnew(1),rynew(1),rznew(1))
c      endif

c     --- grow new chain conformation
      call rosenbluth ( .true., lterm,i,i,imolty,islen,ibox,igrow
     &     ,vdum,lfixnow,cwtorfn,1 )
 
c --- termination of cbmc attempt due to walk termination ---
      if ( lterm ) then 
c        write(6,*) 'termination of growth',i
        return
      endif

      if (llrig) then
         call rigfix(.true.,i,ibox,imolty,lterm,wrig)
         if ( lterm ) return
         weight = weight * wrig
      endif

      if (llplace(imolty).and.lfixnow) then
         call place(.true.,lterm,i,imolty,ibox,islen,wplace)
         if ( lterm ) return
         weight = weight * wplace 
      endif


c     --- grow old chain conformation
      call rosenbluth ( .false.,lterm,i,i,imolty,islen,ibox,igrow
     &     ,vdum,lfixnow,cwtorfo,1)

c     --- termination of old walk due to problems generating orientations
      if ( lterm ) then
         write(6,*) 'CONFIG: old growth rejected'
         return
      endif

      if (llrig) then
         call rigfix(.false.,i,ibox,imolty,lterm,wrig)
         if ( lterm ) then
            write(6,*) 'CONFIG: old rigid fix rejected'
            return
         endif
         weiold = weiold * wrig
      endif

      if (llplace(imolty).and.lfixnow) then
         call place(.false.,lterm,i,imolty,ibox,islen,wplace)
         
         if ( lterm ) then
            write(6,*) 'CONFIG: old hydrogen placement rejected'
            return
         endif
         weiold = weiold * wplace 
      endif

C -----------------------------------------------------------------------------
c     Begin DC-CBMC, Explicit Atom and Ewald-sum Corrections 

      if ( ldual .or. lewald .or. iunit .ne. igrow 
     &     .or. ((.not. lchgall) .and. lelect(imolty)) ) then
c     --- Put on hydrogens for explicit AA model for calculation of COM
c     --- and assign all of the grown new and old beads to rxuion
c     --- with old = 1, new = 2
         do j=1,igrow
            rxuion(j,1)=rxu(i,j)
            ryuion(j,1)=ryu(i,j)
            rzuion(j,1)=rzu(i,j)
            qquion(j,1)=qqu(i,j)
         enddo	 
         do j = 1,igrow
            rxuion(j,2) = rxnew(j)
            ryuion(j,2) = rynew(j)
            rzuion(j,2) = rznew(j)
            qquion(j,2) = qquion(j,1)
         enddo
         nchp1=nchain+1
         nboxi(nchp1) = ibox
         moltyp(nchp1) = imolty
         moltion(1) = imolty
         moltion(2) = imolty

         if ( igrow .ne. iunit ) then
c           -- iii = 1 old conformation
            do j = igrow+1, iunit
               rxuion(j,1) = rxu(i,j)
               ryuion(j,1) = ryu(i,j)
               rzuion(j,1) = rzu(i,j)
               qquion(j,1) = qqu(i,j)
            enddo
c           -- iii = 2 new conformation
            do j=1, igrow
               rxu(nchp1,j) = rxnew(j)
               ryu(nchp1,j) = rynew(j)
               rzu(nchp1,j) = rznew(j)
            enddo
            call explct(nchp1,vtornew,.false.,.false.)
            do j=igrow+1, iunit
               rxuion(j,2) = rxu(nchp1,j)
               ryuion(j,2) = ryu(nchp1,j)
               rzuion(j,2) = rzu(nchp1,j)
               qquion(j,2) = qquion(j,1)
            enddo
         endif
      endif

      if (ldual .or. ((.not. lchgall) .and. lelect(imolty))
     &     .or. (lchgall .and. lewald .and. (.not. ldual))) then
         istt = 1
         iett = igrow

c        -- check new before old
         do iii = 2,1,-1
c          calculate the Full rcut Lennard-Jones energy for the grown beads
c          iii = 1 old conformation
c          iii = 2 new conformation
            
            call energy (i,imolty,v,vintra,vinter,vext,velect
     &           ,vewald,iii,ibox,istt,iett,.true.,ovrlap,
     &           .false.,vdum,.false.,.false.)

            if (ovrlap .and. (iii .eq. 1)) then
c            if (ovrlap) then
               write(6,*) 'disaster: overlap in old conf config',i
               stop
            endif

            if (iii .eq. 2) then
               delen = v - ( vnewinter + vnewext + vnewelect +
     &              vnewewald + vnewintra) 
               weight    = weight*dexp(-(beta*delen))
               vnewt     = vnewt + delen
               vnewinter = vinter
               vnewext   = vext
               vnewelect = velect
               vnewintra = vintra
               vnewewald = vewald
            else
               deleo = v - ( voldinter + voldext + voldelect +
     &              voldewald + voldintra) 
               weiold    = weiold*dexp(-(beta*deleo))
               voldt     = voldt + deleo
               voldinter = vinter
               voldext   = vext
               voldelect = velect
               voldintra = vintra
               voldewald = vewald
            endif
         enddo
         
      endif

      if ( iunit .ne. igrow ) then
         istt = igrow+1
         iett = iunit

c        -- check new before old
         do iii = 2,1,-1
c           calculate the true Lennard-Jones energy for the hydrogens
c           iii=1 new conformation
c           iii=2 old conformation
c           hydrogens were placed and rxuion was assigned above

            if (iii .eq. 1) then
               ltors = .true.
            else
               ltors = .false.
            endif

c Calculate the energy of the non-backbone beads 
            call energy (i,imolty,v,vintra,vinter,vext,velect
     &           ,vewald,iii,ibox,istt,iett, .true.,ovrlap
     &           ,ltors,vtorold,.true.,.false.)

            if (iii .eq. 2) then
               if (ovrlap) return
               delen = v + vtornew
               if ( delen*beta .gt. (2.3d0*softcut) ) then
c                  write(6,*) '##softcut in config caught explicit atoms'
                  return
               endif
               weight = weight*dexp(-(beta*delen))
               vnewt  = vnewt + delen
               vnewintra = vnewintra + vintra 
               vnewinter = vnewinter + vinter 
               vnewext   = vnewext + vext 
               vnewtg = vnewtg + vtornew
               vnewelect = vnewelect + velect
               vnewewald = vnewewald + vewald
            else
               if (ovrlap) then
                  write(6,*) 'ovrlap problem in old confomation -CONFIG'
                  return
               endif
               deleo = v + vtorold
               weiold = weiold*dexp(-(beta*deleo))
               if ( weiold .lt. softlog ) then
                  write(6,*) '##old weight for explicit too low'
               endif
               voldt     = voldt + deleo
               voldintra = voldintra + vintra
               voldinter = voldinter + vinter
               voldext   = voldext + vext
               voldtg    = voldtg + vtorold
               voldelect = voldelect + velect
               voldewald = voldewald + vewald
            endif
         enddo
      endif

      if ( lewald .and. lelect(imolty) ) then
c        --- reciprocal space sum ---
c        --- rxuion: 1= old configuration; 2= new configuration
         call recip(ibox,vrecipn,vrecipo,1)
         delen = vrecipn
         deleo = vrecipo
         weight = weight * dexp(-(beta*vrecipn))
         weiold = weiold * dexp(-(beta*vrecipo))
         vnewelect = vnewelect + vrecipn
         voldelect = voldelect + vrecipo
         vnewt = vnewt + vrecipn
         voldt = voldt + vrecipo
      endif

c     End of DC-CBMC, Explicit Atom and Ewald-sum Corrections

c *** check for acceptance of trial configuration ***
      wnlog = dlog10 ( weight )
      wolog = dlog10 ( weiold )
c      write(6,*) 'weight:',weight
c      write(6,*) 'weiold:',weiold
      wdlog = wnlog - wolog
      if ( wdlog .lt. -softcut ) then
c         write(6,*) 'cbmc softcut',i
         return
      endif
 
      if (lfixnow) then
         wratio = weight * cwtorfo / ( weiold * cwtorfn)
         fbscb(imolty,1,findex-1) = 
     &        fbscb(imolty,1,findex-1) + 1.0d0
      else
         wratio = weight / weiold
         bscb(imolty,1,total) = bscb(imolty,1,total) + 1.0d0
      endif

      if ( random() .le. wratio ) then
c         write(6,*) 'CONFIG accepted',i,ibox
c        --- we can now accept !!!!! ***
         if (lfixnow) then
            fbscb(imolty,2,findex-1) = fbscb(imolty,2,findex-1) 
     &           + 1.0d0
         else
            bscb(imolty,2,total) = bscb(imolty,2,total) + 1.0d0
         endif


         vbox(ibox)    = vbox(ibox)    + ( vnewt - voldt )
         vinterb(ibox) = vinterb(ibox) + (vnewinter - voldinter)
         vintrab(ibox) = vintrab(ibox) + (vnewintra- voldintra)
         vvibb(ibox)   =  vvibb(ibox)  + (vnewbvib- voldbvib)
         vtgb(ibox)    = vtgb(ibox)    + (vnewtg- voldtg)
         vextb(ibox)   = vextb(ibox)   + (vnewext - voldext)
         vbendb(ibox)  = vbendb(ibox)  + (vnewbb - voldbb)
         velectb(ibox) = velectb(ibox) + (vnewelect - voldelect)
     &        + (vnewewald - voldewald)

         do ic = 1, igrow
            rxu(i,ic) = rxnew(ic)
            ryu(i,ic) = rynew(ic)
            rzu(i,ic) = rznew(ic)
         enddo
         do ic = igrow+1, iunit
            rxu(i,ic)  = rxuion(ic,2)
            ryu(i,ic)  = ryuion(ic,2)
            rzu(i,ic)  = rzuion(ic,2)
         enddo

         if ( lewald .and. lelect(imolty) ) then
c *** update reciprocal-space sum
            call recip(ibox,vdum,vdum,2)
         endif

c ***    update center of mass
         call ctrmas(.false.,ibox,i,7)
c *** update linkcell, if applicable
         if ( licell .and. (ibox.eq.boxlink)) then
            call linkcell(2,i,vdum,vdum,vdum,ddum)
         endif

c ***    update the neighbour map ***
         if ( lneigh ) call updnn( i )

         if ( lneighbor ) then
            do ic = 1, neigh_cnt(i)
               j = neighbor(ic,i)
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. i ) then
                     neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                     neigh_cnt(j) = neigh_cnt(j)-1
                     goto 10
                  endif
               enddo
            enddo
 10         neigh_cnt(i) = neigh_icnt
            do ic = 1,neigh_icnt
               j = neighi(ic)
               neighbor(ic,i)=j
               lneighij = .false.
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. i ) then
                     lneighij = .true.
                  endif
               enddo
               if ( .not. lneighij ) then
                  neigh_cnt(j) = neigh_cnt(j)+1
                  neighbor(neigh_cnt(j),j) = i
               endif
            enddo
         endif

       endif

       if (lpresim.or.lfixnow) then
c     --- record bond distances for presimulation and reweighting
          counthist = counthist + 1
          do iw = 1, islen
             do count = 1, grownum(iw)
                k = growlist(iw,count)
                do j = 1, nunit(imolty)
                   if (k.eq.j) goto 100
                   x = rxu(i,j) - rxu(i,k)
                   y = ryu(i,j) - ryu(i,k)
                   z = rzu(i,j) - rzu(i,k)
                   
                   bin = anint(10.0d0*dsqrt(x**2+y**2+z**2))
                   
                   if (bin.gt.maxbin) goto 100
                   
                   hist(j,k,bin) = hist(j,k,bin) + 1
                   hist(k,j,bin) = hist(k,j,bin) + 1
 100               continue
                enddo
             enddo
          enddo
       endif
         
c -----------------------------------------------------------------
c       write(6,*) 'end CONFIG'
       return
       end

