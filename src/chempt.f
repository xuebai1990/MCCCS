      subroutine chempt(bsswap,imolty,ntries,qelect)

c chempt
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
 
c    ********************************************************************
c    ** Ghost insertion of molecule into the boxes to compute chem pot **
c    ** using CBMC insertion techniques.  Works for linear or branched **
c    ** molecules and for anisotropic and Explicit atom                **
c    ** Written M.G. Martin  9-24-97                                   **
c    ********************************************************************
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'system.inc'
      include 'cbmc.inc'
      include 'rosen.inc' 
      include 'boltzmann.inc'
      include 'external.inc'
      include 'connect.inc'
      include 'inputdata.inc'

      logical lctrl
      parameter (lctrl=.false.)

      logical ovrlap,lterm,lbo,ltors

      integer boxins,ichoi,ip,iwalk,ntries,itry,iincre
      integer istt,iett
      integer ibranp(20),nbranp

      integer iu,iutry,iulast,iut,iub,icbu,islen,iins
     &       ,iii,ibr,ibr1,j,jj,jj2,jj3,jj4,invtry,ibox,iunit
     &       ,imolty,jmt,igrow
      double precision v,vintra,vinter,vext,velect,vtornew,vtordum
     &     ,delen,vewald


      double precision bsswap
      double precision random,rdir,qelect,rbf,bsum
      double precision waddnew

      double precision v1insext,v1ins,w1ins
     &                ,v1insint
     &                ,volins,rho,arg,coru,v1inselc
      dimension bsswap(ntmax,4)
      dimension qelect(nntype)

C --------------------------------------------------------------------

c      write(2,*) 'start CHEMP'

c *** store number of units in iunit and igrow ***
      iunit = nunit(imolty)
      igrow = nugrow(imolty)
c *** give i a phony number ***
      iins = nchain + 1
      moltyp(iins) = imolty

c *** give i a phony number ***
      iins = nchain + 1
      moltyp(iins) = imolty
c     give charges to phony number 
      do icbu = 1, iunit
         qqu(iins,icbu) = qelect(ntype(imolty,icbu))
      enddo

      do 500 itry = 1,ntries
c ---    select a box at random
      
         if (random().lt.0.5d00) then
            boxins=1
         else
            boxins=2
         endif

c      X = 2 is # attempts into box 2 X=3 is success into box 1
c      X = 4 is success into box 2
c     bsswap is the same thing but keeps track of successful growths
         
         bsswap(imolty,boxins) = bsswap(imolty,boxins) + 1.0d0

c *** select a position of the first/starting unit at RANDOM ***
c *** and calculate the boltzmann weight                     ***
c *** for the chain to be INSERTED                           ***
         ichoi = nchoi1(imolty)
         do icbu = 1,ichoi
            rxp(1,icbu) = boxlx(boxins) * random()
            ryp(1,icbu) = boxly(boxins) * random()
            if (lpbcz) then
               rzp(1,icbu) = boxlz(boxins) * random()
            else if ( lsami .or. lmuir .or. ljoe ) then
               rzp(1,icbu) = 20*random()-10
            else
               rzp(1,icbu) = 0.0d0
            endif
         enddo
         
c *** select starting unit ***
c --- always using bead 1 as the starting unit
         iutry = 1
         lbo = .true.

c --  insert the first atom
         call boltz( lbo,iins,iins,imolty,igrow,ovrlap,boxins,.true.
     &        ,0.0d0,0,ichoi)
         bnchem(boxins,imolty) = bnchem(boxins,imolty) + 1.0d0
         if ( ovrlap ) goto 500

c *** perform the walk according to the availibility of the choices ***
c *** and calculate the correct weight for the trial walk           ***

         w1ins = 0.0d0
         do ip = 1, ichoi
            w1ins = w1ins + bfac(ip)
         enddo

c --- check for termination of walk ---
         if ( w1ins .lt. softlog ) goto 500

c --- select one position at random ---
         if ( ichoi .gt. 1 ) then
            rbf = w1ins * random()
            bsum = 0.0d0 
            do ip = 1, ichoi
               if ( .not. lovr(ip) ) then
                  bsum = bsum + bfac(ip)
                  if ( rbf .lt. bsum ) then
c                    --- select ip position ---
                     iwalk = ip
                     goto 180
                  endif
               endif
            enddo
         else
            iwalk = 1
         endif

 180     v1ins =  vtry(iwalk)  
         v1insext = vtrext(iwalk)
         v1insint = vtrinter(iwalk)
         v1inselc = vtrelect(iwalk)
         
         rxnew(1) = rxp(1,iwalk)
         rynew(1) = ryp(1,iwalk)
         rznew(1) = rzp(1,iwalk)

c         if ( lelect(imolty) ) then
c            call qqcheck(iins,boxins,rxnew(1),rynew(1),rznew(1))
c         endif

c *** set walk conditions ***
         invtry = invib(imolty,iutry)
         if ( invtry .eq. 0 ) then
c --- Bead 1 is the only bead to be grown ---
            islen = 0
            goto 100
         elseif ( invtry .eq. 1 ) then
c --- Bead 1 is an endpoint ---
            nbranp = 0
            iincre = 1

            iut = ijvib(imolty,iutry,1)
            lexist(iut) = .false.
            wsched(1) = iut
            islen = 1
         else
c --- Bead 1 is a midpoint or branchpoint ---
            nbranp = 0
c --- all branches will be regrown - starting at the lowest/highest ---
            rdir = random()
            if ( rdir .le. 0.5d0 ) then
               iincre = -1
            else
               iincre = 1
            endif
            if ( iincre .eq. -1 ) then
               write(2,*) 'imolty,iutry,invtry',imolty,iutry,invtry
               iut = ijvib(imolty,iutry,invtry)
               do iii = invtry-1, 1, -1
                  nbranp = nbranp + 1
                  ibranp(nbranp) = ijvib(imolty,iutry,iii)
               enddo
            else
               iut = ijvib(imolty,iutry,1)
               do iii = 2, invtry
                  nbranp = nbranp + 1
                  ibranp(nbranp) = ijvib(imolty,iutry,iii)
               enddo
            endif
            lexist(iut) = .false.
            wsched(1) = iut
            islen = 1
         endif

c *** calculate number of trial segments ***
c *** set LEXIST of trial segments to false ***
c *** generate growing schedule ***
         if ( iincre .eq. 1 ) then
c --- find trial segment with lowest index ---
 21         iulast = wsched(islen)
            iut = 0
            iutry = 0
            do iu = 1, invib(imolty,iulast)
               iub = ijvib(imolty,iulast,iu)
               if ( lexist(iub) ) then
                  iut = iub
                  go to 22
               endif
            enddo
 22         if ( iut .ne. 0 ) then
               if ( nbranp .gt. 0 ) then
c --- compare to index of lowest branchpoint ---
                  if ( iut .gt. ibranp(1) ) then
c --- select branchpoint and remove it from list of branchpoints ---
                     iutry = ibranp(1)
                     lexist(iutry) = .false.
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = 0
                     else
                        do ibr = 2, nbranp
                           ibr1 = ibr - 1
                           ibranp(ibr1) = ibranp(ibr)
                        enddo
                        ibranp(nbranp) = 0
                     endif
                     nbranp = nbranp - 1
                  else
c --- select IUT and leave branches ---
                     iutry = iut
                     lexist(iutry) = .false.
                  endif
               else
c --- select IUT ---
                  iutry = iut
                  lexist(iutry) = .false.
               endif
c --- add other connections (including non-selected IUT) to branchpoints ---
 23            do iu = 1, invib(imolty,iulast)
                  iub = ijvib(imolty,iulast,iu)
                  if ( lexist(iub) ) then
                     nbranp = nbranp + 1
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = iub
                     else
                        do ibr = nbranp-1, 1, -1
                           ibr1 = ibr + 1
                           if ( iub .lt. ibranp(ibr) ) then
                              ibranp(ibr1) = ibranp(ibr)
                              if ( ibr .eq. 1 ) ibranp(ibr) = iub
                           else
                              ibranp(ibr1) = iub
                              go to 24
                           endif
                        enddo
                     endif
                  endif
 24               continue
               enddo
            else
               if ( nbranp .gt. 0 ) then
c --- take lowest branchpoint ---
                  iutry = ibranp(1)
                  lexist(iutry) = .false.
                  do ibr = 2, nbranp
                     ibr1 = ibr - 1
                     ibranp(ibr1) = ibranp(ibr)
                  enddo
                  nbranp = nbranp - 1
               else
c --- iulast is the last unit in the walkschedule ---
                  go to 100
               endif
            endif

            islen = islen + 1
            wsched(islen) = iutry
            go to 21
         else

c --- find trial segment with highest index ---
 31         iulast = wsched(islen)
            write(2,*) '*************'
            write(2,*) 'iulast',iulast
            write(2,*) 'ibranp', (ibranp(iu),iu=1,nbranp)
            iut = 0
            iutry = 0
            do iu = invib(imolty,iulast), 1, -1
               iub = ijvib(imolty,iulast,iu)
               if ( lexist(iub) ) then
                  iut = iub
                  go to 32
               endif
            enddo
 32         if ( iut .ne. 0 ) then
               if ( nbranp .gt. 0 ) then
c --- compare to index of highest branchpoint ---
                  if ( iut .lt. ibranp(1) ) then
c --- select branchpoint and remove it from list of branchpoints ---
                     iutry = ibranp(1)
                     lexist(iutry) = .false.
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = 0
                     else
                        do ibr = 2, nbranp
                           ibr1 = ibr - 1
                           ibranp(ibr1) = ibranp(ibr)
                        enddo
                        ibranp(nbranp) = 0
                     endif
                     nbranp = nbranp - 1
                  else
c --- select IUT and leave branches ---
                     iutry = iut
                     lexist(iutry) = .false.
                  endif
               else
c --- select IUT ---
                  iutry = iut
                  lexist(iutry) = .false.
               endif
c --- add other connections (including non-selected IUT) to branchpoints ---
 33            do iu = invib(imolty,iulast), 1, -1
                  iub = ijvib(imolty,iulast,iu)
                  if ( lexist(iub) ) then
                     nbranp = nbranp + 1
                     if ( nbranp .eq. 1 ) then
                        ibranp(1) = iub
                     else
                        do ibr = nbranp-1, 1, -1
                           ibr1 = ibr + 1
                           if ( iub .gt. ibranp(ibr) ) then
                              ibranp(ibr1) = ibranp(ibr)
                              if ( ibr .eq. 1 ) ibranp(ibr) = iub
                           else
                              ibranp(ibr1) = iub
                              go to 34
                           endif
                        enddo
                     endif
                  endif
 34               continue
               enddo
            else
               if ( nbranp .gt. 0 ) then
c --- take highest branchpoint ---
                  iutry = ibranp(1)
                  lexist(iutry) = .false.
                  do ibr = 2, nbranp
                     ibr1 = ibr - 1
                     ibranp(ibr1) = ibranp(ibr)
                  enddo
                  nbranp = nbranp - 1
               else
c --- iulast is the last unit in the walkschedule ---
                  go to 100
               endif
            endif

            islen = islen + 1
            wsched(islen) = iutry
            go to 31

         endif

 100     continue

c - create walk schedule for bonded interactions -
         do iii = 1, igrow
            lexist(iii) = .true.
         enddo
         do iii = 1, islen
            iut = wsched(iii)
            lexist(iut) = .false.
         enddo
         do iii = 1, islen
            iut = wsched(iii)
            do j = 1, invib(imolty,iut)
               jj = ijvib(imolty,iut,j)
               if ( lexist(jj) ) then
                  wschvib(iut,j) = .true.
               else
                  wschvib(iut,j) = .false.
               endif
            enddo
            do j = 1, inben(imolty,iut)
               jj2 = ijben2(imolty,iut,j)
               jj3 = ijben3(imolty,iut,j)
               if ( lexist(jj2) .and. lexist(jj3) ) then
                  wschben(iut,j) = .true.
               else
                  wschben(iut,j) = .false.
               endif
            enddo
            do j = 1, intor(imolty,iut)
               jj2 = ijtor2(imolty,iut,j)
               jj3 = ijtor3(imolty,iut,j)
               jj4 = ijtor4(imolty,iut,j)
               if ( lexist(jj2) .and. lexist(jj3) 
     &              .and. lexist(jj4) ) then
                  wschtor(iut,j) = .true.
               else
                  wschtor(iut,j) = .false.
               endif
            enddo
            lexist(iut) = .true.
         enddo

c ------------------------------------------------------------------
c     --- not currently working 2-11-98
         lterm = .false.
c         call rosnbr ( iins, imolty, islen, boxins, lterm, igrow  )
c --- termination of cbmc attempt due to walk termination ---
         if ( lterm ) goto 500


c --- Begin DC-CBMC Corrections for NEW configuration
         waddnew = 1.0d0
         if (ldual .or. ((.not. lchgall) .and. lelect(imolty))) then 
c     compute the energy of the inserted molecule
c           calculate the full site-site energy
c           iii=1 new conformation
            do j=1,igrow
               rxuion(j,1) = rxnew(j)
               ryuion(j,1) = rynew(j)
               rzuion(j,1) = rznew(j)
               qquion(j,1) = qqu(iins,j)
            enddo

            ibox=boxins
            nboxi(iins) = ibox

            istt=1
            iett = igrow
            call energy (iins,imolty, v, vintra,vinter,vext,velect
     &           ,vewald,1, ibox, istt, iett, .false.,ovrlap,.false.
     &           ,vtordum,.false.,.false.)
            
            if (ovrlap) goto 500
            delen = v - ( vnewintra + vnewinter + vnewext +vnewelect)
     &           - v1ins

            waddnew = waddnew*dexp(-(beta*delen))

            vnewt     = vnewt + delen
            vnewinter = vinter - v1insint
            vnewext   = vext - v1insext
            vnewelect = velect - v1inselc
         endif
c     End DC-CBMC Corrections for NEW configuration

c     Begin Explicit Atom Corrections for NEW configuration
         if ( iunit .ne. igrow ) then
c        calculate the true Lennard-Jones energy for the hydrogens
c        iii=1 new conformation
            do j=1,iunit
               rxu(iins,j)=rxnew(j)
               ryu(iins,j)=rynew(j)
               rzu(iins,j)=rznew(j)
            enddo
            ibox = boxins
            call explct(iins,vtornew,.false.,.false.)
            ltors = .false.
            
            do j=1,iunit
               rxuion(j,1) = rxu(iins,j)
               ryuion(j,1) = ryu(iins,j)
               rzuion(j,1) = rzu(iins,j)
               qquion(j,1) = qqu(iins,j)
            enddo
c Calculate the energy of the non-backbone beads 
            istt = igrow+1
            iett = iunit
            call energy (iins,imolty,v, vintra,vinter,vext,velect
     &           ,vewald,1
     &           ,ibox,istt,iett, .true.,ovrlap,ltors,vtordum,
     +           .true.,.false.)
            if (ovrlap) goto 500
            delen = v - vnewintra + vtornew

            waddnew = waddnew*dexp(-beta*delen)

            vnewt     = vnewt + delen
            vnewintra = vintra
            vnewinter = vnewinter + vinter 
            vnewext   = vnewext + vext
            vnewtg    = vnewtg + vtornew
            vnewelect = vnewelect + velect
         endif

c     End Explicit Atom Corrections for NEW configuration

c     add in the contribution from the first bead
         vnewt     = vnewt + v1ins
         vnewinter = vnewinter + v1insint
         vnewext   = vnewext + v1insext
         vnewelect = vnewelect + v1inselc

         if (lpbcz) then
            volins=boxlx(boxins)*boxly(boxins)*boxlz(boxins)
         else
            volins=boxlx(boxins)*boxly(boxins)
         endif

         arg = w1ins * waddnew * weight * volins 
     &        / dble( ncmt(boxins,imolty)+1 )

         if (ltailc) then
            do jmt = 1, nmolty
               if ( jmt .eq. imolty ) then
                  rho = dble(ncmt(boxins,jmt)+1)/volins
               else
                  rho = dble(ncmt(boxins,jmt))/volins
               endif
               arg=arg*dexp(-(beta*2.0d0*coru(imolty,jmt,rho,boxins)))
            enddo
         endif

         acchem(boxins,imolty) = acchem(boxins,imolty)+arg
         bsswap(imolty,boxins+2) = bsswap(imolty,boxins+2) + 1.0d0

 500  continue


c       write(2,*) 'end CHEMPT'

      return
      end

