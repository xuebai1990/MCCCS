      subroutine safeschedule(igrow,imolty,islen,iutry,findex,movetype)

c safeschedule
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
 
c     *********************************************************************
c     **  Self-Adapting Fixed-Endpoint Configurational-Bias Monte Carlo  **
c     **                    --- SAFE-CBMC ---                            **
c     *********************************************************************
c     **    Determines logic for a CBMC Move Between Fixed End Points    **
c     **           for Linear, Branched, and Cyclic Molecules            **
c     *********************************************************************
c     **       Originally completed by Collin Wick around 1-1-2000       **
c     *********************************************************************


c     *********************************************************************
c     **     Most work is in this subroutine, safecbmc.f, and close.f.   ** 
c     *********************************************************************
c     **     Presently, works for linear molecules with rigid bond       **
c     **     lengths, and can (with little program changes) work for     **
c     **     branched molecules with rigid bonds, as long as it closes   **
c     **     at a binary or tertiary segment.                            **
c     *********************************************************************
c     **     Does work for any branched molecule with flexible bond      **
c     **     lengths.                                                    **
c     *********************************************************************

c     *********************************************************************
c     **                 NEW LOGIC ONLY FOR SAFE-CMBC                    ** 
c     **  -------------------------------------------------------------  **
c     **  iend = beads to grow to                                        **
c     **  ipast = one bead past iend                                     **
c     **  inext = two beads past iend                                    **
c     **  ibef = one bead before iend                                    **
c     **  iwbef = two beads before iend                                  **
c     **  fclose(iu) = beads to calculate interaction with from iu       **
c     **  fcount(iu) = number of fcloses for iu                          **
c     **  COLLIN = MASTER of the known universe                          **
c     *********************************************************************




      implicit none

      include 'control.inc'
      include 'cbmc.inc'
      include 'connect.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'fix.inc'
      include 'inputdata.inc'
      include 'rosen.inc'


      logical lcount,lpick,lterm,lfixed,lfix,lfind

      integer igrow,imolty,count,counta,iw,ivib,iv,iu,ju,iutry
      integer j,ja,kickout,invtry,index,fintnum,fint,k,islen
      integer ffrom,fprev,flist,fnum,fnuma,findex,countb,iv1
      integer movetype,fmaxgrow,ku,kickouta,iufrom,iuprev
      integer num,inum,inuma,max

      parameter(max=10)

      double precision random,rbf

      dimension fint(numax),ffrom(numax,max),fprev(numax,max)
     +     ,flist(numax,max,max),fnum(numax),fnuma(numax,max)
     +     ,lpick(numax),lfix(numax),inum(max),inuma(max),lfind(numax)

c     --------------------------------------------------------------------
      
c     *** safecbmc scheduler ****

c     --- now determine findex excluding 1, it won't do anything

c     --- set begining conditions
      if (movetype.eq.2) then
         if (.not.lring(imolty)) then
            stop 'you can not use safecbmc for swap unless it is a ring'
         endif
         fmaxgrow = nunit(imolty)
      else
         fmaxgrow = maxgrow(imolty) + 1
      endif
      lterm = .false.
      kickout = 0

      kickouta = 0

 100  continue

      
      do j = 1, nunit(imolty)
         lfix(j) = .false.
      enddo
      
      if (kickouta.eq.5) stop

      if (movetype.eq.2) then
         findex = fmaxgrow
      else
         findex = int( random() * dble(fmaxgrow - 1)) + 2
      endif

c     ******************************
c      findex = 7
c     ******************************

      lfixed = .false.

      if (findex.eq.2) then
         lcrank = .true.
      else
         lcrank = .false.
      endif

      do j = 1, igrow
         lexshed(j) = .true.
      enddo

      if (kickout.gt.250) then
         stop 'SAFESCHEDULE KICKED YOU OUT'
      endif
 
c     --- find iutry 
      if (movetype.eq.2) then

         iutry = 1
         invtry = invib(imolty,iutry) - 1

         if (invtry.lt.1) then
            stop 'You need a ring to do this'
         endif

         ffrom(1,1) = iutry
         fprev(1,1) = 0
         do iv = 1, invtry
            flist(1,1,iv) = ijvib(imolty,iutry,iv)
            lexshed(flist(1,1,iv)) = .false.
         enddo
         fnum(1) = 1
         fnuma(1,1) = invtry
      else
c--- JLR 11-11-09
c--- factoring in icbsta for safecbmc 

c --- OLD WAY   
c         iutry = int( random() * dble(iring(imolty)) ) + 1
c --- NEW WAY - picks between icbsta and last unit of chain
         iutry = int( random() * dble(iring(imolty)+icbsta(imolty)))+1
     &        -icbsta(imolty)

c --- need the following of scheduler gets confused!
         if (lrplc(imolty)) then
            if(iutry .eq.7 .or. iutry .eq.8) then 
               iutry=5
            endif
            if (iutry.eq.5.and.lcrank) goto 100
         endif
c --- END JLR 11-11-09

c     *************************
c         iutry = 1
c     *************************

         if (lplace(imolty,iutry)) then
            kickout = kickout + 1
            goto 100
         endif

         invtry = invib(imolty,iutry)
         if (invtry.eq.0) then
            stop 'cant do safecbmc on single bead'
         elseif(invtry.eq.1) then
c     --- At the end point of a molecule

c     --- we will let regular cbmc handle the end points
            kickout = kickout + 1
            goto 100

            fprev(1,1) = 0
            ivib = 0
         else
c     --- at a branch point, decide which way not to grow
c --- JLR 11-11-09
c --- adding statements so we don't get messed up at the branch point in ODS chains 
c            ivib = int(random() * dble(invtry)) + 1
 13         ivib = int(random() * dble(invtry)) + 1
c     ********************************
c            ivib = 2
c     *******************************

            fprev(1,1) = ijvib(imolty,iutry,ivib)
            if (lrplc(imolty)) then
               if (icbdir(imolty).eq.1.and.fprev(1,1).gt.iutry) goto 13
            endif
c --- END JLR 11-11-09


            if (fprev(1,1).gt.iring(imolty)
     &           .or.lplace(imolty,fprev(1,1))) then
               kickout = kickout + 1
               goto 100
            endif
         endif
         ffrom(1,1) = iutry
         count = 0
         do iv = 1, invtry
            if (iv.ne.ivib.and.
     &           .not.lplace(imolty,ijvib(imolty,iutry,iv))) then
               count = count + 1
               flist(1,1,count) = ijvib(imolty,iutry,iv)
               lexshed(flist(1,1,count)) = .false.
            endif
         enddo

         if (count.eq.0) then
            kickout = kickout + 1
            goto 100
         endif
         fnum(1) = 1
         fnuma(1,1) = count
      endif

c     --- find all branches going to maxgrow or end of molecule
      do iw = 2, nunit(imolty)
         count = 0
         do j = 1, fnum(iw-1)
            do ja = 1, fnuma(iw-1,j)
               counta = 0
               lcount = .false.
               iu = flist(iw-1,j,ja)
               if (.not. (lfix(iu).and.lfixed).and.
     &              .not.lrigi(imolty,iu)) then
                  do iv = 1, invib(imolty,iu)
                     ju = ijvib(imolty,iu,iv)
                     if (ju.ne.ffrom(iw-1,j).and.
     &                    .not.(lplace(imolty,ju).and.
     &                    iu.le.iring(imolty))) then
                        if (lfixed) then
                           if (ju.gt.iring(imolty)) then
                              counta = counta + 1
                              flist(iw,count+1,counta) = ju
                              lexshed(ju) = .false.
                              lcount = .true.
                           endif
                        else
                           counta = counta + 1
                           flist(iw,count+1,counta) = ju
                           lexshed(ju) = .false.
                           lcount = .true.
                        endif                     
                     endif
                  enddo
                  if (lcount) then
                     count = count + 1
                     fprev(iw,count) = ffrom(iw-1,j)
                     ffrom(iw,count) = iu
                     fnuma(iw,count) = counta
                  endif
               endif
            enddo
         enddo
         
         if (count.eq.0) then
c     --- we hit the end, lets get out of here
            index = iw - 1
            goto 110
         endif
         fnum(iw) = count
         if (iw.eq.findex) then
            lfixed = .true.
            do j = 1, fnum(findex)
               do ja = 1, fnuma(findex,j)
                  if (flist(findex,j,ja).le.iring(imolty)) then
                     lfix(flist(findex,j,ja)) = .true.
                  endif
               enddo
            enddo
         endif
      enddo
c      index = fmaxgrow
 110  continue

c     --- don't allow 1-bead regrowths, it won't do anything
      if (index.lt.2) then
         kickout = kickout + 1
         goto 100
      endif

      if (index.lt.findex) then
         kickout = kickout + 1
         goto 100
      endif

      lfixed = .false.

c     --- lets set logic so rosenbluth can read it
      count = 0
      do iw = 1, index 

         if (iw.eq.findex) then
            lfixed = .true.
         endif

         do j = 1, fnum(iw)
            kickout = 0

            if (lfixed) then
               if (flist(iw,j,1).le.iring(imolty)) goto 122
            endif
            count = count + 1 
            grownum(count) = fnuma(iw,j)
            growfrom(count) = ffrom(iw,j)
            growprev(count) = fprev(iw,j)
            do ja = 1, fnuma(iw,j)
               lpick(ja) = .false.
            enddo
            do ja = 1, fnuma(iw,j)
 115           continue
               if (kickout.gt.150) then
                  stop 'Randomizer in FECMBC kicked you out'
               endif
c     --- this is the only random part
               counta = int(random() * dble(fnuma(iw,j))) + 1
               if (lpick(counta)) then
                  kickout = kickout + 1
                  goto 115
               endif
               lpick(counta) = .true.
               growlist(count,counta) = flist(iw,j,ja)
            enddo
c     --- reset our logic keep have counta coicide 
            do ja = 1, fnuma(iw,j)
               flist(iw,j,ja) = growlist(count,ja)
            enddo
 122        continue 
         enddo
      enddo
      islen = count

c     --- set all fixed points to true
      do iw = findex, index 
         do count = 1, fnum(iw)
            do counta = 1, fnuma(iw,count)
               iu = flist(iw,count,counta)
               if (iu.le.iring(imolty)) then
                  lexshed(iu) = .true.
               endif
            enddo
         enddo
      enddo

c     --- find ends, ipasts, and inexts
      count = 0
      do j = 1, fnum(findex)
         do ja = 1, fnuma(findex,j)
            if (flist(findex,j,ja).le.iring(imolty)
     &           .and..not.lplace(imolty,flist(findex,j,ja))) then
               count = count + 1
               iend(count) = flist(findex,j,ja)
               counta = 0
               do iv = 1, invib(imolty,iend(count))
                  iu = ijvib(imolty,iend(count),iv)
                  if (iu.ne.ffrom(findex,j)) then
                     counta = counta + 1
                     ipast(iend(count),counta) = iu
                     countb = 0
                     do iv1 = 1, invib(imolty,iu)
                        ju = ijvib(imolty,iu,iv1)
                        if (ju.ne.iend(count)) then
                           countb = countb + 1
                           inext(iu,countb) = ju
                        endif
                     enddo
                     nextnum(iu) = countb
                  endif
               enddo
               pastnum(iend(count)) = counta
            endif
         enddo
         endnum = count
      enddo

c     --- now that we found iends and ipasts, 
c     --- determine which beads to close with each iend
      
      do j = 1, igrow
         fcount(j) = 0
      enddo
      
      do count = 1, endnum
         fintnum = 1
         fint(1) = iend(count)
         do iw = findex, 2, -1
            counta = 0
            do k = 1, fintnum
               do iv = 1, invib(imolty,fint(k))
                  ju = ijvib(imolty,fint(k),iv)
                  do j = 1, fnum(iw)
                     if (ju.eq.ffrom(iw,j)) then
                        fcount(ju) = fcount(ju) + 1
                        fclose(ju,fcount(ju)) = iend(count)
                        counta = counta + 1
                        fint(counta) = ju
                     endif
                  enddo
               enddo
            enddo
         enddo
         fintnum = counta
      enddo

c     --- define iwbef and ibef
      count = 0
      if (lcrank) then
         do j = 1, fnum(findex-1)
            counta = 0
            do ja = 1, fnuma(findex-1,j)
               iu = flist(findex-1,j,ja)
               if (fcount(iu).ne.0) then
                  counta = counta + 1
                  ibef(j,counta) = iu
               endif
            enddo
            befnum(j) = counta
         enddo
         wbefnum = fnum(findex-1)
      else
         do j = 1, fnum(findex - 1)
            iu = ffrom(findex-1,j)
            if (fcount(iu).ne.0) then
               count = count + 1
               iwbef(count) = iu
               counta = 0
               do ja = 1, fnuma(findex-1,j)
                  ju = flist(findex-1,j,ja)
                  if (fcount(ju).ne.0) then
                     counta = counta + 1
                     ibef(count,counta) = ju
                  endif
               enddo
               befnum(count) = counta
            endif
         enddo
         wbefnum = count
      endif
      
c      write(iou,*) '**********************************'
c      write(iou,*) iutry,ivib,findex,fprev(1,1)


c     --- set up place move logic -----------

      do j = 1, nunit(imolty)
         lpnow(j) = .false.
         pnum(j) = 0
      enddo

      nplace = 0
     
      counta = 0

      iw = 1
      do j = 1, fnum(iw)
         iufrom = ffrom(iw,j)
         if (iufrom.le.iring(imolty)) then
            counta = counta + 1
            do ivib = 1, invib(imolty,iufrom)
               iu = ijvib(imolty,iufrom,ivib)
               
               if (lplace(imolty,iu)) then
                  pnum(counta) = pnum(counta) + 1
                  if (pnum(counta).eq.1) then
                     nplace = nplace + 1
                  endif
                  lexshed(iu) = .false.
                  pprev(nplace) = fprev(iw,j)
                  pfrom(nplace) = iufrom
                  iplace(nplace,pnum(counta)) = iu
                  lpnow(iu) = .true. 
               endif
            enddo
            if (pnum(counta).eq.0) counta = counta - 1
         endif
      enddo

      do iw = 1, findex - 1
         do j = 1, fnum(iw)
            iuprev = ffrom(iw,j)
            do ja = 1, fnuma(iw,j)
               iufrom = flist(iw,j,ja)
               if (iufrom.le.iring(imolty)) then
                  counta = counta + 1
                  do ivib = 1, invib(imolty,iufrom)
                     iu = ijvib(imolty,iufrom,ivib)
                     
                     if (lplace(imolty,iu)) then
                        pnum(counta) = pnum(counta) + 1
                        if (pnum(counta).eq.1) then
                           nplace = nplace + 1
                        endif
                        lexshed(iu) = .false.
                        pprev(nplace) = iuprev
                        pfrom(nplace) = iufrom
                        iplace(nplace,pnum(counta)) = iu
                        lpnow(iu) = .true. 
                     endif
                  enddo
                  if (pnum(counta).eq.0) counta = counta - 1
               endif
            enddo
         enddo
      enddo
      
c     ------ end place move logic setup ---------

      

c     ---- begin part for rig logic

      do iu = 1, nunit(imolty)
         lfind(iu) = .false.
      enddo

      counta = 0
      do iw = 1, islen
         iufrom = growfrom(iw)
         do count = 1, grownum(iw)
            iu = growlist(iw,count)

            if (lrigi(imolty,iu)) then
               counta = counta + 1
               rfrom(counta) = iu
               rprev(counta) = iufrom

               lfind(iu) = .true.
               lfind(iufrom) = .true.

               ja = 0
               do ivib = 1, invib(imolty,iu)
                  ju = ijvib(imolty,iu,ivib)
                  if (ju.ne.iufrom.and..not.lpnow(ju)) then
                     ja = ja + 1
                     rlist(counta,ja) = ju
                     lfind(ju) = .true.
                     lexshed(ju) = .false.
                  endif
               enddo

               if (ja.eq.0) then
                  stop 'PROBLEM WITH RIG LOGIC IN SAFESCHEDULE'
               else
                  rnum(counta) = ja
               endif
            endif

         enddo
      enddo
      
      if (counta.eq.0) then
         llrig = .false.
      else
         do iu = 1, islen
            lsave(iw) = .false.
         enddo
         llrig = .true.
         nrigi = counta
      endif

      if (llrig) then
c     --- cycle through the rest of the rigid beads to set lexshed to false

         do iw = 1, nrigi
            do count = 1, rnum(iw)
               iu = rlist(iw,count)
               inum(count) = iu 
            enddo
            
            num = rnum(iw)
 5          continue
            ja = 0
            
            do count = 1, num
               iu = inum(count)

               do iv = 1, invib(imolty,iu)

                  ju = ijvib(imolty,iu,iv)

                  if (.not.lfind(ju)) then

                     ja = ja + 1
                     
                     lfind(ju) = .true.
                     inuma(ja) = ju
                     
                     lexshed(ju) = .false.
                  endif
               enddo
            enddo
            if (ja.gt.10) then
               write(iou,*) 'ja',ja
               stop 'need to set max greater in safeschedule'
            endif

            num = ja

            if (ja.ne.0) then
               do j = 1, ja
                  inum(j) = inuma(j)
               enddo
               goto 5
            endif
         enddo

      endif

      return


c     ****** take out return for diagnostics ************

      do iw = 1, islen
         write(iou,*) growfrom(iw),(growlist(iw,count),count=1
     &        ,grownum(iw))

      enddo

      write(iou,*) '---------------------------------------------'

      do iw = 1, nplace
         write(iou,*) pfrom(iw),(iplace(iw,count),count=1,pnum(iw))
      enddo

      write(iou,*) '--------------------------------------------'

      do iw = 1, nrigi
         write(iou,*) rfrom(iw),(rlist(iw,count),count=1,rnum(iw))
      enddo
      
      stop
            
c     ***************************************************

      return
      end


















