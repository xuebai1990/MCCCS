      subroutine safeschedule(igrow,imolty,islen,iutry,findex,movetype)

!     *********************************************************************
!     **  Self-Adapting Fixed-Endpoint Configurational-Bias Monte Carlo  **
!     **                    --- SAFE-CBMC ---                            **
!     *********************************************************************
!     **    Determines logic for a CBMC Move Between Fixed End Points    **
!     **           for Linear, Branched, and Cyclic Molecules            **
!     *********************************************************************
!     **       Originally completed by Collin Wick around 1-1-2000       **
!     *********************************************************************


!     *********************************************************************
!     **     Most work is in this subroutine, safecbmc.f, and close.f.   ** 
!     *********************************************************************
!     **     Presently, works for linear molecules with rigid bond       **
!     **     lengths, and can (with little program changes) work for     **
!     **     branched molecules with rigid bonds, as long as it closes   **
!     **     at a binary or tertiary segment.                            **
!     *********************************************************************
!     **     Does work for any branched molecule with flexible bond      **
!     **     lengths.                                                    **
!     *********************************************************************

!     *********************************************************************
!     **                 NEW LOGIC ONLY FOR SAFE-CMBC                    ** 
!     **  -------------------------------------------------------------  **
!     **  iend = beads to grow to                                        **
!     **  ipast = one bead past iend                                     **
!     **  inext = two beads past iend                                    **
!     **  ibef = one bead before iend                                    **
!     **  iwbef = two beads before iend                                  **
!     **  fclose(iu) = beads to calculate interaction with from iu       **
!     **  fcount(iu) = number of fcloses for iu                          **
!     **  COLLIN = MASTER of the known universe                          **
!     *********************************************************************

      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_runtime,only:err_exit
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'connect.inc'
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'fix.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'rosen.inc'


      logical::lcount,lpick,lterm,lfixed,lfix,lfind

      integer(KIND=normal_int)::igrow,imolty,count,counta,iw,ivib,iv,iu ,ju,iutry
      integer(KIND=normal_int)::j,ja,kickout,invtry,index,fintnum,fint,k ,islen
      integer(KIND=normal_int)::ffrom,fprev,flist,fnum,fnuma,findex ,countb,iv1
      integer(KIND=normal_int)::movetype,fmaxgrow,kickouta,iufrom ,iuprev
      integer(KIND=normal_int)::num,inum,inuma,max

      parameter(max=10)

      real(KIND=double_precision)::random

      dimension fint(numax),ffrom(numax,max),fprev(numax,max) ,flist(numax,max,max),fnum(numax),fnuma(numax,max) ,lpick(numax) ,lfix(numax),inum(max),inuma(max),lfind(numax)

!     --------------------------------------------------------------------
      
!     *** safecbmc scheduler ****

!     --- now determine findex excluding 1, it won't do anything

!     --- set begining conditions
      if (movetype.eq.2) then
         if (.not.lring(imolty)) then
            call err_exit('you can not use safecbmc for swap unless it is a ring')
         end if
         fmaxgrow = nunit(imolty)
      else
         fmaxgrow = maxgrow(imolty) + 1
      end if
      lterm = .false.
      kickout = 0

      kickouta = 0

 100  continue

      
      do j = 1, nunit(imolty)
         lfix(j) = .false.
      end do
      
      if (kickouta.eq.5) call err_exit('')

      if (movetype.eq.2) then
         findex = fmaxgrow
      else
         findex = int( random() * dble(fmaxgrow - 1)) + 2
      end if

!     ******************************
!      findex = 7
!     ******************************

      lfixed = .false.

      if (findex.eq.2) then
         lcrank = .true.
      else
         lcrank = .false.
      end if

      do j = 1, igrow
         lexshed(j) = .true.
      end do

      if (kickout.gt.250) then
         call err_exit('SAFESCHEDULE KICKED YOU OUT')
      end if
 
!     --- find iutry 
      if (movetype.eq.2) then

         iutry = 1
         invtry = invib(imolty,iutry) - 1

         if (invtry.lt.1) then
            call err_exit('You need a ring to do this')
         end if

         ffrom(1,1) = iutry
         fprev(1,1) = 0
         do iv = 1, invtry
            flist(1,1,iv) = ijvib(imolty,iutry,iv)
            lexshed(flist(1,1,iv)) = .false.
         end do
         fnum(1) = 1
         fnuma(1,1) = invtry
      else
!--- JLR 11-11-09
!--- factoring in icbsta for safecbmc 

! --- OLD WAY   
!         iutry = int( random() * dble(iring(imolty)) ) + 1
! --- NEW WAY - picks between icbsta and last unit of chain
         iutry = int( random() * dble(iring(imolty)+icbsta(imolty)))+1 -icbsta(imolty)

! --- need the following of scheduler gets confused!
         if (lrplc(imolty)) then
            if(iutry .eq.7 .or. iutry .eq.8) then 
               iutry=5
            end if
            if (iutry.eq.5.and.lcrank) goto 100
         end if
! --- END JLR 11-11-09

!     *************************
!         iutry = 1
!     *************************

         if (lplace(imolty,iutry)) then
            kickout = kickout + 1
            goto 100
         end if

         invtry = invib(imolty,iutry)
         if (invtry.eq.0) then
            call err_exit('cant do safecbmc on single bead')
         elseif(invtry.eq.1) then
!     --- At the end point of a molecule

!     --- we will let regular cbmc handle the end points
            kickout = kickout + 1
            goto 100

            fprev(1,1) = 0
            ivib = 0
         else
!     --- at a branch point, decide which way not to grow
! --- JLR 11-11-09
! --- adding statements so we don't get messed up at the branch point in ODS chains 
!            ivib = int(random() * dble(invtry)) + 1
 13         ivib = int(random() * dble(invtry)) + 1
!     ********************************
!            ivib = 2
!     *******************************

            fprev(1,1) = ijvib(imolty,iutry,ivib)
            if (lrplc(imolty)) then
               if (icbdir(imolty).eq.1.and.fprev(1,1).gt.iutry) goto 13
            end if
! --- END JLR 11-11-09


            if (fprev(1,1).gt.iring(imolty) .or.lplace(imolty,fprev(1,1))) then
               kickout = kickout + 1
               goto 100
            end if
         end if
         ffrom(1,1) = iutry
         count = 0
         do iv = 1, invtry
            if (iv.ne.ivib.and. .not.lplace(imolty,ijvib(imolty,iutry,iv))) then
               count = count + 1
               flist(1,1,count) = ijvib(imolty,iutry,iv)
               lexshed(flist(1,1,count)) = .false.
            end if
         end do

         if (count.eq.0) then
            kickout = kickout + 1
            goto 100
         end if
         fnum(1) = 1
         fnuma(1,1) = count
      end if

!     --- find all branches going to maxgrow or end of molecule
      do iw = 2, nunit(imolty)
         count = 0
         do j = 1, fnum(iw-1)
            do ja = 1, fnuma(iw-1,j)
               counta = 0
               lcount = .false.
               iu = flist(iw-1,j,ja)
               if (.not. (lfix(iu).and.lfixed).and. .not.lrigi(imolty,iu)) then
                  do iv = 1, invib(imolty,iu)
                     ju = ijvib(imolty,iu,iv)
                     if (ju.ne.ffrom(iw-1,j).and. .not.(lplace(imolty,ju).and. iu.le.iring(imolty))) then
                        if (lfixed) then
                           if (ju.gt.iring(imolty)) then
                              counta = counta + 1
                              flist(iw,count+1,counta) = ju
                              lexshed(ju) = .false.
                              lcount = .true.
                           end if
                        else
                           counta = counta + 1
                           flist(iw,count+1,counta) = ju
                           lexshed(ju) = .false.
                           lcount = .true.
                        end if                     
                     end if
                  end do
                  if (lcount) then
                     count = count + 1
                     fprev(iw,count) = ffrom(iw-1,j)
                     ffrom(iw,count) = iu
                     fnuma(iw,count) = counta
                  end if
               end if
            end do
         end do
         
         if (count.eq.0) then
!     --- we hit the end, lets get out of here
            index = iw - 1
            goto 110
         end if
         fnum(iw) = count
         if (iw.eq.findex) then
            lfixed = .true.
            do j = 1, fnum(findex)
               do ja = 1, fnuma(findex,j)
                  if (flist(findex,j,ja).le.iring(imolty)) then
                     lfix(flist(findex,j,ja)) = .true.
                  end if
               end do
            end do
         end if
      end do
!      index = fmaxgrow
 110  continue

!     --- don't allow 1-bead regrowths, it won't do anything
      if (index.lt.2) then
         kickout = kickout + 1
         goto 100
      end if

      if (index.lt.findex) then
         kickout = kickout + 1
         goto 100
      end if

      lfixed = .false.

!     --- lets set logic so rosenbluth can read it
      count = 0
      do iw = 1, index 

         if (iw.eq.findex) then
            lfixed = .true.
         end if

         do j = 1, fnum(iw)
            kickout = 0

            if (lfixed) then
               if (flist(iw,j,1).le.iring(imolty)) goto 122
            end if
            count = count + 1 
            grownum(count) = fnuma(iw,j)
            growfrom(count) = ffrom(iw,j)
            growprev(count) = fprev(iw,j)
            do ja = 1, fnuma(iw,j)
               lpick(ja) = .false.
            end do
            do ja = 1, fnuma(iw,j)
 115           continue
               if (kickout.gt.150) then
                  call err_exit('Randomizer in FECMBC kicked you out')
               end if
!     --- this is the only random part
               counta = int(random() * dble(fnuma(iw,j))) + 1
               if (lpick(counta)) then
                  kickout = kickout + 1
                  goto 115
               end if
               lpick(counta) = .true.
               growlist(count,counta) = flist(iw,j,ja)
            end do
!     --- reset our logic keep have counta coicide 
            do ja = 1, fnuma(iw,j)
               flist(iw,j,ja) = growlist(count,ja)
            end do
 122        continue 
         end do
      end do
      islen = count

!     --- set all fixed points to true
      do iw = findex, index 
         do count = 1, fnum(iw)
            do counta = 1, fnuma(iw,count)
               iu = flist(iw,count,counta)
               if (iu.le.iring(imolty)) then
                  lexshed(iu) = .true.
               end if
            end do
         end do
      end do

!     --- find ends, ipasts, and inexts
      count = 0
      do j = 1, fnum(findex)
         do ja = 1, fnuma(findex,j)
            if (flist(findex,j,ja).le.iring(imolty) .and..not.lplace(imolty,flist(findex,j,ja))) then
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
                        end if
                     end do
                     nextnum(iu) = countb
                  end if
               end do
               pastnum(iend(count)) = counta
            end if
         end do
         endnum = count
      end do

!     --- now that we found iends and ipasts, 
!     --- determine which beads to close with each iend
      
      do j = 1, igrow
         fcount(j) = 0
      end do
      
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
                     end if
                  end do
               end do
            end do
         end do
         fintnum = counta
      end do

!     --- define iwbef and ibef
      count = 0
      if (lcrank) then
         do j = 1, fnum(findex-1)
            counta = 0
            do ja = 1, fnuma(findex-1,j)
               iu = flist(findex-1,j,ja)
               if (fcount(iu).ne.0) then
                  counta = counta + 1
                  ibef(j,counta) = iu
               end if
            end do
            befnum(j) = counta
         end do
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
                  end if
               end do
               befnum(count) = counta
            end if
         end do
         wbefnum = count
      end if
      
!      write(iou,*) '**********************************'
!      write(iou,*) iutry,ivib,findex,fprev(1,1)


!     --- set up place move logic -----------

      do j = 1, nunit(imolty)
         lpnow(j) = .false.
         pnum(j) = 0
      end do

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
                  end if
                  lexshed(iu) = .false.
                  pprev(nplace) = fprev(iw,j)
                  pfrom(nplace) = iufrom
                  iplace(nplace,pnum(counta)) = iu
                  lpnow(iu) = .true. 
               end if
            end do
            if (pnum(counta).eq.0) counta = counta - 1
         end if
      end do

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
                        end if
                        lexshed(iu) = .false.
                        pprev(nplace) = iuprev
                        pfrom(nplace) = iufrom
                        iplace(nplace,pnum(counta)) = iu
                        lpnow(iu) = .true. 
                     end if
                  end do
                  if (pnum(counta).eq.0) counta = counta - 1
               end if
            end do
         end do
      end do
      
!     ------ end place move logic setup ---------

      

!     ---- begin part for rig logic

      do iu = 1, nunit(imolty)
         lfind(iu) = .false.
      end do

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
                  end if
               end do

               if (ja.eq.0) then
                  call err_exit('PROBLEM WITH RIG LOGIC IN SAFESCHEDULE')
               else
                  rnum(counta) = ja
               end if
            end if

         end do
      end do
      
      if (counta.eq.0) then
         llrig = .false.
      else
         do iu = 1, islen
            lsave(iw) = .false.
         end do
         llrig = .true.
         nrigi = counta
      end if

      if (llrig) then
!     --- cycle through the rest of the rigid beads to set lexshed to false

         do iw = 1, nrigi
            do count = 1, rnum(iw)
               iu = rlist(iw,count)
               inum(count) = iu 
            end do
            
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
                  end if
               end do
            end do
            if (ja.gt.10) then
               write(iou,*) 'ja',ja
               call err_exit('need to set max greater in safeschedule')
            end if

            num = ja

            if (ja.ne.0) then
               do j = 1, ja
                  inum(j) = inuma(j)
               end do
               goto 5
            end if
         end do

      end if

      return


!     ****** take out return for diagnostics ************

      do iw = 1, islen
         write(iou,*) growfrom(iw),(growlist(iw,count),count=1 ,grownum(iw))

      end do

      write(iou,*) '---------------------------------------------'

      do iw = 1, nplace
         write(iou,*) pfrom(iw),(iplace(iw,count),count=1,pnum(iw))
      end do

      write(iou,*) '--------------------------------------------'

      do iw = 1, nrigi
         write(iou,*) rfrom(iw),(rlist(iw,count),count=1,rnum(iw))
      end do
      
      call err_exit('')
            
!     ***************************************************

      return
      end


















