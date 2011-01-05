      subroutine schedule(igrow,imolty,index,iutry,iprev,movetype)

!     *******************************************************************
!     ** computes the growth shedule for CBMC type moves               **
!     *******************************************************************
!     movetype=1: called from config.f
      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'rosen.inc'
!$$$      include 'connect.inc'
!$$$      include 'inputdata.inc'

      logical::lprint,lfind
      integer(KIND=normal_int)::random_index
      integer(KIND=normal_int)::kickout,icbu,igrow,imolty,iutry,iut ,invtry,iu,ju
      integer(KIND=normal_int)::ibead,count,ivib,idir,movetype,iprev ,itry,i,ib2
      integer(KIND=normal_int)::temp_store,temp_count,izz,outer_sites ,index,outer_num,outer_prev,iufrom,outer_try
      real(KIND=double_precision)::dbgrow,random
      parameter (lprint = .false.)
      dimension temp_store(numax),outer_sites(numax),outer_prev(numax) ,lfind(numax)
! ------------------------------------------------------------------

!      write(iou,*) 'start SCHEDULE'

!     DECODER for the new growth logic
!     lexshed is true if the bead exists at that time of the growth 
!          it is false when a bead has not yet been grown this time
!     growfrom(index) gives the unit number that is to be grown from for
!          the index step
!     growprev(index) gives the bead that exists connected to growfrom(index)
!     grownum(index) gives the number of beads to be grown from growfrom(index)
!     growlist(index,count) gives the unit numbers for the beads grown 
!          from growfrom(index), get count from grownum(index)
     
      kickout = 0

! * dummy loop to set position of 11 for movetype 1
 11   do iu = 1,1
      end do

!     --- initialize temp_count
      temp_count = 0

!     this part is just for config right now
      if ( movetype .eq. 1 ) then 
    
!     --- initialize lexshed so all beads currently exist
         do iu = 1,igrow
            lexshed(iu) = .true.
         end do

         if ( kickout .eq. 50 ) call cleanup('kickout is 50')
 
!        -- select the first bead to grow from
         if (lrigid(imolty)) then
! N.R. Randomly select one of the grow points
            random_index = int(dble(rindex(imolty))*random()+1)
            iutry = riutry(imolty,random_index)
         else if ( lrig(imolty).and.nrig(imolty).gt.0 ) then
            dbgrow = random()*dble(nrig(imolty)) + 1.0d0
            iutry = irig(imolty,int(dbgrow))
         else if ( icbsta(imolty) .gt. 0 ) then
            iutry = icbsta(imolty)
         else if ( icbsta(imolty) .eq. 0 ) then
            dbgrow = dble( igrow )
            iutry = int( dbgrow*random() ) + 1
         else
            dbgrow = dble( igrow + icbsta(imolty) + 1 ) 
            iutry = int( dbgrow*random() ) - icbsta(imolty)       
         end if
        
         if (lrig(imolty).and.(nrig(imolty).eq.0)) then
                        
            ju = int(random() * dble(nrigmax(imolty)  - nrigmin(imolty) + 1)) + nrigmin(imolty)
         end if
         

!         write(iou,*) '********************************************'
!         write(iou,*) 'iutry',iutry

         idir = icbdir(imolty)
         growfrom(1) = iutry
         invtry = invib(imolty,iutry)
         index = 1

         if ( invtry .eq. 0 ) then
!           --- problem, cannot do config move on a 1 bead molecule
            call cleanup('cannot do CBMC on a one grow unit molecule')

         elseif ( invtry .eq. 1 ) then
!           --- regrow entire molecule, check nmaxcbmc
            if ( nmaxcbmc(imolty) .lt. igrow - 1 ) then
               kickout = kickout + 1
               goto 11
            end if
            
            growprev(1) = 0
            
            grownum(1) = invtry
            ivib = invtry
            iut = ijvib(imolty,iutry,ivib)
            if ( idir .eq. 1 .and. iut .lt. iutry ) then
!              --- cannot grow from this bead, try again
               kickout = kickout + 1
               goto 11
            end if
            growlist(1,ivib) = iut
            lexshed(iut) = .false.
         else

!           --- at a branch point, decide how many branches to regrow
!           --- regrow all of legal branches (check idir )
            count = 0
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)               
               if ( idir .eq. 1 .and. iut .lt. iutry ) then
!                 --- this is the one and only previous nongrown bead
                  growprev(1) = iut
               else
!                 --- grow these branches
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               end if
            end do

            do izz = temp_count,1,-1
!              --- choose grow bead randomly from temp_store
               itry = int( dble(izz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(1,count) = iut
               lexshed(iut) = .false.

!              --- update temp_store for next iteration
               temp_store(itry) = temp_store(izz)
            end do
            grownum(1) = count

            if ( count .eq. invtry ) then
               if ( random() .gt. pmall(imolty) ) then
!                 --- not pmall so select a previous bead
 20               ivib = int( random()*dble(invtry) ) + 1

                  iu = growlist(1,ivib)
                  
                  if (lrig(imolty).and.nrig(imolty).gt.0) then
                     if (iu.ne.frig(imolty,int(dbgrow))) goto 20
                  end if
!     --- we should start to determine a random section to keep rigid
                       
                  if (lrigid(imolty)) then
                     if (iu.lt.riutry(imolty,1)) goto 20
                  end if

                  growprev(1) = iu
                  
!                 --- replace this unit with the last unit in growlist
                  growlist(1,ivib) = growlist(1,invtry)

!                 --- reduce numgrow by 1
                  grownum(1) = grownum(1) - 1
!                 --- add this back into lexshed
                  lexshed(iu) = .true.
               else
!              --- we regrew all branches, no previous bead
                  growprev(1) = 0
               end if
            else if ( invtry - count .ne. 1 ) then
!              --- problem in logic, should only be one nongrown bead
               write(iou,*) 'invtry,count',invtry,count
               write(iou,*) 'igrow,imolty',igrow,imolty
               call cleanup('logic problem in schedule')
            end if
         end if   

!     ---end the part that is specific for config
      elseif ( movetype .eq. 2 ) then
!     --- begin the part that is specific for swap
!        --- iutry is the first bead inserted - need to grow its neighbors
         do iu = 1,igrow
            lexshed(iu) = .true.
         end do

         growfrom(1) = iutry
 
         invtry = invib(imolty,iutry)
         growprev(1) = 0
         if ( invtry .eq. 0 ) then
!           --- Bead iutry is the only bead to be grown ---
            index = 0
         else
!           --- grow all of the beads connected to bead iutry
            index = 1
            grownum(1) = invtry
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)
               if (lrigid(imolty)) iut = iutry - 1
!              --- grow these branches
               temp_count = temp_count + 1
               temp_store(temp_count) = iut
            end do

            count = 0
            do izz = temp_count,1,-1
!              --- choose grow bead randomly from temp_store
               itry = int( dble(izz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(1,count) = iut
               lexshed(iut) = .false.
               
!              --- update temp_store for next iteration
               temp_store(itry) = temp_store(izz)
            end do
         end if

!        ---end the part that is specific for swap

      elseif ( movetype .eq. 3 ) then
         do iu = 1,igrow
            lexshed(iu) = .true.
         end do
!        --- begin part that is specific for swatch 
         if ( iutry .eq. 0 ) then
!           --- no beads to be regrown via cbmc
            index = 0
            return
         end if
         growfrom(1) = iutry
         
         invtry = invib(imolty,iutry)
         growprev(1) = iprev

         if ( invtry .eq. 0 ) then
!           --- Bead 1 is the only bead to be grown ---
            index = 0
         else

!           --- grow all (except iprev) of the beads connected to bead 1
            index = 1
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)
               if ( iut .ne. iprev ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               end if
            end do

            count = 0
            do izz = temp_count,1,-1
!              --- choose grow bead randomly from temp_store
               itry = int( dble(izz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(1,count) = iut
               lexshed(iut) = .false.

!              --- update temp_store for next iteration
               temp_store(itry) = temp_store(izz)
            end do

            grownum(1) = count
         end if
!        --- end part that is specific for swatch
         
      elseif ( movetype .eq. 4 ) then
         do iu = 1,igrow
            lexshed(iu) = .true.
         end do
!        --- begin part that is specific for rigid molecules

         if ( rindex(imolty) .eq. 0 ) then
!        --- entire molecule is rigid
            index = 0
         else
            index = rindex(imolty)
           
            do i = 1, rindex(imolty)
!        --- grow all non-rigid beads from the rigid part
               temp_count = 0
               invtry = invib(imolty,riutry(imolty,i)) - 1
               grownum(i) = invtry 
               growfrom(i) = riutry(imolty,i)

!        ---  for rigid molecules, have rigid vib last
               growprev(i)=ijvib(imolty,growfrom(i),invtry+1) 

               do ivib = 1, invtry
                  iut = ijvib(imolty,riutry(imolty,i),ivib)
!         --- grow these branches
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               end do

               count = 0

               do izz = temp_count,1,-1
!                 --- choose grow bead randomly from temp_store
                  itry = int ( dble(izz) * random() ) + 1
                  iut = temp_store(itry)
                  count = count + 1
                  growlist(i,count) = iut
                  lexshed(iut) = .false.

!                 --- update temp_store for next iteration
                  temp_store(itry) = temp_store(izz)
               end do
            end do
         end if


!        --- end part that is specific for rigid molecules
      elseif (movetype .eq. 5) then
         if ( iutry .eq. 0 ) then
!           --- no beads to be regrown via cbmc
            return
         end if
         growfrom(index+1) = iutry
         
         invtry = invib(imolty,iutry)
         growprev(index+1) = iprev

         if ( invtry .eq. 0 ) then
!           --- Bead 1 is the only bead to be grown ---
         else

!           --- grow all (except iprev) of the beads connected to bead 1
            index = index + 1
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)
               if ( iut .ne. iprev ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               end if
            end do

            count = 0
            do izz = temp_count,1,-1
!              --- choose grow bead randomly from temp_store
               itry = int( dble(izz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(index,count) = iut
               lexshed(iut) = .false.
!              --- update temp_store for next iteration
               temp_store(itry) = temp_store(izz)
            end do

            grownum(index) = count
         end if

! * end iswatch stuff *
 
      else

!        --- non-existant move type
         write(iou,*) 'schedule movetype ',movetype
         call cleanup('non-valid move type')
      end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     from here on down config, swap, and swatch are the same


!     OLD METHOD - not fully random - removed 6-13-98

!      if ( .false. ) then 
!         write(iou,*) 'old method'
!      ibead = 0
!      index = 0
! 50   if ( index .lt. islen ) then
!         index = index + 1
!         do icbu = 1,grownum(index)
!            ibead = ibead + 1
!            iu = growlist(index,icbu)
!            invtry = invib(imolty,iu)
!            if ( invtry .ne. 1) then
!              --- this will be the next growing bead
!               islen = islen + 1
!               growprev(islen) = growfrom(index)
!               growfrom(islen) = iu
!               count = 0
!               do ivib = 1,invtry
!                  iut = ijvib(imolty,iu,ivib)
!                  if ( iut .ne. growprev(islen) ) then
!                     count = count + 1
!                     growlist(islen,count) = iut
!                     lexshed(iut) = .false.
!                  end if
!               end do
!               grownum(islen) = count
!            end if
!         end do
!         goto 50
!      end if
      
!     --- set up list of "outer" beads that have further growth sites
!     --- this method implemented 6-13-98

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
! ibead should be set here, otherwise it's not initialized below
      ibead=0

      if ( index .gt. 0 ) then
         ibead = 0
         outer_num = 0
         iufrom = growfrom(index)
         do icbu = 1,grownum(index)
!           --- increment counter ibead for total number of beads grown
            ibead = ibead + 1

!           --- determine whether this bead has any non-grown neighbors
            iu = growlist(index,icbu)
            invtry = invib(imolty,iu)
         
            if ( invtry .gt. 1 ) then
!              --- add one to the stack of outer_sites
               outer_num = outer_num + 1
               outer_sites(outer_num) = iu
               outer_prev(outer_num) = iufrom
            end if
         end do

!        --- begin while loop to grow all outer beads until done
 70      if ( outer_num .gt. 0 ) then

!           --- choose one site randomly from the stack
            outer_try = int( dble(outer_num) * random() ) + 1

!           --- increment index to show this is the next growfrom
            index = index + 1

            iu = outer_sites(outer_try)
            iufrom = outer_prev(outer_try)
!           --- assign growfrom and growprev for this index
            growfrom(index) = iu
            growprev(index) = iufrom

            invtry = invib(imolty,iu)
!           --- assign the grow beads in random order
            temp_count = 0
            do ivib = 1,invtry
               iut = ijvib(imolty,iu,ivib)
               if ( iut .ne. iufrom ) then
!                 --- add to the list of beads to be grown from iu
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               end if
            end do

            count = 0
            do izz = temp_count, 1, -1
               itry = int ( dble(izz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
!              --- assign growlist for current index and count
               growlist(index,count) = iut
               lexshed(iut) = .false.

!              --- update temp_store for next iteration
               temp_store(itry) = temp_store(izz)
            end do
!           --- assign grownum for this index
            grownum(index) = count

!           --- update list of "outer" beads
!           --- remove bead that was just grown from outer list
            outer_sites(outer_try) = outer_sites(outer_num)
            outer_prev(outer_try) = outer_prev(outer_num)
            outer_num = outer_num - 1
         
!           --- add the new beads if they have more to be grown
            iufrom = iu
            do icbu = 1,grownum(index)
!              --- increment counter ibead for total number of beads grown
               ibead = ibead + 1

!              --- determine whether this bead has any non-grown neighbors
               iu = growlist(index,icbu)
               invtry = invib(imolty,iu)

               if ( invtry .gt. 1 ) then
!                 --- add one to the stack of outer_sites
                  outer_num = outer_num + 1
                  outer_sites(outer_num) = iu
                  outer_prev(outer_num) = iufrom
               end if
            end do

!           --- end of while loop 70
            goto 70
         end if
      end if

 100  if ( (movetype .eq. 1) .and. (ibead .gt. nmaxcbmc(imolty)) ) then
         kickout = kickout + 1
         goto 11
      end if

      if ( lprint ) then
         write(iou,*) 'movetype',movetype
         write(iou,*) 'index',index
         do ibead = 1,index
            write(iou,*) 'ibead',ibead
            write(iou,*) 'growfrom(ibead)',growfrom(ibead)
            write(iou,*) 'growprev(ibead)',growprev(ibead)
            write(iou,*) 'grownum(ibead)',grownum(ibead)
            do count = 1,grownum(ibead)
               write(iou,*) 'count,growlist(ibead,count)',count ,growlist(ibead,count)
            end do
         end do
         do iu = 1,igrow
            write(iou,*) 'iu,lexshed(iu)',iu,lexshed(iu)
         end do
      end if

      
      if (lrig(imolty).and.nrig(imolty).eq.0) then
         ib2 = index - ju
         if (ib2.lt.1) then
            llrig = .false.
            return
         else

            do ibead = 1, index
               lsave(ibead) = .false.
            end do

            nrigi = 1
            llrig = .true.
            
            rfrom(1) = growfrom(ib2)
            rprev(1) = growprev(ib2)
            
            rnum(1) = grownum(ib2)
            lsave(ib2) = .true.

            do count = 1, grownum(ib2)
               iu = growlist(ib2,count)

               lfind(iu) = .true.

               rlist(1,count) = iu
               lexshed(iu) = .false.
            end do
            

         end if
         
!     --- cycle through the rest of the rigid beads
         
         do ibead = ib2 + 1, index
            iufrom = growfrom(ibead)

            if (lfind(iufrom)) then
               
               lsave(ibead) = .true.
               
               do count = 1, grownum(ibead)
                  
                  iu = growlist(ibead,count)

                  lfind(iu) = .true.

               end do

            end if
         end do
      else
         llrig = .false.
      end if

      


!      write(iou,*) 'end SCHEDULE'
  
      return
      end
