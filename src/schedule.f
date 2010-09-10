      subroutine schedule(igrow,imolty,index,iutry,iprev,movetype)

c schedule
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
 
c     *******************************************************************
c     ** computes the growth shedule for CBMC type moves               **
c     *******************************************************************
 
      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'cbmc.inc'
      include 'rosen.inc'
      include 'connect.inc'
      include 'inputdata.inc'

      logical::lprint,lfind
      integer::random_index
      integer::kickout,icbu,igrow,imolty,iutry,iut,invtry,iu,ju
      integer::ibead,count,ivib,idir,movetype,iprev,itry,i,ib2
      integer::temp_store,temp_count,zz,outer_sites,index,outer_num
     &     ,outer_prev,iufrom,outer_try
      real(8)::dbgrow,random
      parameter (lprint = .false.)
      dimension temp_store(numax),outer_sites(numax),outer_prev(numax)
     &     ,lfind(numax)
c ------------------------------------------------------------------

c      write(iou,*) 'start SCHEDULE'

c     DECODER for the new growth logic
c     lexshed is true if the bead exists at that time of the growth 
c          it is false when a bead has not yet been grown this time
c     growfrom(index) gives the unit number that is to be grown from for
c          the index step
c     growprev(index) gives the bead that exists connected to growfrom(index)
c     grownum(index) gives the number of beads to be grown from growfrom(index)
c     growlist(index,count) gives the unit numbers for the beads grown 
c          from growfrom(index), get count from grownum(index)
     
      kickout = 0

c * dummy loop to set position of 11 for movetype 1
 11   do iu = 1,1
      enddo

c     --- initialize temp_count
      temp_count = 0

c     this part is just for config right now
      if ( movetype .eq. 1 ) then 
    
c     --- initialize lexshed so all beads currently exist
         do iu = 1,igrow
            lexshed(iu) = .true.
         enddo

         if ( kickout .eq. 50 ) call cleanup('kickout is 50')
 
c        -- select the first bead to grow from
         if (lrigid(imolty)) then
! N.R. Randomly select one of the grow points
            random_index = int(dble(rindex(imolty))*random()+1)
            iutry = riutry(imolty,random_index)
         elseif ( lrig(imolty).and.nrig(imolty).gt.0 ) then
            dbgrow = random()*dble(nrig(imolty)) + 1.0d0
            iutry = irig(imolty,int(dbgrow))
         elseif ( icbsta(imolty) .gt. 0 ) then
            iutry = icbsta(imolty)
         elseif ( icbsta(imolty) .eq. 0 ) then
            dbgrow = dble( igrow )
            iutry = int( dbgrow*random() ) + 1
         else
            dbgrow = dble( igrow + icbsta(imolty) + 1 ) 
            iutry = int( dbgrow*random() ) - icbsta(imolty)       
         endif
        
         if (lrig(imolty).and.(nrig(imolty).eq.0)) then
                        
            ju = int(random() * dble(nrigmax(imolty) 
     &           - nrigmin(imolty) + 1)) + nrigmin(imolty)
         endif
         

c         write(iou,*) '********************************************'
c         write(iou,*) 'iutry',iutry

         idir = icbdir(imolty)
         growfrom(1) = iutry
         invtry = invib(imolty,iutry)
         index = 1

         if ( invtry .eq. 0 ) then
c           --- problem, cannot do config move on a 1 bead molecule
            call cleanup('cannot do CBMC on a one grow unit molecule')

         elseif ( invtry .eq. 1 ) then
c           --- regrow entire molecule, check nmaxcbmc
            if ( nmaxcbmc(imolty) .lt. igrow - 1 ) then
               kickout = kickout + 1
               goto 11
            endif
            
            growprev(1) = 0
            
            grownum(1) = invtry
            ivib = invtry
            iut = ijvib(imolty,iutry,ivib)
            if ( idir .eq. 1 .and. iut .lt. iutry ) then
c              --- cannot grow from this bead, try again
               kickout = kickout + 1
               goto 11
            endif
            growlist(1,ivib) = iut
            lexshed(iut) = .false.
         else

c           --- at a branch point, decide how many branches to regrow
c           --- regrow all of legal branches (check idir )
            count = 0
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)               
               if ( idir .eq. 1 .and. iut .lt. iutry ) then
c                 --- this is the one and only previous nongrown bead
                  growprev(1) = iut
               else
c                 --- grow these branches
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo

            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               itry = int( dble(zz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(1,count) = iut
               lexshed(iut) = .false.

c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo
            grownum(1) = count

            if ( count .eq. invtry ) then
               if ( random() .gt. pmall(imolty) ) then
c                 --- not pmall so select a previous bead
 20               ivib = int( random()*dble(invtry) ) + 1

                  iu = growlist(1,ivib)
                  
                  if (lrig(imolty).and.nrig(imolty).gt.0) then
                     if (iu.ne.frig(imolty,int(dbgrow))) goto 20
                  endif
c     --- we should start to determine a random section to keep rigid
                       
                  if (lrigid(imolty)) then
                     if (iu.lt.riutry(imolty,1)) goto 20
                  endif

                  growprev(1) = iu
                  
c                 --- replace this unit with the last unit in growlist
                  growlist(1,ivib) = growlist(1,invtry)

c                 --- reduce numgrow by 1
                  grownum(1) = grownum(1) - 1
c                 --- add this back into lexshed
                  lexshed(iu) = .true.
               else
c              --- we regrew all branches, no previous bead
                  growprev(1) = 0
               endif
            elseif ( invtry - count .ne. 1 ) then
c              --- problem in logic, should only be one nongrown bead
               write(iou,*) 'invtry,count',invtry,count
               write(iou,*) 'igrow,imolty',igrow,imolty
               call cleanup('logic problem in schedule')
            endif
         endif   

c     ---end the part that is specific for config
      elseif ( movetype .eq. 2 ) then
c     --- begin the part that is specific for swap
c        --- iutry is the first bead inserted - need to grow its neighbors
         do iu = 1,igrow
            lexshed(iu) = .true.
         enddo

         growfrom(1) = iutry
 
         invtry = invib(imolty,iutry)
         growprev(1) = 0
         if ( invtry .eq. 0 ) then
c           --- Bead iutry is the only bead to be grown ---
            index = 0
         else
c           --- grow all of the beads connected to bead iutry
            index = 1
            grownum(1) = invtry
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)
               if (lrigid(imolty)) iut = iutry - 1
c              --- grow these branches
               temp_count = temp_count + 1
               temp_store(temp_count) = iut
            enddo

            count = 0
            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               itry = int( dble(zz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(1,count) = iut
               lexshed(iut) = .false.
               
c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo
         endif

c        ---end the part that is specific for swap

      elseif ( movetype .eq. 3 ) then
         do iu = 1,igrow
            lexshed(iu) = .true.
         enddo
c        --- begin part that is specific for swatch 
         if ( iutry .eq. 0 ) then
c           --- no beads to be regrown via cbmc
            index = 0
            return
         endif
         growfrom(1) = iutry
         
         invtry = invib(imolty,iutry)
         growprev(1) = iprev

         if ( invtry .eq. 0 ) then
c           --- Bead 1 is the only bead to be grown ---
            index = 0
         else

c           --- grow all (except iprev) of the beads connected to bead 1
            index = 1
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)
               if ( iut .ne. iprev ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo

            count = 0
            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               itry = int( dble(zz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(1,count) = iut
               lexshed(iut) = .false.

c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo

            grownum(1) = count
         endif
c        --- end part that is specific for swatch
         
      elseif ( movetype .eq. 4 ) then
         do iu = 1,igrow
            lexshed(iu) = .true.
         enddo
c        --- begin part that is specific for rigid molecules

         if ( rindex(imolty) .eq. 0 ) then
c        --- entire molecule is rigid
            index = 0
         else
            index = rindex(imolty)
           
            do i = 1, rindex(imolty)
c        --- grow all non-rigid beads from the rigid part
               temp_count = 0
               invtry = invib(imolty,riutry(imolty,i)) - 1
               grownum(i) = invtry 
               growfrom(i) = riutry(imolty,i)

c        ---  for rigid molecules, have rigid vib last
               growprev(i)=ijvib(imolty,growfrom(i),invtry+1) 

               do ivib = 1, invtry
                  iut = ijvib(imolty,riutry(imolty,i),ivib)
c         --- grow these branches
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               enddo

               count = 0

               do zz = temp_count,1,-1
c                 --- choose grow bead randomly from temp_store
                  itry = int ( dble(zz) * random() ) + 1
                  iut = temp_store(itry)
                  count = count + 1
                  growlist(i,count) = iut
                  lexshed(iut) = .false.

c                 --- update temp_store for next iteration
                  temp_store(itry) = temp_store(zz)
               enddo
            enddo
         endif


c        --- end part that is specific for rigid molecules
      elseif (movetype .eq. 5) then
         if ( iutry .eq. 0 ) then
c           --- no beads to be regrown via cbmc
            return
         endif
         growfrom(index+1) = iutry
         
         invtry = invib(imolty,iutry)
         growprev(index+1) = iprev

         if ( invtry .eq. 0 ) then
c           --- Bead 1 is the only bead to be grown ---
         else

c           --- grow all (except iprev) of the beads connected to bead 1
            index = index + 1
            do ivib = 1,invtry
               iut = ijvib(imolty,iutry,ivib)
               if ( iut .ne. iprev ) then
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo

            count = 0
            do zz = temp_count,1,-1
c              --- choose grow bead randomly from temp_store
               itry = int( dble(zz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
               growlist(index,count) = iut
               lexshed(iut) = .false.
c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo

            grownum(index) = count
         endif

c * end iswatch stuff *
 
      else

c        --- non-existant move type
         write(iou,*) 'schedule movetype ',movetype
         call cleanup('non-valid move type')
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     from here on down config, swap, and swatch are the same


c     OLD METHOD - not fully random - removed 6-13-98

c      if ( .false. ) then 
c         write(iou,*) 'old method'
c      ibead = 0
c      index = 0
c 50   if ( index .lt. islen ) then
c         index = index + 1
c         do icbu = 1,grownum(index)
c            ibead = ibead + 1
c            iu = growlist(index,icbu)
c            invtry = invib(imolty,iu)
c            if ( invtry .ne. 1) then
c              --- this will be the next growing bead
c               islen = islen + 1
c               growprev(islen) = growfrom(index)
c               growfrom(islen) = iu
c               count = 0
c               do ivib = 1,invtry
c                  iut = ijvib(imolty,iu,ivib)
c                  if ( iut .ne. growprev(islen) ) then
c                     count = count + 1
c                     growlist(islen,count) = iut
c                     lexshed(iut) = .false.
c                  endif
c               enddo
c               grownum(islen) = count
c            endif
c         enddo
c         goto 50
c      endif
      
c     --- set up list of "outer" beads that have further growth sites
c     --- this method implemented 6-13-98

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
c ibead should be set here, otherwise it's not initialized below
      ibead=0

      if ( index .gt. 0 ) then
         ibead = 0
         outer_num = 0
         iufrom = growfrom(index)
         do icbu = 1,grownum(index)
c           --- increment counter ibead for total number of beads grown
            ibead = ibead + 1

c           --- determine whether this bead has any non-grown neighbors
            iu = growlist(index,icbu)
            invtry = invib(imolty,iu)
         
            if ( invtry .gt. 1 ) then
c              --- add one to the stack of outer_sites
               outer_num = outer_num + 1
               outer_sites(outer_num) = iu
               outer_prev(outer_num) = iufrom
            endif
         enddo

c        --- begin while loop to grow all outer beads until done
 70      if ( outer_num .gt. 0 ) then

c           --- choose one site randomly from the stack
            outer_try = int( dble(outer_num) * random() ) + 1

c           --- increment index to show this is the next growfrom
            index = index + 1

            iu = outer_sites(outer_try)
            iufrom = outer_prev(outer_try)
c           --- assign growfrom and growprev for this index
            growfrom(index) = iu
            growprev(index) = iufrom

            invtry = invib(imolty,iu)
c           --- assign the grow beads in random order
            temp_count = 0
            do ivib = 1,invtry
               iut = ijvib(imolty,iu,ivib)
               if ( iut .ne. iufrom ) then
c                 --- add to the list of beads to be grown from iu
                  temp_count = temp_count + 1
                  temp_store(temp_count) = iut
               endif
            enddo

            count = 0
            do zz = temp_count, 1, -1
               itry = int ( dble(zz) * random() ) + 1
               iut = temp_store(itry)
               count = count + 1
c              --- assign growlist for current index and count
               growlist(index,count) = iut
               lexshed(iut) = .false.

c              --- update temp_store for next iteration
               temp_store(itry) = temp_store(zz)
            enddo
c           --- assign grownum for this index
            grownum(index) = count

c           --- update list of "outer" beads
c           --- remove bead that was just grown from outer list
            outer_sites(outer_try) = outer_sites(outer_num)
            outer_prev(outer_try) = outer_prev(outer_num)
            outer_num = outer_num - 1
         
c           --- add the new beads if they have more to be grown
            iufrom = iu
            do icbu = 1,grownum(index)
c              --- increment counter ibead for total number of beads grown
               ibead = ibead + 1

c              --- determine whether this bead has any non-grown neighbors
               iu = growlist(index,icbu)
               invtry = invib(imolty,iu)

               if ( invtry .gt. 1 ) then
c                 --- add one to the stack of outer_sites
                  outer_num = outer_num + 1
                  outer_sites(outer_num) = iu
                  outer_prev(outer_num) = iufrom
               endif
            enddo

c           --- end of while loop 70
            goto 70
         endif
      endif

 100  if ( (movetype .eq. 1) .and. (ibead .gt. nmaxcbmc(imolty)) ) then
         kickout = kickout + 1
         goto 11
      endif

      if ( lprint ) then
         write(iou,*) 'movetype',movetype
         write(iou,*) 'index',index
         do ibead = 1,index
            write(iou,*) 'ibead',ibead
            write(iou,*) 'growfrom(ibead)',growfrom(ibead)
            write(iou,*) 'growprev(ibead)',growprev(ibead)
            write(iou,*) 'grownum(ibead)',grownum(ibead)
            do count = 1,grownum(ibead)
               write(iou,*) 'count,growlist(ibead,count)',count
     &              ,growlist(ibead,count)
            enddo
         enddo
         do iu = 1,igrow
            write(iou,*) 'iu,lexshed(iu)',iu,lexshed(iu)
         enddo
      endif

      
      if (lrig(imolty).and.nrig(imolty).eq.0) then
         ib2 = index - ju
         if (ib2.lt.1) then
            llrig = .false.
            return
         else

            do ibead = 1, index
               lsave(ibead) = .false.
            enddo

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
            enddo
            

         endif
         
c     --- cycle through the rest of the rigid beads
         
         do ibead = ib2 + 1, index
            iufrom = growfrom(ibead)

            if (lfind(iufrom)) then
               
               lsave(ibead) = .true.
               
               do count = 1, grownum(ibead)
                  
                  iu = growlist(ibead,count)

                  lfind(iu) = .true.

               enddo

            endif
         enddo
      else
         llrig = .false.
      endif

      


c      write(iou,*) 'end SCHEDULE'
  
      return
      end
