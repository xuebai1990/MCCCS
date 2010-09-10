      subroutine molsetup(imolty)

c molsetup
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

      integer::i,j,k,n,iunit,imolty,dum,iu,countbend,counttor
     &     ,atype,btype,ctype,dtype,tortype
     &     ,ibend,bendtype,ntor,itor,ju,ku,nu,type,sitelook
     &     ,vibtype,countvib

c     ************************************************************


      call vibcheck(1,dum,dum,dum)
      call bendcheck(1,dum,dum,dum,dum)
      call torcheck(1,dum,dum,dum,dum,dum)


      iunit = nunit(imolty)
      
      masst(imolty) = 0.0d0

      do i = 1, iunit

         read(4,*)
         read(4,*) iu,ntype(imolty,i)
         
         masst(imolty) = masst(imolty) + mass(ntype(imolty,i))

         write(iou,*) 'bead ',iu,' beadtype ',ntype(imolty,i)

         read(4,*) 
         read(4,*) invib(imolty,i),(ijvib(imolty,i,j)
     &        ,j=1,invib(imolty,i))

         if (invib(imolty,i).gt.6) then

            write(iou,*) 'imolty',imolty,'   i',i,'   invib'
     &           ,invib(imolty,i)
            stop 'too many vibrations'
         endif

      enddo
      
      do iu = 1, iunit
         
         countbend = 0
         counttor = 0
         countvib = 0

         atype = ntype(imolty,iu)

         do j = 1, invib(imolty,iu)

            ju = ijvib(imolty,iu,j)

            btype = ntype(imolty,ju)

            call vibcheck(2,atype,btype,vibtype)
            
            if (vibtype.eq.0) then
               write(iou,*) 'atype,btype',atype,btype
               stop 'screwup in vibrations'
            endif
                        
            countvib = countvib + 1

            ijvib(imolty,iu,countvib) = ju
            itvib(imolty,iu,countvib) = vibtype

            do k = 1, invib(imolty,ju)

               ku = ijvib(imolty,ju,k)

               if (ku.ne.iu) then

                  ctype = ntype(imolty,ku)

                  call bendcheck(2,atype,btype,ctype,bendtype)

                  if (bendtype.eq.0) then
                     write(iou,*) 'atype,btype,ctype',atype
     &                    ,btype,ctype
                     stop 'screwup in bending angles'
                  endif

                  countbend = countbend + 1

                  ijben2(imolty,iu,countbend) = ju
                  ijben3(imolty,iu,countbend) = ku
                  itben(imolty,iu,countbend) = bendtype

                  do n = 1, invib(imolty,ku)
                     nu = ijvib(imolty,ku,n)
                     if (nu.ne.ju) then
                        
                        dtype = ntype(imolty,nu)

                        call torcheck(2,atype,btype,ctype,dtype
     &                       ,tortype)

                        if (tortype.eq.0) then
                           write(iou,*) 'atype,btype,ctype,dtype',atype
     &                          ,btype,ctype,dtype
                           stop 'screwup in torsion angles'
                        endif
                        
                        counttor = counttor + 1

                        ijtor2(imolty,iu,counttor) = ju
                        ijtor3(imolty,iu,counttor) = ku
                        ijtor4(imolty,iu,counttor) = nu
                        ittor(imolty,iu,counttor) = tortype

                     endif
                                       
                  enddo
               
               endif   

            enddo
            
         enddo
      
         inben(imolty,iu) = countbend
         intor(imolty,iu) = counttor
         

      enddo


      return
      end



c     *************************************************************

      subroutine vibcheck(iinit,atype,btype,vibtype)

      logical::lfinda,lfindb,lfound

      integer::iinit,atype,btype,vibtype,nsite,isite,vbtype
     &     ,ntvib,n,i,ia,ib

      dimension nsite(20,2),isite(20,2,7),vbtype(20)

      save vbtype,nsite,isite,ntvib
      
      if (iinit.eq.1) then

         read(60,*)
         read(60,*) ntvib

         do n = 1, ntvib

            read(60,*)
            read(60,*) vbtype(n)
            read(60,*) nsite(n,1),(isite(n,1,i),i=1,nsite(n,1))
            read(60,*) nsite(n,2),(isite(n,2,i),i=1,nsite(n,2))
            
         enddo

      else

         vibtype = 0
         lfound = .false.

         do n = 1, ntvib

            lfinda = .false.
            lfindb = .false.

            do i = 1, nsite(n,1)
               if (atype.eq.isite(n,1,i)) then
                  lfinda = .true.
                  ia = i
                  goto 105
               elseif (btype.eq.isite(n,1,i)) then
                  lfindb = .true.
                  ia = i
                  goto 105
               endif
            enddo

 105        continue

            do i = 1, nsite(n,2)
               if (lfinda) then
                  if (btype.eq.isite(n,2,i)) then
                     lfindb = .true.
                     ib = i
                     goto 110
                  endif
               elseif (lfindb) then
                  if (atype.eq.isite(n,2,i)) then
                     lfinda = .true.
                     ib = i
                     goto 110
                  endif
               endif
            enddo
 110        continue
            
            if (lfinda.and.lfindb) then
               if (lfound) then
                  stop 'vibration type not distinguishable'
               endif
               vibtype = vbtype(n)
               lfound = .true.
            endif
         enddo

      endif

      return
      end

c     **************************************************************

      subroutine bendcheck(iinit,atype,btype,ctype,bendtype)

      logical::lfinda,lfindb,lfindc,lfound

      integer::iinit,atype,btype,ctype,bendtype,nsite,isite
     &     ,bntype,ntbend,n,i,ia,ib,ic

      dimension nsite(20,3),isite(20,3,7),bntype(20)
      
      save bntype,nsite,isite,ntbend
      
      
      if (iinit.eq.1) then

         read(61,*) 
         read(61,*) ntbend

         do n = 1, ntbend

            read(61,*)
            read(61,*) bntype(n)
            read(61,*) nsite(n,1),(isite(n,1,i),i=1,nsite(n,1))
            read(61,*) nsite(n,2),(isite(n,2,i),i=1,nsite(n,2))
            read(61,*) nsite(n,3),(isite(n,3,i),i=1,nsite(n,3))

         enddo

      else
         bendtype = 0
         lfound = .false.
         do n = 1, ntbend
                        
            lfinda = .false.
            lfindb = .false.
            lfindc = .false.
            
            do i = 1, nsite(n,1)
               if (atype.eq.isite(n,1,i)) then
                  lfinda = .true.
                  ia = i
                  goto 105
               elseif (ctype.eq.isite(n,1,i)) then
                  lfindc = .true.
                  ia = i
                  goto 105
               endif
            enddo

 105        continue

            do i = 1, nsite(n,2)
               if (btype.eq.isite(n,2,i)) then
                  lfindb = .true.
                  ib = i
                  goto 110
               endif

            enddo

 110        continue

            do i = 1, nsite(n,3)
               if (lfinda) then
                  if (ctype.eq.isite(n,3,i)) then
                     lfindc = .true.
                     ic = i
                     goto 115
                  endif
               elseif (lfindc) then
                  if (atype.eq.isite(n,3,i)) then
                     lfinda = .true.
                     ic = i
                     goto 115
                  endif
               endif
            enddo

 115        continue

            if (lfinda.and.lfindb.and.lfindc) then
               if (lfound) then
                  stop 'bend type not distinguishable'
               endif
               bendtype = bntype(n)
               lfound = .true.

            endif

            

            
         enddo
      endif

      return
      end
      
c     *************************************************************


      subroutine torcheck(iinit,atype,btype,ctype,dtype,tortype)


      logical::lfinda,lfindb,lfindc,lfindd,lfound,lrev

      integer::iinit,atype,btype,ctype,dtype,tortype,isite,nsite
     &     ,n,i,ia,ib,ic,id,trtype,nttor,ir


      dimension trtype(35),nsite(35,7),isite(35,7,7)

      save trtype,nsite,isite,nttor

      if (iinit.eq.1) then
         read(62,*) 
         read(62,*) nttor
         do n = 1, nttor

            read(62,*)
            read(62,*) trtype(n)
            read(62,*) nsite(n,1),(isite(n,1,i),i=1,nsite(n,1))
            read(62,*) nsite(n,2),(isite(n,2,i),i=1,nsite(n,2))
            read(62,*) nsite(n,3),(isite(n,3,i),i=1,nsite(n,3))
            read(62,*) nsite(n,4),(isite(n,4,i),i=1,nsite(n,4))

         enddo
      else
         tortype = 0
         lfound = .false.

         do n = 1,nttor

            lfinda = .false.
            lfindb = .false.
            lfindc = .false.
            lfindd = .false.
            lrev = .false.

            do i = 1, nsite(n,1)
               if (atype.eq.isite(n,1,i)) then
                  lfinda = .true.
                  ia = i
                  goto 104
               elseif (dtype.eq.isite(n,1,i)) then
                  lfindd = .true.
                  ia = i
                  goto 104
               endif
            enddo

 104        continue

            do i = 1, nsite(n,4)
               if (lfinda) then
                  if (dtype.eq.isite(n,4,i)) then
                     lrev = .true.
                     ir = i
                     goto 105
                  endif
               elseif (lfindd) then
                  if (atype.eq.isite(n,4,i)) then
                     lrev = .true.
                     ir = i
                     goto 105
                  endif
               endif
            enddo

 105        continue
            
            do i = 1, nsite(n,2)
               if (lfinda) then
                  if (btype.eq.isite(n,2,i)) then
                     lfindb = .true.
                     ib = i
                     goto 110
                  endif
               elseif (lfindd) then
                  if (ctype.eq.isite(n,2,i)) then
                     lfindc = .true.
                     ib = i
                     goto 110
                  endif
               endif
               
            enddo
            
 110        continue
            
            do i = 1, nsite(n,3)
               if (lfindb) then
                  if (ctype.eq.isite(n,3,i)) then
                     lfindc = .true.
                     ic = i
                     goto 115
                  endif
               elseif (lfindc) then
                  if (btype.eq.isite(n,3,i)) then
                     lfindb = .true.
                     ic = i
                     goto 115
                  endif
               endif
            enddo
            
 115        continue
           
            if (.not.lfindb.or..not.lfindc) goto 120

            do i = 1, nsite(n,4)
               if (lfinda) then
                  if (dtype.eq.isite(n,4,i)) then
                    lfindd = .true.
                    id = i
                    goto 120
                 endif
              elseif (lfindd) then
                 if (atype.eq.isite(n,4,i)) then
                    lfinda = .true.
                    id = i
                    goto 120
                 endif
              endif
           enddo
           
 120       continue
 
           if (lfinda.and.lfindb.and.lfindc.and.lfindd) then
              if (lfound) then
                 write(iou,*) 'a,b,c,d',atype,btype,ctype,dtype
                 stop 'torsion type not distinguishable'
              endif
              
              tortype = trtype(n)
              lfound = .true.
              
           elseif (lrev) then
              lrev = .false.
              if (lfinda) then
                 lfinda = .false.
                 lfindd = .true.
                 id = ir
              elseif (lfindd) then
                 lfindd = .false.
                 lfinda = .true.
                 ia = ir
              endif
              goto 105
           endif
            
        enddo
                
      endif
      
      return
      end
      






