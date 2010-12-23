      subroutine molsetup(imolty)

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
!$$$      include 'connect.inc'

      integer(KIND=normal_int)::i,j,k,n,iunit,imolty,dum,iu,countbend
     & ,counttor,atype,btype,ctype,dtype,tortype,ibend,bendtype,ntor
     & ,itor,ju,ku,nu,type,sitelook,vibtype,countvib

!     ************************************************************


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
            call cleanup('too many vibrations')
         end if

      end do
      
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
               call cleanup('screwup in vibrations')
            end if
                        
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
                     call cleanup('screwup in bending angles')
                  end if

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
                           call cleanup('screwup in torsion angles')
                        end if
                        
                        counttor = counttor + 1

                        ijtor2(imolty,iu,counttor) = ju
                        ijtor3(imolty,iu,counttor) = ku
                        ijtor4(imolty,iu,counttor) = nu
                        ittor(imolty,iu,counttor) = tortype

                     end if
                                       
                  end do
               
               end if   

            end do
            
         end do
      
         inben(imolty,iu) = countbend
         intor(imolty,iu) = counttor
         

      end do


      return
      end



!     *************************************************************

      subroutine vibcheck(iinit,atype,btype,vibtype)

      use var_type
      implicit none

      logical::lfinda,lfindb,lfound

      integer(KIND=normal_int)::iinit,atype,btype,vibtype,nsite,isite
     & ,vbtype,ntvib,n,i,ia,ib

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
            
         end do

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
               end if
            end do

 105        continue

            do i = 1, nsite(n,2)
               if (lfinda) then
                  if (btype.eq.isite(n,2,i)) then
                     lfindb = .true.
                     ib = i
                     goto 110
                  end if
               elseif (lfindb) then
                  if (atype.eq.isite(n,2,i)) then
                     lfinda = .true.
                     ib = i
                     goto 110
                  end if
               end if
            end do
 110        continue
            
            if (lfinda.and.lfindb) then
               if (lfound) then
                  call cleanup('vibration type not distinguishable')
               end if
               vibtype = vbtype(n)
               lfound = .true.
            end if
         end do

      end if

      return
      end

!     **************************************************************

      subroutine bendcheck(iinit,atype,btype,ctype,bendtype)

      use var_type
      implicit none

      logical::lfinda,lfindb,lfindc,lfound

      integer(KIND=normal_int)::iinit,atype,btype,ctype,bendtype,nsite
     & ,isite,bntype,ntbend,n,i,ia,ib,ic

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

         end do

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
               end if
            end do

 105        continue

            do i = 1, nsite(n,2)
               if (btype.eq.isite(n,2,i)) then
                  lfindb = .true.
                  ib = i
                  goto 110
               end if

            end do

 110        continue

            do i = 1, nsite(n,3)
               if (lfinda) then
                  if (ctype.eq.isite(n,3,i)) then
                     lfindc = .true.
                     ic = i
                     goto 115
                  end if
               elseif (lfindc) then
                  if (atype.eq.isite(n,3,i)) then
                     lfinda = .true.
                     ic = i
                     goto 115
                  end if
               end if
            end do

 115        continue

            if (lfinda.and.lfindb.and.lfindc) then
               if (lfound) then
                  call cleanup('bend type not distinguishable')
               end if
               bendtype = bntype(n)
               lfound = .true.

            end if

            

            
         end do
      end if

      return
      end
      
!     *************************************************************


      subroutine torcheck(iinit,atype,btype,ctype,dtype,tortype)

      use global_data,only:iou
      use var_type
      implicit none

      logical::lfinda,lfindb,lfindc,lfindd,lfound,lrev

      integer(KIND=normal_int)::iinit,atype,btype,ctype,dtype,tortype
     & ,isite,nsite,n,i,ia,ib,ic,id,trtype,nttor,ir


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

         end do
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
               end if
            end do

 104        continue

            do i = 1, nsite(n,4)
               if (lfinda) then
                  if (dtype.eq.isite(n,4,i)) then
                     lrev = .true.
                     ir = i
                     goto 105
                  end if
               elseif (lfindd) then
                  if (atype.eq.isite(n,4,i)) then
                     lrev = .true.
                     ir = i
                     goto 105
                  end if
               end if
            end do

 105        continue
            
            do i = 1, nsite(n,2)
               if (lfinda) then
                  if (btype.eq.isite(n,2,i)) then
                     lfindb = .true.
                     ib = i
                     goto 110
                  end if
               elseif (lfindd) then
                  if (ctype.eq.isite(n,2,i)) then
                     lfindc = .true.
                     ib = i
                     goto 110
                  end if
               end if
               
            end do
            
 110        continue
            
            do i = 1, nsite(n,3)
               if (lfindb) then
                  if (ctype.eq.isite(n,3,i)) then
                     lfindc = .true.
                     ic = i
                     goto 115
                  end if
               elseif (lfindc) then
                  if (btype.eq.isite(n,3,i)) then
                     lfindb = .true.
                     ic = i
                     goto 115
                  end if
               end if
            end do
            
 115        continue
           
            if (.not.lfindb.or..not.lfindc) goto 120

            do i = 1, nsite(n,4)
               if (lfinda) then
                  if (dtype.eq.isite(n,4,i)) then
                    lfindd = .true.
                    id = i
                    goto 120
                 end if
              elseif (lfindd) then
                 if (atype.eq.isite(n,4,i)) then
                    lfinda = .true.
                    id = i
                    goto 120
                 end if
              end if
           end do
           
 120       continue
 
           if (lfinda.and.lfindb.and.lfindc.and.lfindd) then
              if (lfound) then
                 write(iou,*) 'a,b,c,d',atype,btype,ctype,dtype
                 call cleanup('torsion type not distinguishable')
              end if
              
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
              end if
              goto 105
           end if
            
        end do
                
      end if
      
      return
      end
      






