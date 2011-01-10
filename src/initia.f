      subroutine initia(file_struct)

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_files
      implicit none
!$$$      include 'mpi.inc'
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'system.inc'
!$$$      include 'poten.inc'
!$$$      include 'conver.inc'
!$$$      include 'inpar.inc'
!$$$      include 'ensemble.inc' 
!$$$      include 'connect.inc'
!$$$      include 'cbmc.inc'

      character(len=*),intent(in)::file_struct
      logical::lhere(nmax)
      integer(KIND=normal_int)::io_struct,jerr,i,j,m,m1,m2,n,nn,ic,jc,kc,it,ip1,ip2,ip3 ,ii,jj,iivib,jjben,jjtor,intemp,imol,ibtype,imolty,ibuild ,rand_id,offset,count_chain

      integer(KIND=normal_int)::iboxst,iboxed,ibox,check,chktot
      integer(KIND=normal_int)::mcmt(ntmax,nbxmax),pcmt(ntmax) ,mcmtma(ntmax,nbxmax)
      integer(KIND=normal_int)::unitc,ntii

      integer(KIND=normal_int)::bmap(numax),imap(numax),zzz,prev,ifrom ,nsave
      logical::lacc(numax),lgrow,lterm,lgrown(ntmax)

      real(KIND=double_precision)::xtemp(numax),ytemp(numax) ,ztemp(numax)
      real(KIND=double_precision)::ddum

      real(KIND=double_precision)::ux,uy,uz,xshift,dic ,xnext,ynext ,znext,xynext,angold,angnew,rot ,x1,y1,z1,d1,x2,y2,z2,d2,bang ,blen,random

      real(KIND=double_precision)::rxui,ryui,rzui,xaa1,yaa1,zaa1,daa1 ,xa1a2,ya1a2,za1a2,da1a2 

      real(KIND=double_precision)::xcc,ycc,zcc,tcc,spltor

      real(KIND=double_precision)::vtorso,vbend,vtg,thetac,theta,dot ,aben,ator

      real(KIND=double_precision)::rxvec(numax,numax),ryvec(numax,numax) ,rzvec(numax,numax),distanceij(numax,numax)

      real(KIND=double_precision)::samx(ntmax,numax),samy(ntmax,numax) ,samz(ntmax,numax)

      real(KIND=double_precision)::vdummy

      dimension ux(nbxmax),uy(nbxmax),uz(nbxmax)
      dimension check(ntmax)

! --------------------------------------------------------------------
      if (myid.eq.0) then
         write(iou,*) 
         write(iou,*) 'subroutine initia'
         write(iou,*) 
      end if

!     --- initialize nchbox ---
      do i=1,nbox
         nchbox(i) = 0
      end do

      iboxst = 1
      iboxed = nbox

         chktot = 0

         do i = 1,nmolty
            check(i) = 0
            do j = 1, nbox
               nchbox(j) = nchbox(j) + ininch(i,j)
               check(i)= check(i) + ininch(i,j) 
            end do

            chktot = chktot + check(i)
         end do

         if ( chktot .ne. nchain ) then
            write(iou,*) 'inconsistant number of chains in INITIA'
            do j = 1,nbox
               write(iou,*) 'ininch',j,(ininch(i,j),i=1,nmolty)
            end do
            write(iou,*) 'nchain',nchain
            call cleanup('')
         end if
         
         do i = 1, nmolty
            if ( temtyp(i) .ne. check(i) ) then
               write(iou,*) 'inconsistant number of chains in INITIA'
               write(iou,*) 'moltyp',i,(ininch(i,j),j=1,nbox)
               write(iou,*) 'temtyp:',temtyp(i)
               call cleanup('')
            end if            
         end do

         do i = iboxst,iboxed
            unitc = inix(i)*iniy(i)*iniz(i)
            if ( nchbox(i) .gt. unitc ) then
               write(iou,*) 'unit cell too small in box',i
               call cleanup('')
            end if
         end do
         
         
! -----------------------------------------------------------------------------
 
! *** calculation of unit cell dimensions ***
      do i = 1,nbox
         ux(i) = boxlx(i) / dble(inix(i)) 
         uy(i) = boxly(i) / dble(iniy(i))
         uz(i) = boxlz(i) / dble(iniz(i))
         if (myid.eq.0) then
            write(iou,*) 'box',i
            write(iou,*) 'ini',inix(i),iniy(i),iniz(i)
            write(iou,*) 'box',boxlx(i),boxly(i),boxlz(i)
            write(iou,*) 'uni',ux(i),uy(i),uz(i)
         end if
      end do

! - count number of molecules of each type -
      do i = 1,nmolty
         mcmt(i,1) = 0
      end do

      do i=1,nchain
         imolty = moltyp(i)
         mcmt(imolty,1) = mcmt(imolty,1) + 1
         lhere(i) = .false.
      end do

      do i = 1,nmolty
         if ( mcmt(i,1) .ne. check(i) ) then
            write(iou,*) 'inconsistant number of type in INITIA'
            write(iou,*) 'mcmt(i,total),check(i)',mcmt(i,1),check(i)
            call cleanup('')
         end if
      end do

      do i = 1,nmolty
         do j = 1, nbox
            mcmt(i,j) = ininch(i,j)
            mcmtma(i,j) = 0
         end do
      end do
      do ibox = 1,nbox
         mcmtma(1,ibox) = mcmt(1,ibox)
         do i = 2, nmolty
            mcmtma(i,ibox) = mcmtma(i-1,ibox) + mcmt(i,ibox)
         end do
      end do

      if (myid.eq.0) then
         write(iou,*) 'nmolty',nmolty
         write(iou,*) '   mcmt',((mcmt(i,ibox),i=1,nmolty),ibox=1,nbox)
      end if

! *****************************
! *** calculate coordinates ***

!     read sample structure from unit 78 -
      io_struct=get_iounit()
      open(unit=io_struct,access='sequential',action='readwrite',file=file_struct,form='formatted',iostat=jerr,status='unknown')
      if (jerr.ne.0) then
         call cleanup('cannot open input file')
      end if

      do i = 1, nmolty         
         lgrown(i) = .false.
         if ( lbranch(i) ) then
            read(io_struct,*)
            do m = 1, nunit(i)
               read(io_struct,*) samx(i,m), samy(i,m), samz(i,m)
            end do
         else if ( .not. lbranch(i)) then
! * if lbranch is false but the molecule is not linear attempt
! * to grow it with cbmc
            lgrow = .false.
            do m = 1,nunit(i)
               if (invib(i,m) .gt. 2) then
                  lgrow = .true.
               end if
            end do
            if (lgrow) then

               if (myid.eq.0) then
                  write(iou,*) 'growing a sample structure with CBMC'
               end if

               if (nunit(i) .ne. nugrow(i)) then
                  write(iou,*) 'Cant grow molecule.  Please', ' provide a structure via fort.78'
                  call cleanup('')
               end if
! * put the first bead at the origin
               rxnew(1) = 0.0d0
               rynew(1) = 0.0d0
               rznew(1) = 0.0d0

! * determine the growth schedule
               call schedule(nugrow(i),i,ifrom,1,0,2)

! * actually grow the structure
               nsave = nchain

               moltyp(1) = i
               do m = 1,nunit(i)
                  qqu(1,m) = qelect(ntype(i,m))
               end do

               nchain = 1

               call rosenbluth( .true.,lterm,1,1,i,ifrom ,1,nugrow(i),ddum,.false.,ddum,2 )

               if (lterm) then
                  write(iou,*) 'error in initia growing molecule'
                  write(iou,*) 'maybe increasing nchoi would help?'
                  call cleanup('')
               end if

! * return the value of nchain
               nchain = nsave
! * assign the coordinates
               do m = 1,nunit(i)
                  samx(i,m) = rxnew(m)
                  samy(i,m) = rynew(m)
                  samz(i,m) = rznew(m)
               end do

               lgrown(i) = .true.

            end if
         end if
      end do

      close(io_struct)

      do i = 1,nmolty
         if (lgrown(i)) then
            lbranch(i) = .true.
         end if
      end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     - inimix = 0 : take molecules at random
!     - inimix > 0 : take molecules in order (first type I etc.)
!     - inimix < 0 : take molecules in alternating order
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

      count_chain = 0 
      offset = 0

      do ibox = iboxst,iboxed
         if(nmolty .gt. 1) then
           if(inimix(ibox).eq.0) then
              if (ibox.eq.1) then
                 offset = 0
              else
                 offset = offset+nchbox(ibox-1)
              end if
 18           rand_id = idint(dble(nchbox(ibox))*random())+ 1 + offset
              if (.not.lhere(rand_id)) then
                 count_chain = count_chain + 1
                 lhere(rand_id) = .true.
                 do imolty = 1,nmolty  
                    if (count_chain.le.(mcmtma(imolty,ibox)+offset)) then
                       moltyp(rand_id) = imolty     
                       goto 20 
                    end if 
                 end do
 20              continue
!             write(iou,*) count_chain, rand_id,moltyp(rand_id) 
              else
                 goto 18
              end if
              if (count_chain.lt.(nchbox(ibox)+offset)) then
                 goto 18
              end if 
           end if
         end if
      end do   


      nn = 0

      do 102 ibox = iboxst,iboxed
         do imol = 1,nmolty
            pcmt(imol) = 0
         end do

         n = 0

         do 101 kc = 0, iniz(ibox)-1
            if ( mod(kc,2) .eq. 0) then
               xshift = 0.0d0
            else
               xshift = dshift(ibox)
            end if
            
            do 100 ic = 0, inix(ibox)-1
               
               do 99 jc = 0, iniy(ibox)-1

                  if ( mod(jc,2) .eq. 0) then
                     dic = 0.0d0
                  else
                     dic = 0.5d0
                  end if

                  n=n+1

                  if (n .le. nchbox(ibox) ) then
                     nn=nn+1
                     rxu(nn,1) = ( dble(ic) + dic )*ux(ibox)+ xshift
                     ryu(nn,1) = dble(jc) * uy(ibox) 
                     rzu(nn,1) = dble(kc) * uz(ibox) + zshift(ibox)
                     nboxi(nn) = ibox
                  else
                     goto 102
                  end if
                  
!                  write(iou,*) 'nn',nn
!                  write(iou,*) 'ic',ic,'   jc',jc,'   kc',kc
                  
!     - inimix > 0 : take molecules in order (first type I etc.)
!     - inimix < 0 : take molecules in alternating order

                  if ( nmolty .gt. 1 ) then
                     if ( inimix(ibox) .gt. 0 ) then
                        do imol = 1, nmolty
                           if ( n .le. mcmtma(imol,ibox) ) then
                              intemp = imol
                              exit
                           end if
                        end do
                     elseif ( inimix(ibox) .lt. 0 ) then
!     This is not right. It doesn't consider the number of molecules of each type in the box.
!                        do imol = 1, nmolty
!                           nt = n - imol
!                           if ( mod( nt, nmolty ) .eq. 0 ) then
!                              intemp = imol
!                           end if
!                        end do
                        intemp = mod(n,nmolty)
                     end if
                  else
                     intemp = 1
                  end if

                  if (inimix(ibox).eq.0) then
                      intemp = moltyp(nn)
                  else
                       moltyp(nn) = intemp
                  end if 

                  ncmt(ibox,intemp) = ncmt(ibox,intemp) + 1
                  
!                  write(iou,*) 'intemp', intemp     

                  if ( lbranch(intemp) ) then
                     ibuild = nunit(intemp)
                  else
                     ibuild = nugrow(intemp)
                  end if

                  if ( .not. lbranch(intemp)) then
! *************************************************
! *** start determination of linear chain positions
! *** allowing for numbering out of order
! *** note: doesn't exactly create equilibrium structure with
! *** respect to bond angles or torsions, but that will shake out
! *** with CBMC anyways.  Should at least take away the overlaps of
! *** the previous method
! *************************************************
! * first need to determine re-mapped bead order- search through connectivity
!
!   call the results map(i) where i=1 is one chain end, and its
!   value is equal to the bead number of that end
! 
!   for example, methanol oxygen first, then hydrogen, then CH3
!   
!                        H---O--CH3
!
!        bead numbers:   2 - 1 - 3
!
!        bmap(1) = 2
!        bmap(2) = 1
!        bmap(3) = 3
!
!        the inverse map is just the opposite:
!
!        imap(1) = 2
!        imap(2) = 1
!        imap(3) = 3
!
! * initialize accounted for variable
                     do m = 1,ibuild
                        lacc(m) = .false.
                     end do

! * first find the end with the lowest number
                     zzz = ibuild
                     do m = 1,ibuild
                        if (invib(intemp,m) .le. 1) then
                           if (m .le. zzz) then
                              zzz = m
                           end if
                        elseif (invib(intemp,m) .gt. 2) then
                           write(iou,*) 'initia only works for linear', ' molecules!  Maybe you should make', ' a fort.78 file and use lbranch?'
                           call cleanup('')
                        end if
                     end do

                     bmap(1) = zzz
                     imap(zzz) = 1
                     lacc(zzz) = .true.

! * now determine the rest
                     do m = 2,ibuild
                        prev = bmap(m-1)
                        do zzz = 1,ibuild
                           if (ijvib(intemp,zzz,1) .eq. prev .or. ijvib(intemp,zzz,2) .eq. prev) then
                              if (.not. lacc(zzz)) then
                                 bmap(m) = zzz
                                 imap(zzz) = m
                                 lacc(zzz) = .true.
                              end if
                           end if
                        end do
                     end do

! * now use old method with re-mapped numbers:
! * put first end at origin:
                     xtemp(1) = 0.0d0
                     ytemp(1) = 0.0d0
                     ztemp(1) = 0.0d0

! * now we need to loop over all the other beads:
                     do m = 2, ibuild
                     
                        m1 = m - 1
                        m2 = m - 2

                        if ( inirot(ibox) .eq. 0 ) then
                           rot = random() * 360.0d0 * degrad
                        elseif ( inirot(ibox) .gt. 0 ) then
                           rot = dble(inirot(ibox)) * degrad
                        else
                           if ( mod(jc,2) .eq. 0 ) then
                              rot = dble(inirot(ibox)) * degrad
                           else
                              rot = -(dble(inirot(ibox)) * degrad)
                           end if
                        end if

                        if ( inben(intemp,bmap(m1)) .gt. 0 ) then
                           ibtype = itben(intemp,bmap(m1),1)
                           angold = brben(ibtype) / 2.0d0

                           if ( m .eq. 2 ) then
                              ibtype = itben(intemp,bmap(m1),1)
                              angnew = brben(ibtype) - angold
                           else
                              ibtype = itben(intemp,bmap(m2),1)
                              angnew = brben(ibtype) - angold
                           end if 
!     write(iou,*) 'angold',angold*raddeg,
!     +                       '   angnew',angnew*raddeg
                           angold = angnew
                           
! * need to search for proper bond length
                           do zzz = 1,invib(intemp,bmap(m))
                              if (ijvib(intemp,bmap(m),zzz)  .eq. bmap(m1)) then
                                 ibtype = itvib(intemp,bmap(m),zzz)
                              end if
                           end do

                           ztemp(m) = dsin(angnew) * brvib(ibtype)
                           xynext = dcos(angnew) * brvib(ibtype)
!                           write(iou,*) 'znext',znext,'   xynext',xynext
                        else
! * need to search for proper bond length
                           do zzz = 1,invib(intemp,bmap(m))
                              if (ijvib(intemp,bmap(m),zzz)  .eq. bmap(m1)) then
                                 ibtype = itvib(intemp,bmap(m),zzz)
                              end if
                           end do

!                           ztemp(m) = dsin(angnew) * brvib(ibtype)
!                           xynext = dcos(angnew) * brvib(ibtype)

                           ztemp(m) = brvib(ibtype)
                           xynext = 0.0d0
                        end if
                        
                        if ( mod(m,2) .eq. 0 ) then
                           xtemp(m) = dcos(rot) * xynext
                           ytemp(m) = dsin(rot) * xynext
                        else
                           xtemp(m) = -(dcos(rot) * xynext)
                           ytemp(m) = -(dsin(rot) * xynext)
                        end if

                        xtemp(m) = xtemp(m1) + xtemp(m)
                        ytemp(m) = ytemp(m1) + ytemp(m)
                        ztemp(m) = ztemp(m1) + ztemp(m)


                     end do

! * translate so that first bead number is at origin
                     do m = 1,ibuild
                        if (m .ne. imap(1)) then
                          xtemp(m) = xtemp(m) - xtemp(imap(1)) 
                          ytemp(m) = ytemp(m) - ytemp(imap(1)) 
                          ztemp(m) = ztemp(m) - ztemp(imap(1)) 
                        end if
                     end do
                     xtemp(imap(1)) = 0.0d0
                     ytemp(imap(1)) = 0.0d0
                     ztemp(imap(1)) = 0.0d0

                  end if
! *** end linear determination
! ****************************
!                  write(iou,*) 'ibuild',ibuild
                  do 98 m = 2, ibuild
                     
                     m1 = m - 1
                     m2 = m - 2
!                     write(iou,*) 'intemp',intemp 
                     if ( lbranch(intemp) ) then
!     - branched molecule with sample structure -
                        xnext = samx(intemp,m) -samx(intemp,m1)
                        ynext = samy(intemp,m) -samy(intemp,m1)
                        znext = samz(intemp,m) -samz(intemp,m1)
                     else
! * linear molecule determined above- replacing old code that follows.
                        xnext = xtemp(bmap(m)) - xtemp(bmap(m1))
                        ynext = ytemp(bmap(m)) - ytemp(bmap(m1))
                        znext = ztemp(bmap(m)) - ztemp(bmap(m1))
!$$$c     - alkane type molecule -
!$$$                        if ( inirot(ibox) .eq. 0 ) then
!$$$                           rot = random() * 360.0d0 * degrad
!$$$                        elseif ( inirot(ibox) .gt. 0 ) then
!$$$                           rot = dble(inirot(ibox)) * degrad
!$$$                        else
!$$$                           if ( mod(jc,2) .eq. 0 ) then
!$$$                              rot = dble(inirot(ibox)) * degrad
!$$$                           else
!$$$                              rot = -(dble(inirot(ibox)) * degrad)
!$$$                           end if
!$$$                        end if
!$$$                        
!$$$                        if ( inben(intemp,m1) .gt. 0 ) then
!$$$                           ibtype = itben(intemp,1,1)
!$$$                           angold = brben(ibtype) / 2.0d0
!$$$                           if ( m .eq. 2 ) then
!$$$                              ibtype = itben(intemp,m1,1)
!$$$                              angnew = brben(ibtype) - angold
!$$$                           else
!$$$                              ibtype = itben(intemp,m2,1)
!$$$                              angnew = brben(ibtype) - angold
!$$$                           end if 
!$$$c     write(iou,*) 'angold',angold*raddeg,
!$$$c     +                       '   angnew',angnew*raddeg
!$$$                           angold = angnew
!$$$                           
!$$$                           ibtype = itvib(intemp,m,1)
!$$$                           znext = dsin(angnew) * brvib(ibtype)
!$$$                           xynext = dcos(angnew) * brvib(ibtype)
!$$$c                           write(iou,*) 'znext',znext,'   xynext',xynext
!$$$                        else
!$$$                           ibtype = itvib(intemp,m,1)
!$$$                           znext = brvib(ibtype)
!$$$                           xynext = 0.0d0
!$$$                        end if
!$$$                        
!$$$                        if ( mod(m,2) .eq. 0 ) then
!$$$                           xnext = dcos(rot) * xynext
!$$$                           ynext = dsin(rot) * xynext
!$$$                        else
!$$$                           xnext = -(dcos(rot) * xynext)
!$$$                           ynext = -(dsin(rot) * xynext)
!$$$                        end if
                     end if

                     if (n.le.nchbox(ibox)) then
                        rxu(nn,m) = rxu(nn,m1) + xnext
                        ryu(nn,m) = ryu(nn,m1) + ynext
                        rzu(nn,m) = rzu(nn,m1) + znext
                     end if
                
 98               continue

 99            continue
 100        continue
 101     continue
 102  continue
! -----------------------------------------------------

! *** check initial structure ***

      nn = nchain
      if (lgrand) nn=nchain

      aben = 0.0d0
      ator = 0.0d0

      do n = 1, nn

         imolty = moltyp(n)
!         write(iou,*) 'n',n,'   imolty',imolty

            
         if ( lbranch(imolty) ) then
! - branched molecule with connectivity table -
! - go through entire chain -
! - calculate all bonds vectors and lengths
! - calculate all stretching, bending, and torsional potentials
! - that have an end-bead with an index smaller than the current bead
            do ii = 1, nunit(imolty)
               rxui=rxu(n,ii)
               ryui=ryu(n,ii)
               rzui=rzu(n,ii)

               if ( n .eq. 1 .or. m .eq. 1 ) write(iou," ('coord. unit:   ',2i4,3f9.3,i6)") n,ii,rxui,ryui,rzui ,nboxi(n)
               do iivib = 1, invib(1,ii)
                  jj = ijvib(1,ii,iivib)
                  rxvec(ii,jj) = rxu(n,jj) - rxui
                  ryvec(ii,jj) = ryu(n,jj) - ryui
                  rzvec(ii,jj) = rzu(n,jj) - rzui
                  distanceij(ii,jj) = dsqrt( rxvec(ii,jj)**2 + ryvec(ii,jj)**2 + rzvec(ii,jj)**2 )
                  if ( nunit(imolty) .ne. nugrow(imolty) )then
!                 --- account for explct atoms in opposite direction
                     rxvec(jj,ii)   = -rxvec(ii,jj)
                     ryvec(jj,ii)   = -ryvec(ii,jj)
                     rzvec(jj,ii)   = -rzvec(ii,jj)
                     distanceij(jj,ii) = distanceij(ii,jj)
                  end if
               end do
            end do
 
            do j = 1, nunit(imolty)

! - vibrations -
               do iivib = 1, invib(1,j)
                  jj = ijvib(1,j,iivib)
                  if ( n .eq. 1 ) write(iou,"('bond with units:',2i3 ,'   length:',f9.4)") j,jj,distanceij(j,jj)
               end do

! - bending -
               do jjben = 1, inben(imolty,j)
                  ip2 = ijben3(imolty,j,jjben)
                  ip1 = ijben2(imolty,j,jjben)
                  it  = itben(imolty,j,jjben)
                  thetac = ( rxvec(ip1,j)*rxvec(ip1,ip2) + ryvec(ip1,j)*ryvec(ip1,ip2) + rzvec(ip1,j)*rzvec(ip1,ip2) ) / ( distanceij(ip1,j)*distanceij(ip1,ip2) )
                  theta = dacos(thetac)
                  vbend = brbenk(it) * (theta-brben(it))**2
                  aben = aben + vbend
!                  if ( n .eq. 1 ) then
!                  write(iou,*) 'theta',theta,'vbend',vbend
!                  write(iou,*) 'brben',brben(it),'brbenk',brbenk(it)
!                  end if
                  if ( n .eq. 1 ) write(iou,"('bend with units:',3i3 ,'   type:',i3,'   angle:',f9.4,f9.2)") j,ip1,ip2,it ,theta*raddeg,vbend
               end do

! - torsions -
               do jjtor = 1, intor(imolty,j)
                  ip3 = ijtor4(imolty,j,jjtor)
                  ip1 = ijtor2(imolty,j,jjtor)
                  ip2 = ijtor3(imolty,j,jjtor)
                  it  = ittor(imolty,j,jjtor)
! *** calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                  xaa1 = ryvec(ip1,j) * rzvec(ip2,ip1) + rzvec(ip1,j) * ryvec(ip1,ip2)
                  yaa1 = rzvec(ip1,j) * rxvec(ip2,ip1) + rxvec(ip1,j) * rzvec(ip1,ip2)
                  zaa1 = rxvec(ip1,j) * ryvec(ip2,ip1) + ryvec(ip1,j) * rxvec(ip1,ip2)
                  xa1a2 = ryvec(ip1,ip2) * rzvec(ip2,ip3) + rzvec(ip1,ip2) * ryvec(ip3,ip2)
                  ya1a2 = rzvec(ip1,ip2) * rxvec(ip2,ip3) + rxvec(ip1,ip2) * rzvec(ip3,ip2)
                  za1a2 = rxvec(ip1,ip2) * ryvec(ip2,ip3) + ryvec(ip1,ip2) * rxvec(ip3,ip2)
! *** calculate lengths of cross products ***
                  daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                  da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
! *** calculate dot product of cross products ***
                  dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                  thetac = -(dot / ( daa1 * da1a2 ))
                  if (thetac.gt.1.0d0) thetac=1.0d0
                  if (thetac.lt.-1.0d0) thetac=-1.0d0

!     KEA -- added for extending range to +/- 180
                  if (it .ge. 500) then
!     *** calculate cross product of cross products ***
                     xcc = yaa1*za1a2 - zaa1*ya1a2
                     ycc = zaa1*xa1a2 - xaa1*za1a2
                     zcc = xaa1*ya1a2 - yaa1*xa1a2
!     *** calculate scalar triple product ***
                     tcc = xcc*rxvec(ip1,ip2) + ycc*ryvec(ip1,ip2) + zcc*rzvec(ip1,ip2)
                     theta=dacos(thetac)

                     if (tcc .lt. 0.0d0) theta = -theta
                     if (it.le.600) then
                        call splint(theta,spltor,it)
                     else
                        call lininter(theta,spltor,it)
                     end if

                     vtg=spltor
                  else
                     vtg = vtorso (thetac,it)
                  end if

                  ator = ator + vtg
!                  if ( n .eq. 1 ) write(iou,*) 'thetac',thetac,'vtg',vtg
                  if ( n .eq. 1 ) write(iou,"('tors with units:',4i3 ,'   type:',i3,'   angle:',f9.4,f9.2)") j,ip1,ip2,ip3 ,it,dacos(thetac)*raddeg,vtg
               end do
            end do

         else
! ---                  
            ibuild = nugrow(imolty)
            do m = 1, ibuild

               m1 = m - 1
               m2 = m - 2

               if ( m1 .gt. 0 ) then
                  x1 = rxu(n,m) - rxu(n,m1)
                  y1 = ryu(n,m) - ryu(n,m1)
                  z1 = rzu(n,m) - rzu(n,m1)
                  d1 = dsqrt( x1**2 + y1**2 + z1**2 ) 
                  if ( m2 .gt. 0 ) then
                     x2 = rxu(n,m2) - rxu(n,m1)
                     y2 = ryu(n,m2) - ryu(n,m1)
                     z2 = rzu(n,m2) - rzu(n,m1)
                     d2 = dsqrt( x2**2 + y2**2 + z2**2 ) 
                     bang = dacos((x1*x2+y1*y2+ z1*z2)/(d1*d2))*raddeg
                  else
                     bang  = 0.0d0
                  end if
                  blen = d1
               else
                  blen = 0.0d0
                  bang = 0.0d0
               end if

!               if ( n .eq. 1 .or. m .eq. 1 )
!     &              write(iou,'(2i4,5f9.3,i6)') n,m,rxu(n,m),ryu(n,m),rzu(n,m),
!     &                            blen,bang,nboxi(n)

!               write(iou,'(2i4,5f9.3,i6)') n,m,rxu(n,m),ryu(n,m),rzu(n,m),
!     &                            blen,bang,nboxi(n)

            end do

         end if
      end do

      do i=1,nchain
         imolty = moltyp(i)
         if ( nugrow(imolty) .ne. nunit(imolty) ) then
            call explct(i,vdummy,.true.,.false.)
         end if
      end do
!     --- set up intial charges on the atoms
      do i = 1,nchain
         imolty = moltyp(i)
         do ii = 1,nunit(imolty)
            ntii = ntype(imolty,ii)
            qqu(i,ii) = qelect(ntii)
         end do
      end do

      if (myid.eq.0) then
         write(iou,*) 'aben',aben/2.0d0,'ator',ator/2.0d0
      end if

      return
      end
