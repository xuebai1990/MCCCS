      subroutine initia(qelect)

c initia
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

c *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'poten.inc'
      include 'conver.inc'
      include 'inpar.inc'
      include 'ensemble.inc' 
      include 'connect.inc'
      include 'cbmc.inc'


      logical lhere(nmax)       
      integer i,j,m,m1,m2,n,nn,ic,jc,kc,it,ip1,ip2,ip3,ii,jj
     &       ,iivib,jjben,jjtor,intemp,imol,nt
     &       ,ibtype,imolty,ibuild,rand_id,offset,count_chain

      integer iboxst,iboxed,ibox,pct,check,chktot
      integer mcmt(ntmax,nbxmax),pcmt(ntmax),mcmtma(ntmax,nbxmax)
      integer unitc,ntii,ichain

      integer bmap(numax),imap(numax),zzz,prev,ifrom,nsave
      logical lacc(numax),lgrow,lterm,lgrown(ntmax)

      double precision xtemp(numax),ytemp(numax),ztemp(numax)
      double precision ddum

      double precision ux,uy,uz,xshift,dic
     &     ,xnext,ynext,znext,xynext,angold,angnew,rot
     &     ,x1,y1,z1,d1,x2,y2,z2,d2,bang,blen,random

      double precision rxui,ryui,rzui,xaa1,yaa1,zaa1,daa1
     &                ,xa1a2,ya1a2,za1a2,da1a2 

      double precision xcc,ycc,zcc,tcc,spltor

      double precision vtorso,vbend,vtg,thetac,theta,dot,aben,ator

      double precision xvec(numax,numax),yvec(numax,numax)
     &                ,zvec(numax,numax),distij(numax,numax)

      double precision samx(ntmax,numax),samy(ntmax,numax)
     &     ,samz(ntmax,numax)

      double precision vdummy
      double precision qelect
      dimension qelect(nntype)

      dimension ux(nbxmax),uy(nbxmax),uz(nbxmax)
      dimension check(ntmax)

C --------------------------------------------------------------------

      write(iou,*) 
      write(iou,*) 'subroutine initia'
      write(iou,*) 


c     --- initialize nchbox ---
      do i=1,nbox
         nchbox(i) = 0
      enddo

      iboxst = 1
      iboxed = nbox

         chktot = 0

         do i = 1,nmolty
            check(i) = 0
            do j = 1, nbox
               nchbox(j) = nchbox(j) + ininch(i,j)
               check(i)= check(i) + ininch(i,j) 
            enddo

            chktot = chktot + check(i)
         enddo

         if ( chktot .ne. nchain ) then
            write(iou,*) 'inconsistant number of chains in INITIA'
            do j = 1,nbox
               write(iou,*) 'ininch',j,(ininch(i,j),i=1,nmolty)
            enddo
            write(iou,*) 'nchain',nchain
            stop
         endif
         
         do i = 1, nmolty
            if ( temtyp(i) .ne. check(i) ) then
               write(iou,*) 'inconsistant number of chains in INITIA'
               write(iou,*) 'moltyp',i,(ininch(i,j),j=1,nbox)
               write(iou,*) 'temtyp:',temtyp(i)
               stop
            endif            
         enddo

         do i = iboxst,iboxed
            unitc = inix(i)*iniy(i)*iniz(i)
            if ( nchbox(i) .gt. unitc ) then
               write(iou,*) 'unit cell too small in box',i
               stop
            endif
         enddo

C -----------------------------------------------------------------------------
 
c *** calculation of unit cell dimensions ***
      do i = 1,nbox
         ux(i) = boxlx(i) / dble(inix(i)) 
         uy(i) = boxly(i) / dble(iniy(i))
         uz(i) = boxlz(i) / dble(iniz(i))
         write(iou,*) 'box',i
         write(iou,*) 'ini',inix(i),iniy(i),iniz(i)
         write(iou,*) 'box',boxlx(i),boxly(i),boxlz(i)
         write(iou,*) 'uni',ux(i),uy(i),uz(i)
      enddo

c - count number of molecules of each type -
      do i = 1,nmolty
         mcmt(i,1) = 0
      enddo

      do i=1,nchain
         imolty = moltyp(i)
         mcmt(imolty,1) = mcmt(imolty,1) + 1
         lhere(i) = .false.
      enddo

      do i = 1,nmolty
         if ( mcmt(i,1) .ne. check(i) ) then
            write(iou,*) 'inconsistant number of type in INITIA'
            write(iou,*) 'mcmt(i,total),check(i)',mcmt(i,1),check(i)
            stop
         endif
      enddo

      do i = 1,nmolty
         do j = 1, nbox
            mcmt(i,j) = ininch(i,j)
            mcmtma(i,j) = 0
         enddo
      enddo
      do ibox = 1,nbox
         mcmtma(1,ibox) = mcmt(1,ibox)
         do i = 2, nmolty
            mcmtma(i,ibox) = mcmtma(i-1,ibox) + mcmt(i,ibox)
         enddo
      enddo

      write(iou,*) 'nmolty',nmolty
      write(iou,*) '   mcmt',((mcmt(i,ibox),i=1,nmolty),ibox=1,nbox)

c *****************************
c *** calculate coordinates ***

c     read sample structure from unit 78 -
      open(unit=78,FILE='input_struc.xyz',status="unknown")

      do i = 1, nmolty         
         lgrown(i) = .false.
         if ( lbranch(i) ) then
            read(78,*)
            do m = 1, nunit(i)
               read(78,*) samx(i,m), samy(i,m), samz(i,m)
            enddo
         elseif ( .not. lbranch(i)) then
c * if lbranch is false but the molecule is not linear attempt
c * to grow it with cbmc
            lgrow = .false.
            do m = 1,nunit(i)
               if (invib(i,m) .gt. 2) then
                  lgrow = .true.
               endif
            enddo
            if (lgrow) then

               write(iou,*) 'growing a sample structure with CBMC'

               if (nunit(i) .ne. nugrow(i)) then
                  write(iou,*) 'Cant grow molecule.  Please',
     &                 ' provide a structure via fort.78'
                  stop
               endif
c * put the first bead at the origin
               rxnew(1) = 0.0d0
               rynew(1) = 0.0d0
               rznew(1) = 0.0d0

c * determine the growth schedule
               call schedule(nugrow(i),i,ifrom,1,0,2)

c * actually grow the structure
               nsave = nchain

               moltyp(1) = i
               do m = 1,nunit(i)
                  qqu(1,m) = qelect(ntype(i,m))
               enddo

               nchain = 1

               call rosenbluth( .true.,lterm,1,1,i,ifrom
     &              ,1,nugrow(i),ddum,.false.,ddum,2 )

               if (lterm) then
                  write(iou,*) 'error in initia growing molecule'
                  write(iou,*) 'maybe increasing nchoi would help?'
                  stop
               endif

c * return the value of nchain
               nchain = nsave
c * assign the coordinates
               do m = 1,nunit(i)
                  samx(i,m) = rxnew(m)
                  samy(i,m) = rynew(m)
                  samz(i,m) = rznew(m)
               enddo

               lgrown(i) = .true.

            endif
         endif
      enddo

      close(unit=78)

      do i = 1,nmolty
         if (lgrown(i)) then
            lbranch(i) = .true.
         endif
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     - inimix = 0 : take molecules at random
c     - inimix > 0 : take molecules in order (first type I etc.)
c     - inimix < 0 : take molecules in alternating order
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

      count_chain = 0 
      offset = 0

      do ibox = iboxst,iboxed
         if(nmolty .gt. 1) then
           if(inimix(ibox).eq.0) then
              if (ibox.eq.1) then
                 offset = 0
              else
                 offset = offset+nchbox(ibox-1)
              endif
 18           rand_id = idint(dble(nchbox(ibox))*random())+ 1 + offset
              if (.not.lhere(rand_id)) then
                 count_chain = count_chain + 1
                 lhere(rand_id) = .true.
                 do imolty = 1,nmolty  
                    if (count_chain.le.(mcmtma(imolty,ibox)+offset)) 
     &                                                          then
                       moltyp(rand_id) = imolty     
                       goto 20 
                    endif 
                 enddo
 20              continue
!             write(6,*) count_chain, rand_id,moltyp(rand_id) 
              else
                 goto 18
              endif
              if (count_chain.lt.(nchbox(ibox)+offset)) then
                 goto 18
              endif 
           endif
         endif
      enddo   


      nn = 0

      do 102 ibox = iboxst,iboxed
         do imol = 1,nmolty
            pcmt(imol) = 0
         enddo

         n = 0

         do 101 kc = 0, iniz(ibox)-1
            if ( mod(kc,2) .eq. 0) then
               xshift = 0.0d0
            else
               xshift = dshift(ibox)
            endif
            
            do 100 ic = 0, inix(ibox)-1
               
               do 99 jc = 0, iniy(ibox)-1

                  if ( mod(jc,2) .eq. 0) then
                     dic = 0.0d0
                  else
                     dic = 0.5d0
                  endif

                  n=n+1

                  if (n .le. nchbox(ibox) ) then
                     nn=nn+1
                     rxu(nn,1) = ( dble(ic) + dic )*ux(ibox)+ xshift
                     ryu(nn,1) = dble(jc) * uy(ibox) 
                     rzu(nn,1) = dble(kc) * uz(ibox) + zshift(ibox)
                     nboxi(nn) = ibox
                  else
                     goto 102
                  endif
                  
c                  write(iou,*) 'nn',nn
c                  write(iou,*) 'ic',ic,'   jc',jc,'   kc',kc
                  
c     - inimix > 0 : take molecules in order (first type I etc.)
c     - inimix < 0 : take molecules in alternating order

                  if ( nmolty .gt. 1 ) then
                     if ( inimix(ibox) .gt. 0 ) then
                        do imol = 1, nmolty
                           if ( n .le. mcmtma(imol,ibox) ) then
                              intemp = imol
                              goto 19
                           endif
                        enddo
 19                     continue
                     elseif ( inimix(ibox) .lt. 0 ) then
                        do imol = 1, nmolty
                           nt = n - imol
                           if ( mod( nt, nmolty ) .eq. 0 ) then
                              intemp = imol
                           endif
                        enddo
                     endif
                  else
                     intemp = 1
                  endif

                  if (inimix(ibox).eq.0) then
                      intemp = moltyp(nn)
                  else
                       moltyp(nn) = intemp
                  endif 

                  ncmt(ibox,intemp) = ncmt(ibox,intemp) + 1
                  
c                  write(iou,*) 'intemp', intemp     

                  if ( lbranch(intemp) ) then
                     ibuild = nunit(intemp)
                  else
                     ibuild = nugrow(intemp)
                  endif

                  if ( .not. lbranch(intemp)) then
c *************************************************
c *** start determination of linear chain positions
c *** allowing for numbering out of order
c *** note: doesn't exactly create equilibrium structure with
c *** respect to bond angles or torsions, but that will shake out
c *** with CBMC anyways.  Should at least take away the overlaps of
c *** the previous method
c *************************************************
c * first need to determine re-mapped bead order- search through connectivity
c
c   call the results map(i) where i=1 is one chain end, and its
c   value is equal to the bead number of that end
c 
c   for example, methanol oxygen first, then hydrogen, then CH3
c   
c                        H---O--CH3
c
c        bead numbers:   2 - 1 - 3
c
c        bmap(1) = 2
c        bmap(2) = 1
c        bmap(3) = 3
c
c        the inverse map is just the opposite:
c
c        imap(1) = 2
c        imap(2) = 1
c        imap(3) = 3
c
c * initialize accounted for variable
                     do m = 1,ibuild
                        lacc(m) = .false.
                     enddo

c * first find the end with the lowest number
                     zzz = ibuild
                     do m = 1,ibuild
                        if (invib(intemp,m) .le. 1) then
                           if (m .le. zzz) then
                              zzz = m
                           endif
                        elseif (invib(intemp,m) .gt. 2) then
                           write(iou,*) 'initia only works for linear',
     &                          ' molecules!  Maybe you should make',
     &                          ' a fort.78 file and use lbranch?'
                           stop
                        endif
                     enddo

                     bmap(1) = zzz
                     imap(zzz) = 1
                     lacc(zzz) = .true.

c * now determine the rest
                     do m = 2,ibuild
                        prev = bmap(m-1)
                        do zzz = 1,ibuild
                           if (ijvib(intemp,zzz,1) .eq. prev
     &                          .or. ijvib(intemp,zzz,2) .eq. prev) then
                              if (.not. lacc(zzz)) then
                                 bmap(m) = zzz
                                 imap(zzz) = m
                                 lacc(zzz) = .true.
                              endif
                           endif
                        enddo
                     enddo

c * now use old method with re-mapped numbers:
c * put first end at origin:
                     xtemp(1) = 0.0d0
                     ytemp(1) = 0.0d0
                     ztemp(1) = 0.0d0

c * now we need to loop over all the other beads:
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
                           endif
                        endif

                        if ( inben(intemp,bmap(m1)) .gt. 0 ) then
                           ibtype = itben(intemp,bmap(m1),1)
                           angold = brben(ibtype) / 2.0d0

                           if ( m .eq. 2 ) then
                              ibtype = itben(intemp,bmap(m1),1)
                              angnew = brben(ibtype) - angold
                           else
                              ibtype = itben(intemp,bmap(m2),1)
                              angnew = brben(ibtype) - angold
                           endif 
c     write(iou,*) 'angold',angold*raddeg,
c     +                       '   angnew',angnew*raddeg
                           angold = angnew
                           
c * need to search for proper bond length
                           do zzz = 1,invib(intemp,bmap(m))
                              if (ijvib(intemp,bmap(m),zzz) 
     &                             .eq. bmap(m1)) then
                                 ibtype = itvib(intemp,bmap(m),zzz)
                              endif
                           enddo

                           ztemp(m) = dsin(angnew) * brvib(ibtype)
                           xynext = dcos(angnew) * brvib(ibtype)
c                           write(iou,*) 'znext',znext,'   xynext',xynext
                        else
c * need to search for proper bond length
                           do zzz = 1,invib(intemp,bmap(m))
                              if (ijvib(intemp,bmap(m),zzz) 
     &                             .eq. bmap(m1)) then
                                 ibtype = itvib(intemp,bmap(m),zzz)
                              endif
                           enddo

c                           ztemp(m) = dsin(angnew) * brvib(ibtype)
c                           xynext = dcos(angnew) * brvib(ibtype)

                           ztemp(m) = brvib(ibtype)
                           xynext = 0.0d0
                        endif
                        
                        if ( mod(m,2) .eq. 0 ) then
                           xtemp(m) = dcos(rot) * xynext
                           ytemp(m) = dsin(rot) * xynext
                        else
                           xtemp(m) = -(dcos(rot) * xynext)
                           ytemp(m) = -(dsin(rot) * xynext)
                        endif

                        xtemp(m) = xtemp(m1) + xtemp(m)
                        ytemp(m) = ytemp(m1) + ytemp(m)
                        ztemp(m) = ztemp(m1) + ztemp(m)


                     enddo

c * translate so that first bead number is at origin
                     do m = 1,ibuild
                        if (m .ne. imap(1)) then
                          xtemp(m) = xtemp(m) - xtemp(imap(1)) 
                          ytemp(m) = ytemp(m) - ytemp(imap(1)) 
                          ztemp(m) = ztemp(m) - ztemp(imap(1)) 
                        endif
                     enddo
                     xtemp(imap(1)) = 0.0d0
                     ytemp(imap(1)) = 0.0d0
                     ztemp(imap(1)) = 0.0d0

                  endif
c *** end linear determination
c ****************************
c                  write(iou,*) 'ibuild',ibuild
                  do 98 m = 2, ibuild
                     
                     m1 = m - 1
                     m2 = m - 2
c                     write(iou,*) 'intemp',intemp 
                     if ( lbranch(intemp) ) then
c     - branched molecule with sample structure -
                        xnext = samx(intemp,m) -samx(intemp,m1)
                        ynext = samy(intemp,m) -samy(intemp,m1)
                        znext = samz(intemp,m) -samz(intemp,m1)
                     else
c * linear molecule determined above- replacing old code that follows.
                        xnext = xtemp(bmap(m)) - xtemp(bmap(m1))
                        ynext = ytemp(bmap(m)) - ytemp(bmap(m1))
                        znext = ztemp(bmap(m)) - ztemp(bmap(m1))
c$$$c     - alkane type molecule -
c$$$                        if ( inirot(ibox) .eq. 0 ) then
c$$$                           rot = random() * 360.0d0 * degrad
c$$$                        elseif ( inirot(ibox) .gt. 0 ) then
c$$$                           rot = dble(inirot(ibox)) * degrad
c$$$                        else
c$$$                           if ( mod(jc,2) .eq. 0 ) then
c$$$                              rot = dble(inirot(ibox)) * degrad
c$$$                           else
c$$$                              rot = -(dble(inirot(ibox)) * degrad)
c$$$                           endif
c$$$                        endif
c$$$                        
c$$$                        if ( inben(intemp,m1) .gt. 0 ) then
c$$$                           ibtype = itben(intemp,1,1)
c$$$                           angold = brben(ibtype) / 2.0d0
c$$$                           if ( m .eq. 2 ) then
c$$$                              ibtype = itben(intemp,m1,1)
c$$$                              angnew = brben(ibtype) - angold
c$$$                           else
c$$$                              ibtype = itben(intemp,m2,1)
c$$$                              angnew = brben(ibtype) - angold
c$$$                           endif 
c$$$c     write(iou,*) 'angold',angold*raddeg,
c$$$c     +                       '   angnew',angnew*raddeg
c$$$                           angold = angnew
c$$$                           
c$$$                           ibtype = itvib(intemp,m,1)
c$$$                           znext = dsin(angnew) * brvib(ibtype)
c$$$                           xynext = dcos(angnew) * brvib(ibtype)
c$$$c                           write(iou,*) 'znext',znext,'   xynext',xynext
c$$$                        else
c$$$                           ibtype = itvib(intemp,m,1)
c$$$                           znext = brvib(ibtype)
c$$$                           xynext = 0.0d0
c$$$                        endif
c$$$                        
c$$$                        if ( mod(m,2) .eq. 0 ) then
c$$$                           xnext = dcos(rot) * xynext
c$$$                           ynext = dsin(rot) * xynext
c$$$                        else
c$$$                           xnext = -(dcos(rot) * xynext)
c$$$                           ynext = -(dsin(rot) * xynext)
c$$$                        endif
                     endif

                     if (n.le.nchbox(ibox)) then
                        rxu(nn,m) = rxu(nn,m1) + xnext
                        ryu(nn,m) = ryu(nn,m1) + ynext
                        rzu(nn,m) = rzu(nn,m1) + znext
                     endif
                
 98               continue

 99            continue
 100        continue
 101     continue
 102  continue
c -----------------------------------------------------

c *** check initial structure ***

      nn = nchain
      if (lgrand) nn=nchain

      aben = 0.0d0
      ator = 0.0d0

      do 200 n = 1, nn

         imolty = moltyp(n)
c         write(iou,*) 'n',n,'   imolty',imolty

            
         if ( lbranch(imolty) ) then
c - branched molecule with connectivity table -
c - go through entire chain -
c - calculate all bonds vectors and lengths
c - calculate all stretching, bending, and torsional potentials
c - that have an end-bead with an index smaller than the current bead
            do ii = 1, nunit(imolty)
               rxui=rxu(n,ii)
               ryui=ryu(n,ii)
               rzui=rzu(n,ii)

               if ( n .eq. 1 .or. m .eq. 1 )
     &              write(iou,1002) n,ii,rxui,ryui,rzui,nboxi(n)
               do iivib = 1, invib(1,ii)
                  jj = ijvib(1,ii,iivib)
                  xvec(ii,jj) = rxu(n,jj) - rxui
                  yvec(ii,jj) = ryu(n,jj) - ryui
                  zvec(ii,jj) = rzu(n,jj) - rzui
                  distij(ii,jj) = dsqrt( xvec(ii,jj)**2
     +                 + yvec(ii,jj)**2 + zvec(ii,jj)**2 )
                  if ( nunit(imolty) .ne. nugrow(imolty) )then
c                 --- account for explct atoms in opposite direction
                     xvec(jj,ii)   = -xvec(ii,jj)
                     yvec(jj,ii)   = -yvec(ii,jj)
                     zvec(jj,ii)   = -zvec(ii,jj)
                     distij(jj,ii) = distij(ii,jj)
                  endif
               enddo
            enddo
 
            do j = 1, nunit(imolty)

c - vibrations -
               do iivib = 1, invib(1,j)
                  jj = ijvib(1,j,iivib)
                  if ( n .eq. 1 ) write(iou,1003) j,jj,distij(j,jj)
               enddo

c - bending -
               do jjben = 1, inben(imolty,j)
                  ip2 = ijben3(imolty,j,jjben)
                  ip1 = ijben2(imolty,j,jjben)
                  it  = itben(imolty,j,jjben)
                  thetac = ( xvec(ip1,j)*xvec(ip1,ip2) +
     +                 yvec(ip1,j)*yvec(ip1,ip2) +
     +                 zvec(ip1,j)*zvec(ip1,ip2) ) /
     +                 ( distij(ip1,j)*distij(ip1,ip2) )
                  theta = dacos(thetac)
                  vbend = brbenk(it) * (theta-brben(it))**2
                  aben = aben + vbend
c                  if ( n .eq. 1 ) then
c                  write(iou,*) 'theta',theta,'vbend',vbend
c                  write(iou,*) 'brben',brben(it),'brbenk',brbenk(it)
c                  endif
                  if ( n .eq. 1 ) write(iou,1004) j,ip1,ip2,it,
     &                 theta*raddeg,vbend
               enddo

c - torsions -
               do jjtor = 1, intor(imolty,j)
                  ip3 = ijtor4(imolty,j,jjtor)
                  ip1 = ijtor2(imolty,j,jjtor)
                  ip2 = ijtor3(imolty,j,jjtor)
                  it  = ittor(imolty,j,jjtor)
c *** calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                  xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     +                 zvec(ip1,j) * yvec(ip1,ip2)
                  yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     +                 xvec(ip1,j) * zvec(ip1,ip2)
                  zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     +                 yvec(ip1,j) * xvec(ip1,ip2)
                  xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     +                 zvec(ip1,ip2) * yvec(ip3,ip2)
                  ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     +                 xvec(ip1,ip2) * zvec(ip3,ip2)
                  za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     +                 yvec(ip1,ip2) * xvec(ip3,ip2)
c *** calculate lengths of cross products ***
                  daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                  da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
c *** calculate dot product of cross products ***
                  dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                  thetac = -(dot / ( daa1 * da1a2 ))
                  if (thetac.gt.1.0d0) thetac=1.0d0
                  if (thetac.lt.-1.0d0) thetac=-1.0d0

c     KEA -- added for extending range to +/- 180
                  if (it .ge. 500) then
c     *** calculate cross product of cross products ***
                     xcc = yaa1*za1a2 - zaa1*ya1a2
                     ycc = zaa1*xa1a2 - xaa1*za1a2
                     zcc = xaa1*ya1a2 - yaa1*xa1a2
c     *** calculate scalar triple product ***
                     tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2)
     &                    + zcc*zvec(ip1,ip2)
                     theta=dacos(thetac)

                     if (tcc .lt. 0.0d0) theta = -theta
                     if (it.le.600) then
                        call splint(theta,spltor,it)
                     else
                        call lininter(theta,spltor,it)
                     endif

                     vtg=spltor
                  else
                     vtg = vtorso (thetac,it)
                  endif

                  ator = ator + vtg
c                  if ( n .eq. 1 ) write(iou,*) 'thetac',thetac,'vtg',vtg
                  if ( n .eq. 1 ) write(iou,1005) j,ip1,ip2,ip3,it,
     &                 dacos(thetac)*raddeg,vtg
               enddo
            enddo

         else
c ---                  
            ibuild = nugrow(imolty)
            do 199 m = 1, ibuild

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
                  endif
                  blen = d1
               else
                  blen = 0.0d0
                  bang = 0.0d0
               endif

c               if ( n .eq. 1 .or. m .eq. 1 )
c     &              write(iou,1001) n,m,rxu(n,m),ryu(n,m),rzu(n,m),
c     &                            blen,bang,nboxi(n)

c               write(iou,1001) n,m,rxu(n,m),ryu(n,m),rzu(n,m),
c     &                            blen,bang,nboxi(n)

 199        continue

         endif

 200  continue

      do i=1,nchain
         imolty = moltyp(i)
         if ( nugrow(imolty) .ne. nunit(imolty) ) then
            call explct(i,vdummy,.true.,.false.)
         endif
      enddo
c     --- set up intial charges on the atoms
      do i = 1,nchain
         imolty = moltyp(i)
         do ii = 1,nunit(imolty)
            ntii = ntype(imolty,ii)
            qqu(i,ii) = qelect(ntii)
         enddo
      enddo


      write(iou,*) 'aben',aben/2.0d0,'ator',ator/2.0d0

 1001 format(2i4,5f9.3,i6)
 1002 format('coord. unit:   ',2i4,3f9.3,i6)
 1003 format('bond with units:',2i3,'   length:',f9.4)
 1004 format('bend with units:',3i3,'   type:',i3,'   angle:',f9.4,f9.2)
 1005 format('tors with units:',4i3,'   type:',i3,'   angle:',f9.4,f9.2)

      return
      end
