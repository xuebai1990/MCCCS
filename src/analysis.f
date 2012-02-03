      subroutine analysis(switch)

c     *******************************************************************
c     *** Modified to perform analysis on the fly based on anal10.f   ***
c     *** [Marcus Martin] by Neeraj Rai 07/14/04                      ***
c     *******************************************************************

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'inputdata.inc'
      include 'connect.inc'
      include 'system.inc'
      include 'cell.inc' 
      include 'gor.inc'

      integer switch,ichain
      logical lskip
      integer bin,k,kk,box,z
      integer tempbx,dummy,xx,yy,g,gg
      integer chnum,imolty,i,ii,j,jj
     &     ,jmolty,binadj,ntij,ntii,ntjj,jstart,istart,ntji
      integer iivib,ip1,ip2,ip3,gaudef,uu,dum,units,zzz,zz1

      double precision avolume,rho,binstep,vec_hist
     &     ,rxui,ryui,rzui,xxideal
     &     ,rxuij,ryuij,rzuij,bx,hbx,by,hby,bz,hbz,ruijsq,ruij
     &     ,numxx,numyy,count,const,rlower,rupper,nideal,analhist
      double precision comanalhist
      double precision boxmin,shlsumx,shlsumy,rcutsq

C   NOW DEFINED IN CONTROL.INC

Cc     ntmax = number of moltyps, ntdifmx = max number of diff beads in sim
C      parameter (nbinmx=200,nbxmax=1,ntmax=5,nmax=1600,numax=18
C     & ,ntdifmx=8)

      integer bend,iv,iuvib,iuv,iutest
      double precision onepi,ang_bin_size,
     &    value,total,degree
      parameter (onepi = 3.141592654d0)
      integer torsion, tor_code,itor,iutor,patt,bthree
     &     ,decimal,power
      double precision xcc,ycc,zcc,tcc,fplus,fminus
     &     ,ftrans
      double precision xvec,yvec,zvec,distij,xaa1,yaa1,zaa1,xa1a2
     &     ,ya1a2,za1a2,daa1,da1a2,dot,thetac,theta
     &     ,ratio
     &     ,psum,tempzcm,tempmasst,slab_vol

      dimension patt(tor_max),tempzcm(nbxmax),tempmasst(nbxmax)
      dimension tor_code(numax,numax)
      dimension vec_hist(nbxmax,ntmax,nbinmax_ete)  

********* Taking out charge part*************************

CCc     --- variables used in the charge parts
CC     logical qhere
CC      integer qbin,qbinmax,qbins,qqcode,imol,iunit
CC      integer qbin,qbins,qqcode,imol,iunit
CC    NOW DEFINED IN CONTROL.INC
CC      parameter (qbinmax=1000)

CC      double precision qdisp
CC     double precision qmin,qmax,qdiff,qstep,qdummy
CC      dimension qanalhist(numax*ntdifmx+numax,nbxmax,qbinmax)
CC      dimension qcount(numax*ntdifmx+numax,nbxmax)

CCc     --- variables used in the dipole parts
CC     logical ldipole
CC     integer qblock,block,nblock
CC      double precision dipole,dicount,dipx,dipy,dipz,diconv,stddev
CC     &     ,avera,diprev,dcprev,dipblk
CC     dimension dipole(ntdifmx,nbxmax),dicount(ntdifmx,nbxmax)
CC     dimension diprev(ntdifmx,nbxmax),dcprev(ntdifmx,nbxmax)
CC      dimension dipblk(ntdifmx,nbxmax,20)
********** end charge part ********************************

      dimension analhist(ntdifmx*ntdifmx*ntmax*ntmax,nbinmx)
      dimension comanalhist(ntmax*ntmax,nbinmx)
      dimension avolume(nbxmax)
      dimension rho(ntmax*ntmax*ntdifmx)
      dimension count(nbxmax,ntdifmx*ntmax)
      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)
     &     ,distij(numax,numax)
      dimension boxmin(nbxmax)


      if(switch.eq.0) then

      if(nhere.gt.ntdifmx) then
        write(6,*) 'nhere greater than ntdifmax', nhere, ntdifmx
        stop 'choose a larger ntdifmx in control.inc'
      endif
 
      if(nbin.gt.nbinmx) then
         write(6,*) 'number of bins "nbin" .gt. nbinmx', nbin, nbinmx
         stop 'choose a larger nbinmx in control.inc'
       endif  
      
      do kk=1,nbxmax
         max_boxlz(kk)=boxlz(kk)
      enddo
      

      if (lrdf) then
         if ( lintra ) then
            nskip = 0.5d0
         else
            nskip = 1.5d0
         endif

      endif

      lstretch = .false.

*************************************************************************
C  THIS PART SHOULD GO TO READDAT

!      if (lcharge) then
!         write(6,*) 'input minimum charge for dist (qmin)'
!         read(5,*) qmin
!         write(6,*) 'input maximum charge for dist (qmax)'
!         read(5,*) qmax
!         write(6,*) 'input qbins must be less than',qbinmax
!         read(5,*) qbins
!         if ( qmax .lt. qmin ) stop 'qmax cannont be less than qmin'
!         qdiff = qmax - qmin
!         qstep = qdiff/dble(qbins)
!         write(6,*) 'input additional disp for charges (0.0 for none)'
!         read(5,*) qdisp
!      endif
!      write(6,*) 'do you want the average dipole moments?'
!      read(5,*) ldipole
!      if ( ldipole ) then
!         write(6,*) 'number of blocks for dipole?'
!         read(5,*) qblock
!      endif
*************************************************************************



      do kk = 1,nbox
         boxmin(kk) = 1000.0d0
      enddo

C  end to end vector probability distribution

      if(lete) then
         do imolty = i,nmolty
           max_length(imolty)=0.0d0
           do i = 1,nunit(imolty)
               do j = 1,invib(imolty,i)
                  if(i<ijvib(imolty,i,j)) then
                     max_length(imolty)=max_length(imolty)
     &                                  +brvib(itvib(imolty,i,j))
                  endif 
               enddo
           enddo
         enddo  
     
         do imolty=1,nmolty
           i = (dint((max_length(imolty))/bin_width)+1) 
           if(i.gt.nbinmax_ete) then
	      write(6,*) 'number of bins greater than nbinmax_ete',i,
     &                nbinmax_ete
              stop 'choose larger nbinmax_ete in control.inc'
           endif 
          enddo

c      initializing the arrays

         do kk= 1,nbxmax
          do i=1,nmolty
            do j= 1,nbinmax_ete
               end_to_end(kk,i,j)=0.0d0
            enddo
          enddo   
         enddo

      endif 

!      print*, '*************'
!      print*, '*************'
!      print*,  max_length(1)
!      print*,  max_length(2)
!      print*,  max_length(3)
!      print*, '*************'
!      print*, '*************'


      if(lrhoz) then


         do kk=1,nbox
            i=(dint(boxlz(kk)/bin_width)+1)
            if(i.gt.nbinmax_ete) then
               write(6,*) 'number of bins greater than nbinmax_ete',i,
     &                nbinmax_ete
              stop 'choose larger nbinmax_ete in control.inc'
           endif
         enddo 
        
         do kk=1,nbox
          do i = 1,nmolty 
            do j = 1, nbinmax_ete
               bigrhoz(kk,i,j)=0.0d0
            enddo
            do j=-nbinmax_ete,nbinmax_ete
               bigboxcom_rhoz(kk,i,j)=0.0d0
            enddo
          enddo
         enddo
      endif 




      if ( lbend ) then
c     --- compute ang_bin_size
         ang_bin_size = onepi/dble(ang_bin_max)
c        --- figure out how many angles exist and assign them numbers
         do imolty = 1,nmolty
            bend = 0
            do ii = 1,nunit(imolty)
               do iv = 1,invib(imolty,ii)
                  iuvib = ijvib(imolty,ii,iv)
c                 --- determine whether the beads connected to this unit
c                 --- are of higher index than ii
                  do iuv = 1,invib(imolty,iuvib)
                     iutest = ijvib(imolty,iuvib,iuv)
                     if ( iutest .gt. ii ) then
                        bend = bend + 1
                        angle_1(imolty,bend) = ii
                        angle_2(imolty,bend) = iuvib
                        angle_3(imolty,bend) = iutest
                     endif
                  enddo
               enddo
            enddo
            angle_num(imolty) = bend
            if(bend.gt.angle_max) then
               write(6,*) 'number of bends greater than angle_max',
     &               'molecule type', imolty,' bends', bend   
               stop 'choose a larger angle_max in control.inc'
             endif
         enddo
      endif


      if ( lgvst ) then
c        --- compute tor_bin_size
         tor_bin_size = 360.0d0/dble(tor_bin_max)
c        --- figure out how many torsions exist and assign them numbers
         do imolty = 1,nmolty
            torsion = 0
            do ii = 1,nunit(imolty)
               do itor = 1,intor(imolty,ii)
                  iutor = ijtor4(imolty,ii,itor)
c                 --- determine whether final bead connected to this unit
c                 --- is of higher index than ii
                  if ( iutor .gt. ii ) then
                     torsion = torsion + 1
                     tor_1(imolty,torsion) = ii
                     tor_2(imolty,torsion) = ijtor2(imolty,ii,itor)
                     tor_3(imolty,torsion) = ijtor3(imolty,ii,itor)
                     tor_4(imolty,torsion) = iutor
                     tor_code(ii,iutor) = torsion
                  endif
               enddo
            enddo
            tor_num(imolty) = torsion
            if(torsion.gt.tor_max) then
               write(6,*) 'number of torsions greater than tor_max',
     &               'molecule type', imolty,'torsions', torsion
               stop 'choose a larger tor_max in control.inc*** requires
     &            memory 3**torsion so choose it equal to torsions'
             endif
         enddo
      endif


c --- initialize the arrays

      if ( lbend ) then
         do kk = 1,nbox
            do imolty = 1,nmolty
               do bend = 1,angle_num(imolty)
                  do bin = 1,ang_bin_max
                     angle_bin(kk,imolty,bend,bin) = 0.0d0
                     angle_tot(kk,imolty,bend) = 0.0d0
                  enddo
               enddo
            enddo
         enddo
      endif

      if ( lgvst ) then
         do kk = 1,nbox
            do imolty = 1,nmolty
               do torsion = 1,tor_num(imolty)
                  do bin = 1,tor_bin_max
                     tor_bin(kk,imolty,torsion,bin) = 0.0d0
                     tor_tot(kk,imolty,torsion) = 0.0d0
                  enddo

               enddo
            enddo
         enddo

         do xx = 1,ntmax
            do zzz = 1,nbox
               do uu = 1,tor_max
                  gdefect(zzz,xx,uu) = 0.0d0
                  btrans(zzz,xx,uu) = 0.0d0
                  bg_plus(zzz,xx,uu) = 0.0d0
                  bg_minus(zzz,xx,uu) = 0.0d0
               enddo
               gdefect(zzz,xx,uu+1) = 0.0d0
            enddo
         enddo

      endif

      if ( lrdf ) then
         dummy = ntdifmx*ntdifmx*ntmax*ntmax
         do ntij = 1,dummy
            do kk=1,nbox
               nnone(kk,ntij) = 0.0d0
               binadj = (kk-1)*dummy + ntij
               do bin = 1,nbin
                  biganalhist(binadj,bin) = 0.0d0
                  shell(binadj,bin,1) = 0.0d0
                  shell(binadj,bin,2) = 0.0d0
               enddo
            enddo
         enddo

         dummy = ntmax*ntmax
         do ntij = 1,dummy
            do kk=1,nbox
               comnone(kk,ntij) = 0.0d0
               binadj = (kk-1)*dummy + ntij
               do bin = 1,nbin
                  combiganalhist(binadj,bin) = 0.0d0
                  comshell(binadj,bin,1) = 0.0d0
                  comshell(binadj,bin,1) = 0.0d0 
               enddo
            enddo
         enddo
      endif

!      if ( lcharge ) then
!         do kk = 1, nbox
!            do imol = 1,ntdifmx
!               do iunit = 1,numax
!                  qqcode = numax*(imol-1) + iunit
!                  do qbin = 1,qbins
!                     qanalhist(qqcode,kk,qbin) = 0.0d0
!                     qcount(qqcode,kk) = 0.0d0
!                  enddo
!               enddo
!            enddo
!         enddo
!      endif
      
!      if ( ldipole ) then
!         do kk = 1,nbox
!            do imol = 1,nmolty
!               dipole(imol,kk) = 0.0d0
!               dicount(imol,kk) = 0.0d0
!               diprev(imol,kk) = 0.0d0
!               dcprev(imol,kk) = 0.0d0
!            enddo
!         enddo
!      endif
      return
      endif 




c     cycle through all of the frames and calculate the radial dist func

      if(switch.eq.1) then
         
         do kk=1,nbox
             if(boxlz(kk).gt.max_boxlz(kk)) then
                max_boxlz(kk)=boxlz(kk)
             endif
         enddo

         nframe = nframe + 1.0d0
         
         do kk = 1,nbox
           if(lsolid(kk).and.(.not.(lrect(kk)))) then
                 avolume(kk) = hmat(kk,1)*hmat(kk,5)*hmat(kk,9)
	   else      
                 avolume(kk) = boxlx(kk)*boxly(kk)*boxlz(kk)
 	  endif 
         enddo

c     initialize frame specific arrays
         do ii = 1,nmolty
            do gg = 1,nbox
               cmolec(ii,gg) = 0.0d0
               do g = 1,nhere
                  dummy = nhere*(ii-1) + g
                  count(gg,dummy) = 0.0d0
               enddo
            enddo
         enddo
         
         do ichain = 1,nchain
            imolty= moltyp(ichain)
            tempbx = nboxi(ichain)
            cmolec(imolty,tempbx) = cmolec(imolty,tempbx) + 1.0d0
            do z = 1, nunit(imolty)
               dummy = nhere*(imolty-1)+decode(ntype(imolty,z))
               count(tempbx,dummy) = count(tempbx,dummy) + 1.0d0
            enddo
         enddo

!         if (lcharge) then
!            do kk = 1,nbox
!               do i = 1,nchain
!                  if ( nboxi(i) .eq. kk ) then
!                     imolty = moltyp(i)
!                     do ii = 1,nunit(imolty)
!                        qdummy = qqu(i,ii)
!                        if ( qdummy .lt. qmin .or. qdummy .gt. qmax) then
!                           write(6,*) 'q value out of range'
!                           write(6,*) 'frame,i,ii,q',k,i,ii,qdummy
!                           stop
!                        endif
!                        qdummy = qdummy - qmin
!                        qbin = dint(qdummy/qstep) + 1
!                        qqcode = numax*(imolty-1) + ii
!                        qanalhist(qqcode,kk,qbin) = 
!     &                       qanalhist(qqcode,kk,qbin) + 1.0d0
!                        qcount(qqcode,kk) = qcount(qqcode,kk) + 1.0d0
!                     enddo
!                  endif
!               enddo
!            enddo
!         endif

!         if ( ldipole ) then
!            do i = 1,nchain
!               kk = nboxi(i)
!               imolty = moltyp(i)
!               dipx = 0.0d0
!               dipy = 0.0d0
!               dipz = 0.0d0
!               do iunit = 1,nunit(imolty)
!                  qdummy = qqu(i,iunit)
!                  dipx = dipx + qdummy*rxu(i,iunit)
!                  dipy = dipy + qdummy*ryu(i,iunit)
!                  dipz = dipz + qdummy*rzu(i,iunit)
!               enddo
!               dipole(imolty,kk) = dipole(imolty,kk) 
!     &              + dsqrt(dipx*dipx + dipy*dipy + dipz*dipz)
!               dicount(imolty,kk) = dicount(imolty,kk) + 1.0d0
!            enddo
!
!            if ( mod(k,block) .eq. 0 ) then
!               nblock = nblock + 1
!               do kk = 1,nbox
!                  do imol = 1,nmolty
!                     if ( dicount(imolty,kk) .gt. 0.5d0 ) then
!                        dipblk(imol,kk,nblock) = 
!     &                       (dipole(imolty,kk) - diprev(imolty,kk))/
!     &                       (dicount(imolty,kk) - dcprev(imolty,kk))
!                     endif
!                     diprev(imolty,kk) = dipole(imolty,kk)
!                     dcprev(imolty,kk) = dicount(imolty,kk)
!                  enddo
!               enddo
!            endif
!
!         endif
        
!   Compute end to end vector distribution


        if(lete) then
      
           do kk=1,nbox
              do imolty=1,nmolty
                  do i=1,nbinmax_ete
                      vec_hist(kk,imolty,i) =0.0d0
                  enddo
               enddo
            enddo      

           do i=1,nchain
              kk = nboxi(i)
              imolty = moltyp(i)
 
              rxuij = rxu(i,1)- rxu(i,nunit(imolty))
              ryuij = ryu(i,1)- ryu(i,nunit(imolty))
              rzuij = rzu(i,1)- rzu(i,nunit(imolty))

              ruijsq = rxuij*rxuij + ryuij*ryuij+rzuij*rzuij
              
              ruij   = sqrt(ruijsq)
            
              vec_hist(kk,imolty,(dint(ruij/bin_width)+1))=
     &         vec_hist(kk,imolty,(dint(ruij/bin_width)+1))+1

           enddo
 
           do kk=1,nbox
             do imolty=1,nmolty
                 xx =  (dint(max_length(imolty)/bin_width)+1)
                do bin=1,xx 
                  vec_hist(kk,imolty,bin)=vec_hist(kk,imolty,bin)/
     &                             ncmt(kk,imolty)
                
                  end_to_end(kk,imolty,bin)=end_to_end(kk,imolty,bin)+
     &                     vec_hist(kk,imolty,bin)
                enddo    
             enddo
           enddo 

        endif 

        if(lrhoz) then

         do kk=1,nbox
           do i = 1,nmolty
             do j = 1, nbinmax_ete
                 rhoz(kk,i,j)=0.0d0
             enddo
             do j=-nbinmax_ete,nbinmax_ete
                 boxcom_rhoz(kk,i,j)=0.0d0
             enddo
           enddo
         enddo
           
!! CALCULATE CENTER OF MASS OF EACH BOX

         do kk=1,nbox
            tempzcm(kk)=0.0d0 
            tempmasst(kk) = 0.0d0
            do i =1,nchain
               if(nboxi(i).eq.kk) then
                   tempzcm(kk)   = tempzcm(kk)+mass(moltyp(i))*zcm(i) 
                   tempmasst(kk) = tempmasst(kk)+mass(moltyp(i))
               endif
            enddo
            tempzcm(kk)=tempzcm(kk)/tempmasst(kk)  
         enddo
         
         do i=1,nchain
            kk=nboxi(i)
            imolty=moltyp(i) 
            rhoz(kk,imolty,dint(zcm(i)/bin_width)+1)=rhoz(kk,imolty,
     &      dint(zcm(i)/bin_width)+1) +1

            rzuij=zcm(i)-tempzcm(kk)
            boxcom_rhoz(kk,imolty,floor(rzuij/bin_width))=boxcom_rhoz( 
     &         kk,imolty,floor(rzuij/bin_width))+1          
         enddo  

         do kk=1,nbox
            do i=1,nmolty
               xx = (dint(max_boxlz(kk)/bin_width))+1
               do bin=1,xx
                 slab_vol = hmat(kk,1)*hmat(kk,5)*bin_width
                 rhoz(kk,imolty,bin)=rhoz(kk,imolty,bin)/slab_vol
                 bigrhoz(kk,imolty,bin)=bigrhoz(kk,imolty,bin)+
     &                        rhoz(kk,imolty,bin)
               enddo
              
               do bin=-xx,xx
                 boxcom_rhoz(kk,imolty,bin)=boxcom_rhoz(kk,imolty,bin)/
     &                        slab_vol
                 bigboxcom_rhoz(kk,imolty,bin)=bigboxcom_rhoz(kk,imolty,
     &           bin)+boxcom_rhoz(kk,imolty,bin)
               enddo 
            enddo
         enddo 
        endif

 
       if (lrdf) then
c        --- compute the radial distribution analhistograms for this frame
         do kk = 1,nbox
            bx = boxlx(kk)
            hbx = bx/2.0d0
            by = boxly(kk)
            hby = by/2.0d0
            bz = boxlz(kk)
            hbz = bz/2.0d0
            if ( hbx .lt. boxmin(kk) ) boxmin(kk) = hbx
            rcutsq = rcut*rcut
            binstep = rcut/dble(nbin)
c           initiallize analhist
            do g = 1,nmolty*nmolty*nhere*nhere
               do gg = 1,nbin
                  analhist(g,gg) = 0.0d0
               enddo
            enddo
c           initiallize comanalhist
            do g = 1,nmolty*nmolty
               do gg = 1,nbin
                  comanalhist(g,gg) = 0.0d0
               enddo
            enddo

            do i = 1, nchain
 
c  check if i is in relevant box 
               if ( nboxi(i) .eq. kk ) then
                  imolty = moltyp(i)
c     --- loop over all beads ii of chain i 
                  do ii = 1, nunit(imolty)
                     ntii = nhere*(imolty-1) + decode(ntype(imolty,ii))
                     rxui = rxu(i,ii)
                     ryui = ryu(i,ii)
                     rzui = rzu(i,ii)

c     --- loop over all chains j with j>=i 
                     if ( lintra ) then
                        istart = i
                     else
                        istart = i+1
                     endif
                     do j = istart, nchain
c                       ### check for simulation box ###
                        if ( nboxi(j) .eq. kk ) then
                           jmolty = moltyp(j)
c                       --- loop over all beads jj of chain j 
                           if (i .eq. j) then 
                              jstart = ii+1
                           else
                              jstart = 1
                           endif
                           do jj = jstart, nunit(jmolty)
                        
                              ntjj = nhere*(jmolty-1) + 
     &                             decode( ntype(jmolty,jj) )
                              if ( ntii .gt. ntjj ) then
                                 ntij = (ntjj-1)*nhere*nmolty + ntii
                              else
                                 ntij = (ntii-1)*nhere*nmolty + ntjj
                              endif

                              rxuij = rxui - rxu(j,jj)
                              ryuij = ryui - ryu(j,jj)
                              rzuij = rzui - rzu(j,jj)

c *** minimum image the pair separations ***
                              if ( rxuij .gt. hbx ) then
                                 rxuij=rxuij-bx
                              else
                                 if (rxuij.lt.-hbx) rxuij=rxuij+bx
                              endif

                              if ( ryuij .gt. hby ) then
                                 ryuij=ryuij-by
                              else
                                 if (ryuij.lt.-hby) ryuij=ryuij+by
                              endif

                              if (rzuij.gt.hbz) then
                                 rzuij=rzuij-bz
                              else
                                 if (rzuij.lt.-hbz) rzuij=rzuij+bz
                              endif

                              ruijsq = rxuij*rxuij + ryuij*ryuij 
     &                             + rzuij*rzuij

                              if (ruijsq .lt. rcutsq) then
                                 ruij = dsqrt(ruijsq)

                                 bin = dint(ruij/binstep) + 1 
                                 analhist(ntij,bin) = analhist(ntij,
     &                                                 bin) + 1.0d0
                              endif

                           enddo
                        endif
                     enddo
                  enddo
               endif
            enddo

c     normalize the analhistogram and add it to the big analhistogram
            
            do ii = 1,nmolty
               do xx  = 1,nhere
                  ntii = nhere*(ii-1) + xx
                  rho(ntii) = (4.0d0 * 3.1415d0 * count(kk,ntii)) 
     &                 / (3.0d0 * avolume(kk) )
               enddo
            enddo
            
            do ii = 1,nmolty
               do xx = 1,nhere
                  ntii = nhere*(ii-1)+xx
                  numxx = count(kk,ntii)
                  
                  do jj = 1,nmolty
                     do yy = 1,nhere
                        ntjj = nhere*(jj-1)+yy
                        numyy = count(kk,ntjj)
                        ntij = (ntii-1)*nhere*nmolty + ntjj

                        ntji = (ntjj-1)*nhere*nmolty + ntii

                        binadj = (kk-1)*nhere*nhere*nmolty*nmolty 
     &                       + ntij
                        if ( ntii .eq. ntjj ) then
                           const = rho(ntjj)/2.0d0
                        else
                           const = rho(ntjj)
                        endif
                        
                        shlsumx = 0.0d0
                        shlsumy = 0.0d0  
c     check to see if there are enough molecules to have any interactions
c     if not then lskip is set to true

                        lskip = .false.
                        if ( ntii .eq. ntjj ) then
                           if ( cmolec(ii,kk) .lt. nskip ) then
                              lskip = .true.
                           elseif ( numxx .lt. 0.5d0 ) then
                              lskip =.true.
                           endif
                        else
                           if ( cmolec(ii,kk)*cmolec(jj,kk) 
     &                          .lt. 0.5d0 .or. numxx*numyy 
     &                          .lt. 0.5d0) then
                              lskip = .true.
                           endif
                        endif
                        
                        if ( .not. lskip  ) then

                           do bin = 1,nbin

                              if ( ntii .eq. ntjj ) then
                                    shlsumx = shlsumx
     &                                   + 2.0d0*analhist(ntij,bin)
     &                                   /numxx
                                    shell(binadj,bin,1) =
     &                                   shell(binadj,bin,1)
     &                                   + shlsumx
                                 else
                                    shlsumx = shlsumx
     &                                   + analhist(ntij,bin)/numxx
                                    shlsumy = shlsumy
     &                                   + analhist(ntij,bin)/numyy
                                    shell(binadj,bin,1) =
     &                                   shell(binadj,bin,1)
     &                                   + shlsumx
                                    shell(binadj,bin,2) =
     &                                   shell(binadj,bin,2)
     &                                   + shlsumy
                                 endif

                              
                              rlower = dble(bin-1)*binstep
                              rupper = rlower + binstep
                              nideal = const*(rupper**3 -rlower**3)
                              xxideal = numxx*nideal

                              biganalhist(binadj,bin) = 
     &                             biganalhist(binadj,bin) 
     &                             + analhist(ntij,bin)/xxideal


                           enddo
                           
                        else
                           nnone(kk,ntij) = nnone(kk,ntij) + 1.0d0
                        endif
                     enddo
                  enddo
               enddo
            enddo

c           --- Center of Mass Radial Distribution Functions
            do i = 1, nchain
 
c  check if i is in relevant box 
               if ( nboxi(i) .eq. kk ) then
                  imolty = moltyp(i)

                  rxui = xcm(i)
                  ryui = ycm(i)
                  rzui = zcm(i)

c     --- loop over all chains j with j>=i 
                  istart = i+1
                  do j = istart, nchain
c                       ### check for simulation box ###
                     if ( nboxi(j) .eq. kk ) then
                        jmolty = moltyp(j)

                        if ( imolty .gt. jmolty ) then
                           ntij = (jmolty-1)*nmolty + imolty
                        else
                           ntij = (imolty-1)*nmolty + jmolty
                        endif

                        rxuij = rxui - xcm(j)
                        ryuij = ryui - ycm(j)
                        rzuij = rzui - zcm(j)

c *** minimum image the pair separations ***
                        if ( rxuij .gt. hbx ) then
                           rxuij=rxuij-bx
                        else
                           if (rxuij.lt.-hbx) rxuij=rxuij+bx
                        endif

                        if ( ryuij .gt. hby ) then
                           ryuij=ryuij-by
                        else
                           if (ryuij.lt.-hby) ryuij=ryuij+by
                        endif
                        
                        if (rzuij.gt.hbz) then
                           rzuij=rzuij-bz
                        else
                           if (rzuij.lt.-hbz) rzuij=rzuij+bz
                        endif
                        
                        ruijsq = rxuij*rxuij + ryuij*ryuij 
     &                       + rzuij*rzuij
                        
                        if (ruijsq .lt. rcutsq) then
                           ruij = dsqrt(ruijsq)
                           
                           bin = dint(ruij/binstep) + 1 
                           comanalhist(ntij,bin) = comanalhist(ntij,bin) 
     &                          + 1.0d0
                        endif

                     endif
                  enddo
               endif
            enddo

c     normalize the COM analhistogram and add it to the big analhistogram
            
            do ii = 1,nmolty
               rho(ii) = (4.0d0 * 3.1415d0 * cmolec(ii,kk)) 
     &                 / (3.0d0 * avolume(kk) )
            enddo
            
            do ii = 1,nmolty
               ntii = ii
               numxx = cmolec(ntii,kk)
               do jj = 1,nmolty

                  ntjj = jj
                  numyy = cmolec(ntjj,kk)
                  
                  ntij = (ntii-1)*nmolty + ntjj
                  ntji = (ntjj-1)*nmolty + ntii

                  binadj = (kk-1)*nmolty*nmolty + ntij
                  
                  if ( ntii .eq. ntjj ) then
                     const = rho(ntjj)/2.0d0
                  else
                     const = rho(ntjj)
                  endif
                  
                  shlsumx = 0.0d0
                  shlsumy = 0.0d0

c     check to see if there are enough molecules to have any interactions
c     if not then lskip is set to true

                  lskip = .false.
                  
                  if ( ntii .eq. ntjj ) then
                     if ( cmolec(ii,kk) .lt. nskip ) then
                        lskip = .true.
                     endif
                  elseif ( cmolec(ii,kk)*cmolec(jj,kk) 
     &                    .lt. 0.5d0 .or. numxx*numyy 
     &                    .lt. 0.5d0) then
                     lskip = .true.
                  endif
                           
                  if ( .not. lskip  ) then

                     do bin = 1,nbin

	               if ( ntii .eq. ntjj ) then
                           shlsumx = shlsumx
     &                          + 2.0d0*comanalhist(ntij,bin)
     &                          /numxx
                           comshell(binadj,bin,1) =
     &                          comshell(binadj,bin,1)+ shlsumx
                        else
                           shlsumx = shlsumx
     &                          + comanalhist(ntij,bin)/numxx
                           shlsumy = shlsumy
     &                          + comanalhist(ntij,bin)/numyy
                           comshell(binadj,bin,1) =
     &                          comshell(binadj,bin,1)+ shlsumx
                           comshell(binadj,bin,2) =
     &                          comshell(binadj,bin,2)+ shlsumy
                        endif


                        rlower = dble(bin-1)*binstep
                        rupper = rlower + binstep
                        nideal = const*(rupper**3 -rlower**3)
                        xxideal = numxx*nideal

                        if (jj .ge. ii) then

                           combiganalhist(binadj,bin) =  
     &                          combiganalhist(binadj,bin) + 
     &                          comanalhist(ntij,bin)/xxideal

                        endif

                     enddo
                              
                  else
                     comnone(kk,ntij) = comnone(kk,ntij) + 1.0d0
                  endif
               enddo
            enddo

         enddo
        endif

cccccccccccccccccccccccccc

c     analyse the torsional angles

        if (lstretch .or. lbend .or. lgvst ) then
           do kk=1,nbox
              do i = 1, nchain
          
               imolty = moltyp(i)

c  check if i is in relevant box 
               if ( nboxi(i) .eq. kk ) then

c                 --- calculate all bonds vectors and lengths

                  do ii = 1, nunit(imolty)
                     rxui=rxu(i,ii)
                     ryui=ryu(i,ii)
                     rzui=rzu(i,ii)
                     do iivib = 1, invib(imolty,ii)
                        jj = ijvib(imolty,ii,iivib)
                        xvec(ii,jj) = rxu(i,jj) - rxui
                        yvec(ii,jj) = ryu(i,jj) - ryui
                        zvec(ii,jj) = rzu(i,jj) - rzui
                        distij(ii,jj) = dsqrt( xvec(ii,jj)**2
     +                       + yvec(ii,jj)**2 + zvec(ii,jj)**2 )
                     enddo
                  enddo
                  
c - stretching -
                  if ( lstretch ) then
                  endif

c - bending -
                  if ( lbend ) then
                     do bend = 1,angle_num(imolty)
                        j = angle_1(imolty,bend)
                        ip1 = angle_2(imolty,bend)
                        ip2 = angle_3(imolty,bend)

                        thetac = ( xvec(ip1,j)*xvec(ip1,ip2) +
     +                       yvec(ip1,j)*yvec(ip1,ip2) +
     +                       zvec(ip1,j)*zvec(ip1,ip2) ) /
     +                       ( distij(ip1,j)*distij(ip1,ip2) )
                        theta = dacos(thetac)
                        angle_avg(kk,imolty,bend) = 
     &                       angle_avg(kk,imolty,bend) + theta
                        bin = dint(theta/ang_bin_size)+1
                        angle_bin(kk,imolty,bend,bin) = 
     &                       angle_bin(kk,imolty,bend,bin) + 1.0d0
                        angle_tot(kk,imolty,bend) = 
     &                       angle_tot(kk,imolty,bend)+ 1.0d0
                     enddo
                  endif

c - torsions -
                  if ( lgvst ) then
c  molecule with dihedral potenials 
                     gaudef = 1
                     dum = 0
                     do torsion = 1, tor_num(imolty)
                        j = tor_1(imolty,torsion)
                        ip1 = tor_2(imolty,torsion)
                        ip2 = tor_3(imolty,torsion)
                        ip3 = tor_4(imolty,torsion)

c              ***  calculate cross product d_a x d_a-1  ***
                        xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     +                       zvec(ip1,j) * yvec(ip1,ip2)
                        yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     +                       xvec(ip1,j) * zvec(ip1,ip2)
                        zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     +                       yvec(ip1,j) * xvec(ip1,ip2)
c              ***  calculate cross product d_a-1 x d_a-2 ***
                        xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     +                       zvec(ip1,ip2) * yvec(ip3,ip2)
                        ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     +                       xvec(ip1,ip2) * zvec(ip3,ip2)
                        za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     +                       yvec(ip1,ip2) * xvec(ip3,ip2)
                              
c     *** calculate lengths of cross products ***
                        daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                        da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
c              *** calculate dot product of cross products ***
                        dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                        thetac = - dot / ( daa1 * da1a2 )
                        theta = dacos(thetac)

c              *** calculate cross product of cross products ***
                        xcc = yaa1*za1a2 - zaa1*ya1a2
                        ycc = zaa1*xa1a2 - xaa1*za1a2
                        zcc = xaa1*ya1a2 - yaa1*xa1a2
c              *** calculate scalar triple product ***
                        tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2) 
     &                       + zcc*zvec(ip1,ip2) 
                        if ( tcc .lt. 0.0d0 ) theta = - theta

c              *** convert theta to degrees ***
                        theta = theta*(180.0d0/onepi)

c                       --- bin the torsion --
                        bin = dint( (theta+180.0d0)/tor_bin_size)+1
                        tor_bin(kk,imolty,torsion,bin) = 
     &                       tor_bin(kk,imolty,torsion,bin)+1.0d0
                        tor_tot(kk,imolty,torsion) = 
     &                       tor_tot(kk,imolty,torsion) + 1.0d0

                        if ( theta .lt. -60.0d0 ) then
c                          --- gauch minus
                           g_minus(kk,imolty)=g_minus(kk,imolty) + 1.0d0
                           gaudef = gaudef + 1
                           bg_minus(kk,imolty,torsion) = 
     &                          bg_minus(kk,imolty,torsion) + 1.0d0
                           dum = dum + 0 * 3**(torsion-1)
                        elseif ( theta .lt. 60.0d0 ) then
c                          --- trans
                           trans(kk,imolty) = trans(kk,imolty)+1.0d0
                           btrans(kk,imolty,torsion) = 
     &                          btrans(kk,imolty,torsion) + 1.0d0
                           dum = dum + 1 * 3**(torsion-1)
                        else
c                          --- gauch plus
                           g_plus(kk,imolty)=g_plus(kk,imolty) + 1.0d0
                           gaudef = gaudef + 1
                           bg_plus(kk,imolty,torsion) = 
     &                          bg_plus(kk,imolty,torsion) + 1.0d0
                           dum = dum + 2 * 3**(torsion-1)
                        endif

                     enddo

                     gdefect(kk,imolty,gaudef) = 
     &                    gdefect(kk,imolty,gaudef) + 1.0d0 
                     gdefect(kk,imolty,tor_max+1) = 
     &                    gdefect(kk,imolty,tor_max+1) + 1.0d0
                     pattern(kk,imolty,dum) = pattern(kk,imolty,dum) 
     &                    + 1.0d0
                  endif
               endif
              enddo
           enddo
        endif

      return
      endif

      if(switch.eq.2) then

      if(lete) then

         open (unit=141,FILE="end2end_box1",STATUS="unknown")
         open (unit=142,FILE="end2end_box2",STATUS="unknown")
         open (unit=143,FILE="end2end_box3",STATUS="unknown")

         do kk=1,nbox
           do imolty=1,nmolty
              xx =  (dint(max_length(imolty)/bin_width)+1)
              write(140+kk,*) 'molecule type', imolty
              write(140+kk,*) 
              do bin=1,xx
                  rxuij=bin_width*(dble(bin)-0.5d0)
                  write(140+kk,*) rxuij,end_to_end(kk,imolty,bin)/
     &                                    nframe
              enddo
           enddo
         enddo              
        
         close(141)
         close(142)
         close(143)

      endif


      if(lrhoz) then

         open (unit=151,FILE="rhoz_box1",STATUS="unknown")
         open (unit=152,FILE="rhoz_box2",STATUS="unknown")
         open (unit=153,FILE="rhoz_box3",STATUS="unknown")
  
         open (unit=154,FILE="comrhoz_box1",STATUS="unknown")
         open (unit=155,FILE="comrhoz_box2",STATUS="unknown")
         open (unit=156,FILE="comrhoz_box3",STATUS="unknown")

         do kk=1,nbox
           do imolty=1,nmolty
              xx =  (dint(max_boxlz(kk)/bin_width)+1)
              write(150+kk,*) 'molecule type',imolty
              write(153+kk,*) 'molecule type',imolty
              write(153+kk,*) 
              write(150+kk,*)
              do bin=1,xx
                  rxuij=bin_width*(dble(bin)-0.5d0)
                  write(150+kk,*) rxuij,bigrhoz(kk,imolty,bin)/
     &                                    nframe
              enddo
              do bin=-xx,xx
                 if (bin.eq.0) then
                     rxuij=0.0
                 else
                     rxuij=bin_width*(dble(bin)+0.5d0) 
                 endif
                 write(153+kk,*) rxuij,bigboxcom_rhoz(kk,imolty,bin)/
     &                           nframe
              enddo 
           enddo
         enddo

         close(151)
         close(152)
         close(153)
         close(154)
         close(155)
         close(156)

      endif


c     whew, now that is done we just have to divide out the number of
c     frames in biganalhist


      if (lrdf) then
         open (unit=101,FILE="beadrdf_box1",STATUS="unknown")
         open (unit=102,FILE="beadrdf_box2",STATUS="unknown")
         open (unit=103,FILE="beadrdf_box3",STATUS="unknown")

         open (unit=104,FILE="comrdf_box1",STATUS="unknown")
         open (unit=105,FILE="comrdf_box2",STATUS="unknown")
         open (unit=106,FILE="comrdf_box3",STATUS="unknown")

         open (unit=111,FILE="beadnum_box1",STATUS="unknown")
         open (unit=112,FILE="beadnum_box2",STATUS="unknown")
         open (unit=113,FILE="beadnum_box3",STATUS="unknown")

         open (unit=114,FILE="comnum_box1",STATUS="unknown")
         open (unit=115,FILE="comnum_box2",STATUS="unknown")
         open (unit=116,FILE="comnum_box3",STATUS="unknown")

c         write(6,*)
c         write(6,*) 'bead bead radial distribution functions in *gor*'
c         write(6,*)
c         write(6,*) 'avg. number of beads vs. shell size in *num*'
         
!         print*, 'nhere', nhere
         
      do kk = 1,nbox
        binstep = rcut/dble(nbin)
         do ii = 1,nmolty
            do xx = 1,nhere
               ntii = nhere*(ii-1)+xx
               do jj = 1,nmolty
                  do yy = 1,nhere
                   ntjj = nhere*(jj-1)+yy
                   if ( ntii .le. ntjj ) then 
                      ntij = (ntii-1)*nhere*nmolty + ntjj
                   
                      binadj = (kk-1)*nhere*nhere*nmolty*nmolty + ntij

                      aframe = nframe-nnone(kk,ntij)

                      if ( aframe .gt. 0.5d0 ) then

                         write(100+kk,50)0.0d0, 0.0d0,ii,beadtyp(xx)
     &                        ,jj,beadtyp(yy) 

                         write(110+kk,50) 0.0d0,0.0d0,ii,beadtyp(xx)
     &                     ,jj,beadtyp(yy) 
 50                   format(2f7.2,4i5)
                      do bin = 1,nbin

                            rxuij =  binstep*(dble(bin)-0.5d0)
 
                            biganalhist(binadj,bin) = biganalhist(binadj
     &                                        ,bin)/aframe

                            write(100+kk,*) rxuij,biganalhist(binadj,
     &                                                     bin)

                             shell(binadj,bin,1) = shell(binadj,bin,1)
     &                        /aframe
                         write(110+kk,*) rxuij,shell(binadj,bin,1)
                         
                      enddo
                      write(100+kk,*)
	              write(110+kk,*)
 
                      if (ntii .ne. ntjj) then
                          write(110+kk,50) 0.0d0,0.0d0,jj,beadtyp(yy)
     &                           ,ii,beadtyp(xx)
                          do bin = 1,nbin
                              rxuij =  binstep*(dble(bin)-0.5d0)
                              shell(binadj,bin,2) =
     &                             shell(binadj,bin,2) /aframe
                              write(110+kk,*) rxuij,shell(binadj,bin,2)
                            enddo
                            write(110+kk,*)
                      endif
                      
                   elseif( aframe .lt. -0.5d0 ) then
                      write(6,*) 'aframe',aframe
                      write(6,*) 'nframe',nframe
                      write(6,*) 'nnone(kk,ntij),kk,ntij',
     &                     nnone(kk,ntij),kk,ntij
                      stop 'srewup aframe'
                     endif
                   endif
                 enddo
               enddo
            enddo
         enddo
       enddo

c     --- same thing for the COM rdf
       do kk = 1,nbox
         binstep = rcut/dble(nbin)
         do ii = 1,nmolty
            ntii = ii
            do jj = 1,nmolty

               ntjj = jj
               ntij = (ntii-1)*nmolty + ntjj

               binadj = (kk-1)*nmolty*nmolty + ntij
               aframe = nframe-comnone(kk,ntij)

               if ( aframe .gt. 0.5d0 ) then
                  write(103+kk,51) 0.0d0,ii,jj
                  write(113+kk,51) 0.0d0,0.0d0,ii,jj
 51               format(2f7.2,2i5)
                  do bin = 1,nbin
                     rxuij =  binstep*(dble(bin)-0.5d0)
                     
                     combiganalhist(binadj,bin) = combiganalhist(bina
     &                        dj,bin)/aframe
                     write(103+kk,*) rxuij,combiganalhist(binadj,bin)

                     comshell(binadj,bin,1) = comshell(binadj,bin,1)
     &                    /aframe
                     write(113+kk,*) rxuij,comshell(binadj,bin,1)

                  enddo

                  write(103+kk,*)
                  write(113+kk,*)
                  
                  if(ntii.ne.ntjj) then  
                    write(113+kk,51) 0.0d0,0.0d0,jj,ii
                    do bin = 1,nbin
                       rxuij =  binstep*(dble(bin)-0.5d0)
                       comshell(binadj,bin,2) =
     &                       comshell(binadj,bin,2) /aframe
                       write(113+kk,*) rxuij,comshell(binadj,bin,2)
	            enddo
                    write(113+kk,*)
                  endif    
               elseif( aframe .lt. -0.5d0 ) then
                  write(6,*) 'aframe',aframe
                  write(6,*) 'nframe',nframe
                  write(6,*) 'comnone(kk,ntij),kk,ntij',
     &                 comnone(kk,ntij),kk,ntij
                  stop 'srewup aframe'
               endif
            enddo
         enddo
       enddo

       close(101)
       close(102)
       close(103)
       close(104)
       close(105)
       close(106)
       close(111)
       close(112)
       close(113)
       close(114)
       close(115)
       close(116)
      endif 

!      if ( lcharge ) then
!         write(6,*) 'charge distributions in fort.35+box'
!         do kk = 1, nbox
!            do imol = 1,nmolty
!               do iunit = 1,nunit(imolty)
!                  qqcode = numax*(imol-1) + iunit
!                  qhere = .false.
!                  do qbin = 1,qbins
!                     if ( qcount(qqcode,kk) .gt. 0.5d0 ) then
!                        qhere = .true.
!                        qanalhist(qqcode,kk,qbin) = 
!     &                       qanalhist(qqcode,kk,qbin)/qcount(qqcode,kk)
!                     endif
!                  enddo
!                  if (qhere) then
!                     qdummy = ( (dble(1)-0.5d0)*qstep)+qmin
!                     write(35+kk,*) qdummy,qanalhist(qqcode,kk,1)+qdisp
!     &                    ,imol,iunit
!                     
!                     do qbin = 2,qbins
!                        qdummy = ( (dble(qbin)-0.5d0)*qstep)+qmin
!                        write(35+kk,*) qdummy
!     &                       ,qanalhist(qqcode,kk,qbin)+qdisp
!                     enddo
!                     write(35+kk,*)
!                  endif
!               enddo
!            enddo
!         enddo
!      endif

!      if ( ldipole ) then
!         diconv = (1.602d-19)/((1d10)*(3.336d-30))
!         do kk = 1, nbox
!            write(6,*) 'Dipoles [e A], [D] in Box',kk
!            do imol = 1,nmolty
!               if ( dicount(imol,kk) .gt. 0.5d0 ) then
!                  dipole(imol,kk) = dipole(imol,kk)/dicount(imol,kk)
!               endif
!               write(6,*) 'Moltyp,dipole',imolty
!     &              ,dipole(imol,kk),dipole(imol,kk)*diconv
!            enddo
!         enddo
!
!         if ( block .gt. 1 ) then
!            write(6,*) 'Number of blocks input and used',block,nblock
!            do kk = 1,nbox
!               write(6,*) 'Box     Moltyp   Dipole [D]   Std dev [D]'
!               do imol = 1,nmolty
!                  avera = 0.0d0
!                  do k=1,nblock
!                     avera = avera + dipblk(imol,kk,k)
!                  enddo
!                  avera = avera/nblock
!                  stddev = 0.0d0
!                  do k = 1,nblock
!                     stddev = stddev + (dipblk(imol,kk,k)-avera)**2
!                  enddo
!                  stddev = dsqrt( stddev/dble(nblock-1) )
!                  stddev = stddev*diconv
!                  avera = avera*diconv
!                  write(6,'(i3,4x,i3,4x,2f10.4)') kk,imol,avera,stddev
!               enddo
!            enddo
!         endif

!      endif
               
      if ( lbend ) then
         open (UNIT=121, FILE="bendang_dist_box1",status="unknown")
         open (UNIT=122, FILE="bendang_dist_box2",status="unknown")
         open (UNIT=123, FILE="bendang_dist_box3",status="unknown")
c     --- output the bending angle distributions for each angle in the molecule
c     --- write to 37+box

         do kk = 1,nbox
            do imolty = 1,nmolty
               do bend = 1,angle_num(imolty)
                  total = angle_tot(kk,imolty,bend)
                  if ( total .gt. 0.5d0 ) then
                     value = ang_bin_size/2.0d0
                     degree = value*180.0d0/onepi
                     write(120+kk,*) degree
     &                       ,angle_bin(kk,imolty,bend,bin)/total
     &                    ,'Moltyp ',imolty,' angle '
     &                    ,angle_1(imolty,bend),angle_2(imolty,bend)
     &                    ,angle_3(imolty,bend)
                     do bin = 2,ang_bin_max
                        value = value + ang_bin_size
                        degree = value*180.0d0/onepi
                        write(120+kk,*) degree
     &                       ,angle_bin(kk,imolty,bend,bin)/total
                     enddo
                     write(120+kk,*)
c                    --- output the average angle to the screen
                     value = angle_avg(kk,imolty,bend)/total
                     degree = value*180.0d0/onepi
                     write(6,*) 
     &                    'Moltyp ',imolty,' angle '
     &                    ,angle_1(imolty,bend),angle_2(imolty,bend)
     &                    ,angle_3(imolty,bend),' average ',degree
                     
                  endif
               enddo
            enddo
         enddo                
         close(121)
         close(122)
         close(123)
      endif

      if ( lgvst ) then


         open (unit=131,FILE="tors_frac_box1",status="unknown")
         open (unit=132,FILE="tors_frac_box2",status="unknown")
         open (unit=133,FILE="tors_frac_box3",status="unknown")

         open (unit=134,FILE="torsprob_box1",status="unknown")
         open (unit=135,FILE="torsprob_box2",status="unknown")
         open (unit=136,FILE="torsprob_box3",status="unknown")


c     output the gauch versus trans data for each moltype
         write(6,*) 'moltyp  box  trans      g+        g-       g frac'
         do kk  = 1,nbox
            do imolty = 1,nmolty
               write(130+kk,*) 'Moltype ',imolty
               write(130+kk,41)
               gaudef = g_plus(kk,imolty)+g_minus(kk,imolty)
               ratio = gaudef/(trans(kk,imolty)+gaudef)
               write(6,40)imolty,kk ,trans(kk,imolty),g_plus(kk,imolty)
     &              ,g_minus(kk,imolty),ratio

               do torsion = 1,tor_num(imolty)
                  total = bg_plus(kk,imolty,torsion) + 
     &                 bg_minus(kk,imolty,torsion) +
     &                 btrans(kk,imolty,torsion)

                  if ( total .gt. 0.5d0) then
c                    --- output torsion type fractions 
                     fplus = bg_plus(kk,imolty,torsion)/total
                     fminus = bg_minus(kk,imolty,torsion)/total
                     ftrans = btrans(kk,imolty,torsion)/total
                     write(130+kk,42) tor_1(imolty,torsion)
     &                    ,tor_2(imolty,torsion),tor_3(imolty,torsion)
     &                    ,tor_4(imolty,torsion),fplus,fminus,ftrans

C                    --- output the torsion probablility distributions
                     value = -180.0d0 + (0.5d0*tor_bin_size)
                     write(133+kk,*) value,tor_bin(kk,imolty,torsion,1)
     &                    ,imolty,tor_1(imolty,torsion)
     &                    ,tor_2(imolty,torsion),tor_3(imolty,torsion)
     &                    ,tor_4(imolty,torsion)

                     do bin = 2,tor_bin_max
                        value = value + tor_bin_size
                        write(133+kk,*) value
     &                       ,tor_bin(kk,imolty,torsion,bin)
                     enddo
                     write(133+kk,*)
                  endif
               enddo
            enddo
         enddo
 40      format(i3,5x,i3,3(2x,e8.2),f8.4)
 41      format('Units',8x,'Frac g+',1x,'Frac g-',1x,'Frac t')
 42      format(4(i2,1x),1x,3(f5.3,3x))

         close(131)
         close(132)
         close(133)
         close(134)
         close(135)
         close(136)

!         write(6,*)
!         write(6,*) 'gauch fraction vs. torsion in fort.90+box'
!         write(6,*)
!         write(6,*) 'torsion angle distribution in fort.95+box'
!         write(6,*)
!         write(6,*) 'probabiltiy vs. number of defects per chain'
!         write(6,*) 'shown in fort.5*'
!         write(6,*)
         
         open(unit=137,FILE="Gauchedefects_box1",status="unknown")
         open(unit=138,FILE="Gauchedefects_box2",status="unknown")
         open(unit=139,FILE="Gauchedefects_box3",status="unknown")

         open(unit=124,FILE="decoder",status="unknown")
         open(unit=125,FILE="pattern_box1",status="unknown")
         open(unit=126,FILE="pattern_box2",status="unknown")
         open(unit=127,FILE="pattern_box3",status="unknown")



c     analyse the number of gauch defects per chain
         do kk = 1,nbox
            do imolty = 1,nmolty
               total = gdefect(kk,imolty,tor_max+1)
               if ( total .gt. 0.01d0 ) then
                  write(136+kk,*) 'molecule type ',imolty
                  do torsion = 1,tor_num(imolty)
                     gdefect(kk,imolty,torsion) = 
     &                    gdefect(kk,imolty,torsion)/total
                     write(136+kk,*) torsion,gdefect(kk,imolty,torsion)
                  enddo
               endif

               psum = 0.0d0
               if ( tor_num(imolty) .gt. 0 ) then
c                 --- write out a decoder for the torsions
                  if ( kk .eq. 1 ) then
                     write(124,*) 'Moltyp ',imolty,' torsions'
                     do torsion = 1,tor_num(imolty)
                        write(124,*) torsion,'   ',tor_1(imolty,torsion)
     &                       ,tor_2(imolty,torsion)
     &                       ,tor_3(imolty,torsion)
     &                       ,tor_4(imolty,torsion)
                     enddo
                  endif
                  units = 3**( tor_num(imolty) )
                  do uu = 0,units
                     psum = psum + pattern(kk,imolty,uu)
                  enddo
                  write(124+kk,*) 'Code     Prob.     Torsion pattern'
                  do uu = 0,units-1
                     if ( pattern(kk,imolty,uu) .gt. 0.5d0 ) then
                        decimal = uu
                        power = units
                        torsion = tor_num(imolty)
                        do dummy = torsion,1,-1
                           power = power/3
                           bthree = decimal/power 
                           decimal = decimal - bthree*power
                           patt(dummy) = bthree-1
                        enddo

                        write(124+kk,52) uu,
     &                       pattern(kk,imolty,uu)/psum
     &                       ,(patt(dummy),dummy=1,torsion)
 52                     format(i8,2x,f5.3,2x,20(i2,1x))
                     endif
                  enddo 
               endif

            enddo
         enddo
         close(137)
         close(138)
         close(139)
         close(124)
         close(125)
         close(126)
         close(127)

!         write(6,*) 'patterns of gauch defects in fort.2*'
!         write(6,*) 'transform first number into base 3 to get the'
!         write(6,*) 'pattern of gauch defects where -1=g-, 0=trans 1=g+'
!         write(6,*)
         write(6,*) 'the maximum rcut value that could be used is'
         write(6,*) boxmin
         write(6,*)
         

      endif
      return 
      endif
       
      end
