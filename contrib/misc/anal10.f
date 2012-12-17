      program anal10

c     *******************************************************************
c     *** this program analyses the fort.10 'movie file' and outputs  ***
c     *** the bead-bead radial distribution functions, and some trans ***
c     *** and gauche defect information into other fort.* files       ***
c     *** written by Marcus Martin last modified 2-23-98              ***
c     *******************************************************************

      implicit none

      integer bin,nbin,k,kk,box,nframe,nchain,nhere,z,beadtyp
      integer nmolty,nboxmx,nbox,ntype,ntypmx,decode,tempbx
     &     ,temp,xx,yy,nmax,numax,fstep,ncmt,ntmax,g,gg,nmolec
      integer chnum,imolty,nunit,nboxi,i,ii,moltyp,j,jj
     &     ,jmolty,binadj,ntij,ntii,ntjj,jstart,istart
      integer invib,ijvib,intor,ijtor2,ijtor3,ijtor4,iivib,jjtor
     &     ,ip1,ip2,ip3,gaudef,uu,dum,units,nbinmx,zzz,zz1

      logical lintra,lstretch,lskip,lcharge,lrdf

      double precision volume,rho,rcut,rcutsq,binstep
     &     ,rxu,ryu,rzu,boxlx,boxly,boxlz,rxui,ryui,rzui,xxideal
     &     ,rxuij,ryuij,rzuij,bx,hbx,by,hby,bz,hbz,ruijsq,ruij
     &     ,numxx,numyy,count,const,rlower,rupper,nideal,hist,bighist
      double precision nnone,comnone,aframe,nskip,cmolec,addist
     &     ,comhist,combighist
      double precision xcm,ycm,zcm,qqu,boxmin,shlsumx,shlsumy
     &     ,shell,comshell

c     ntypmx = number of moltyps, ntmax = max number of diff beads in sim
      parameter (nbinmx=100,nboxmx=1,ntypmx=6,nmax=1600,numax=30
     & ,ntmax=8)

c     --- variables used in the bending parts
      logical lbend
      integer angle_max,ang_bin_max,bend,iv,iuvib,iuv,iutest
      double precision angle_bin,angle_tot,angle_num,angle_1,angle_2
     &     ,angle_3,onepi,ang_bin_size,value,total,degree,angle_avg
      parameter (onepi = 3.141592654d0,angle_max = 50,ang_bin_max=360)
      dimension angle_num(ntypmx)
      dimension angle_1(ntypmx,angle_max),angle_2(ntypmx,angle_max)
     &     ,angle_3(ntypmx,angle_max),angle_tot(nboxmx,ntypmx,angle_max)
      dimension angle_avg(nboxmx,ntypmx,angle_max)
      dimension angle_bin(nboxmx,ntypmx,angle_max,ang_bin_max)

c     --- variables used in the torsion parts
      logical lgvst
      integer tor_bin_max,torsion,tor_1,tor_2,tor_3
     &     ,tor_4,tor_code,tor_num,tor_max,itor,iutor,patt,bthree
     &     ,decimal,power
      parameter (tor_bin_max=360,tor_max=15)
      double precision tor_bin,tor_tot,xcc,ycc,zcc,tcc,fplus,fminus
     &     ,ftrans,tor_bin_size
      double precision xvec,yvec,zvec,distij,xaa1,yaa1,zaa1,xa1a2
     &     ,ya1a2,za1a2,daa1,da1a2,dot,thetac,theta
     &     ,ratio,gdefect,btrans,bg_plus,bg_minus,g_plus,g_minus
     &     ,trans,psum,pattern

      dimension patt(tor_max)
      dimension tor_1(ntypmx,tor_max),tor_2(ntypmx,tor_max)
     &     ,tor_3(ntypmx,tor_max),tor_4(ntypmx,tor_max)
      dimension tor_num(ntypmx)
      dimension tor_code(numax,numax)
      dimension tor_bin(nboxmx,ntypmx,tor_max,tor_bin_max)
      dimension tor_tot(nboxmx,ntypmx,tor_max)
      dimension trans(nboxmx,ntypmx),g_minus(nboxmx,ntypmx)
     &     ,g_plus(nboxmx,ntypmx)
      dimension btrans(nboxmx,ntypmx,tor_max)
     & ,bg_plus(nboxmx,ntypmx,tor_max),bg_minus(nboxmx,ntypmx,tor_max)
      dimension gdefect(nboxmx,ntypmx,tor_max+1)
      dimension pattern(nboxmx,ntypmx,0:3**tor_max)

c     --- variables used in the charge parts
      logical qhere
      integer qbin,qbinmax,qbins,qqcode,imol,iunit
      parameter (qbinmax=1000)
      double precision qhist,qcount,qdisp
      double precision qmin,qmax,qdiff,qstep,qtemp
      dimension qhist(numax*ntmax+numax,nboxmx,qbinmax)
      dimension qcount(numax*ntmax+numax,nboxmx)

c     --- variables used in the dipole parts
      logical ldipole
      integer qblock,block,nblock
      double precision dipole,dicount,dipx,dipy,dipz,diconv,stddev
     &     ,avera,diprev,dcprev,dipblk
      dimension dipole(ntmax,nboxmx),dicount(ntmax,nboxmx)
      dimension diprev(ntmax,nboxmx),dcprev(ntmax,nboxmx)
      dimension dipblk(ntmax,nboxmx,20)

      dimension xcm(nmax),ycm(nmax),zcm(nmax)
      dimension hist(ntmax*ntmax*ntypmx*ntypmx,nbinmx)
      dimension bighist(nboxmx*ntypmx*ntmax*ntypmx*ntmax,nbinmx)
      dimension shell(2*ntypmx*ntypmx*ntmax*ntmax,nbinmx,2)
      dimension comhist(ntypmx*ntypmx,nbinmx)
      dimension combighist(nboxmx*ntypmx*ntypmx,nbinmx)
      dimension comshell(2*ntypmx*ntypmx,nbinmx,2)
      dimension volume(nboxmx),boxlx(nboxmx),boxly(nboxmx),boxlz(nboxmx)
      dimension rho(ntypmx*ntypmx*ntmax)
      dimension beadtyp(ntmax),nunit(ntypmx),decode(100)
      dimension rxu(nmax,numax),ryu(nmax,numax),rzu(nmax,numax)
     &     ,ntype(nmax,numax),qqu(nmax,numax)
      dimension ncmt(nboxmx,ntmax)
      dimension count(nboxmx,ntmax*ntypmx)
      dimension moltyp(nmax),nboxi(nmax)
      dimension invib(ntypmx,numax),ijvib(ntypmx,numax,4)
      dimension intor(ntypmx,numax),ijtor2(ntypmx,numax,4)
     &     ,ijtor3(ntypmx,numax,4),ijtor4(ntypmx,numax,4)
      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)
     &     ,distij(numax,numax)
      dimension rcut(nboxmx),boxmin(nboxmx)
      dimension nnone(nboxmx,ntypmx*ntypmx*ntmax*ntmax)
      dimension comnone(nboxmx,ntypmx*ntypmx)
      dimension cmolec(ntypmx,nboxmx)


      write(6,*) 'Please answer the questions .true. or .false.'

      write(6,*) 'do you want radial distribution functions?'
      read(5,*) lrdf
      if (lrdf) then
         write(6,*) 'do you want intramolecular rdf .true. or .false. ?'
         read(5,*) lintra

         if ( lintra ) then
            nskip = 0.5d0
         else
            nskip = 1.5d0
         endif
         write(6,*) 'input number of bins for rdf (less than',nbinmx,')'
         read(5,*) nbin
         
         write(6,*) 'amount of additional displacement (0.0 for none)?'
         read(5,*) addist

      endif

      lstretch = .false.

      write(6,*) 'do you want to compute the bond angle distribution?'
      read(5,*) lbend

      write(6,*) 'do you want to analyze torsion angles ?'
      read(5,*) lgvst

      write(6,*) 'do you want charge distribution?'
      read(5,*) lcharge

      if (lcharge) then
         write(6,*) 'input minimum charge for dist (qmin)'
         read(5,*) qmin
         write(6,*) 'input maximum charge for dist (qmax)'
         read(5,*) qmax
         write(6,*) 'input qbins must be less than',qbinmax
         read(5,*) qbins
         if ( qmax .lt. qmin ) stop 'qmax cannont be less than qmin'
         qdiff = qmax - qmin
         qstep = qdiff/dble(qbins)
         write(6,*) 'input additional disp for charges (0.0 for none)'
         read(5,*) qdisp
      endif
      write(6,*) 'do you want the average dipole moments?'
      read(5,*) ldipole
      if ( ldipole ) then
         write(6,*) 'number of blocks for dipole?'
         read(5,*) qblock
      endif

c     read initial information from fort.10
      read(10,*) nframe,nchain,nmolty,(rcut(z),z=1,nboxmx)
      read(10,*) nhere,(beadtyp(z),z=1,nhere)

c     --- safety checks
      if ( nmolty .gt. ntypmx ) then
         write(6,*) 'nmolty of ',nmolty,' gt ntypmx of ',ntypmx
         stop
      endif

      if ( nchain .gt. nmax ) then
         write(6,*) 'nchain of ',nchain,' gt nmax of ',nmax
         stop
      endif

      if ( ldipole ) then
         block = nframe/qblock
         write(6,*) 'block size',block
         nblock = 0
      endif

      if ( nhere .gt. ntmax ) then
         write(6,*) 'nhere =',nhere,' is greater than ntmax =',ntmax
         stop
      endif

      do kk = 1,nhere
         temp = beadtyp(kk)
         decode(temp) = kk
      enddo

      do imolty = 1,nmolty
         read(10,*) nunit(imolty)
         if ( nunit(imolty) .gt. numax ) stop 'nunit .gt. numax'
c        --- read bond connectivity information
         do ii = 1,nunit(imolty)
            read(10,*) invib(imolty,ii)
     &           ,(ijvib(imolty,ii,z),z=1,invib(imolty,ii))
         enddo

c        --- read torsional connectivity information
         do j = 1,nunit(imolty)
            read(10,*) intor(imolty,j)
            do ii = 1,intor(imolty,j)
               read(10,*) ijtor2(imolty,j,ii)
     &              ,ijtor3(imolty,j,ii),ijtor4(imolty,j,ii)
            enddo
         enddo


      enddo
      
c     place where someday could specify the number of boxes used
c     in the simulation
      nbox = nboxmx

      do kk = 1,nbox
         boxmin(kk) = 1000.0d0
      enddo

      if ( lbend ) then
c        --- compute ang_bin_size
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

         do xx = 1,ntypmx
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
         temp = ntmax*ntmax*ntypmx*ntypmx
         do ntij = 1,temp
            do kk=1,nbox
               nnone(kk,ntij) = 0.0d0
               binadj = (kk-1)*temp + ntij
               do bin = 1,nbin
                  bighist(binadj,bin) = 0.0d0
                  shell(binadj,bin,1) = 0.0d0
                  shell(binadj,bin,2) = 0.0d0
               enddo
            enddo
         enddo

         temp = ntypmx*ntypmx
         do ntij = 1,temp
            do kk=1,nbox
               comnone(kk,ntij) = 0.0d0
               binadj = (kk-1)*temp + ntij
               do bin = 1,nbin
                  combighist(binadj,bin) = 0.0d0
                  comshell(binadj,bin,1) = 0.0d0
                  comshell(binadj,bin,2) = 0.0d0
               enddo
            enddo
         enddo
      endif

      if ( lcharge ) then
         do kk = 1, nbox
            do imol = 1,ntmax
               do iunit = 1,numax
                  qqcode = numax*(imol-1) + iunit
                  do qbin = 1,qbins
                     qhist(qqcode,kk,qbin) = 0.0d0
                     qcount(qqcode,kk) = 0.0d0
                  enddo
               enddo
            enddo
         enddo
      endif
      
      if ( ldipole ) then
         do kk = 1,nbox
            do imol = 1,nmolty
               dipole(imol,kk) = 0.0d0
               dicount(imol,kk) = 0.0d0
               diprev(imol,kk) = 0.0d0
               dcprev(imol,kk) = 0.0d0
            enddo
         enddo
      endif

c     cycle through all of the frames and calculate the radial dist func

      do k = 1, nframe
c     read in frame step number
         read(10,*) fstep

c     read in ncmt for the boxes
         do kk = 1,nbox
            read(10,*) (ncmt(kk,z),z=1,nmolty)

            read(10,*) boxlx(kk),boxly(kk),boxlz(kk)

            volume(kk) = boxlx(kk)*boxly(kk)*boxlz(kk)
         enddo

c     initialize frame specific arrays
         do ii = 1,nmolty
            do gg = 1,nbox
               cmolec(ii,gg) = 0.0d0
               do g = 1,nhere
                  temp = nhere*(ii-1) + g
                  count(gg,temp) = 0.0d0
               enddo
            enddo
         enddo

         do nmolec = 1,nchain
            read(10,*) chnum,imolty,nunit(imolty),nboxi(chnum)
     &           ,xcm(chnum),ycm(chnum),zcm(chnum)
            moltyp(chnum) = imolty
            if (nmolec .ne. chnum) then
               write(6,*) 'nmolec,chnum',nmolec,chnum
               write(6,*) 'nunit,nboxi',nunit(imolty),nboxi(chnum)
               write(6,*) '*cm',xcm(chnum),ycm(chnum),zcm(chnum)
               stop 'nmolec ne chnum'
            endif
            tempbx = nboxi(chnum)

            cmolec(imolty,tempbx) = cmolec(imolty,tempbx) + 1.0d0

            do z = 1, nunit(imolty)
               read(10,*) rxu(chnum,z),ryu(chnum,z),rzu(chnum,z)
     &              ,qqu(chnum,z),ntype(chnum,z)
               temp = nhere*(imolty-1)+decode(ntype(chnum,z))
               count(tempbx,temp) = count(tempbx,temp) + 1.0d0
            enddo
         enddo

         if (lcharge) then
            do kk = 1,nbox
               do i = 1,nchain
                  if ( nboxi(i) .eq. kk ) then
                     imolty = moltyp(i)
                     do ii = 1,nunit(imolty)
                        qtemp = qqu(i,ii)
                        if ( qtemp .lt. qmin .or. qtemp .gt. qmax) then
                           write(6,*) 'q value out of range'
                           write(6,*) 'frame,i,ii,q',k,i,ii,qtemp
                           stop
                        endif
                        qtemp = qtemp - qmin
                        qbin = dint(qtemp/qstep) + 1
                        qqcode = numax*(imolty-1) + ii
                        qhist(qqcode,kk,qbin) = 
     &                       qhist(qqcode,kk,qbin) + 1.0d0
                        qcount(qqcode,kk) = qcount(qqcode,kk) + 1.0d0
                     enddo
                  endif
               enddo
            enddo
         endif

         if ( ldipole ) then
            do i = 1,nchain
               kk = nboxi(i)
               imolty = moltyp(i)
               dipx = 0.0d0
               dipy = 0.0d0
               dipz = 0.0d0
               do iunit = 1,nunit(imolty)
                  qtemp = qqu(i,iunit)
                  dipx = dipx + qtemp*rxu(i,iunit)
                  dipy = dipy + qtemp*ryu(i,iunit)
                  dipz = dipz + qtemp*rzu(i,iunit)
               enddo
               dipole(imolty,kk) = dipole(imolty,kk) 
     &              + dsqrt(dipx*dipx + dipy*dipy + dipz*dipz)
               dicount(imolty,kk) = dicount(imolty,kk) + 1.0d0
            enddo

            if ( mod(k,block) .eq. 0 ) then
               nblock = nblock + 1
               do kk = 1,nbox
                  do imol = 1,nmolty
                     if ( dicount(imolty,kk) .gt. 0.5d0 ) then
                        dipblk(imol,kk,nblock) = 
     &                       (dipole(imolty,kk) - diprev(imolty,kk))/
     &                       (dicount(imolty,kk) - dcprev(imolty,kk))
                     endif
                     diprev(imolty,kk) = dipole(imolty,kk)
                     dcprev(imolty,kk) = dicount(imolty,kk)
                  enddo
               enddo
            endif

         endif

        if (lrdf) then
c        --- compute the radial distribution histograms for this frame
         do kk = 1,nbox
            bx = boxlx(kk)
            hbx = bx/2.0d0
            by = boxly(kk)
            hby = by/2.0d0
            bz = boxlz(kk)
            hbz = bz/2.0d0
            if ( hbx .lt. boxmin(kk) ) boxmin(kk) = hbx
            rcutsq = rcut(kk)*rcut(kk)
            binstep = rcut(kk)/dble(nbin)
c           initiallize hist
            do g = 1,nmolty*nmolty*nhere*nhere
               do gg = 1,nbin
                  hist(g,gg) = 0.0d0
               enddo
            enddo
c           initiallize comhist
            do g = 1,nmolty*nmolty
               do gg = 1,nbin
                  comhist(g,gg) = 0.0d0
               enddo
            enddo

            do i = 1, nchain
 
c ### check if i is in relevant box ###
               if ( nboxi(i) .eq. kk ) then
                  imolty = moltyp(i)
c     --- loop over all beads ii of chain i 
                  do ii = 1, nunit(imolty)
                     ntii = nhere*(imolty-1) + decode(ntype(i,ii) )
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
     &                             decode( ntype(j,jj) )
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
                                 hist(ntij,bin) = hist(ntij,bin) + 1.0d0
                              endif

                           enddo
                        endif
                     enddo
                  enddo
               endif
            enddo

c     normalize the histogram and add it to the big histogram
            
            do ii = 1,nmolty
               do xx  = 1,nhere
                  ntii = nhere*(ii-1) + xx
                  rho(ntii) = (4.0d0 * 3.1415d0 * count(kk,ntii)) 
     &                 / (3.0d0 * volume(kk) )
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
                        if ( ntii .le. ntjj ) then 
                           ntij = (ntii-1)*nhere*nmolty + ntjj
                           binadj = (kk-1)*nhere*nhere*nmolty*nmolty 
     &                          + ntij
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
     &                             .lt. 0.5d0 .or. numxx*numyy 
     &                             .lt. 0.5d0) then
                                 lskip = .true.
                              endif
                           endif
                           
                           if ( .not. lskip  ) then
                              do bin = 1,nbin
                                 if ( ntii .eq. ntjj ) then
                                    shlsumx = shlsumx 
     &                                   + 2.0d0*hist(ntij,bin)
     &                                   /numxx
                                    shell(binadj,bin,1) = 
     &                                   shell(binadj,bin,1) 
     &                                   + shlsumx
                                 else
                                    shlsumx = shlsumx 
     &                                   + hist(ntij,bin)/numxx
                                    shlsumy = shlsumy
     &                                   + hist(ntij,bin)/numyy
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
                                 hist(ntij,bin) = hist(ntij,bin)
     &                                /xxideal
                                 bighist(binadj,bin) = 
     &                                bighist(binadj,bin) 
     &                                + hist(ntij,bin)
                              enddo
                              
                           else
                              nnone(kk,ntij) = nnone(kk,ntij) + 1.0d0
                           endif
                        endif
                     enddo
                  enddo
               enddo
            enddo

c           --- Center of Mass Radial Distribution Functions
            do i = 1, nchain
 
c ### check if i is in relevant box ###
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
                           comhist(ntij,bin) = comhist(ntij,bin) 
     &                          + 1.0d0
                        endif

                     endif
                  enddo
               endif
            enddo

c     normalize the COM histogram and add it to the big histogram
            
            do ii = 1,nmolty
               rho(ii) = (4.0d0 * 3.1415d0 * cmolec(ii,kk)) 
     &                 / (3.0d0 * volume(kk) )
            enddo
            
            do ii = 1,nmolty
               ntii = ii
               numxx = cmolec(ntii,kk)
                  
               do jj = ii,nmolty
                  ntjj = jj
                  numyy = cmolec(ntjj,kk)
                  
                  ntij = (ntii-1)*nmolty + ntjj
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
     &                          + 2.0d0*comhist(ntij,bin)
     &                          /numxx
                           comshell(binadj,bin,1) = 
     &                          comshell(binadj,bin,1)+ shlsumx
                        else
                           shlsumx = shlsumx 
     &                          + comhist(ntij,bin)/numxx
                           shlsumy = shlsumy 
     &                          + comhist(ntij,bin)/numyy
                           comshell(binadj,bin,1) = 
     &                          comshell(binadj,bin,1)+ shlsumx
                           comshell(binadj,bin,2) = 
     &                          comshell(binadj,bin,2)+ shlsumy
                        endif
                        
                        rlower = dble(bin-1)*binstep
                        rupper = rlower + binstep
                        nideal = const*(rupper**3 -rlower**3)
                        xxideal = numxx*nideal
                        comhist(ntij,bin) = comhist(ntij,bin)
     &                       /xxideal
                        combighist(binadj,bin) = combighist(binadj,bin) 
     &                       + comhist(ntij,bin)
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

               
c ### check if i is in relevant box ###
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
c     ### molecule with dihedral potenials ###
                     gaudef = 1
                     dum = 0
                     
                     do torsion = 1, tor_num(imolty)
                        
                        j = tor_1(imolty,torsion)
                        ip1 = tor_2(imolty,torsion)
                        ip2 = tor_3(imolty,torsion)
                        ip3 = tor_4(imolty,torsion)

c     ***  calculate cross product d_a x d_a-1  ***
                        xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     +                       zvec(ip1,j) * yvec(ip1,ip2)
                        yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     +                       xvec(ip1,j) * zvec(ip1,ip2)
                        zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     +                       yvec(ip1,j) * xvec(ip1,ip2)
c     ***  calculate cross product d_a-1 x d_a-2 ***
                        xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     +                       zvec(ip1,ip2) * yvec(ip3,ip2)
                        ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     +                       xvec(ip1,ip2) * zvec(ip3,ip2)
                        za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     +                       yvec(ip1,ip2) * xvec(ip3,ip2)
                        
c     *** calculate lengths of cross products ***
                        daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                        da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
c     *** calculate dot product of cross products ***
                        dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                        thetac = - dot / ( daa1 * da1a2 )
                        theta = dacos(thetac)
                        
c     *** calculate cross product of cross products ***
                        xcc = yaa1*za1a2 - zaa1*ya1a2
                        ycc = zaa1*xa1a2 - xaa1*za1a2
                        zcc = xaa1*ya1a2 - yaa1*xa1a2
c              *** calculate scalar triple product ***
                        tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2) 
     &                       + zcc*zvec(ip1,ip2) 
                        if ( tcc .lt. 0.0d0 ) theta = - theta
                        
c              *** convert theta to degrees ***
                        theta = theta*(180.0d0/onepi)
                        
c     --- bin the torsion --
                        bin = dint( (theta+180.0d0)/tor_bin_size)+1
                        tor_bin(kk,imolty,torsion,bin) = 
     &                       tor_bin(kk,imolty,torsion,bin)+1.0d0
                        tor_tot(kk,imolty,torsion) = 
     &                       tor_tot(kk,imolty,torsion) + 1.0d0
                        
                        if ( theta .lt. -60.0d0 ) then
c     --- gauch minus
                           g_minus(kk,imolty)=g_minus(kk,imolty) + 1.0d0
                           gaudef = gaudef + 1
                           bg_minus(kk,imolty,torsion) = 
     &                          bg_minus(kk,imolty,torsion) + 1.0d0
                           dum = dum + 0 * 3**(torsion-1)
                        elseif ( theta .lt. 60.0d0 ) then
c     --- trans
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
                 pattern(kk,imolty,dum) = pattern(kk,imolty,dum) + 1.0d0

                  endif
               endif
              enddo
           enddo
        endif

      enddo


c     whew, now that is done we just have to divide out the number of
c     frames in bighist

      if (lrdf) then
       write(6,*)
       write(6,*) 'bead bead radial distribution functions in fort.3*'
       write(6,*)
       write(6,*) 'avg. number of beads vs. shell size in fort.6*'

       do kk = 1,nbox
         binstep = rcut(kk)/dble(nbin)
         do ii = 1,nmolty
            do xx = 1,nhere
               ntii = nhere*(ii-1)+xx
               do jj = 1,nmolty
                  do yy = 1,nhere
                   ntjj = nhere*(jj-1)+yy
                   if ( ntii .le. ntjj ) then
                      ntij = (ntii-1)*nhere*nmolty + ntjj

                      binadj = (kk-1)*nhere*nhere*nmolty*nmolty + ntij
               
                      aframe = dble(nframe)-nnone(kk,ntij)
                      if ( aframe .gt. 0.5d0 ) then
                         write(30+kk,50) 0.0d0,addist,ii,beadtyp(xx)
     &                        ,jj,beadtyp(yy) 
                         write(60+kk,50) 0.0d0,0.0d0,ii,beadtyp(xx)
     &                        ,jj,beadtyp(yy) 
 50                      format(2f7.2,4i5)
                         do bin = 1,nbin
                            bighist(binadj,bin) = bighist(binadj,bin)
     &                           /aframe
                            
                            rxuij =  binstep*(dble(bin)-0.5d0)
                            write(30+kk,*) rxuij,bighist(binadj,bin)
     &                           +addist
                            shell(binadj,bin,1) = shell(binadj,bin,1)
     &                           /aframe
                            write(60+kk,*) rxuij,shell(binadj,bin,1)
                     
                         enddo
                         write(30+kk,*)
                         write(60+kk,*)

                         if ( ntii .ne. ntjj ) then
c                     --- now write the corresponding 60 file for 2
                            write(60+kk,50) 0.0d0,0.0d0,jj,beadtyp(yy) 
     &                           ,ii,beadtyp(xx)
                            do bin = 1,nbin
                               rxuij =  binstep*(dble(bin)-0.5d0)
                               shell(binadj,bin,2) = 
     &                              shell(binadj,bin,2) /aframe
                               write(60+kk,*) rxuij,shell(binadj,bin,2)
                            enddo
                            write(60+kk,*)
                         endif
                      elseif( aframe .lt. -0.5d0 ) then
                         write(6,*) 'aframe',aframe
                         write(6,*) 'nframe',nframe
                         write(6,*) 'nnone(kk,ntij),kk,ntij',
     &                        nnone(kk,ntij),kk,ntij
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
         binstep = rcut(kk)/dble(nbin)
         do ii = 1,nmolty
            ntii = ii
            do jj = ii,nmolty
               ntjj = jj
               ntij = (ntii-1)*nmolty + ntjj
               binadj = (kk-1)*nmolty*nmolty + ntij
               
               aframe = dble(nframe)-comnone(kk,ntij)
               if ( aframe .gt. 0.5d0 ) then
                  write(32+kk,51) 0.0d0,addist,ii,jj
                  write(62+kk,51) 0.0d0,0.0d0,ii,jj
 51               format(2f7.2,2i5)
                  do bin = 1,nbin
                     combighist(binadj,bin) = combighist(binadj,bin)
     &                    /aframe
                     
                     rxuij =  binstep*(dble(bin)-0.5d0)
                     write(32+kk,*) rxuij,combighist(binadj,bin)
     &                    +addist
                     comshell(binadj,bin,1) = 
     &                    comshell(binadj,bin,1) /aframe
                     write(62+kk,*) rxuij,comshell(binadj,bin,1)
                     
                  enddo
                  write(32+kk,*)
                  write(62+kk,*)

c                 --- now output 60 for 2
                  if ( ntii .ne. ntjj ) then
                     write(62+kk,51) 0.0d0,0.0d0,jj,ii
                     do bin = 1,nbin
                        rxuij =  binstep*(dble(bin)-0.5d0)
                        comshell(binadj,bin,2) = 
     &                       comshell(binadj,bin,2) /aframe
                        write(62+kk,*) rxuij,comshell(binadj,bin,2)
                     
                     enddo
                     write(62+kk,*)
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

      endif 

      if ( lcharge ) then
         write(6,*) 'charge distributions in fort.35+box'
         do kk = 1, nbox
            do imol = 1,nmolty
               do iunit = 1,nunit(imolty)
                  qqcode = numax*(imol-1) + iunit
                  qhere = .false.
                  do qbin = 1,qbins
                     if ( qcount(qqcode,kk) .gt. 0.5d0 ) then
                        qhere = .true.
                        qhist(qqcode,kk,qbin) = 
     &                       qhist(qqcode,kk,qbin)/qcount(qqcode,kk)
                     endif
                  enddo
                  if (qhere) then
                     qtemp = ( (dble(1)-0.5d0)*qstep)+qmin
                     write(35+kk,*) qtemp,qhist(qqcode,kk,1)+qdisp
     &                    ,imol,iunit
                     
                     do qbin = 2,qbins
                        qtemp = ( (dble(qbin)-0.5d0)*qstep)+qmin
                        write(35+kk,*) qtemp
     &                       ,qhist(qqcode,kk,qbin)+qdisp
                     enddo
                     write(35+kk,*)
                  endif
               enddo
            enddo
         enddo
      endif

      if ( ldipole ) then
         diconv = (1.602d-19)/((1d10)*(3.336d-30))
         do kk = 1, nbox
            write(6,*) 'Dipoles [e A], [D] in Box',kk
            do imol = 1,nmolty
               if ( dicount(imol,kk) .gt. 0.5d0 ) then
                  dipole(imol,kk) = dipole(imol,kk)/dicount(imol,kk)
               endif
               write(6,*) 'Moltyp,dipole',imolty
     &              ,dipole(imol,kk),dipole(imol,kk)*diconv
            enddo
         enddo

         if ( block .gt. 1 ) then
            write(6,*) 'Number of blocks input and used',block,nblock
            do kk = 1,nbox
               write(6,*) 'Box     Moltyp   Dipole [D]   Std dev [D]'
               do imol = 1,nmolty
                  avera = 0.0d0
                  do k=1,nblock
                     avera = avera + dipblk(imol,kk,k)
                  enddo
                  avera = avera/nblock
                  stddev = 0.0d0
                  do k = 1,nblock
                     stddev = stddev + (dipblk(imol,kk,k)-avera)**2
                  enddo
                  stddev = dsqrt( stddev/dble(nblock-1) )
                  stddev = stddev*diconv
                  avera = avera*diconv
                  write(6,'(i3,4x,i3,4x,2f10.4)') kk,imol,avera,stddev
               enddo
            enddo
         endif

      endif
               
      if ( lbend ) then
c     --- output the bending angle distributions for each angle in the molecule
c     --- write to 37+box

         do kk = 1,nbox
            do imolty = 1,nmolty
               do bend = 1,angle_num(imolty)
                  total = angle_tot(kk,imolty,bend)
                  if ( total .gt. 0.5d0 ) then
                     value = ang_bin_size/2.0d0
                     degree = value*180.0d0/onepi
                     write(37+kk,*) degree
     &                       ,angle_bin(kk,imolty,bend,bin)/total
     &                    ,'Moltyp ',imolty,' angle '
     &                    ,angle_1(imolty,bend),angle_2(imolty,bend)
     &                    ,angle_3(imolty,bend)
                     do bin = 2,ang_bin_max
                        value = value + ang_bin_size
                        degree = value*180.0d0/onepi
                        write(37+kk,*) degree
     &                       ,angle_bin(kk,imolty,bend,bin)/total
                     enddo
                     write(37+kk,*)
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
      endif

      if ( lgvst ) then
c     output the gauch versus trans data for each moltype
         write(6,*) 'moltyp  box  trans      g+        g-       g frac'
         do kk  = 1,nbox
            do imolty = 1,nmolty
               write(90+kk,*) 'Moltype ',imolty
               write(90+kk,41)
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
                     write(90+kk,42) tor_1(imolty,torsion)
     &                    ,tor_2(imolty,torsion),tor_3(imolty,torsion)
     &                    ,tor_4(imolty,torsion),fplus,fminus,ftrans

C                    --- output the torsion probablility distributions
                     value = -180.0d0 + (0.5d0*tor_bin_size)
                     write(95+kk,*) value,tor_bin(kk,imolty,torsion,1)
     &                    ,imolty,tor_1(imolty,torsion)
     &                    ,tor_2(imolty,torsion),tor_3(imolty,torsion)
     &                    ,tor_4(imolty,torsion)

                     do bin = 2,tor_bin_max
                        value = value + tor_bin_size
                        write(95+kk,*) value
     &                       ,tor_bin(kk,imolty,torsion,bin)
                     enddo
                     write(95+kk,*)
                  endif
               enddo
            enddo
         enddo
 40      format(i3,5x,i3,3(2x,e8.2),f8.4)
 41      format('Units',8x,'Frac g+',1x,'Frac g-',1x,'Frac t')
 42      format(4(i2,1x),1x,3(f5.3,3x))
         write(6,*)
         write(6,*) 'gauch fraction vs. torsion in fort.90+box'
         write(6,*)
         write(6,*) 'torsion angle distribution in fort.95+box'
         write(6,*)
         write(6,*) 'probabiltiy vs. number of defects per chain'
         write(6,*) 'shown in fort.5*'
         write(6,*)

c     analyse the number of gauch defects per chain
         do kk = 1,nbox
            do imolty = 1,nmolty
               total = gdefect(kk,imolty,tor_max+1)
               if ( total .gt. 0.01d0 ) then
                  write(50+kk,*) 'molecule type ',imolty
                  do torsion = 1,tor_num(imolty)
                     gdefect(kk,imolty,torsion) = 
     &                    gdefect(kk,imolty,torsion)/total
                     write(50+kk,*) torsion,gdefect(kk,imolty,torsion)
                  enddo
               endif

               psum = 0.0d0
               if ( tor_num(imolty) .gt. 0 ) then
c                 --- write out a decoder for the torsions
                  if ( kk .eq. 1 ) then
                     write(20,*) 'Moltyp ',imolty,' torsions'
                     do torsion = 1,tor_num(imolty)
                        write(20,*) torsion,'   ',tor_1(imolty,torsion)
     &                       ,tor_2(imolty,torsion)
     &                       ,tor_3(imolty,torsion)
     &                       ,tor_4(imolty,torsion)
                     enddo
                  endif
                  units = 3**( tor_num(imolty) )
                  do uu = 0,units
                     psum = psum + pattern(kk,imolty,uu)
                  enddo
                  write(20+kk,*) 'Code     Prob.     Torsion pattern'
                  do uu = 0,units-1
                     if ( pattern(kk,imolty,uu) .gt. 0.5d0 ) then
                        decimal = uu
                        power = units
                        torsion = tor_num(imolty)
                        do temp = torsion,1,-1
                           power = power/3
                           bthree = decimal/power 
                           decimal = decimal - bthree*power
                           patt(temp) = bthree-1
                        enddo

                        write(20+kk,52) uu,
     &                       pattern(kk,imolty,uu)/psum
     &                       ,(patt(temp),temp=1,torsion)
 52                     format(i8,2x,f5.3,2x,20(i2,1x))
                     endif
                  enddo 
               endif

            enddo
         enddo

         write(6,*) 'patterns of gauch defects in fort.2*'
         write(6,*) 'transform first number into base 3 to get the'
         write(6,*) 'pattern of gauch defects where -1=g-, 0=trans 1=g+'
         write(6,*)
         write(6,*) 'the maximum rcut value that could be used is'
         write(6,*) boxmin
         write(6,*)
         

      endif

      end
