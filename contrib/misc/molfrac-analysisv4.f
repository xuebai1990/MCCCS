      program NTanalysis

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   modified version of Becky's RDF/NUM code which computes
ccc   COM number integrals and rdf's
ccc   also calculates the local mole fraction enhancement
ccc   also calculates the distribution of dipole vector angles
ccc   in theory, works for more than one box, but has not been tested
ccc   05/11/07
ccc   modified 06/26/07 to calculate angles for the C-N vector
ccc   which is much closer to the actual dipole vector
ccc   unformatted read statement to handle movie file from new code
ccc   05/27/09 calculates orientational order parameter S
ccc   ldipole needs to be true for zangle and order parameter
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer bmax,nmax,numax,ntmax,nbxmax,nntype,nvmax,fmax
      integer ncmax,ios
      parameter (bmax=50000,nmax=2050,numax=20,ntmax=12,nbxmax=3)
      parameter (nntype=220,nvmax=30,fmax=50, ncmax=8000)
      integer i,ii,j,jj,k,kk,l,ll,nbox,btype1,btype2,nconfig,nchain
      integer nmolty,nbeadty,beadlist(nntype),nunit(ntmax)
      integer nvib(ntmax,numax),step,ijvib(ntmax,numax,nvmax)
      integer ntor(ntmax,numax),ijtor2(ntmax,numax,nvmax)
      integer ijtor3(ntmax,numax,nvmax),ijtor4(ntmax,numax,nvmax)
      integer ncmt(nbxmax,ntmax),nc,idum,imolty(nmax),nboxi(nmax)
      integer boxi,boxj,bin(nbxmax),moltyj,maxbin(nbxmax),minbin(nbxmax)
      integer mtype1,mtype2,moltyi
      integer nfile,ilist,tconfig,beadty(ntmax,numax)
      integer len1, len2, len3, count, btype3
      double precision pi,rcut(nbxmax),boxlx(nbxmax)
      double precision boxlz(nbxmax),rlower,const1,const2,boxly(nbxmax)
      double precision xcm(nmax),ycm(nmax),zcm(nmax),xcoord(nmax,numax)
      double precision xnn, ynn, znn, rnn
      double precision ycoord(nmax,numax),zcoord(nmax,numax),ddum,rupper
      double precision q(nmax, numax)
      double precision maxR(nbxmax),xij,yij,zij,rij,vol,dens1, dens2
      double precision nideal1, nideal2, binwidth(nbxmax),rx,
     &     nct1(nbxmax,ncmax), nct2(nbxmax,ncmax), 
     &     nct1avg(nbxmax), nct2avg(nbxmax) 
      double precision binwidth2, maxbin2, minbin2
      double precision n11(nbxmax,bmax), n12(nbxmax,bmax),
     &     n21(nbxmax,bmax), n22(nbxmax,bmax)
      double precision g11(nbxmax,bmax), g12(nbxmax,bmax), 
     &     g21(nbxmax,bmax),g22(nbxmax,bmax),  
     &     nint11(nbxmax,bmax), nint12(nbxmax,bmax), 
     &     nint21(nbxmax,bmax), nint22(nbxmax,bmax)
      double precision gavg11(nbxmax,bmax), gavg12(nbxmax,bmax), 
     &     gavg21(nbxmax,bmax), gavg22(nbxmax,bmax)
      double precision numavg11(nbxmax,bmax), numavg12(nbxmax,bmax), 
     &     numavg21(nbxmax,bmax), numavg22(nbxmax,bmax)
      double precision num11(nbxmax), num12(nbxmax), num21(nbxmax),
     &     num22(nbxmax)
      integer boxcount(nbxmax, ncmax), boxcountavg(nbxmax)
      double precision molefrac(nbxmax,nntype)
      double precision local11(nbxmax,bmax),local12(nbxmax,bmax), 
     &     local21(nbxmax,bmax),local22(nbxmax,bmax),
     &     totalnum1(nbxmax,bmax),totalnum2(nbxmax,bmax),
     &     localavg11(nbxmax,bmax), localavg12(nbxmax,bmax),
     &     localavg21(nbxmax,bmax), localavg22(nbxmax,bmax)
      integer  jzero11(nbxmax,bmax),jzero12(nbxmax,bmax),
     &     jzero21(nbxmax,bmax), jzero22(nbxmax,bmax)

      integer diptyp, bin2(nbxmax), bin3(nbxmax)
      double precision theta1(nbxmax,bmax*2), theta2(nbxmax,bmax*2),
     &     theta3(nbxmax,bmax*2), theta4(nbxmax,bmax*2),
     &     avgtheta1(nbxmax,bmax*2), avgtheta2(nbxmax,bmax*2),
     &     avgtheta3(nbxmax,bmax*2), avgtheta4(nbxmax,bmax*2)
      double precision rcm, rdip1, rdip2,rcmx,rcmy,rcmz,costheta,
     &     rdipx1, rdipy1, rdipz1, rdipx2, rdipy2, rdipz2
      integer count1(nbxmax), count2(nbxmax), count3(nbxmax), 
     &     count4(nbxmax), bincount(nbxmax)
      double precision avg1(nbxmax), avg2(nbxmax), avg3(nbxmax),
     &     avg4(nbxmax)
      double precision jzero11avg(nbxmax,bmax), jzero12avg(nbxmax,bmax),
     &     jzero21avg(nbxmax,bmax), jzero22avg(nbxmax,bmax)
      double precision interval1, interval2
      double precision dipx, dipy, dipz, dip, avgdip
      double precision rxdip, rydip, rzdip
      integer dipcount, dipmolecule, zcount, avgzcount
      double precision zangle(bmax), avgzangle(bmax), avgz 
      double precision s, avgs
      character*40 filename1,filename2,filelist(fmax)
      character*40 filename3, filename4, filename5
      character*40 filename6, filename7, filename8
      character*40 filename9, filename10, filename11,filename12
      character*40 filename13,filename14,filename15,filename16
      character*6 atom(nntype),char1,char2,char3
      character*25 fileformat, fileformat2, fileformat3, fileformat4
      character*25 fileformat5
      logical lskip(nbxmax), langle, ldipole
      
      pi=4.0d0*atan(1.0d0)

c --  atom types from suijtab.f

c      atom(3)='CH4 '
c      atom(4)='CH3 '
      atom(5)='CH2  '
      atom(57)='chg  '
      atom(58)='chg  '
      atom(59)='C    '
      atom(61)='H    '
      atom(62)='O    '
      atom(63)='CH3  '
      atom(64)='CH2  '
      atom(92)='H    '
      atom(99)='CH3  '
      atom(100)='C    '
      atom(101)='O    '
      atom(102)='O    '
      atom(103)='H    '
      atom(107)='O    '
      atom(108)='H    '
      atom(109)='O    '
      atom(110)='H    '
      atom(111)='O    '
      atom(112)='H    '
      atom(113)='M    '
      atom(114)='O    '
      atom(115)='H    '
      atom(116)='M    '
      atom(117)='O    '
      atom(118)='H    '
      atom(119)='LP   '
      atom(124)='HE   '
      atom(129)='C    '
      atom(130)='O    '
      atom(175)='O    '
      atom(176)='H    '
      atom(1) = 'O    '
      atom(2) = 'H    '
      atom(3) = 'BC   '
      atom(4) = 'LP   '
      atom(200)='BEN  '
      atom(201)='LJC12'
      atom(202)='LJNB '
      atom(203)='LJmNT'
      atom(204)='LJoNT'
      atom(205)='LJC10'
      atom(206)='mNT  '
      atom(207)='oNT  '
      atom(208)='C10  '
      atom(209)='BUT'
      atom(210)='ACR'
      atom(211)='MMA'
      atom(212)='HEA'
      atom(213)='HEM'
      atom(214)='IOA'
      atom(215)='EHA'
      atom(216)='ETH'
      atom(217)='HEX'
      atom(218)='WAT'
      atom(219)='MET'

      atom(6) = 'C6 '
      atom(36) = 'C36'

      open(4,file='analinput')
      read(4,*)
      read(4,*)
      read(4,*) nfile
        read(4,*)
      do i=1,nfile
         read(4,*) filelist(i)
      enddo
      read(4,*)
      read(4,*) nbox
        read(4,*)
      do i=1,nbox
         read(4,*) binwidth(i)
        enddo
      read(4,*)
      read(4,*)
      read(4,*) btype1,btype2
      write(6,*) btype1, btype2
      read(4,*)
      read(4,*) langle
c      if (langle) then
         read(4,*)
         read(4,*) diptyp
         read(4,*)
         read(4,*) btype3
         read(4,*)
         read(4,*) interval1
         read(4,*) 
         read(4,*) interval2
         read(4,*)
         read(4,*) binwidth2
c      endif
      read(4,*)
      read(4,*) ldipole
c      if (ldipole) then
         read(4,*)
         read(4,*) dipmolecule
c      endif
c --- dipole moment calculation assumes 1 box
c --- as does orienation with z vector

c      write(6,*) 'langle ', langle

C     set histogram counters and logicals
      tconfig=0

c     keep track of nct1 and nct2 for each box for each config
c     in case of fluctuating numbers of particles (i.e. swaps)
      
      do i=1,nbox
         do k=1,ncmax
            nct1(i,k)=0.0d0
            nct2(i,k)=0.0d0
            boxcount(i,k)=0
         enddo
         lskip(i)=.false.
         nct1avg(i)=0.0d0
         nct2avg(i)=0.0d0
         boxcountavg(i) = 0

      enddo

      count = 0
      dipcount = 0
      zcount = 0
      avgzcount = 0
      dipx = 0.0d0
      dipy = 0.0d0
      dipz = 0.0d0
      avgdip = 0.0d0
      avgz = 0.0d0
      avgs = 0.0d0
      
      do i=1,nbox
         do j=1,bmax
            gavg11(i,j)=0.0d0
            gavg12(i,j)=0.0d0
            gavg21(i,j)=0.0d0
            gavg22(i,j)=0.0d0
            
            numavg11(i,j)=0.0d0
            numavg12(i,j)=0.0d0
            numavg21(i,j)=0.0d0
            numavg22(i,j)=0.0d0

            localavg11(i,j)=0.0d0
            localavg12(i,j)=0.0d0
            localavg21(i,j)=0.0d0
            localavg22(i,j)=0.0d0

            jzero11(i,j)=0
            jzero12(i,j)=0
            jzero21(i,j)=0
            jzero22(i,j)=0
 
            jzero11avg(i,j) = 0.0d0
            jzero12avg(i,j) = 0.0d0
            jzero21avg(i,j) = 0.0d0
            jzero22avg(i,j) = 0.0d0

            avgzangle(j) = 0.0d0

         enddo
      enddo
      
c      open(80)

      do ilist=1,nfile

         open(10,file=filelist(ilist))

C     Read header
         read(10,*) nconfig,nchain,nmolty,(rcut(i),i=1,nbox)
c         write(6,*) nconfig, nchain, nmolty, (rcut(i), i=1,nbox)
         read(10,*) nbeadty,(beadlist(ii),ii=1,nbeadty)
         tconfig=tconfig+nconfig
         
         do i=1,nmolty
            read(10,*) nunit(i)
            do ii=1,nunit(i)
               read(10,*) nvib(i,ii),(ijvib(i,ii,j),j=1,nvib(i,ii))
            enddo
            
            do ii=1,nunit(i)
               read(10,*) ntor(i,ii)
               do j=1,ntor(i,ii)
                  read(10,*) ijtor2(i,ii,j),ijtor3(i,ii,j),
     &                 ijtor4(i,ii,j)
               enddo
            enddo
         enddo
         write(6,*)'nconfig',nconfig
c         write(6,*)'nchain', nchain
c         write(6,*)'nmolty', nmolty
c         write(6,*)'nbeadty',nbeadty

C     Loop over configurations
         do nc=1,nconfig
            count = count + 1

C     Read bead coordinates
            read(10,*) step
            
            do i=1,nbox
               read(10,*) (ncmt(i,j), j=1,nmolty)
               read(10,*) boxlx(i),boxly(i),boxlz(i)
c               write(6,*)'boxlength', boxlx(i), boxly(i), boxlz(i)

C              set-up bins
               maxR(i)=boxlz(i)/2.0d0
c               maxR(i) = boxlz(i)
               maxbin(i)=int(maxR(i)/binwidth(i))
               minbin(i)=int(-maxR(i)/binwidth(i))
               maxbin2=int(maxR(i)/binwidth2)
               minbin2=int(-maxR(i)/binwidth2)
c               write(6,*) 'maxR(i)', maxR(i)
c               write(6,*) 'binwidth(i) ', binwidth(i)
c               write(6,*) 'maxbin', maxbin(i)

c --- initialize accumulators to zero
               do j=1, maxbin(i)
                  n11(i,j)=0.0d0
                  n12(i,j)=0.0d0
                  n21(i,j)=0.0d0
                  n22(i,j)=0.0d0
                  g11(i,j)=0.0d0
                  g12(i,j)=0.0d0
                  g21(i,j)=0.0d0
                  g22(i,j)=0.0d0
                  nint11(i,j)=0.0d0
                  nint12(i,j)=0.0d0
                  nint21(i,j)=0.0d0
                  nint22(i,j)=0.0d0
                  local11(i,j)=0.0d0
                  local12(i,j)=0.0d0
                  local21(i,j)=0.0d0
                  local22(i,j)=0.0d0
               enddo
               num11(i)=0.0d0
               num12(i)=0.0d0
               num21(i)=0.0d0
               num22(i)=0.0d0
               zcount=0

               if (langle) then
                  count1(i) = 0
                  count2(i) = 0
                  count3(i) = 0
                  count4(i) = 0
                  avg1(i)=0.0d0
                  avg2(i)=0.0d0
                  avg3(i)=0.0d0
                  avg4(i)=0.0d0
                  bincount(i)=0.0d0
                  do jj=minbin(i), maxbin(i)
                     theta1(i,jj)=0.0d0
                     theta2(i,jj)=0.0d0
                     theta3(i,jj)=0.0d0
                     theta4(i,jj)=0.0d0
                     zangle(jj) = 0.0d0
                   enddo
                  if (count .eq. 1) then
                     avgtheta1(i,jj)=0.0d0
                     avgtheta2(i,jj)=0.0d0
                     avgtheta3(i,jj)=0.0d0
                     avgtheta4(i,jj)=0.0d0
                  endif
               endif
            enddo
            do i=1,nchain
ccc --- new movie file has different format than old
ccc --- use unformatted read so that it works for both
               read(10,*) idum,imolty(i),
     &              nunit(imolty(i)),nboxi(i),xcm(i),ycm(i),zcm(i)
               if (idum.ne.i) stop 
               do ii=1,nunit(imolty(i))
                  read(10,*) xcoord(i,ii),ycoord(i,ii),
     &                 zcoord(i,ii),q(i,ii),beadty(imolty(i),ii)
c$$$               read(10,'(4(1x,I5),3(1X,F10.6))') idum,imolty(i),
c$$$     &              nunit(imolty(i)),nboxi(i),xcm(i),ycm(i),zcm(i)
c$$$               if (idum.ne.i) stop 
c$$$               do ii=1,nunit(imolty(i))
c$$$                  read(10,'(4(1x,F10.6),I5)') xcoord(i,ii),ycoord(i,ii),
c$$$     &                 zcoord(i,ii),ddum,beadty(imolty(i),ii)
c                  write(6,'(4(1x,F10.6),I5)') xcoord(i,ii),ycoord(i,ii),
c     &                 zcoord(i,ii),ddum,beadty(imolty(i),ii)

               enddo
            enddo


c$$$C      if a molecule is not in the box, set lskip to true 
c$$$            do i=1,nbox
c$$$               if (ncmt(i,1).eq.0.or.ncmt(i,2).eq.0) then
c$$$                  lskip(i)=.true.
c$$$                  write(6,*)i,lskip(i)
c$$$               endif
c$$$            enddo




cc - lskip made program a bit more efficient, but threw off the averaging at the end
cc - modified how lskip works


C     Calculate specific pair distances

C           Loop over chains - type 1 becomes the origin bead
C            rdf is directional type 1 to type 2, 
C            although type 2 to type 1 should be same
C            number integrals will be different depending on direction

            do l=1,nchain

               moltyi=imolty(l)
               boxi=nboxi(l)
               boxcount(boxi, count) = boxcount(boxi,count) + 1
c               write(6,*) 'nboxi(',l,')', nboxi(l)               

C              Count number of molecules of type 1 that could interact
               if (moltyi.eq.1) then
                  nct1(boxi,count)=nct1(boxi,count)+1.0d0
               endif
C              Count number of molecules of type 2 that could interact
               if (moltyi.eq.2) then
                  nct2(boxi,count)=nct2(boxi,count)+1.0d0
               endif

               if (ldipole .and. moltyi.eq.dipmolecule) then
                  dipcount = dipcount + 1
                  dipx = 0.0
                  dipy = 0.0
                  dipz = 0.0
                  do i=1,nunit(moltyi)
                     rxdip = xcoord(l,i) - xcm(l)
c                     write(6,*) 'molecule ', l, ' unit ', i
c                     write(6,*) 'xcoord ', xcoord(l,i), 
c     &                    ' xcm ', xcm(l), ' rxdip ', rxdip
                     rydip = ycoord(l,i) - ycm(l)
                     rzdip = zcoord(l,i) - zcm(l)

                     dipx = dipx + rxdip*q(l,i)

                     dipy = dipy + rydip*q(l,i)
                     dipz = dipz + rzdip*q(l,i)
c                     write(6,*) 'q ', q(l,i), ' dipx ', dipx,
c     &                    '  dipy ', dipy, ' dipz ', dipz
                  enddo

                  dip = dsqrt(dipx**2 + dipy**2 + dipz**2)
c                  write(6,*) '  DIPOLE MOMENT molecule ', l, moltyi,
c     &                 dip*4.802d0
                  
                  avgdip = avgdip + dip


c ---- calculate angle between dipole vector and z-direction (direction of field)

c                          calculate vector between carbon (unit 4)  
c                           and nitrogen atom (unit 3)
c                         (not actual dipole vector, but very close)
                  rdipx1=xcoord(l,3)-xcoord(l,4)
                  rdipy1=ycoord(l,3)-ycoord(l,4)
                  rdipz1=zcoord(l,3)-zcoord(l,4)
                                 
                  rdipx1=rdipx1-boxlx(boxi)*
     &                 anint(rdipx1/boxlx(boxi))
                  rdipy1=rdipy1-boxly(boxi)*
     &                 anint(rdipy1/boxly(boxi))
                  rdipz1=rdipz1-boxlz(boxi)*
     &                 anint(rdipz1/boxlz(boxi))
                        

                  rdip1=dsqrt(rdipx1**2+rdipy1**2
     &                 +rdipz1**2)

c                                 write(6,*) 'rdip1 ', rdipx1, rdipy1, 
c     &                                rdipz1, rdip1
                  rdipx2=0.0d0
                  rdipy2=0.0d0
                  rdipz2=-1.0d0
                  rdip1=dsqrt(rdipx2**2+rdipy2**2
     &                 +rdipz2**2)
                                 
c                             calculate the cos of the angle between the vectors

                  costheta=(rdipx1*rdipx2+rdipy1*rdipy2+
     &                 rdipz1*rdipz2)/(rdip1*rdip2)


                  if (costheta.ge.1.0d0) costheta=1.0d0
                  if (costheta.le.-1.0d0) costheta=-1.0d0

c                  write(6,*) 'costheta ', costheta
                  
                  avgz = avgz + costheta
                  avgzcount = avgzcount + 1
c                            add 1 to cosines so that all are positive
c                            will subtract this after binning

                  costheta=costheta+1.0d0

c                 determine bin
                  bin3(boxi)=int(costheta/
     &                 binwidth2)+1
                  if (bin3(boxi).le.maxbin2) then
                     zangle(bin3(boxi))=
     &                    zangle(bin3(boxi))+
     &                    1.0d0
                     zcount = zcount + 1  
c                     write(6,*) 'zcount ', zcount
c                     write(6,*) 'zangle ', zangle(bin3(boxi)), 
c     &                    bin3(boxi)
                  endif

c ---   determine order parameter S = 0.5<3cos**2(theta)-1>
                  costheta = costheta-1.0d0
                  s = 3.0d0*costheta**2 - 1.0d0
c                  write(6,*) 's ', s
                  avgs = avgs + s

               endif   ! end if ldipole

C              Loop over interacting chains - type 1 and type 2
               do ll=1,nchain
c                 write(6,*)'ll', ll
                  moltyj=imolty(ll)
                  boxj=nboxi(ll)

                  if (l .ne. ll) then
C              if molecule is type 1
                     if (moltyi.eq.1 ) then
c ---  calculate 1-1 and 1-2 nint

C                       Check that the interacting molecule is type 1 
C                         and is in the same box
c ----  calculate 1-1 nint

                        if (moltyj.eq. 1 .and.boxi.eq.boxj) then

                                 
C                          Calculate COM distances
                           xij=xcm(l)-xcm(ll)
                           yij=ycm(l)-ycm(ll)
                           zij=zcm(l)-zcm(ll)
                         
C                       Minimum image
                           xij=xij-boxlx(boxi)*anint
     &                          (xij/boxlx(boxi))
                           yij=yij-boxly(boxi)*
     &                          anint(yij/boxly(boxi))
                           zij=zij-boxlz(boxi)*
     &                          anint(zij/boxlz(boxi))
                                 
c                                    write(6,*)'x',l,ll,'=', xij
c                                    write(6,*)'y',l,ll,'=', yij
c                                    write(6,*)'z',l,ll,'=', zij
C                          Bin distances
                           rij=dsqrt(xij**2+yij**2+zij**2)
                           bin(boxi)=int(rij/binwidth(boxi))+1

c$$$                        write(6,*) l, xcm(l), ycm(l), zcm(l)
c$$$                        write(6,*) ll, xcm(ll), ycm(ll), zcm(ll) 
c$$$                        write(6,*) xij, yij, zij

c                        write(6,*)'rij ', l, ll, rij
c                                    write(6,*)'bin(',boxi,')', bin(boxi)
                        
                           if (bin(boxi).le.maxbin(boxi)) then
                              
                              n11(boxi,bin(boxi))=n11(boxi,
     &                             bin(boxi))+1.0d0
c                              if (boxi .eq. 1) then
c                                 write(6,*)'n11(',boxi,bin(boxi), 
c     &                                ')', n11(boxi,bin(boxi))
c                              endif
                           endif
                        endif 

C                       Check that the interacting molecule is type 2 
C                         and is in the same box
c ----  calculate 1-2 nint

                        if (moltyj.eq. 2.and.boxi.eq.boxj)then

                                 
C                          Calculate COM distances
                           xij=xcm(l)-xcm(ll)
                           yij=ycm(l)-ycm(ll)
                           zij=zcm(l)-zcm(ll)
                        
C                          Minimum image
                           xij=xij-boxlx(boxi)*anint
     &                          (xij/boxlx(boxi))
                           yij=yij-boxly(boxi)*
     &                          anint(yij/boxly(boxi))
                           zij=zij-boxlz(boxi)*
     &                          anint(zij/boxlz(boxi))
                        
c                                    write(6,*)'x',l,ll,'=', xij
c                                    write(6,*)'y',l,ll,'=', yij
c                                    write(6,*)'z',l,ll,'=', zij
C                          Bin distances
                           rij=dsqrt(xij**2+yij**2+zij**2)
                           bin(boxi)=int(rij/binwidth(boxi))+1
 
c$$$                        write(6,*) l, xcm(l), ycm(l), zcm(l)
c$$$                        write(6,*) ll, xcm(ll), ycm(ll), zcm(ll)   
c$$$                        write(6,*) xij, yij, zij
                        

c                                    write(6,*)'bin(',boxi,')', bin(boxi)
                        
                           if (bin(boxi).le.maxbin(boxi)) then
                              
                              n12(boxi,bin(boxi))=n12(boxi,
     &                             bin(boxi))+1.0d0
c                               write(6,*)'n12(',boxi,bin(boxi), 
c     &                                   ')', n12(boxi,bin(boxi))
                           endif
                        endif 
                     endif


C              if molecule is type 2
c ---  calculate 2-1 and 2-2 nint

                     if (moltyi.eq.2 .and. .not.lskip(boxi)) then

C                       Check that the interacting molecule is type 1 
C                         and is in the same box
c ----  calculates 2-1 nint
                        if (moltyj.eq. 1 .and.boxi.eq.boxj) then

                                 
C                          Calculate COM distances
                           xij=xcm(l)-xcm(ll)
                           yij=ycm(l)-ycm(ll)
                           zij=zcm(l)-zcm(ll)
                        
C                          Minimum image
                           xij=xij-boxlx(boxi)*anint
     &                          (xij/boxlx(boxi))
                           yij=yij-boxly(boxi)*
     &                          anint(yij/boxly(boxi))
                           zij=zij-boxlz(boxi)*
     &                          anint(zij/boxlz(boxi))
                                 
c                                    write(6,*)'x',l,ll,'=', xij
c                                    write(6,*)'y',l,ll,'=', yij
c                                    write(6,*)'z',l,ll,'=', zij
C                          Bin distances
                           rij=dsqrt(xij**2+yij**2+zij**2)
                           bin(boxi)=int(rij/binwidth(boxi))+1

c$$$                        write(6,*) l, xcm(l), ycm(l), zcm(l)
c$$$                        write(6,*) ll, xcm(ll), ycm(ll), zcm(ll)
c$$$                        write(6,*) xij, yij, zij

c                        write(6,*)'rij', l, ll, rij
c                                    write(6,*)'bin(',boxi,')', bin(boxi)
                        
                           if (bin(boxi).le.maxbin(boxi)) then
                           
                              n21(boxi,bin(boxi))=n21(boxi,
     &                             bin(boxi))+1.0d0
c                               write(6,*)'n21(',boxi,bin(boxi), 
c     &                                   ')', n21(boxi,bin(boxi))
                           endif
                        endif 

C                       Check that the interacting molecule is type 2 
C                         and is in the same box
c ----  calculates 2-2 nint
                        if (moltyj.eq. 2.and.boxi.eq.boxj) then

                                 
C                          Calculate COM distances
                           xij=xcm(l)-xcm(ll)
                           yij=ycm(l)-ycm(ll)
                           zij=zcm(l)-zcm(ll)
                        
C                       Minimum image
                           xij=xij-boxlx(boxi)*anint
     &                          (xij/boxlx(boxi))
                           yij=yij-boxly(boxi)*
     &                          anint(yij/boxly(boxi))
                           zij=zij-boxlz(boxi)*
     &                          anint(zij/boxlz(boxi))
                           
c                                    write(6,*)'x',l,ll,'=', xij
c                                    write(6,*)'y',l,ll,'=', yij
c                                    write(6,*)'z',l,ll,'=', zij
C                          Bin distances
                           rij=dsqrt(xij**2+yij**2+zij**2)
                           bin(boxi)=int(rij/binwidth(boxi))+1
c$$$                        write(6,*) l, xcm(l), ycm(l), zcm(l)
c$$$                        write(6,*) ll, xcm(ll), ycm(ll), zcm(ll)
c$$$                        write(6,*) xij, yij, zij
                                    
c                                    write(6,*)'rij', l, ll, rij
c                                    write(6,*)'bin(',boxi,')', bin(boxi)
                        
                           if (bin(boxi).le.maxbin(boxi)) then
                              
                              n22(boxi,bin(boxi))=n22(boxi,
     &                             bin(boxi))+1.0d0
c                               write(6,*)'n22(',boxi,bin(boxi), 
c     &                                   ')', n22(boxi,bin(boxi))
                           endif
                        endif 
                     endif




c ---  calculates angle between dipole vectors of moltyp diptyp
                     if (langle) then
                        if(boxi.eq.boxj) then
                           if(moltyj.eq.diptyp.and.moltyi.eq.diptyp)then
c                          calculate  COM distance and NN distance
                              rcmx=xcm(l)-xcm(ll)
                              rcmy=ycm(l)-ycm(ll)
                              rcmz=zcm(l)-zcm(ll)
                              xnn = xcoord(l,3) - xcoord(ll,3)
                              ynn = ycoord(l,3) - ycoord(ll,3)
                              znn = zcoord(l,3) - zcoord(ll,3)
c                          minimum image
                              rcmx=rcmx-boxlx(boxi)*
     &                             anint(rcmx/boxlx(boxi))
                              rcmy=rcmy-boxly(boxi)*
     &                             anint(rcmy/boxly(boxi))
                              rcmz=rcmz-boxlz(boxi)*
     &                             anint(rcmz/boxlz(boxi))
                              xnn = xnn - boxlx(boxi)*
     &                             anint(xnn/boxlx(boxi))
                              ynn = ynn - boxly(boxi)*
     &                             anint(ynn/boxly(boxi))
                              znn = znn - boxlz(boxi)*
     &                             anint(znn/boxlz(boxi))
                              
                              rcm=dsqrt(rcmx**2+rcmy**2+rcmz**2)
                              rnn=dsqrt(xnn**2+ynn**2+znn**2)
c                              write(6,*) 'rcm ', rcm

c                           switch from com (rcm) to N-N (rnn)
                              if (rnn .lt. interval2) then
c                          calculate vector between carbon (unit 4)  
c                           and nitrogen atom (unit 3)
                                 rdipx1=xcoord(l,3)-xcoord(l,4)
                                 rdipy1=ycoord(l,3)-ycoord(l,4)
                                 rdipz1=zcoord(l,3)-zcoord(l,4)
                                 
                                 rdipx1=rdipx1-boxlx(boxi)*
     &                                anint(rdipx1/boxlx(boxi))
                                 rdipy1=rdipy1-boxly(boxi)*
     &                                anint(rdipy1/boxly(boxi))
                                 rdipz1=rdipz1-boxlz(boxi)*
     &                                anint(rdipz1/boxlz(boxi))
                                 

                                 rdip1=dsqrt(rdipx1**2+rdipy1**2
     &                                +rdipz1**2)

c                                 write(6,*) 'rdip1 ', rdipx1, rdipy1, 
c     &                                rdipz1, rdip1
                                 rdipx2=xcoord(ll,3)-xcoord(ll,4)
                                 rdipy2=ycoord(ll,3)-ycoord(ll,4)
                                 rdipz2=zcoord(ll,3)-zcoord(ll,4)
                                 
                                 rdipx2=rdipx2-boxlx(boxi)*
     &                                anint(rdipx2/boxlx(boxi))
                                 rdipy2=rdipy2-boxly(boxi)*
     &                                anint(rdipy2/boxly(boxi))
                                 rdipz2=rdipz2-boxlz(boxi)*
     &                                anint(rdipz2/boxlz(boxi))                

                                 rdip2=dsqrt(rdipx2**2+rdipy2**2+
     &                                rdipz2**2)
                                 
c                                 write(6,*) 'rdip2 ', rdipx2, rdipy2, 
c     &                                rdipz2, rdip2
c                             calculate the cos of the angle between the vectors

                                 costheta=(rdipx1*rdipx2+rdipy1*rdipy2+
     &                                rdipz1*rdipz2)/(rdip1*rdip2)


                                 if (costheta.ge.1.0d0) costheta=1.0d0
                                 if (costheta.le.-1.0d0) costheta=-1.0d0

c                                 write(6,*) 'costheta ', costheta

c                            add 1 to cosines so that all are positive
c                            will subtract this after binning

                                 costheta=costheta+1.0d0

c                             determine proper histogram
                                 if (rnn.le.interval1) then
c                               determine bin
                                    bin2(boxi)=int(costheta/
     &                                   binwidth2)+1
                                    if (bin2(boxi).le.maxbin2) then
                                       theta1(boxi,bin2(boxi))=
     &                                      theta1(boxi,bin2(boxi))+
     &                                      1.0d0
                                       count1(boxi)=count1(boxi)+1
          
                                    endif
                                    
                                 else
c                             determine bin
                                    bin2(boxi)=int(costheta/
     &                                   binwidth2)+1
                                    if (bin2(boxi).le.maxbin2) then
                                       theta2(boxi,bin2(boxi)) = 
     &                                  theta2(boxi,bin2(boxi))+
     &                                      1.0d0
                                       count2(boxi)=count2(boxi)+1
c                                       write(6,*) 'theta2 ', 
c     &                                      theta2(boxi,bin2(boxi))
                                    endif   
                                 endif
                              endif
                           endif
                        endif
                     endif      ! end if (langle)

c$$$                     if (rij .lt. rcut(boxi)) then
c$$$                        write(80,'(A3, 3X, I3, 3X, I3, 3X, F7.4)')'rij', 
c$$$     &                       l, ll, rij
c$$$                     endif

                  endif         ! end if (l .ne. ll)
               enddo            ! end loop over molecules j (interacting molecules ll)
            enddo               ! end loop over molecules i (l)


C     Normalize to ideal gas

            do boxi=1, nbox  ! loop over number of boxes
            
c               write(6,*) 'outside loop j ', j
c               write(6,*) 'maxbin(boxi) ', maxbin(boxi)
c               write(6,*) 'boxi ', boxi
               
               vol=boxlx(boxi)*boxly(boxi)*boxlz(boxi)

               if (nct1(boxi,count).gt.0.0) then
                  dens1=nct1(boxi,count)/vol
               else
                  dens1 = 0.0d0
               endif
               if (nct2(boxi,count) .gt. 0.0) then
                  dens2=nct2(boxi,count)/vol
               else
                  dens2 = 0.0d0
               endif

               const1=4.0d0*pi*dens1/3.0d0
               const2=4.0d0*pi*dens2/3.0d0
               
               j=1
               
c            write(6,*)
c            write(6,*) 'count = ', count

c --- determine mole fraction of each component
c --- should not divide by nchain to get molfrac in each box
c --- in 2 box simulation!
c --- need to divide by total number in the box                
               molefrac(boxi,1) = nct1(boxi,count)/
     &              dble(boxcount(boxi,count))
               molefrac(boxi,2) = nct2(boxi,count)/
     &              dble(boxcount(boxi,count))
c               write(6,*) 'boxcount box ', boxi, boxcount(boxi,count)
c               write(6,*) 'mole fraction of type 1 box ',boxi, 
c     &              molefrac(boxi,1)
c               write(6,*) 'mole fraction of type 2 box ', boxi,
c     &              molefrac(boxi,2)

               do j=1,maxbin(boxi)
c               write(6,*) 'inside loop j ', j

                  rx=binwidth(boxi)*(dble(j)-0.5d0)
                  rlower=binwidth(boxi)*dble(j-1)
                  rupper=rlower+binwidth(boxi)

c               write(6,*) 'rx ', rx
            
                  nideal1=const1*(rupper**3-rlower**3)
                  nideal2=const2*(rupper**3-rlower**3)
 
c               write(6,*) 'vol ', vol
c               write(6,*) 'dens ', dens
c               write(6,*)'const ', const
c               write(6,*)'rupper ', rupper
c               write(6,*)'rlower ', rlower
c               write(6,*) 'nideal ', nideal

c               write(6,*) 'nct1 ', nct1(boxi,count)
c               write(6,*) 'nct2 ', nct2(boxi,count)

                  if (nct1(boxi,count) .gt. 0.0d0) then
                     num11(boxi)=num11(boxi)+n11(boxi,j)/
     &                    nct1(boxi,count)
                     num12(boxi)=num12(boxi)+n12(boxi,j)/
     &                    nct1(boxi,count)
                  endif
                  
                  nint11(boxi,j)=num11(boxi)
                  
                  nint12(boxi,j) = num12(boxi)
                  
                  if (nct2(boxi,count) .gt. 0.0d0) then                 
                     num21(boxi)=num21(boxi)+n21(boxi,j)/
     &                    nct2(boxi,count)
                     
                     num22(boxi)=num22(boxi)+n22(boxi,j)/
     &                    nct2(boxi,count)
                  endif

                  nint21(boxi,j) = num21(boxi)
                  
                  nint22(boxi,j) = num22(boxi)

               enddo

               do j=1,maxbin(boxi)
c ---- rdf           
                  rx=binwidth(boxi)*(dble(j)-0.5d0)
                  rlower=binwidth(boxi)*dble(j-1)
                  rupper=rlower+binwidth(boxi)
                  
c               write(6,*) 'rx ', rx
            
                  nideal1=const1*(rupper**3-rlower**3)
                  nideal2=const2*(rupper**3-rlower**3)
    
                  if (nct1(boxi,count) .gt. 0.0 .and. 
     &                 nct2(boxi,count) .gt. 0.0) then
                     g11(boxi,j)=n11(boxi,j)/nct1(boxi,count)/nideal1
                     g12(boxi,j)=n12(boxi,j)/nct1(boxi,count)/nideal2
                     g21(boxi,j)=n21(boxi,j)/nct2(boxi,count)/nideal1
                     g22(boxi,j)=n22(boxi,j)/nct2(boxi,count)/nideal2
c                     write(6,*) 'g11 ', g11(boxi,j)
                  elseif (nct1(boxi,count) .eq. 0.0) then
                     g11(boxi,j)=0.0d0
c                     write(6,*) 'g11 is zero'
                     g12(boxi,j)=0.0d0
                     g21(boxi,j)=0.0d0
                     g22(boxi,j)=n22(boxi,j)/nct2(boxi,count)/nideal2
                  elseif (nct2 (boxi,count) .eq. 0.0) then
                     g11(boxi,j)=n11(boxi,j)/nct1(boxi,count)/nideal1
                     g12(boxi,count) = 0.0d0
                     g22(boxi,count) = 0.0d0
                     g21(boxi,count) = 0.0d0
                  endif
                                 

c ---- local mole fraction
                  totalnum1(boxi,j)=nint11(boxi,j)+nint12(boxi,j)
                  totalnum2(boxi,j)=nint21(boxi,j)+nint22(boxi,j)

                  if (totalnum1(boxi,j) .ne. 0.00) then
ccc  switch from 0.05 to 1.0 for solubility parameter paper
                     if (nint11(boxi,j) .gt. 0.5d0) then
                        local11(boxi,j)= nint11(boxi,j)/
     &                    totalnum1(boxi,j)/molefrac(boxi,1)
                        jzero11(boxi,j)=1
                     else
                        local11(boxi,j)=0.0d0
                        jzero11(boxi,j)=0
c                        if (boxi.eq.1) then
c                           write(6,*) j, jzero11(boxi,j)
c                        endif
                     endif

c                     if (boxi.eq.1.and.j.eq.300) then
c                        write(6,*) totalnum1(boxi,j), nint11(boxi,j)
c                        write(6,*) local11(boxi,j), jzero11(boxi,j)
c                     endif

c                     if (boxi.eq.1) then
c                        write(6,*) 'local11 ', local11(boxi,j)
c                     endif
                     
ccc  switch from 0.05 to 1.0 for solubility parameter paper
                     if (nint12(boxi,j) .gt. 0.5d0) then
                        local12(boxi,j)= nint12(boxi,j)/
     &                       totalnum1(boxi,j)/molefrac(boxi,2)
                        jzero12(boxi,j)=1
                     else
                        local12(boxi,j)=0.0d0
                        jzero12(boxi,j)=0
                     endif
                  endif

                  if (totalnum2(boxi,j) .ne. 0.00) then
ccc  switch from 0.05 to 1.0 for solubility parameter paper
                     if (nint21(boxi,j) .gt. 0.5d0) then
                        local21(boxi,j)= nint21(boxi,j)/
     &                       totalnum2(boxi,j)/molefrac(boxi,1)
                        jzero21(boxi,j)=1
                     else
                        local21(boxi,j)=0.0d0
                        jzero21(boxi,j)=0
                     endif
         
ccc  switch from 0.05 to 1.0 for solubility parameter paper            
                     if (nint22(boxi,j) .gt. 0.5d0) then
                        local22(boxi,j)= nint22(boxi,j)/
     &                       totalnum2(boxi,j)/molefrac(boxi,2)
                        jzero22(boxi,j)=1
                     else
                        local22(boxi,j)=0.0d0
                        jzero22(boxi,j)=0
                     endif
                  endif


c ---- accumulate averages
                  gavg11(boxi,j)=gavg11(boxi,j)+g11(boxi,j)
                  gavg12(boxi,j)=gavg12(boxi,j)+g12(boxi,j)
                  gavg21(boxi,j)=gavg21(boxi,j)+g21(boxi,j)
                  gavg22(boxi,j)=gavg22(boxi,j)+g22(boxi,j)

c                  if (boxi.eq.1) then
c                     write(6,*) 'gavg21 ', j, gavg21(boxi,j)
c                  endif

                  numavg11(boxi,j)=numavg11(boxi,j)+nint11(boxi,j)
                  numavg12(boxi,j)=numavg12(boxi,j)+nint12(boxi,j)
                  numavg21(boxi,j)=numavg21(boxi,j)+nint21(boxi,j)
                  numavg22(boxi,j)=numavg22(boxi,j)+nint22(boxi,j)
                  
c               write(6,*) 'numavg11 ', numavg11(boxi,j)

                  localavg11(boxi,j)=localavg11(boxi,j)+local11(boxi,j)
                  localavg12(boxi,j)=localavg12(boxi,j)+local12(boxi,j)
                  localavg21(boxi,j)=localavg21(boxi,j)+local21(boxi,j)
                  localavg22(boxi,j)=localavg22(boxi,j)+local22(boxi,j)
                  
                  jzero11avg(boxi,j)=jzero11avg(boxi,j)+jzero11(boxi,j)
                  jzero12avg(boxi,j)=jzero12avg(boxi,j)+jzero12(boxi,j)
                  jzero21avg(boxi,j)=jzero21avg(boxi,j)+jzero21(boxi,j)
                  jzero22avg(boxi,j)=jzero22avg(boxi,j)+jzero22(boxi,j)

c                  if (boxi.eq.1) then 
c                     write(6,*) 'localavg11 ', localavg11(boxi,j)
c                  endif

               enddo            ! end loop j=1,maxbin(boxi)
    
               if (langle) then
c              cos(theta) goes from -1 to 1, so # of bins is:
                  bincount(boxi)=2.0d0/dble(binwidth2)
                  
                  do jj=minbin2,maxbin2
c     if (theta1(boxi,jj) .ne. 0.0) then
c     write(6,*) 'theta1 ', theta1(boxi,jj)
c     write(6,*) 'count1 ', count1(boxi)
c     write(6,*) 'bincount ', bincount(boxi)
c     write(6,*) 'count ', count, 'jj ', jj
c     endif
                     if (count1(boxi) .gt. 0) then
                        theta1(boxi,jj)=(theta1(boxi,jj)/
     &                       dble(count1(boxi)))*bincount(boxi)
                     endif
c                  write(6,*) 'afterwards, theta1 ', theta1(boxi,jj)
                     if (count2(boxi) .gt. 0) then
                        theta2(boxi,jj)=(theta2(boxi,jj)/
     &                       dble(count2(boxi)))*bincount(boxi)
c                        write(6,*) 'theta2 ', theta2(boxi,jj)
                     endif
c                     theta3(boxi,jj)=(theta3(boxi,jj)/
c     &                    dble(count3(boxi)))*bincount(boxi)
c                     theta4(boxi,jj)=(theta4(boxi,jj)/
c     &                    dble(count4(boxi)))*bincount(boxi)

                     avgtheta1(boxi,jj)=avgtheta1(boxi,jj)+
     &                    theta1(boxi,jj)
                     avgtheta2(boxi,jj)=avgtheta2(boxi,jj)+
     &                    theta2(boxi,jj)
c                     write(6,*) 'avgtheta2 ', avgtheta2
c                     avgtheta3(boxi,jj)=avgtheta3(boxi,jj)+
c     &                    theta3(boxi,jj)
c                     avgtheta4(boxi,jj)=avgtheta4(boxi,jj)+
c     &                    theta4(boxi,jj)
c                  if (avgtheta1(boxi,jj) .ne. 0.0 .and. jj.eq.50) then
c                     write(6,*) 'theta1 ', theta1(boxi,jj)
c                     write(6,*)'jj ',jj,' avgtheta1 ',avgtheta1(boxi,jj)
c                  endif
                  enddo

               endif            !end if (langle)

               if (ldipole) then
c              cos(theta) goes from -1 to 1, so # of bins is:
                  bincount(boxi)=2.0d0/dble(binwidth2)
c                  write(6,*) 'bincount ', bincount(boxi)
c                  write(6,*) 'minbin2 ', minbin2, ' maxbin2 ', maxbin2
                  do jj=minbin2, maxbin2
                     if (zcount .gt. 0) then
c                        write(6,*) 'zangle ', zangle(jj), jj
                        zangle(jj) = (zangle(jj)/dble(zcount))
     *                       *bincount(boxi)
                        if (zangle(jj) .gt. 1000) then
                           write(6,*) 'zangle ', zangle(jj)
                        endif
                     endif

                     avgzangle(jj) = avgzangle(jj) + zangle(jj)
c                     write(6,*) 'avgzangle ', avgzangle(jj)
                  enddo
               endif            ! end if ldipole

            enddo   ! end loop over boxes
               
c         if (count .eq. 1) then
c            do j=1, maxbin(boxi)
c               open(50)
c               open(51)
c               open(52)
c               open(53)
c               rx=binwidth(boxi)*(dble(j)-0.5d0)
c               write(50,*) rx, g11(boxi,j)
c               write(51,*) rx, g12(boxi,j)
c               write(52,*) rx, g21(boxi,j)
c               write(53,*) rx, g22(boxi,j)
c            enddo
c         endif

         enddo                  ! end loop over configurations

c$$$c -----    test COMNUM for fluctuations
c$$$               if (rx .gt. 7.54 .and. rx .lt. 7.56) then
c$$$                  if (mod(count,1).eq. 0) then
c$$$c                  write(6,*) 'writing to fort.50'
c$$$                     write(50,*) count, n11(boxi,j)
c$$$                  endif
c$$$               endif
c$$$         close(10)

      enddo                     ! end loop over files

c     only tested for one box

      do boxi=1,nbox
         
         write(6,*)'BOX',boxi
      
c        calculate average number of each type in each box
c        and the average total number in each box

         do k=1,count
            nct1avg(boxi)=nct1avg(boxi)+nct1(boxi,k)
            nct2avg(boxi)=nct2avg(boxi)+nct2(boxi,k)
c            if (boxi.eq.1) then
c               write(6,*) 'nct1 ', nct1(boxi,k), ' nct2 ', nct2(boxi,k)
c            endif
            boxcountavg(boxi) = boxcountavg(boxi) + boxcount(boxi,k)
c            write(6,*) 'boxcountavg ', k, boxcountavg(boxi)
         enddo

c         write(6,*) 'tconfig ', tconfig

         nct1avg(boxi)=nct1avg(boxi)/dble(tconfig)
         nct2avg(boxi)=nct2avg(boxi)/dble(tconfig)
         boxcountavg(boxi)=boxcountavg(boxi)/dble(tconfig)
         
c         write(6,*)'nctavg box',boxi, nct1avg(boxi),nct2avg(boxi),
c     &        boxcountavg(boxi)

cc         write(6,*) 'lskip(',boxi,')', lskip(boxi)
cc         if (.not.lskip(boxi)) then
            
            char1=atom(btype1)
            
            char2=atom(btype2)

            len1 = index(char1,' ')
            len2 = index(char2, ' ')
            len1 = len1 - 1
            len2 = len2 - 1

            if (langle) then
               if (btype3 .eq. btype1) then
                  char3=char1
                  len3=len1
               else
                  char3=char2
                  len3=len2
               endif
            endif

C     Create output filenames
c     first determine the length of the character strings
c     to get the right format statement
            if (len1 .eq. 1) then
               fileformat='(A10,I1,A1,A1,A1,A1,A4)'
               if (len2 .eq. 1) then
                  fileformat2='(A10,I1,A1,A1,A1,A1,A4)'
                  fileformat3='(A10,I1,A1,A1,A1,A1,A4)'
                  fileformat4='(A10,I1,A1,A1,A1,A1,A4)'
               elseif(len2 .eq.2) then
                  fileformat2='(A10,I1,A1,A1,A1,A2,A4)'
                  fileformat3='(A10,I1,A1,A2,A1,A1,A4)'
                  fileformat4='(A10,I1,A1,A2,A1,A2,A4)'
               elseif(len2.eq.3) then
                  fileformat2='(A10,I1,A1,A1,A1,A3,A4)'
                  fileformat3='(A10,I1,A1,A3,A1,A1,A4)'
                  fileformat4='(A10,I1,A1,A3,A1,A3,A4)'
               elseif(len2.eq.4) then
                  fileformat2='(A10,I1,A1,A1,A1,A4,A4)'
                  fileformat3='(A10,I1,A1,A4,A1,A1,A4)'
                  fileformat4='(A10,I1,A1,A4,A1,A4,A4)'
               elseif(len2.eq.5) then
                  fileformat2='(A10,I1,A1,A1,A1,A5,A4)'
                  fileformat3='(A10,I1,A1,A5,A1,A1,A4)'
                  fileformat4='(A10,I1,A1,A5,A1,A5,A4)'
               endif
               
            elseif (len1 .eq. 2) then
               fileformat='(A10,I1,A1,A2,A1,A2,A4)'
               if (len2 .eq. 1) then
                  fileformat2='(A10,I1,A1,A2,A1,A1,A4)'
                  fileformat3='(A10,I1,A1,A1,A1,A2,A4)'
                  fileformat4='(A10,I1,A1,A1,A1,A1,A4)'
               elseif(len2 .eq.2) then
                  fileformat2='(A10,I1,A1,A2,A1,A2,A4)'
                  fileformat3='(A10,I1,A1,A2,A1,A2,A4)'
                  fileformat4='(A10,I1,A1,A2,A1,A2,A4)'
               elseif(len2.eq.3) then
                  fileformat2='(A10,I1,A1,A2,A1,A3,A4)'
                  fileformat3='(A10,I1,A1,A3,A1,A2,A4)'
                  fileformat4='(A10,I1,A1,A3,A1,A3,A4)'
               elseif(len2.eq.4) then
                  fileformat2='(A10,I1,A1,A2,A1,A4,A4)'
                  fileformat3='(A10,I1,A1,A4,A1,A2,A4)'
                  fileformat4='(A10,I1,A1,A4,A1,A4,A4)'
               elseif(len2.eq.5) then
                  fileformat2='(A10,I1,A1,A2,A1,A5,A4)'
                  fileformat3='(A10,I1,A1,A5,A1,A2,A4)'
                  fileformat4='(A10,I1,A1,A5,A1,A5,A4)'
               endif
               
            elseif (len1 .eq. 3) then
               fileformat='(A10,I1,A1,A3,A1,A3,A4)'
               if (len2 .eq. 1) then
                  fileformat2='(A10,I1,A1,A3,A1,A1,A4)'
                  fileformat3='(A10,I1,A1,A1,A1,A3,A4)'
                  fileformat4='(A10,I1,A1,A1,A1,A1,A4)'
               elseif(len2 .eq.2) then
                  fileformat2='(A10,I1,A1,A3,A1,A2,A4)'
                  fileformat3='(A10,I1,A1,A2,A1,A3,A4)'
                  fileformat4='(A10,I1,A1,A2,A1,A2,A4)'
               elseif(len2.eq.3) then
                  fileformat2='(A10,I1,A1,A3,A1,A3,A4)'
                  fileformat3='(A10,I1,A1,A3,A1,A3,A4)'
                  fileformat4='(A10,I1,A1,A3,A1,A3,A4)'
               elseif(len2.eq.4) then
                  fileformat2='(A10,I1,A1,A3,A1,A4,A4)'
                  fileformat3='(A10,I1,A1,A4,A1,A3,A4)'
                  fileformat4='(A10,I1,A1,A4,A1,A4,A4)'
               elseif(len2.eq.5) then
                  fileformat2='(A10,I1,A1,A3,A1,A5,A4)'
                  fileformat3='(A10,I1,A1,A5,A1,A3,A4)'
                  fileformat4='(A10,I1,A1,A5,A1,A5,A4)'
               endif
               
            elseif (len1 .eq. 4) then
               fileformat='(A10,I1,A1,A4,A1,A4,A4)'
               if (len2 .eq. 1) then
                  fileformat2='(A10,I1,A1,A4,A1,A1,A4)'
                  fileformat3='(A10,I1,A1,A1,A1,A4,A4)'
                  fileformat4='(A10,I1,A1,A1,A1,A1,A4)'
               elseif(len2 .eq.2) then
                  fileformat2='(A10,A4,A1,A2,A4)'
                  fileformat3='(A10,I1,A1,A2,A1,A4,A4)'
                  fileformat4='(A10,I1,A1,A2,A1,A2,A4)'
               elseif(len2.eq.3) then
                  fileformat2='(A10,I1,A1,A4,A1,A3,A4)'
                  fileformat3='(A10,I1,A1,A3,A1,A4,A4)'
                  fileformat4='(A10,I1,A1,A3,A1,A3,A4)'
               elseif(len2.eq.4) then
                  fileformat2='(A10,I1,A1,A4,A1,A4,A4)'
                  fileformat3='(A10,I1,A1,A4,A1,A4,A4)'
                  fileformat4='(A10,I1,A1,A4,A1,A4,A4)'
               elseif(len2.eq.5) then
                  fileformat2='(A10,I1,A1,A4,A1,A5,A4)'
                  fileformat3='(A10,I1,A1,A5,A1,A4,A4)'
                  fileformat4='(A10,I1,A1,A5,A1,A5,A4)'
               endif
               
            elseif (len1 .eq. 5) then
               fileformat='(A10,I1,A1,A5,A1,A5,A4)'
               if (len2 .eq. 1) then
                  fileformat2='(A10,I1,A1,A5,A1,A1,A4)'
                  fileformat3='(A10,I1,A1,A1,A1,A5,A4)'
                  fileformat4='(A10,I1,A1,A1,A1,A1,A4)'
               elseif(len2 .eq.2) then
                  fileformat2='(A10,I1,A1,A5,A1,A2,A4)'
                  fileformat3='(A10,I1,A1,A2,A1,A5,A4)'
                  fileformat4='(A10,I1,A1,A2,A1,A2,A4)'
               elseif(len2.eq.3) then
                  fileformat2='(A10,I1,A1,A5,A1,A3,A4)'
                  fileformat3='(A10,I1,A1,A3,A1,A5,A4)'
                  fileformat4='(A10,I1,A1,A3,A1,A3,A4)'
               elseif(len2.eq.4) then
                  fileformat2='(A10,I1,A1,A5,A1,A4,A4)'
                  fileformat3='(A10,I1,A1,A4,A1,A5,A4)'
                  fileformat4='(A10,I1,A1,A4,A1,A4,A4)'
               elseif(len2.eq.5) then
                  fileformat2='(A10,I1,A1,A5,A1,A5,A4)'
                  fileformat3='(A10,I1,A1,A5,A1,A5,A4)'
                  fileformat4='(A10,I1,A1,A5,A1,A5,A4)'
               endif
            endif
            
            write(filename1,fileformat) 
     &           'COMRDF-box',boxi,'-',char1,'-', 
     &           char1,'.dat'   
c            write(6,*) filename1
            
            write(filename2,fileformat2) 
     &           'COMRDF-box',boxi,'-',char1,'-', 
     &           char2,'.dat' 
c         write(6,*) filename2

            write(filename3,fileformat3) 
     &           'COMRDF-box',boxi,'-',char2,'-', 
     &           char1,'.dat' 
c         write(6,*) filename3

            write(filename4,fileformat4) 
     &           'COMRDF-box',boxi,'-',char2,'-', 
     &           char2,'.dat'   
c         write(6,*) filename4

            write(filename5,fileformat) 
     &           'COMNUM-box',boxi,'-',char1,'-', 
     &           char1,'.dat'
c         write(6,*) filename5

            write(filename6,fileformat2) 
     &           'COMNUM-box',boxi,'-',char1,'-', 
     &           char2,'.dat'
c         write(6,*) filename6

            write(filename7,fileformat3) 
     &           'COMNUM-box',boxi,'-',char2,'-', 
     &           char1,'.dat'  
c         write(6,*) filename7

            write(filename8,fileformat4) 
     &           'COMNUM-box',boxi,'-',char2,'-', 
     &           char2,'.dat'
c         write(6,*) filename8
 
            write(filename9,fileformat) 
     &           'molfrac-bx',boxi,'-',char1,'-', char1,'.dat'            
c         write(6,*) filename9
      
            write(filename10,fileformat2) 
     &           'molfrac-bx',boxi,'-',char1,'-', char2,'.dat' 
c         write(6,*) filename10
      
            write(filename11,fileformat3) 
     &           'molfrac-bx',boxi,'-',char2,'-', char1,'.dat' 
c         write(6,*) filename11
      
            write(filename12,fileformat4) 
     &           'molfrac-bx',boxi,'-',char2,'-', char2,'.dat'   
c         write(6,*) filename12
           
            open(unit=21, file=filename1)
            open(unit=22, file=filename2)
            open(unit=23, file=filename3)
            open(unit=24, file=filename4)
            open(unit=25, file=filename5,status="unknown",iostat=ios)

            open(unit=26, file=filename6)
            open(unit=27, file=filename7)
            open(unit=28, file=filename8)

            if (ios.ne.0) stop  

            open(unit=12, file=filename9)
            open(unit=13, file=filename10)
            open(unit=14, file=filename11)
            open(unit=15, file=filename12)

 15         format(F6.2,F12.6)

            write(6,*) 'tconfig ', tconfig
            write(6,*) 'count ', count
            write(6,*) 'avg molefrac type 1 2', nct1avg(boxi)/
     &           dble(boxcountavg(boxi)), nct2avg(boxi)/
     &           dble(boxcountavg(boxi))
            if (tconfig .ne. count ) stop 'count .ne. tconfig'

c     calculate average


            do j=1,maxbin(boxi)
               rx=binwidth(boxi)*(dble(j)-0.5d0)
c     write(6,*) 'rx ', rx
               gavg11(boxi,j) = gavg11(boxi,j)/dble(tconfig)
               gavg12(boxi,j) = gavg12(boxi,j)/dble(tconfig)
               gavg21(boxi,j) = gavg21(boxi,j)/dble(tconfig)
               gavg22(boxi,j) = gavg22(boxi,j)/dble(tconfig)
               
               numavg11(boxi,j) = numavg11(boxi,j)/dble(tconfig)
               numavg12(boxi,j) = numavg12(boxi,j)/dble(tconfig)
               numavg21(boxi,j) = numavg21(boxi,j)/dble(tconfig)
               numavg22(boxi,j) = numavg22(boxi,j)/dble(tconfig)

               localavg11(boxi,j)=localavg11(boxi,j)/dble(tconfig)
               localavg12(boxi,j)=localavg12(boxi,j)/dble(tconfig)
               localavg21(boxi,j)=localavg21(boxi,j)/dble(tconfig)
               localavg22(boxi,j)=localavg22(boxi,j)/dble(tconfig)

               jzero11avg(boxi,j)=dble(jzero11avg(boxi,j))/dble(tconfig)
               jzero12avg(boxi,j)=dble(jzero12avg(boxi,j))/dble(tconfig)
               jzero21avg(boxi,j)=dble(jzero21avg(boxi,j))/dble(tconfig)
               jzero22avg(boxi,j)=dble(jzero22avg(boxi,j))/dble(tconfig)               


cc --- jzeroavg keeps track of how many times a molecule is at the appropriate
cc --- distance (i.e. jzeroavg% of the time a molecule is actually there
cc --- (we don't want to include a spot in the result if there's only a 
cc ---   molecule there 1% of the time)
cc --- this can be more easily accounted for using the number integral.

c               if (boxi.eq.1 .and. rx .lt. rcut(boxi)+1.0d0) then
c                  write(6,*) '  jzero11avg(boxi,j) ', jzero11avg(boxi,j)
c                  write(6,*) 'localavg11(boxi,j) ', rx, 
c     &                 localavg11(boxi,j)
c               endif

c --- cut off number integrals and rdf's at rcut + 1
               if (rx .lt. rcut(boxi)+5.0d0 ) then
c --- make sure nint and rdf not zero
                  if (gavg11(boxi,j).gt.0.0) then
                     write(21,*) rx,gavg11(boxi,j)
c                     write(6,*) boxi, rx, gavg11(boxi,j)
                  endif
                  if (gavg12(boxi,j).gt.0.0) then  
                     write(22,*) rx,gavg12(boxi,j)
                  endif
                  if (gavg21(boxi,j).gt.0.0) then
                     write(23,*) rx,gavg21(boxi,j)
c                     write(6,*) boxi, rx, gavg21(boxi,j)
                  endif
                  if (gavg22(boxi,j).gt.0.0) then
                     write(24,*) rx,gavg22(boxi,j)
                  endif
                  if (numavg11(boxi,j).ne.0.0) then
c                     write(6,*) rx, numavg11(boxi,j)
                     write(25,*) rx,numavg11(boxi,j)
                  endif
                  if (numavg12(boxi,j) .gt. 0.0) then
                     write(26,*) rx,numavg12(boxi,j)
                  endif
                  if (numavg21(boxi,j) .gt. 0.0) then
                     write(27,*) rx,numavg21(boxi,j)
                  endif
                  if (numavg22(boxi,j) .gt. 0.0) then
                     write(28,*) rx,numavg22(boxi,j)
                  endif
                  

cc --  i don't remember what jzero is supposed to do!
cc --  its getting set to zero for the one config with no
cc -- acrylate molecules in box 1, then it never prints out
cc -- local mole frac info

c                  if (jzero11(boxi,j).eq. 1 .and. 
c     &                 localavg11(boxi,j).ne.0.00) then
c                     write(12,*) rx,localavg11(boxi,j)
                  if (boxi.eq.1) then
c                     write(6,*) rx,jzero21avg(boxi,j),localavg21(boxi,j)
                  endif
                  if (localavg11(boxi,j).ne.0.00 .and. 
     &                 jzero11avg(boxi,j) .ge. 0.989) then
                     write(12,*) rx,localavg11(boxi,j)
                     if (boxi.eq.1 ) then
c                        write(6,*) 'localavg11 ', rx, localavg11(boxi,j)
                     endif
                  endif
                  if (localavg12(boxi,j).ne.0.00 .and. 
     &                 jzero12avg(boxi,j) .ge. 0.989) then
                     write(13,*) rx,localavg12(boxi,j)
                  endif
                  if (localavg21(boxi,j).ne.0.00 .and. 
     &                 jzero21avg(boxi,j) .ge. 0.989) then
                     write(14,*) rx,localavg21(boxi,j)
                  endif
c                  if (jzero22(boxi,j) .eq. 1 .and.
c     &                 localavg22(boxi,j) .ne. 0.00) then
                  if (localavg22(boxi,j) .ne. 0.00 .and. 
     &                 jzero22avg(boxi,j) .ge. 0.989) then
                     write(15,*) rx,localavg22(boxi,j)

                  endif

               endif
            enddo


c$$$            if (langle) then

c$$$               if (len3.eq.1) then
c$$$                  fileformat5='(A1,A6,I1,A2,I1,A4,I1,A4)'
c$$$               elseif (len3.eq.2) then
c$$$                  fileformat5='(A2,A6,I1,A2,I1,A4,I1,A4)'
c$$$               elseif (len3.eq.3) then
c$$$                  fileformat5='(A3,A6,I1,A2,I1,A4,I1,A4)'
c$$$               elseif (len3.eq.4) then
c$$$                  fileformat5='(A4,A6,I1,A2,I1,A4,I1,A4)'
c$$$               elseif (len3.eq.5) then
c$$$                  fileformat5='(A5,A6,I1,A2,I1,A4,I1,A4)'
c$$$               endif
c$$$               
c$$$c               write(6,*) 'writing filenames '
c$$$               write(filename13,fileformat5)
c$$$     &              char3,'angle_',3,'to',4,'_box',boxi,'.dat'
c$$$               write(filename14,fileformat5)
c$$$     &              char3,'angle_',4,'to',5,'_box',boxi,'.dat'
c$$$               write(filename15,fileformat5)
c$$$     &              char3,'angle_',5,'to',6,'_box',boxi,'.dat'
c$$$               write(filename16,fileformat5)
c$$$     &              char3,'angle_',6,'to',7,'_box',boxi,'.dat'   
c$$$               
c$$$               open(unit=30,file=filename13)
c$$$               open(unit=31,file=filename14)
c$$$               open(unit=32,file=filename15)
c$$$               open(unit=33,file=filename16)
               
c               write(6,*) 'opened files '
c               write(6,*) filename13, filename14, filename15, filename16


            if (langle) then


               open(unit=30,file='angle-1.dat')
               open(unit=31,file='angle-2.dat')


               do jj=minbin2,maxbin2
c                  write(6,*) 'avgtheta1 ', avgtheta1(boxi,jj)
c                  write(6,*) 'tconfig ', tconfig
                  avgtheta1(boxi,jj)=avgtheta1(boxi,jj)/
     &                 dble(tconfig)
c                  write(6,*) 'after ', avgtheta1(boxi,jj)
                  avgtheta2(boxi,jj)=avgtheta2(boxi,jj)/
     &                 dble(tconfig)
c                  avgtheta3(boxi,jj)=avgtheta3(boxi,jj)/
c     &                 dble(tconfig)
c                  avgtheta4(boxi,jj)=avgtheta4(boxi,jj)/
c     &                 dble(tconfig)
                  
                  rx=binwidth2*(dble(jj)-0.5d0)
                  
                  if (avgtheta1(boxi,jj) .gt. 0.0) then
c                     write(6,*) 'writing to files '
c                    subtract 1 so that costheta again goes from -1 to 1
                     write(30,*) rx-1.0d0, avgtheta1(boxi,jj)
                     write(31,*) rx-1.0d0, avgtheta2(boxi,jj)
c                     write(6,*) rx-1.0d0, avgtheta2(boxi,jj)
c                     write(32,*) rx-1.0d0, avgtheta3(boxi,jj)
c                     write(33,*) rx-1.0d0, avgtheta4(boxi,jj)
                  endif
               enddo
            endif               !end if (langle)
            
cc         endif                  !end if not lskip(i)

            if (ldipole) then
cc - write out average dipole moment
cc - convert to Debye: 1 e*A = 1.602e-29 C*m = 4.802 D
               avgdip = avgdip/dble(dipcount)
               avgz = avgz/dble(avgzcount)
               avgs = 0.50d0*(avgs/dble(avgzcount))

               write(6,*) 'average dipole moment: ', avgdip*4.802d0
               write(6,*) 'cos average z angle: ', avgz 
               write(6,*) 'avg z angle (deg): ', dacos(avgz)*
     &              (180.0d0/3.14159d0)
               write(6,*) 'avg order parameter: ', avgs

               open(32,file='zangle.dat')
               write(6,*) 'tconfig ', tconfig
               do jj=minbin2,maxbin2

                  avgzangle(jj) = avgzangle(jj)/dble(tconfig)
                  rx=binwidth2*dble(jj-0.5d0)
                  if (avgzangle(jj).gt.1d-10) then
                     write(32,*) rx-1.0d0,avgzangle(jj)
c                     write(6,*) rx-1.0d0, avgzangle(jj)
                  endif
               enddo

            endif

c     write a blank line at the end of each file so that xmgrace displays all points
      
         write(21,*)
         write(22,*)
         write(23,*)
         write(24,*)
         write(25,*)
         write(26,*)
         write(27,*)
         write(28,*)
         write(12,*)
         write(13,*)
         write(14,*)
         write(15,*)
         write(30,*)
         write(31,*)
         write(32,*)
         write(33,*)
  
         close(21)
         close(22)
         close(23)
         close(24)
         close(25)
         close(26)
         close(27)
         close(28)
         close(12)
         close(13)
         close(14)
         close(15)
         close(30)
         close(31)
         close(32)
         close(33)
         
      enddo                     ! end second loop over boxes


      
      end
