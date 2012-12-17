      program RDF

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   modified version of Becky's RDF code
ccc   modifed on 01/31/07
ccc   removed slab stuff
ccc   not normalizing properly for NpT-Gibbs acrylate mixtures
ccc   removed lskip to try and fix this
ccc   10/09/08
ccc   correctly averages the number integrals and RDFs 10/16/08 KM
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer bmax,nmax,numax,ntmax,nbxmax,nntype,nvmax,fmax
      parameter (bmax=20000,nmax=1000,numax=20,ntmax=12,nbxmax=3)
      parameter (nntype=450,nvmax=30,fmax=100)
      integer i,ii,j,jj,k,kk,l,ll,nbox,btype1,btype2,nconfig,nchain
      integer nmolty,nbeadty,beadlist(nntype),nunit(ntmax),im
      integer nvib(ntmax,numax),step,ijvib(ntmax,numax,nvmax)
      integer ntor(ntmax,numax),ijtor2(ntmax,numax,nvmax)
      integer ijtor3(ntmax,numax,nvmax),ijtor4(ntmax,numax,nvmax)
      integer ncmt(nbxmax,ntmax),nc,idum,imolty(nmax),nboxi(nmax)
      integer boxi,boxj,bin(nbxmax),moltyj,maxbin(nbxmax),utype1,utype2
      integer beadcount(ntmax,nntype),zct,mtype1,mtype2,moltyi
      integer nfile,ilist,tconfig,x,xx,blnk,blnk2,beadty(ntmax,numax)
      double precision pi,rcut(nbxmax),boxlx(nbxmax),prev,nct1(nbxmax)
      double precision boxlz(nbxmax),rlower,const,boxly(nbxmax)
      double precision xcm(nmax),ycm(nmax),zcm(nmax),xcoord(nmax,numax)
      double precision ycoord(nmax,numax),zcoord(nmax,numax),ddum,rupper
      double precision maxR(nbxmax),xij,yij,zij,rij,vol,dens,binlist
      double precision nideal,g(nbxmax,bmax),binwidth(nbxmax),rx,zcom,
     &     ni,zsum,zmax(nbxmax),zmin(nbxmax),n(nbxmax,bmax),nct2(nbxmax)
      double precision nint(nbxmax,bmax), volavg(nbxmax), MAX(nbxmax)
      double precision num(nbxmax), gavg(nbxmax,bmax), navg(nbxmax,bmax)
      double precision nct1avg(nbxmax), nct2avg(nbxmax)
      double precision numskip(nbxmax)
      character*40 filename,filename2,filelist(fmax),sfilename,
     &     bfilename,sfilename2,bfilename2
      character*4 atom(nntype),astr1,astr2
      logical lsymm,lcalc,lcalc2,lslab,lskip(nbxmax),lset1,lset2
      parameter (lsymm=.false., lslab=.false.)
      
      logical meth(ntmax,numax)

      pi=4.0d0*atan(1.0d0)
      lcalc = .true.
      lcalc2 = .true.

      atom(3)='CH4 '
      atom(407) = 'CH3'
      atom(395) = 'CH3'
 
      atom(57)='chg '
      atom(58)='chg '
      atom(59)='C   '
      atom(61)='H   '
      atom(62)='O   '
      atom(63)='CH3 '
      atom(64)='CH2 '
      atom(71)='O   '
      atom(72)='CH3 '
      atom(92)='H   '
      atom(99)='CH3 '
      atom(100)='C   '
      atom(101)='O   '
      atom(102)='O   '
      atom(103)='H   '
      atom(107)='O   '
      atom(108)='H   '
      atom(109)='O   '
      atom(110)='H   '
      atom(111)='O   '
      atom(112)='H   '
      atom(113)='M   '
      atom(114)='O   '
      atom(115)='H   '
      atom(116)='M   '
      atom(117)='O   '
      atom(118)='H   '
      atom(119)='LP  '
      atom(124)='HE  '
      atom(129)='C   '
      atom(130)='O   '
      atom(175)='O   '
      atom(176)='H   '
c      atom(1) = 'O   '
c      atom(2) = 'H   '
c      atom(3) = 'BC  '
c      atom(4) = 'LP  '
      atom(177)='BC  '
      atom(178)='BC  '
      atom(179)='BC  '
      atom(180)='BC  '
      atom(181)='BC  '
      atom(182)='BC  '
      atom(183)='BC  '
      atom(184)='BC  '
      atom(185)='BC  '
      atom(186)='LP  '
      atom(187)='LP  '
      atom(188)='LP  '
      atom(189)='LP  '
      atom(190)='LP  '
      atom(191)='LP  '
      atom(192)='LP  '
      atom(193)='LP  '
      atom(194)='LP  '
      atom(195)='O  '
      atom(196)='H  '
      atom(197)='N '
      atom(198)='O '

      atom(4)='CH3'
      atom(5)='CH2'
 
      atom(210)='CH3'
      atom(211)='O1 '
      atom(212)='C  '
      atom(213)='O2 '
      atom(214)='CH '
      atom(215)='CH2'
      atom(216)='CH3'
      atom(217)='CH '
      atom(218)='CH2'

c      atom(176)='Cl  '
c      atom(179)='Na  '
c      atom(182)='Aa  '
c      atom(183)='Bb  '
c      atom(184)='Cc  '
c      atom(185)='Dd  '
c      atom(186)='Ee  '
c      atom(187)='Ff  '
c      atom(188)='Gg  '
c      atom(189)='Hh  '
c      atom(190)='Ii  '
c      atom(191)='CIN '

      open(4,file='rdfinput')
c      write(6,*)'Enter number of files'
      read(4,*)
      read(4,*) nfile
      write(6,*) nfile, ' input files'
      read(4,*)
      do i=1,nfile
         read(4,*) filelist(i)
         write(6,*) filelist(i)
      enddo
c      write(6,*)'Enter number of simulation boxes'
      read(4,*)
      read(4,*) nbox
      read(4,*)
      write(6,*) nbox, ' boxes'
      do i=1,nbox
c         write(6,*)'Enter binwidth for box',i
         read(4,*) binwidth(i)
      enddo
      write(6,*)'Enter type of RDF - atom1,atom2'
      write(6,*)'(use beadtype numbers from suijtab.f)'
      read(5,*) btype1,btype2

C     set histogram counters and logicals
      lset1=.false.
      lset2=.false.
      tconfig=0

      do i=1,nbox
         nct1(i)=0.0d0
         nct2(i)=0.0d0
         nct1avg(i)=0.0d0
         nct2avg(i)=0.0d0
         MAX(i)=0.0d0
         num(i)=0.0d0
         numskip(i)=0.0d0
         do k=1,bmax
            n(i,k)=0.0d0
            g(i,k)=0.0d0
            nint(i,k)=0.0d0
            gavg(i,k)=0.0d0
            navg(i,k)=0.0d0
         enddo

      lskip(i)=.false.

      enddo 

      do ilist=1,nfile

         open(10,file=filelist(ilist))

C     Read header
         read(10,*) nconfig,nchain,nmolty,(rcut(i),i=1,nbox)
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


            
C     Read bead coordinates
            read(10,*) step

            do i=1,nbox
               read(10,*) (ncmt(i,j), j=1,nmolty)
               read(10,*) boxlx(i),boxly(i),boxlz(i)
c               write(6,*)'boxlength', boxlx(i), boxly(i), boxlz(i)


C              set-up bins
               maxR(i)=boxlz(i)
               maxbin(i)=int(maxR(i)/binwidth(i))
c               write(6,*) 'maxbin', maxbin(i)
               if (MAX(i) .lt. maxbin(i)) then
                  MAX(i) = maxbin(i)
               endif

cccc  set to zero each time
               nct1(i)=0.0d0
               nct2(i)=0.0d0
               num(i)=0.0d0
               do k=1,maxbin(i)
                  n(i,k)=0.0d0
                  nint(i,k)=0.0d0
                  g(i,k)=0.0d0
               enddo

            enddo  ! end loop over boxes

            do i=1,nchain
               read(10,*) idum,imolty(i),
     &              nunit(imolty(i)),nboxi(i),xcm(i),ycm(i),zcm(i)
               if (idum.ne.i) stop 
               do ii=1,nunit(imolty(i))
                  read(10,*) xcoord(i,ii),ycoord(i,ii),
     &                 zcoord(i,ii),ddum,beadty(imolty(i),ii)
c                  write(6,*) xcoord(i,ii),ycoord(i,ii),
c     &                 zcoord(i,ii),ddum,beadty(imolty(i),ii)
c     write(6,'(4(1x,F10.6),I5)') xcoord(i,ii),ycoord(i,ii),
c     &                 zcoord(i,ii),ddum,beadty(imolty(i),ii)
c     this was included to prevent including too many methyl groups in
c     methacrylate RDFs
c                  if (ii .eq. 8) then
c                     meth(imolty(i),ii) = .true.
c                  else
c                     meth(imolty(i),ii) = .false.
c                  endif

C          Determine the types of molecules interacting 
C            based on specified bead types
                  if (.not.lset1) then
c     write(6,*) 'beadty()', beadty(imolty(i),ii)
c                     write(6,*) 'btype1', btype1
                     if (beadty(imolty(i),ii).eq.btype1)then
c     &                    .and. .not. meth(imolty(i),ii)) then
                        mtype1=imolty(i)
                        utype1=ii
                        lset1=.true.
c                        write(6,*) 'btype1', btype1
c                        write(6,*) 'lset1', lset1
                     endif
                  endif
                  if (.not.lset2) then
c                     write(6,*)'beadty()', beadty(imolty(i),ii)
c                     write(6,*)'btype2', btype2
                     if (beadty(imolty(i),ii).eq.btype2)then 
c     &                    .and. .not. meth(imolty(i),ii)) then
                        mtype2=imolty(i)
                        utype2=ii
                        lset2=.true.
c                        write(6,*)'btype2',btype2
c                        write(6,*) 'lset2', lset2
                     endif
                  endif

               enddo  ! end loop over units (ii)
            enddo     ! end loop over molecules (i)
            
c            write(6,*)'mtypes',mtype1,mtype2
c            write(6,*)'utypes',utype1,utype2
            
            if (.not.lset1 .or. .not.lset2) stop

C       Determine which boxes/molecule types can be skipped
C       i.e., neither mtype1 or mtype2 are currently in this box, 
C       or this type does not contain the specified beads
c            do i=1,nbox
c               write(6,*)'i',i,' ncmt1',ncmt(i,mtype1),' ncmt2',
c     &              ncmt(i,mtype2)
c               if (ncmt(i,mtype1).eq.0.or.ncmt(i,mtype2).eq.0) then
c                  lskip(i)=.true.
c                  write(6,*)i,lskip(i)
c               endif
c            enddo

C     For first configuration, determine number of beads of each type
C     for each molecule type (SPC/E has 2 beads of type 108, etc)
            if (ilist.eq.1.and.nc.eq.1) then
               if (mtype1.eq.mtype2) then
                  im=mtype1
c     write(6,*)'im', im, 'mtype1', mtype1
                  do j=1,nunit(im)
                     do k=1,nbeadty
                        if (beadty(im,j).eq.beadlist(k)) then
                           beadcount(im,beadlist(k))=
     &                          beadcount(im,beadlist(k))+1
                        endif
                     enddo
                  enddo
               else
                  do i=1,2
                     if (i.eq.1) im=mtype1
                     if (i.eq.2) im=mtype2
                     do j=1,nunit(im)
                        do k=1,nbeadty
                           if (beadty(im,j).eq.beadlist(k)) then
                              beadcount(im,beadlist(k))=
     &                             beadcount(im,beadlist(k))+1
                           endif
                        enddo
                     enddo
                  enddo
               endif
c$$$               write(6,*)'beadcount'
c$$$               do i=1,2
c$$$                  if (i.eq.1) im=mtype1
c$$$                  if (i.eq.2) im=mtype2
c$$$                  do k=1,nbeadty
c$$$                     write(6,*) im, beadlist(k), 
c$$$     &                    beadcount(im,beadlist(k))
c$$$                  enddo
c$$$               enddo
            endif

C     Calculate specific pair distances

C           Loop over chains - btype1 becomes the origin bead
C            rdf is directional btype1 to btype2, 
C            although btype2 to btype1 should be same
C            number integrals will be different depending on direction
            do l=1,nchain

               moltyi=imolty(l)
               boxi=nboxi(l)

C              Count number of molecules that could interact (mtype2)
               if (moltyi.eq.mtype2) then
                  nct2(boxi)=nct2(boxi)+1.0d0
               endif
C              Count number of molecules that could interact (mtype1)
               if (moltyi.eq.mtype1) then
                  nct1(boxi)=nct1(boxi)+1.0d0
               endif


C              First find molecule type 1, then set btype1
c               if (moltyi.eq.mtype1 .and. .not.lskip(boxi)) then
               if (moltyi.eq.mtype1 ) then
                  do  k=1,nunit(moltyi)

                     if (beadty(moltyi,k).eq.btype1)then
c     &                    .and..not.meth(imolty(l),k)) then

                        
C                    Loop over interacting chains - only btype2 can interact
                        do ll=1,nchain
c                           write(6,*)'ll', ll
                           moltyj=imolty(ll)
                           boxj=nboxi(ll)

C                       Check that the interacting molecule is mtype2 
C                         and is in the same box

                           if (moltyj.eq.mtype2.and.boxi.eq.boxj
     &                          .and.l.ne.ll) then


                             
                              do kk=1,nunit(moltyj)
c                              write(6,*) 'beadty ', beadty(moltyj,kk),
c     &                             ' btype2 ', btype2 
                                 if(beadty(moltyj,kk).eq.btype2)then
c     &                                .and..not.meth(imolty(ll),kk))then
c                                    write(6,*) 'check interactions'

C                          Calculate distances
                                    xij=xcoord(l,k)-xcoord(ll,kk)
                                    yij=ycoord(l,k)-ycoord(ll,kk)
                                    zij=zcoord(l,k)-zcoord(ll,kk)

C                          Minimum image
                                    xij=xij-boxlx(boxi)*anint
     &                                   (xij/boxlx(boxi))
                                    yij=yij-boxly(boxi)*
     &                                   anint(yij/boxly(boxi))
                                    zij=zij-boxlz(boxi)*
     &                                   anint(zij/boxlz(boxi))

C                          Bin distances
                                    rij=dsqrt(xij**2+yij**2+zij**2)
                                    bin(boxi)=int(rij/binwidth(boxi))+1
                                    
                                    if (bin(boxi).le.maxbin(boxi)) then

                                       n(boxi,bin(boxi))=n(boxi,
     &                                      bin(boxi))+1.0d0

                                    endif

                                 endif
                              enddo ! end loop over units kk
                           endif 
                        enddo   ! end loop over molecules ll
                     endif
                  enddo         ! end loop over units k
               endif
            enddo               ! end loop over molecules l
 

C     Normalize to ideal gas
            do boxi=1,nbox

               
ccc   calculate number integral            
               do j=1,maxbin(boxi)
                  
                  if (nct1(boxi) .gt. 0.0d0) then

                     num(boxi)=num(boxi)+n(boxi,j)/nct1(boxi)

                  endif
                  
                  nint(boxi,j)=num(boxi)
                  
               enddo
               
ccc   calculate rdf
               vol=boxlx(boxi)*boxly(boxi)*boxlz(boxi)
               
               dens=nct2(boxi)*dble(beadcount(mtype2,btype2))/
     &              vol
               const=4.0d0*pi*dens/3.0d0

               do j=1,maxbin(boxi)
                  
                  rx=binwidth(boxi)*(dble(j)-0.5d0)
c     binlist=rx
                  rlower=binwidth(boxi)*dble(j-1)
                  rupper=rlower+binwidth(boxi)
                  
ccc  divide nideal by 2 for methacrylates
ccc  due to overcounting number of beads
                  nideal=const*(rupper**3-rlower**3)
                  
                  if (nct1(boxi).gt.0.00.and.nct2(boxi).gt.0.00) then

                     g(boxi,j)=n(boxi,j)/
     &                    dble(beadcount(mtype1,btype1))/
     &                    nct1(boxi)/nideal

ccc   don't divide by beadcount if methacrylates
ccc   due to overcounting number of beads

c                     g(boxi,j)=n(boxi,j)/
c     &                    nct1(boxi)/nideal
                  else
                     g(boxi,j)=0.0d0
c                     write(6,*) 'skip!'
                  endif

               enddo
               
ccc  accumulate averages

               do j=1,maxbin(boxi)

                  navg(boxi,j)=navg(boxi,j)+nint(boxi,j)
                  gavg(boxi,j)=gavg(boxi,j)+g(boxi,j)

               enddo
            
               nct1avg(boxi)=nct1avg(boxi)+nct1(boxi)
               nct2avg(boxi)=nct2avg(boxi)+nct2(boxi)
               if (nct1(boxi).eq.0.0.or.nct2(boxi).eq.0.0) then
                  numskip(boxi)=numskip(boxi)+1
c                  write(6,*) 'numskip(',boxi,')', numskip(boxi)
               endif

            enddo    ! end loop over boxes
                            
         enddo                  ! end loop over configurations

         close(10)
      enddo                     ! end loop over files
      
      do i=1,nbox

         

         nct2avg(i)=nct2avg(i)/dble(tconfig)
         nct1avg(i)=nct1avg(i)/dble(tconfig)
         write(6,*)'nctavg',nct1avg(i),nct2avg(i)

            
c         if (.not.lskip(i)) then

         astr1=atom(btype1)
         x=index(astr1,' ')
         blnk=x-1
         
         astr2=atom(btype2)
         xx=index(astr2,' ')
         blnk2=xx-1
         
c            write(6,*) 'astr1 ', astr1, ' x ', x, 
c     &           ' blnk ', blnk
c            write(6,*) 'astr2 ', astr2, 'xx ', xx,
c     &           ' blnk2 ', blnk2

C     Create output filenames
         open(9)

         if (blnk.eq.1) then
            if (blnk2.eq.1) then
               write(9,'(A7,I1,A1,A1,A1,A1,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            elseif (blnk2.eq.2) then
               write(9,'(A7,I1,A1,A1,A1,A2,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            else
               write(9,'(A7,I1,A1,A1,A1,A3,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            endif
         elseif (blnk.eq.2) then
            if (blnk2.eq.1) then
               write(9,'(A7,I1,A1,A2,A1,A1,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            elseif (blnk2.eq.2) then
               write(9,'(A7,I1,A1,A2,A1,A2,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            else
               write(9,'(A7,I1,A1,A2,A1,A3,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            endif
         else
            if (blnk2.eq.1) then
               write(9,'(A7,I1,A1,A3,A1,A1,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &                 astr2(:blnk2),'.dat'
            elseif (blnk2.eq.2) then
               write(9,'(A7,I1,A1,A3,A1,A2,A4)') 
     &                 'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            else
               write(9,'(A7,I1,A1,A3,A1,A3,A4)') 
     &              'rdf-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            endif
         endif
         close(9)
         
         open(9)
         read(9,'(A)') filename
         close(9)
         
         open(9)
         if (blnk.eq.1) then
            if (blnk2.eq.1) then
               write(9,'(A8,I1,A1,A1,A1,A1,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            elseif (blnk2.eq.2) then
               write(9,'(A8,I1,A1,A1,A1,A2,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            else
               write(9,'(A8,I1,A1,A1,A1,A3,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            endif
         elseif (blnk.eq.2) then
            if (blnk2.eq.1) then
               write(9,'(A8,I1,A1,A2,A1,A1,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            elseif (blnk2.eq.2) then
               write(9,'(A8,I1,A1,A2,A1,A2,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            else
               write(9,'(A8,I1,A1,A2,A1,A3,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            endif
         else
            if (blnk2.eq.1) then
               write(9,'(A8,I1,A1,A3,A1,A1,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            elseif (blnk2.eq.2) then
               write(9,'(A8,I1,A1,A3,A1,A2,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            else
               write(9,'(A8,I1,A1,A3,A1,A3,A4)') 
     &              'nint-box',i,'-',astr1(:blnk),'-',
     &              astr2(:blnk2),'.dat'
            endif
         endif
         close(9)
         
         open(9)
         read(9,'(A)') filename2
         close(9)
         
         open(12,file=filename)
         open(11,file=filename2)
         
c            write(6,*)'filename12 ', filename
c            write(6,*)'filename11 ', filename2
            
         
c            open(12,file='rdf.dat')
c            open(11,file='nint.dat')

         do j=1, maxbin(i)
            rx=binwidth(i)*(dble(j)-0.5d0)
            write(12,*) rx,gavg(i,j)/dble(tconfig)
            write(11,*) rx, navg(i,j)/dble(tconfig)
c            write(6,*) rx, navg(i,j)/dble(tconfig)

c            if (j.eq.MAX(i)) then
c               write(6,*)'Number integral'
c               write(6,*) nct1(i),nct2(i),navg(i,j)/dble(tconfig)
c               write(6,*) 'beadcount ', beadcount(mtype1,btype1)
c               write(6,*) 'tconfig ', tconfig
c               write(6,*) 'rx ', rx
c            endif 
         enddo
         
c     endif
      enddo                     ! end loop over boxes
      
      end



