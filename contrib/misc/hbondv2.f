      program hbondanal
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  modified version of Kelly's benzoic acid h-bond and cluster analysis   cccc
ccc  assumes 1 molecule type
ccc  slightly more sophisticated than v1
ccc  attempts to work for multiple molecule types, but not tested
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer bmax,nmax,numax,ntmax,nbxmax,nntype,nvmax,fmax, nclusmax
      parameter (bmax=20000,nmax=3000,numax=20,ntmax=5,nbxmax=3)
      parameter (nntype=200,nvmax=50,fmax=120,nclusmax=500)
      integer i,ii,j,jj,nfile,nconfig,nchain,nmolty,nbox,nbeadty
      integer beadlist(nntype),nunit(ntmax),nvib(ntmax,numax),ibox,jbox
      integer ijvib(ntmax,numax,nvmax),ntor(ntmax,numax),nothers
      integer ijtor2(ntmax,numax,nvmax),ijtor3(ntmax,numax,nvmax)
      integer nc,step,ncmt(nbxmax,ntmax),idum,moltyp(nmax),imolty,jmolty
      integer nboxi(nmax),beadty(ntmax,numax),ijtor4(ntmax,numax,nvmax)
      integer hbond(nbxmax,nmax,numax,numax),
     &     hneigh(nmax,numax,numax,100)
      integer hnbd(nmax,4,nmax),iii,jjj,total_meth,
     &     total_solute(nbxmax,numax)
      integer hbond_0(nbxmax),hbond_1(nbxmax),hbond_2(nbxmax),
     &     hbond_3(nbxmax),hbond_4(nbxmax)
      integer hbond_sol(nbxmax),hbond_m3(ntmax),tot_mol,tot_agg
      integer hbond_m0(ntmax),hbond_m1(ntmax),hbond_m2(ntmax)

      integer k,kk,ij,count,cnta,cntb,sum_hbond(nbxmax,nmax),
     &     shbond(nmax,nvmax)
      integer icluster(nmax),nclus,nj(nvmax),nk(500),cluster(nclusmax)
      integer tally,size(nmax),iloop,icyclic(500),ilinear(500),cnt1
      integer cnt2,cnt3,ibranch(nclusmax),cl_mol(nclusmax,nclusmax),
     &     maxtally
      integer cl_shape(nclusmax),big_cluster(nclusmax),
     &     big_arch(nclusmax,3)
      integer sol_cluster(nclusmax),sol_arch(nclusmax,3),sol_agg,sol_mol
      integer solu_cnt(50,50),solv_cnt(50,50),tot

      logical lclus(nmax,nvmax,nclusmax),lcount(nmax),ldone
      logical lneigh,lsolcnt(nclusmax),lhb(nbxmax,ntmax,ntmax),loo
      logical lwait(nbxmax,ntmax,ntmax)

      double precision rcut(nbxmax),boxlx(nbxmax),boxly(nbxmax)
      double precision boxlz(nbxmax),xcm(nmax),ycm(nmax),zcm(nmax)
      double precision xcoord(nmax,numax),ycoord(nmax,numax),ddum
      double precision zcoord(nmax,numax),rxij,ryij,rzij,rij,OOsq
      double precision hdist(nmax,4,nmax),rij_h1(ntmax,ntmax),
     &     rij_h2(ntmax,ntmax),OHsq
      double precision rxij_h,ryij_h,rzij_h,thetac_h1a,thetac_h2a
      double precision avg_meth,avg_sol,aa,bb,rij_h1a,rij_h2a
      double precision thetac_h1(ntmax,ntmax),thetac_h2(ntmax,ntmax),
     &     numdens(nvmax), banumdens(nvmax,nvmax)

      integer a_0(nbxmax),a_1(nbxmax),a_2(nbxmax),a_3(nbxmax),
     &     a_4(nbxmax),atotal,a_m0,a_m1,a_m2,a_m3,cnt_solu

      character*40 filelist(fmax)
      character*2 char(10)

      integer nH, nO, ilist, tconfig
      double precision hu(numax), ou(numax),htyp(numax), otyp(numax)

      OOsq = 3.3d0 ! defining H bond as OO dist < 3.3 Ang
      OHsq = 2.5d0 ! defining H bond as OH dist < 2.5 Ang

      atotal = 0
      tot_mol = 0
      tconfig=0



      open (4,file='hbondinput')
c      write(6,*)'Enter number of files'
      read(4,*)
      read(4,*) nfile
      read(4,*)
      do i=1,nfile
         read(4,*) filelist(i)
      enddo
      read(4,*)
      read(4,*) nbox
c      write(6,*) 'enter number of hydrogens '
      read(4,*)
      read(4,*) nH
      read(4,*)
      do i=1,nH
         read(4,*) htyp(i), hu(i)
      enddo
c      write(6,*) 'enter number of oxygens'
      read(4,*)
      read(4,*) nO
      read(4,*)
      do i=1,nO
         read(4,*) otyp(i), ou(i)
      enddo
      

      do i=1,nbox
         hbond_sol(i) = 0
         hbond_0(i) = 0
         hbond_1(i) = 0
         hbond_2(i) = 0
         hbond_3(i) = 0
         hbond_4(i) = 0         
         a_0(i) = 0
         a_1(i) = 0
         a_2(i) = 0
         a_3(i) = 0
         a_4(i) = 0
         do j=1,nmolty
            total_solute(i,j) = 0
         enddo
      enddo

c      do j=1,nbox
c         do i=1,nchain
c            ibox = nboxi(i)
c            sum_hbond(ibox,i) = 0
c         enddo
c      enddo
      
      do ilist=1,nfile

         open(10,file=filelist(ilist))
         write(6,*) 'file ', filelist(ilist)

C     Read header
         read(10,*) nconfig,nchain,nmolty,(rcut(i),i=1,nbox)
         tconfig = tconfig + nconfig
         read(10,*) nbeadty,(beadlist(ii),ii=1,nbeadty)
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


C        Loop over configurations
         do 100 nc=1,nconfig
C     Read bead coordinates
            read(10,*) step
c         write(6,*) 'step',step
         
            do i=1,nbox
               read(10,*) (ncmt(i,j), j=1,nmolty)
               read(10,*) boxlx(i),boxly(i),boxlz(i)
            enddo
            
            do i=1,nchain
               read(10,*) idum,moltyp(i),
     &              nunit(moltyp(i)),nboxi(i),xcm(i),ycm(i),zcm(i)
c            write(6,*) idum, moltyp(i),
c     &           nunit(moltyp(i)),nboxi(i),xcm(i),ycm(i),zcm(i)

               if (idum.ne.i) stop 
               do ii=1,nunit(moltyp(i))
                  read(10,*) xcoord(i,ii),ycoord(i,ii),
     &                 zcoord(i,ii),ddum,beadty(moltyp(i),ii)
               enddo
            enddo

c    set to zero each time
            do j=1,nbox
               do i=1,nchain
                  do ii=1,nO
                     do jj=1,nmolty
                        hbond(j,i,ou(ii),jj) = 0
                     enddo
                  enddo
               enddo
            enddo
            

c     sum number of solutes to get average number at end
            do j=1,nbox
               do jj=1,nmolty
                  total_solute(j,jj) = total_solute(j,jj)+ncmt(j,jj)
c               write(6,*) 'total_solute(',i,')', total_solute(i,j)
               enddo
            enddo


c     loop over all particles and evaluate O-O distances
c         write(66,*)
            do 23 i=1,nchain-1
               imolty = moltyp(i)
               ibox = nboxi(i)

               do 21 j=i+1,nchain
                  jmolty = moltyp(j)
                  jbox = nboxi(j)

                  if(jbox.ne.ibox) goto 21
c have no idea what this is for
c initialize to large distance to prevent overcounting??
               do ii=1,nO
                 do jj=1,nO
                    rij_h1(ou(ii),ou(jj)) = 10.0d0
                    rij_h2(ou(ii),ou(jj)) = 10.0d0
                  enddo
               enddo

               do k=1,nbox
                  do ii=1,nO
                     do jj=1,nO
                        lhb(k,ou(ii),ou(jj)) = .false.
                        lwait(k,ou(ii),ou(jj))=.false.
                     enddo
                  enddo
               enddo

               loo=.false.
               do ii = 1,nO
c                  if (moltyp(ou(ii)).eq.imolty) then
                     do jj = 1,nO
c                        if (moltyp(ou(jj)).eq.jmolty) then
                           rxij = xcoord(i,ou(ii))-xcoord(j,ou(jj))
                           ryij = ycoord(i,ou(ii))-ycoord(j,ou(jj))
                           rzij = zcoord(i,ou(ii))-zcoord(j,ou(jj))
C     mimage
                           rxij=rxij-boxlx(ibox)*anint(rxij/boxlx(ibox))
                           ryij=ryij-boxly(ibox)*anint(ryij/boxly(ibox))
                           rzij=rzij-boxlz(ibox)*anint(rzij/boxlz(ibox))
                     
C     distance
                           rij=dsqrt(rxij**2+ryij**2+rzij**2)

                        
                           if(rij.lt.OOsq) then
c     now want hydrogen to be between also
                              do iii=1,nH
c     ----                     oxygen chain j/hydrogen chain i
                                 rxij_h = xcoord(j,ou(jj))-
     &                                xcoord(i,hu(iii))
                                 ryij_h = ycoord(j,ou(jj))-
     &                                ycoord(i,hu(iii))
                                 rzij_h = zcoord(j,ou(jj))-
     &                                zcoord(i,hu(iii))
C     mimage
                                 rxij_h=rxij_h-boxlx(ibox)*
     &                                anint(rxij_h/boxlx(ibox))
                                 ryij_h=ryij_h-boxly(ibox)*
     &                                anint(ryij_h/boxly(ibox))
                                 rzij_h=rzij_h-boxlz(ibox)*
     &                                anint(rzij_h/boxlz(ibox))
c                              write(6,*) 'rxij_h ', rxij_h
C     distance
                                 rij_h1(ou(ii),ou(jj))=sqrt(rxij_h**2
     &                                +ryij_h**2+rzij_h**2)

                                    
                                 thetac_h1(ou(ii),ou(jj)) = 
     &                                (rxij_h*rxij+ryij_h*ryij+
     &                                rzij_h*rzij)
     &                                /(rij_h1(ou(ii),ou(jj))*rij)

c     ----                         hydrogen on chain j/oxygen chain i
                                 rxij_h = xcoord(j,hu(iii))
     &                                -xcoord(i,ou(ii))
                                 ryij_h = ycoord(j,hu(iii))
     &                                -ycoord(i,ou(ii))
                                 rzij_h = zcoord(j,hu(iii))
     &                                -zcoord(i,ou(ii))
C     mimage
                                 rxij_h=rxij_h-boxlx(ibox)*
     &                                anint(rxij_h/boxlx(ibox))
                                 ryij_h=ryij_h-boxly(ibox)*
     &                                anint(ryij_h/boxly(ibox))
                                 rzij_h=rzij_h-boxlz(ibox)*
     &                                anint(rzij_h/boxlz(ibox))
                                 
C     distance
                                 rij_h2(ou(ii),ou(jj))=sqrt(rxij_h**2
     &                                +ryij_h**2+rzij_h**2)

                                 
                                 thetac_h2(ou(ii),ou(jj))=(rxij_h*
     &                                rxij+ryij_h*ryij+rzij_h*rzij)/
     &                                (rij*rij_h2(ou(ii),ou(jj)))


                                 loo = .true.

                                 if ((rij_h1(ou(ii),ou(jj)).lt.OHsq
     &                                .and.thetac_h1(ou(ii),ou(jj))
     &                                .lt.-0.1d0) .or.(rij_h2(ou(ii),
     &                                ou(jj)).lt.OHsq.and.thetac_h2(
     &                                ou(ii),ou(jj)).lt. -0.1d0)) then
                                    

                                    lhb(ibox,ou(ii),ou(jj)) =.true.
                                    hbond(ibox,i,ou(ii),jmolty) = 
     &                                   hbond(ibox,i,ou(ii),jmolty) + 1
                                    hbond(ibox,j,ou(jj),imolty) = 
     &                                   hbond(ibox,j,ou(jj),imolty) + 1
                                    hbond_sol(ibox) = hbond_sol(ibox)+2

                                 endif

                              enddo  ! end loop iii=1,nH
c                           write(6,*) 'loo true'
                           endif  ! end if rij .lt. OOsq
c                        endif

                     enddo     ! end loop jj=1,nO
c                  endif
               enddo           ! end loop ii=1,nO

c$$$c have computed ii,jj combinations; now see where Hbonds exist

c$$$               if(loo) then
c$$$                  do ii=1,nO
c$$$                     do jj=1,nO
c$$$                        if((rij_h1(ou(ii),ou(jj)).lt.OHsq.and.
c$$$     &                       thetac_h1(ou(ii),ou(jj)).lt.-0.4d0)
c$$$     &                       .or.(rij_h2(ou(ii),ou(jj)).lt.OHsq
c$$$     &                       .and.thetac_h1(ou(ii),ou(jj)).lt.
c$$$     &                       -0.4d0)) then
c$$$                           lhb(ibox,ou(ii),ou(jj)) =.true.
c$$$c -- this is what I would have done (counts all)
c$$$c -- could move this up into the jj loop
c$$$                           hbond(ibox,i,ou(ii),jmolty) = 
c$$$     &                          hbond(ibox,i,ou(ii),jmolty) + 1
c$$$                           hbond(ibox,j,ou(jj),imolty) = 
c$$$     &                          hbond(ibox,j,ou(jj),imolty) + 1
c$$$                           hbond_sol(ibox) = hbond_sol(ibox)+2             
c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$               endif  ! end if loo (remove if count Kelly's way)


c --- this agrees with Kelly's program, but I don't know why
c --- the ii,ii and jj,jj interactions are exluded if
c --- there is an H-bond with ii,jj or jj,ii

c$$$                  do ii=1,nO
c$$$                     do jj=1,nO
c$$$                        if (ou(ii).ne.ou(jj)) then
c$$$                           if (lhb(ibox,ou(ii),ou(jj))) then
c$$$                              hbond(ibox,i,ou(ii),jmolty) = 
c$$$     &                             hbond(ibox,i,ou(ii),jmolty) + 1
c$$$                              hbond(ibox,j,ou(jj),imolty) = 
c$$$     &                             hbond(ibox,j,ou(jj),imolty) + 1
c$$$                              hbond_sol(ibox) = hbond_sol(ibox)+2
c$$$                              lwait(ibox,ou(ii),ou(jj))=.true.
c$$$                              lwait(ibox,ou(ii),ou(ii)) = .true.
c$$$                              lwait(ibox,ou(jj),ou(jj)) = .true.
c$$$c                              write(9,*) 'hbond ', ou(ii),ou(jj)
c$$$                           endif
c$$$                        endif     
c$$$                     enddo
c$$$                  enddo
c$$$  
c$$$
c$$$                  do ii=1,nO
c$$$                     do jj=1,nO
c$$$                        if (lhb(ibox,ou(ii),ou(jj)).and..not.
c$$$     &                       lwait(ibox,ou(ii),ou(jj))) then
c$$$c                           write(6,*) 'lwait hbond'
c$$$                           hbond(ibox,i,ou(ii),jmolty) = 
c$$$     &                          hbond(ibox,i,ou(ii),jmolty) + 1
c$$$                           hbond(ibox,j,ou(jj),imolty) = 
c$$$     &                          hbond(ibox,j,ou(jj),imolty) + 1
c$$$                           hbond_sol(ibox) = hbond_sol(ibox)+2
c$$$                        endif
c$$$                     enddo
c$$$                  enddo
c$$$               endif     ! end if loo         




 21         enddo               ! end loop j
c     write(6,*) 'chain i ', i

 23      enddo                  ! end loop i                  

         do i=1,nchain
            ibox = nboxi(i)
            imolty = moltyp(i)
            sum_hbond(ibox,i) = 0
            do ii=1,nO
c               do jj=1,nO
                  sum_hbond(ibox,i)= sum_hbond(ibox,i)+
     &                 hbond(ibox,i,ou(ii),imolty)
c               write(6,*) 'sum_hbond(',ibox,i,')', sum_hbond(ibox,i)
c                  k=0
c$$$                  do kk=1,hbond(ibox,i,ou(ii),ou(jj))
c$$$                     k=k+1
c$$$                     shbond(i,k) = hneigh(i,ou(ii),ou(jj),kk)
c$$$                  enddo
c               enddo
            enddo
         enddo

c     tabulate non vs 1 vs 2 hbonded solutes
c     write(32,*) nc
         do 30 i=1,nchain
            atotal = 0
            ibox = nboxi(i)
            imolty = moltyp(i)
            do ii=1,nO
               atotal = atotal + hbond(ibox,i,ou(ii),imolty)
c                  write(6,*) 'hbond ', hbond(ibox,i,ou(ii),imolty)
c               write(6,*) 'atotal ', atotal
            enddo
c            write(6,*) 'atotal ', atotal
c            atotal = hbond(ibox,i,2,2) + hbond(ibox,i,2,3)
c     &            + hbond(ibox,i,3,2) + hbond(ibox,i,3,3)
c            write(3,*) 'atotal ', atotal
            if(atotal.gt.3) then
               a_4(ibox) = a_4(ibox)+1
            elseif(atotal.gt.2) then
               a_3(ibox) = a_3(ibox)+1
            elseif(atotal.gt.1) then
               a_2(ibox) = a_2(ibox)+1
            elseif(atotal.gt.0) then
               a_1(ibox) = a_1(ibox)+1
            elseif(atotal.eq.0) then
               a_0(ibox) = a_0(ibox) + 1
            else
               write(6,*) 'What??',i,' solute'
            endif
 30      continue
c         write(32,*)
c         write(33,*)

c     tabulate non vs 1 vs 2 hbonded solutes
         do 60 i=1,nchain
            ibox=nboxi(i)
            ii=sum_hbond(ibox,i)
            jj=sum_hbond(ibox,i)

c            solu_cnt(ii+1,jj+1)=solu_cnt(ii+1,jj+1)+1 

            if(sum_hbond(ibox,i).gt.3) then
               hbond_4(ibox) = hbond_4(ibox)+1
            elseif(sum_hbond(ibox,i).gt.2) then
               hbond_3(ibox) = hbond_3(ibox)+1
            elseif(sum_hbond(ibox,i).gt.1) then
               hbond_2(ibox) = hbond_2(ibox)+1
            elseif(sum_hbond(ibox,i).gt.0) then
               hbond_1(ibox) = hbond_1(ibox)+1
            elseif(sum_hbond(ibox,i).eq.0) then
               hbond_0(ibox) = hbond_0(ibox) + 1
            else
               write(6,*) 'What??',i, sum_hbond(ibox,i)
            endif
 60      continue

 100  enddo                     ! end loop over configurations
      close(10)
      enddo            ! end loop over files

               
c determine average number of hydrogen bonds per molecule type
      do ibox=1,nbox
         write(6,*) 'BOX ', ibox
c -- assume 1 mol type
         tot = dble(total_solute(ibox,1))
         tot_mol = tot_mol + tot
         avg_sol = dble(hbond_sol(ibox))/tot
         write(6,*)
         write(6,*) 'Solute',hbond_sol(ibox),total_solute(ibox,1)
         write(6,1000) 'Average number of hbonds per solute',
     &        avg_sol
         write(6,*)

         write(6,*)
         write(6,*) 'Solute breakdown',a_0(ibox),a_1(ibox),a_2(ibox),
     &        a_3(ibox),a_4(ibox)
         write(6,1000) 'Fraction of solute with no hbonds',
     &        dble(a_0(ibox))/tot
         write(6,1000) 'Fraction of solute with one hbond',
     &        dble(a_1(ibox))/tot
         write(6,1000) 'Fraction of solute with two hbonds',
     &        dble(a_2(ibox))/tot
         write(6,1000) 'Fraction of solute with three hbonds',
     &        dble(a_3(ibox))/tot
         write(6,1000) 'Fraction of solute with >three hbonds',
     &        dble(a_4(ibox))/tot
         write(6,*)
         write(6,*) 'Solute-solute breakdown',hbond_0(ibox),
     &        hbond_1(ibox),hbond_2(ibox),hbond_3(ibox),
     &        hbond_4(ibox)
         write(6,1000) 'Fraction of solute with no hbonds',
     &        dble(hbond_0(ibox))/dble(total_solute(ibox,1))
         write(6,1000) 'Fraction of solute with one hbond',
     &        dble(hbond_1(ibox))/dble(total_solute(ibox,1))
         write(6,1000) 'Fraction of solute with two hbonds',
     &        dble(hbond_2(ibox))/dble(total_solute(ibox,1))
         write(6,1000) 'Fraction of solute with three hbonds',
     &        dble(hbond_3(ibox))/dble(total_solute(ibox,1))
         write(6,1000) 'Fraction of solute with >three hbonds',
     &        dble(hbond_4(ibox))/dble(total_solute(ibox,1))
         write(6,*)

      enddo
      
      write(6,*) 'TOTAL MOL ', tot_mol
 300  format(A40,I7,F10.5)
 1000 format(A40,F10.5)

      stop
      end

