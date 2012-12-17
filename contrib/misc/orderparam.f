      program order

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  calculates the tetrahedral order parameter
ccc  as defined in Errington and Debenedetti, Nature, 409, 318-321 (2001)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      integer bmax,nmax,numax,ntmax,nbxmax,nntype,nvmax,fmax
      integer ncmax,ios
      parameter (bmax=50000,nmax=2050,numax=20,ntmax=12,nbxmax=3)
      parameter (nntype=220,nvmax=30,fmax=50, ncmax=8000)
      integer i,ii,j,jj,k,kk,l,ll,nbox,nconfig,nchain
      integer nmolty,nbeadty,beadlist(nntype),nunit(ntmax)
      integer nvib(ntmax,numax),step,ijvib(ntmax,numax,nvmax)
      integer ntor(ntmax,numax),ijtor2(ntmax,numax,nvmax)
      integer ijtor3(ntmax,numax,nvmax),ijtor4(ntmax,numax,nvmax)
      integer ncmt(nbxmax,ntmax),nc,idum,imolty(nmax),nboxi(nmax)
      integer nfile,ilist,tconfig,beadty(ntmax,numax)
      integer count, boxi, boxj, moltyi, moltyj
      double precision rcut(nbxmax),boxlx(nbxmax),boxly(nbxmax),
     &     boxlz(nbxmax)
      double precision xcm(nmax),ycm(nmax),zcm(nmax),xcoord(nmax,numax)
      double precision ycoord(nmax,numax),zcoord(nmax,numax), ddum
      double precision rxij, ryij, rzij, rij
      double precision costheta, angleterm, max_dist
      double precision rx1, ry1, rz1, r1, rx2, ry2, rz2, r2
      double precision avgq(nbxmax), q(nbxmax)
      integer otyp, ou, ineigh
      integer neigh(nmax,4), num(nbxmax)
      logical neighbor(nbxmax,nmax,nmax)
      character *40 filelist(fmax)

c   read in fort.10

      open(4,file='orderinput')
      read(4,*)
      read(4,*) nfile
        read(4,*)
      do i=1,nfile
         read(4,*) filelist(i)
      enddo
      read(4,*)
      read(4,*) nbox
      read(4,*)  !molecule type, bead unit number for oxygen (assumes 1)
      read(4,*) otyp, ou

c     initialize accumulators

      do i=1,nbox
         avgq(i) = 0.0d0
         num(i) = 0
      enddo
      
      count = 0
      tconfig = 0

      do ilist=1,nfile

         open(10,file=filelist(ilist))
         write(6,*) 'file ', filelist(ilist)

C     Read header
         read(10,*) nconfig,nchain,nmolty,(rcut(i),i=1,nbox)
c         write(6,*) 'nconfig ', nconfig, ' nchain ', nchain
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
         do nc=1,nconfig
            count = count + 1

C     Read bead coordinates
            read(10,*) step
c         write(6,*) 'step',step
         
            do i=1,nbox
               read(10,*) (ncmt(i,j), j=1,nmolty)
               read(10,*) boxlx(i),boxly(i),boxlz(i)
            enddo
            
            do i=1,nchain
               read(10,*) idum,imolty(i),
     &              nunit(imolty(i)),nboxi(i),xcm(i),ycm(i),zcm(i)

               if (idum.ne.i) stop 
               do ii=1,nunit(imolty(i))
                  read(10,*) xcoord(i,ii),ycoord(i,ii),
     &                 zcoord(i,ii),ddum,beadty(imolty(i),ii)
               enddo
            enddo

c     set to zero at beginning of each config
            do i=1,nbox
               num(i) = 0
               q(i) = 0.0d0
            enddo

 23         do l=1,nchain

               moltyi = imolty(l)
               boxi = nboxi(l)
               
               if (moltyi.ne.otyp) goto 23
               
               num(boxi) = num(boxi) + 1

c           initialize to zero/false
               angleterm = 0.0d0
               do ineigh=1,4
                  neigh(l,ineigh) = 0
               enddo
               do ll=1,nchain
                  neighbor(boxi,l,ll)=.false.
               enddo
              

               do ineigh=1,4

c                  write(6,*) 'checking neighbor ', ineigh, ' of ', l
                  max_dist = 100.0d0
                     
 21               do ll=1,nchain
                     if (l.ne.ll) then

                        moltyj = imolty(ll)
                        
                        if (moltyj.ne.otyp) goto 21
                        boxj = nboxi(ll)

                     
                        if (boxi.eq.boxj) then
c                     write(6,*) 'll, boxj ', ll, boxj
                     
                           if (.not. neighbor(boxi,l,ll)) then

c                           write(6,*) 'checking ', l, ll

c                    calculate o-o distance
                              rxij = xcoord(l,ou) - xcoord(ll,ou)
                              ryij = ycoord(l,ou) - ycoord(ll,ou)
                              rzij = zcoord(l,ou) - zcoord(ll,ou)
                           
c                    mimage
                              rxij =rxij- boxlx(boxi)*
     &                             anint(rxij/boxlx(boxi))
                              ryij =ryij- boxly(boxi)
     &                             *anint(ryij/boxly(boxi))
                              rzij =rzij- boxlz(boxi)*
     &                             anint(rzij/boxlz(boxi))
                              
                              rij = dsqrt(rxij**2+ryij**2+rzij**2)
                     
c   find nearest oxygen neighbors
c   need four closest oxygens
c                           write(6,*) 'll, rij ', ll, rij

                              if (rij.le.max_dist) then
                                 max_dist = rij
                                 neigh(l,ineigh) = ll
c                                 write(6,*) ll, ' a possible neighbor',
c     &                                rij
                              endif
                     
                           endif ! end if not neighbor
                        endif   ! end if boxi.eq.boxj
                     endif      ! end if l.ne.ll
                  enddo         ! end loop over chains ll
                     
                  neighbor(boxi,l,neigh(l,ineigh)) = .true.
c                  write(6,*) 'oxygens ', l, neigh(l,ineigh), 
c     &                 ' are neighbors'

               enddo            ! end loop over ineigh

c    compute O-O-O angles
                  
               angleterm = 0.0d0

               do i=1,3
                  do j=i+1,4
                        
c     calculate distances
                     rx1 = xcoord(l,ou) - xcoord(neigh(l,i),ou)
                     ry1 = ycoord(l,ou) - ycoord(neigh(l,i),ou)
                     rz1 = zcoord(l,ou) - zcoord(neigh(l,i),ou)

                     rx2 = xcoord(l,ou) - xcoord(neigh(l,j),ou)
                     ry2 = ycoord(l,ou) - ycoord(neigh(l,j),ou)
                     rz2 = zcoord(l,ou) - zcoord(neigh(l,j),ou)

c     mimage
                     rx1 = rx1 - boxlx(boxi)*anint(rx1/boxlx(boxi))
                     ry1 = ry1 - boxly(boxi)*anint(ry1/boxly(boxi))
                     rz1 = rz1 - boxlz(boxi)*anint(rz1/boxlz(boxi))
                     
                     rx2 = rx2 - boxlx(boxi)*anint(rx2/boxlx(boxi))
                     ry2 = ry2 - boxly(boxi)*anint(ry2/boxly(boxi))
                     rz2 = rz2 - boxlz(boxi)*anint(rz2/boxlz(boxi))
                     

                     r1 = dsqrt(rx1**2 + ry1**2 + rz1**2)
                     r2 = dsqrt(rx2**2 + ry2**2 + rz2**2)

                     costheta = (rx1*rx2+ry1*ry2+rz1*rz2)/(r1*r2)

                     if (boxi.eq.1) then
c                        write(6,*) 'costheta ', costheta
                     endif
                     
                     if (costheta.ge.1.0d0) costheta=1.0d0
                     if (costheta.le.-1.0d0) costheta=-1.0d0
                     
                     angleterm = angleterm + 
     &                    (costheta+(1.0d0/3.0d0))**2

                  enddo
               enddo
               
               q(boxi) = q(boxi) + 1.0d0 - (3.0d0/8.0d0)*angleterm
                        
            enddo               ! end loop over chains l

c     average q over configuration 
            do i = 1,nbox
               avgq(i) = avgq(i) + q(i)/num(i)
c               if (i.eq.1) then
c                  write(6,*) 'q  for config ', nc, q(i)/num(i)
c               endif
            enddo
                     
         enddo                  ! end loop over configurations (nc)
      enddo                     ! end loop over files (ilist)

      if (count.ne.tconfig) stop 'count.ne.tconfig'

      write(6,*) 'tconfig ', tconfig
      do i=1,nbox
         avgq(i) = avgq(i)/dble(tconfig)
         
         write(6,*)
         write(6,*) 'BOX ', i
         write(6,*) 'average q ', avgq(i)
         write(6,*)
         
      enddo
         
      end
