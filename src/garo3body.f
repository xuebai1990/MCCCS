      subroutine vthreebody(vthree)


      implicit none

      include "control.inc"
      include "coord.inc"
      include 'garofalini.inc'

      integer ntang,nta,ntb,tri,i,j,k
      double precision vthree,vthreea,thetac,g,p
      logical lwrite

      vthreea = 0.0d0
      vthree = 0.0d0
      nta = 0
      ntb = 0
      lwrite=.false.

      do 10 tri = 1,ntr
         i = itr1(tri)
         j = itr2(tri)
         k = itr3(tri)
         
c     skip if i (central atom) is hydrogen
         if(ntype(moltyp(i),1).eq.3) goto 10
c     determine type
         if(ntype(moltyp(j),1).eq.ntype(moltyp(k),1)) then 
            ntang = ntype(moltyp(j),1)
            if(lwrite) then
c             write(30,*) 'ntang determination:',j,i,k,ntype(moltyp(j),1)
              write(30,*) tri,'a:',j,i,k,'(',ntang,')',dij(tri),dik(tri)
            endif
            nta = 1
            ntb = 1
            if ((dij(tri).gt.grij(ntang,1)).or.
     &           (dik(tri).gt.grij(ntang,1))) then
               if(lwrite)
     &         write(30,*) '   catch distance Si-O-Si',dij(tri),dik(tri)
               goto 10
            endif
         else
            ntang = 4
            if(lwrite) then
c             write(30,*) 'ntang determination:',j,i,k,ntype(moltyp(j),1)
              write(30,*) tri,'b:',j,i,k,'(',ntang,')',dij(tri),dik(tri)
            endif
            if(ntype(moltyp(j),1).eq.1) then
               nta = 1
               ntb = 2
               if((dij(tri).gt.grij(4,nta)).or.
     &              (dik(tri).gt.grij(4,ntb))) then
                  if(lwrite)
     &            write(30,*)'  catch distance Si-O-H',dij(tri),dik(tri)
                  goto 10
               endif
            else
               nta = 2
               ntb = 1
               if((dij(tri).gt.grij(4,nta)).or.
     &              (dik(tri).gt.grij(4,ntb))) then
                  if(lwrite)
     &           write(30,*) '  catch distance H-O-Si',dij(tri),dik(tri)
                  goto 10
               endif
            endif
         endif
         
c         dij = dsqrt(dijsq(tri))
c         dik = dsqrt(diksq(tri))
         
         thetac = (dxij(tri)*dxik(tri)+dyij(tri)*dyik(tri)+
     &        dzij(tri)*dzik(tri))/(dij(tri)*dik(tri))
         if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
         if ( thetac .le. -1.0d0 ) thetac = -1.0d0
         
         p = (thetac-gtheta(ntang))**2
         g = dexp((ggamma(ntang,nta)/(dij(tri)-grij(ntang,nta))) +
     &        (ggamma(ntang,ntb)/(dik(tri)-grij(ntang,ntb))))

         vthreea = glambda(ntang)*p*g
         if(lwrite) then
c            write(6,*) 'details:',ntang,' out of vthreebody',vthreea
c            write(6,*) 'thetac',thetac,' gtheta',gtheta(ntang)
               write(69,1201) 'vthree jik',j,i,k,vthreea
              write(69,1202) 'ij',dxij(tri),dyij(tri),dzij(tri),dij(tri)
              write(69,1202) 'ik',dxik(tri),dyik(tri),dzik(tri),dik(tri)
c               write(6,*) 'dij',dij(tri),' dik',dik(tri)
c            write(6,*) 'g',g,' theta',p
c            write(6,*)
         endif
c         write(6,*) j,i,k,'out of vthreebody',vthreea
         
         vthree = vthree + vthreea
 10   continue

c clear accumlators for the next triad call
      do i=1,1000
         itr1(i) = 0
         itr2(i) = 0
         itr3(i) = 0
      enddo
      ntr = 0

 1200 format(A9,I3,I3,I3,A7,F7.3,A7,F7.3,A7,F7.3,A7,I2,I2,I2,A1,L1,A1,L1
     &     ,A1,L1)
 1201 format(A11,I10,I10,I10,F15.6)
 1202 format(A11,3F10.6,F25.6)
c      write(6,*) 'end Vthreebody'

      return
      end


c**************************************************************************
      subroutine triad

      implicit none

      include "control.inc"
      include "coord.inc"
      include "garofalini.inc"
      include "neigh.inc"

      integer i,j,k,imolty,ptr,ptr2

c      write(6,*) 'starting triad'

      ntr = 0

      do  10 i=1,nchain
c         imolty = moltyp(i)
c         do 20 ii=1,nunit(imolty)
c         write(6,*) 'from triad, neigh_cnt(i)',i,neigh_cnt(i)
            do 30 ptr = 1,neigh_cnt(i)-1

               j = neighbor(ptr,i)

               do 40 ptr2 = ptr+1,neigh_cnt(i)
                  k = neighbor(ptr2,i)
                  ntr = ntr + 1
c                  write(65,*) 'triad',ntr,':',j,i,k
c                  write(6,*) 'ndists:',ndij(ptr,i),ndij(ptr2,i)

                  itr1(ntr) = i
                  itr2(ntr) = j
                  itr3(ntr) = k

                  dij(ntr) = ndij(ptr,i)
                  dik(ntr) = ndij(ptr2,i)

                  dxij(ntr) = nxij(ptr,i)
                  dyij(ntr) = nyij(ptr,i)
                  dzij(ntr) = nzij(ptr,i)
                  
                  dxik(ntr) = nxij(ptr2,i)
                  dyik(ntr) = nyij(ptr2,i)
                  dzik(ntr) = nzij(ptr2,i)
 40            continue
 30         continue
c 20      continue
 10   continue

c      write(6,*) 'end triad'

      return
      end

*******************************************************************
      subroutine triad_en(i,vthree,cnt,ni,nrij,nxi,nyi,nzi,lupdate)

      implicit none

      include "control.inc"
      include "coord.inc"
      include "garofalini.inc"
      include "neigh.inc"

      integer i,j,k,m,imolty,jmolty,kmolty,mmolty,atomj,atomk,atomm
      integer nta,ntb,ntang,cnt,ni(nmax),temp_cnt,temp,number(nmax)
      integer temp_nei(nmax),nnn,itype,jtype,ktype,mtype
      logical ltemp(nmax),lupdate,lwrite
      double precision vthree,vthreea,thetac,p,g,nrij(nmax),
     &     nxi(nmax),nyi(nmax),nzi(nmax)
      double precision temp_dist(nmax),temp_x(nmax),temp_y(nmax),
     &     temp_z(nmax)

      lwrite=.false.


c this is to determine correct 3-body interactions for single particles
c     particle moved: i
c     number of pairs within cutoff: cnt
c     identity of pair molecule: ni(cnt)   ---- set in energy

      vthree = 0.0d0
      temp_cnt = 0
      do temp = 1,nmax
         ltemp(temp) = .false.
         temp_dist(temp) = 0.0d0
         temp_x(temp) = 0.0d0
         temp_y(temp) = 0.0d0
         temp_z(temp) = 0.0d0
      enddo

      imolty = moltyp(i)
      itype = ntype(imolty,1)
c      if(i.eq.18) then
c         write(6,*) 
c         write(6,*) 'neigh_icnt',cnt,':',(ni(j),j=1,cnt)
c      endif
      do 10 j = 1,cnt
         atomj = ni(j)
         jmolty = moltyp(atomj)
         jtype = ntype(jmolty,1)

c     skip if i=H; go straight to l-j-i loop
         if (itype.ne.3) then

c     skip if Si-Si or O-O
            if(jtype.eq.itype) goto 10

c     skip if Si-H
            if(jtype.eq.3.and.itype.eq.1) goto 10

c     loop over other pairs with i as central atom
            do 20 k = j+1,cnt
               atomk = ni(k)
               kmolty = moltyp(atomk)
               ktype = ntype(kmolty,1)

               if(ktype.eq.itype) goto 20
               if(ktype.eq.3.and.itype.eq.1) goto 20
            
c     determine type
               if(jtype.eq.ktype) then 
                  ntang = jtype
c                  if(lwrite) write(6,*) 'ntang determination a:',atomj,
c     &                 i,atomk,ntang,' (',i,')'
                  nta = 1
                  ntb = 1
                  if((nrij(j).gt.grij(ntang,1)).or.
     &                 (nrij(k).gt.grij(ntang,1))) then
c                      if(lwrite) write(6,*) 'skipped above'
                      goto 20
                   endif
               else
                  ntang = 4
c                  if(lwrite) write(6,*) 'ntang determination a:',atomj,
c     &                 i,atomk,ntang,' (',i,')'
                  if(jtype.eq.1) then
                     nta = 1
                     ntb = 2
                     if((nrij(j).gt.grij(4,nta)).or.
     &                    (nrij(k).gt.grij(4,ntb))) then
c                        if(lwrite) write(6,*) 'skipped above'
c                        if(i.eq.18) 
c     &                       write(6,*) 'catch distance Si-O-H',
c     &                       atomj,i,atomk,nrij(j),nrij(k)
                        goto 20
                     endif
                  else
                     nta = 2
                     ntb = 1
                     if((nrij(j).gt.grij(4,nta)).or.
     &                    (nrij(k).gt.grij(4,ntb))) then
c                        if(i.eq.18)
c     &                       write(6,*) 'catch distance H-O-Si',
c     &                       atomj,i,atomk,nrij(j),nrij(k)
                        goto 20
                     endif
                  endif
               endif

               thetac = (nxi(j)*nxi(k)+nyi(j)*nyi(k)+
     &              nzi(j)*nzi(k))/(nrij(j)*nrij(k))
               if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
               if ( thetac .le. -1.0d0 ) thetac = -1.0d0
               
               p = (thetac-gtheta(ntang))**2
              g = dexp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta))
     &             )+(ggamma(ntang,ntb)/(nrij(k)-grij(ntang,ntb))))
               
               vthreea = glambda(ntang)*p*g

               if(lwrite) then
c                  write(6,*) 'details:',lupdate,' out of triad_en 1'
c$$$                  write(6,*) 'thetac',thetac
                  write(6,*) 'vthree jik',atomj,i,atomk,vthreea
c$$$                  write(6,*) nxi(j),nyi(j),nzi(j)
c$$$                  write(6,*) nxi(k),nyi(k),nzi(k)
c$$$                  write(6,*) nrij(j),nrij(k)
c                  write(6,*) g,thetac,p
c$$$                  write(6,*) 
               endif
         
               vthree = vthree + vthreea
c actual neighbor counter
               if(.not.ltemp(j).and.lupdate) then
                  temp_cnt = temp_cnt+1
                  number(temp_cnt) = j
                  temp_nei(temp_cnt) = neighi(j)
                  temp_dist(temp_cnt) = ndiji(j)
                  temp_x(temp_cnt) = nxiji(j)
                  temp_y(temp_cnt) = nyiji(j)
                  temp_z(temp_cnt) = nziji(j)
                  if(lwrite) then
                     write(6,*) 'temps a',temp_cnt,temp_nei(temp_cnt)
           write(6,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(6,*) 
                  endif
                  ltemp(j) = .true.
               endif
               if(.not.ltemp(k).and.lupdate) then
                  temp_cnt = temp_cnt+1
                  number(temp_cnt) = k
                  temp_nei(temp_cnt) = neighi(k)
                  temp_dist(temp_cnt) = ndiji(k)
                  temp_x(temp_cnt) = nxiji(k)
                  temp_y(temp_cnt) = nyiji(k)
                  temp_z(temp_cnt) = nziji(k)
                  if(lwrite) then
                     write(6,*) 'temps b',temp_cnt,temp_nei(temp_cnt)
           write(6,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(6,*) 
                  endif
                  ltemp(k) = .true.
               endif
 20         continue
         endif

c     now starting loop to check i-j-m pairs
         if(jtype.ne.3.and..not.(itype.eq.3.and.jtype.eq.1)) then
            if(lwrite) then
               write(6,*) 
            write(6,*) 'j',atomj,' neigh_cnt(j)',neigh_cnt(atomj),':',
     &           (neighbor(m,atomj),m=1,neigh_cnt(atomj))
            endif
            do 30 m = 1,neigh_cnt(atomj)
               atomm = neighbor(m,atomj)
               mmolty = moltyp(atomm)
               mtype = ntype(mmolty,1)

               if(mtype.eq.jtype) goto 30
               if(mtype.eq.3.and.jtype.eq.1) goto 30
               if(atomm.eq.i) goto 30

c     determine type
               if(itype.eq.mtype) then 
                  ntang = mtype
c                  if(lwrite) write(6,*) 'ntang determination b:',i,
c     &                 atomj,atomm,ntang,' (',i,')'
                  nta = 1
                  ntb = 1
                  if((nrij(j).gt.grij(ntang,1)) .or.
     &                 (ndij(m,atomj).gt.grij(ntang,1))) goto 30
               else
                  ntang = 4
c                  if(lwrite) write(6,*) 'ntang determination b:',i,
c     &                 atomj,atomm,ntang,' (',i,')'
                  if(itype.eq.1) then
                     nta = 1
                     ntb = 2
                     if((nrij(j).gt.grij(4,nta)).or.
     &                    (ndij(m,atomj).gt.grij(4,ntb))) then
c                        if(lwrite) write(6,*) 'skipped above'
c     & write(6,*) 'catch distance Si-O-H',nrij(j),ndij(m,atomj)
                        goto 30
                     endif
                  else
                     nta = 2
                     ntb = 1
                     if((nrij(j).gt.grij(4,nta)).or.
     &                    (ndij(m,atomj).gt.grij(4,ntb))) then
c                        if(lwrite) write(6,*) 'skipped above'
c     & write(6,*) 'catch distance H-O-Si',nrij(j),ndij(m,atomj)
                        goto 30
                     endif
                  endif
               endif
c     because the values saved are i-j not j-i must multiply all by -1 
               thetac = ((-nxi(j))*nxij(m,atomj)+(-nyi(j))*
     &              nyij(m,atomj)+(-nzi(j))*nzij(m,atomj))/
     &              (nrij(j)*ndij(m,atomj))
               if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
               if ( thetac .le. -1.0d0 ) thetac = -1.0d0
               
               p = (thetac-gtheta(ntang))**2
              g = dexp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta))
     &             )+(ggamma(ntang,ntb)/(ndij(m,atomj)-
     &              grij(ntang,ntb))))
               
               vthreea = glambda(ntang)*p*g

               if(lwrite) then
c                  write(6,*) 'details:',lupdate,' triad_en 2'
c$$$                  write(6,*) 'thetac',thetac,' ntang',ntang
                  write(6,1201) 'vthree ijm',i,atomj,atomm,vthreea
c$$$                  write(6,*) nxi(j),nyi(j),nzi(j)
c$$$                  write(6,*) nxij(m,atomj),nyij(m,atomj),nzij(m,atomj)
c$$$                  write(6,*) nrij(j),ndij(m,atomj)
c                  write(6,*) g,p
c$$$                  write(6,*)
               endif
              vthree = vthree + vthreea
c actual neighbor counter
               if(.not.ltemp(j).and.lupdate) then
                  temp_cnt = temp_cnt+1
                  number(temp_cnt) = j
                  temp_nei(temp_cnt) = neighi(j)
                  temp_dist(temp_cnt) = ndiji(j)
                  temp_x(temp_cnt) = nxiji(j)
                  temp_y(temp_cnt) = nyiji(j)
                  temp_z(temp_cnt) = nziji(j)
                  ltemp(j) = .true.
                  if(lwrite) then
                  write(6,*) 'temps ijm',temp_cnt,temp_nei(temp_cnt)
           write(6,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(6,*) 
                endif
               endif

 30         continue
         endif
 10   continue

c     updating actual neighbors
c      if(lupdate) then
c         if(i.eq.18) then
c            write(6,*) 'updating actual neighbors',temp_cnt
c            do temp = 1,temp_cnt
c               write(6,*) 'temp_cnt number',temp,number(temp),
c     &              neighi(number(temp))
c            enddo
c         endif
c         neigh_icnt=temp_cnt
c         do temp=1,neigh_icnt
c            neighi(temp) = temp_nei(temp)
c            ndiji(temp) = temp_dist(temp)
c            nxiji(temp) = temp_x(temp)
c            nyiji(temp) = temp_y(temp)
c            nziji(temp) = temp_z(temp)
c         enddo
c         do temp = 1,neigh_icnt
c            write(6,*) 'neigh_icnt number',temp,neighi(temp)
c            write(6,*) 'dists',nxiji(temp),nyiji(temp),nziji(temp)
c         enddo
c      endif
      if(lwrite) write(6,*) 'triad_en vthree',vthree

 1201 format(A13,I5,I5,I5,F15.7)
      return
      end
