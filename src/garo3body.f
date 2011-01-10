      subroutine vthreebody(vthree)


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

!$$$      include "control.inc"
!$$$      include "coord.inc"
!$$$      include 'garofalini.inc'

      integer(KIND=normal_int)::ntang,nta,ntb,tri,i,j,k
      real(KIND=double_precision)::vthree,vthreea,thetac,g,p
      logical::lwrite

      vthreea = 0.0d0
      vthree = 0.0d0
      nta = 0
      ntb = 0
      lwrite=.false.

      do 10 tri = 1,ntr
         i = itr1(tri)
         j = itr2(tri)
         k = itr3(tri)
         
!     skip if i (central atom) is hydrogen
         if(ntype(moltyp(i),1).eq.3) goto 10
!     determine type
         if(ntype(moltyp(j),1).eq.ntype(moltyp(k),1)) then 
            ntang = ntype(moltyp(j),1)
            if(lwrite) then
!             write(30,*) 'ntang determination:',j,i,k,ntype(moltyp(j),1)
              write(30,*) tri,'a:',j,i,k,'(',ntang,')',dij(tri),dik(tri)
            end if
            nta = 1
            ntb = 1
            if ((dij(tri).gt.grij(ntang,1)).or. (dik(tri).gt.grij(ntang,1))) then
               if(lwrite) write(30,*) '   catch distance Si-O-Si',dij(tri),dik(tri)
               goto 10
            end if
         else
            ntang = 4
            if(lwrite) then
!             write(30,*) 'ntang determination:',j,i,k,ntype(moltyp(j),1)
              write(30,*) tri,'b:',j,i,k,'(',ntang,')',dij(tri),dik(tri)
            end if
            if(ntype(moltyp(j),1).eq.1) then
               nta = 1
               ntb = 2
               if((dij(tri).gt.grij(4,nta)).or. (dik(tri).gt.grij(4,ntb))) then
                  if(lwrite) write(30,*)'  catch distance Si-O-H',dij(tri),dik(tri)
                  goto 10
               end if
            else
               nta = 2
               ntb = 1
               if((dij(tri).gt.grij(4,nta)).or. (dik(tri).gt.grij(4,ntb))) then
                  if(lwrite) write(30,*) '  catch distance H-O-Si',dij(tri),dik(tri)
                  goto 10
               end if
            end if
         end if
         
!         dij = dsqrt(dijsq(tri))
!         dik = dsqrt(diksq(tri))
         
         thetac = (dxij(tri)*dxik(tri)+dyij(tri)*dyik(tri)+ dzij(tri)*dzik(tri))/(dij(tri)*dik(tri))
         if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
         if ( thetac .le. -1.0d0 ) thetac = -1.0d0
         
         p = (thetac-gtheta(ntang))**2
         g = dexp((ggamma(ntang,nta)/(dij(tri)-grij(ntang,nta))) + (ggamma(ntang,ntb)/(dik(tri)-grij(ntang,ntb))))

         vthreea = glambda(ntang)*p*g
         if(lwrite) then
!            write(iou,*) 'details:',ntang,' out of vthreebody',vthreea
!            write(iou,*) 'thetac',thetac,' gtheta',gtheta(ntang)
            write(69,'(A11,I10,I10,I10,F15.6)') 'vthree jik',j,i,k ,vthreea
            write(69,'(A11,3F10.6,F25.6)') 'ij',dxij(tri),dyij(tri) ,dzij(tri),dij(tri)
            write(69,'(A11,3F10.6,F25.6)') 'ik',dxik(tri),dyik(tri) ,dzik(tri),dik(tri)
!               write(iou,*) 'dij',dij(tri),' dik',dik(tri)
!            write(iou,*) 'g',g,' theta',p
!            write(iou,*)
         end if
!         write(iou,*) j,i,k,'out of vthreebody',vthreea
         
         vthree = vthree + vthreea
 10   continue

! clear accumlators for the next triad call
      do i=1,1000
         itr1(i) = 0
         itr2(i) = 0
         itr3(i) = 0
      end do
      ntr = 0

!      write(iou,*) 'end Vthreebody'

      return
      end


!**************************************************************************
      subroutine triad

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

!$$$      include "control.inc"
!$$$      include "coord.inc"
!$$$      include "garofalini.inc"
!$$$      include "neigh.inc"

      integer(KIND=normal_int)::i,j,k,ptr,ptr2

!      write(iou,*) 'starting triad'

      ntr = 0

      do  10 i=1,nchain
!         imolty = moltyp(i)
!         do 20 ii=1,nunit(imolty)
!         write(iou,*) 'from triad, neigh_cnt(i)',i,neigh_cnt(i)
            do 30 ptr = 1,neigh_cnt(i)-1

               j = neighbor(ptr,i)

               do 40 ptr2 = ptr+1,neigh_cnt(i)
                  k = neighbor(ptr2,i)
                  ntr = ntr + 1
!                  write(65,*) 'triad',ntr,':',j,i,k
!                  write(iou,*) 'ndists:',ndij(ptr,i),ndij(ptr2,i)

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
! 20      continue
 10   continue

!      write(iou,*) 'end triad'

      return
      end

!******************************************************************
      subroutine triad_en(i,vthree,cnt,ni,nrij,nxi,nyi,nzi,lupdate)

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

!$$$      include "control.inc"
!$$$      include "coord.inc"
!$$$      include "garofalini.inc"
!$$$      include "neigh.inc"

      integer(KIND=normal_int)::i,j,k,m,imolty,jmolty,kmolty,mmolty ,atomj,atomk,atomm
      integer(KIND=normal_int)::nta,ntb,ntang,cnt,ni(nmax),temp_cnt ,itemp,number(nmax)
      integer(KIND=normal_int)::temp_nei(nmax),itype,jtype,ktype ,mtype
      logical::ltemp(nmax),lupdate,lwrite
      real(KIND=double_precision)::vthree,vthreea,thetac,p,g,nrij(nmax), nxi(nmax),nyi(nmax),nzi(nmax)
      real(KIND=double_precision)::temp_dist(nmax),temp_x(nmax) ,temp_y(nmax),temp_z(nmax)

      lwrite=.false.


! this is to determine correct 3-body interactions for single particles
!     particle moved: i
!     number of pairs within cutoff: cnt
!     identity of pair molecule: ni(cnt)   ---- set in energy

      vthree = 0.0d0
      temp_cnt = 0
      do itemp = 1,nmax
         ltemp(itemp) = .false.
         temp_dist(itemp) = 0.0d0
         temp_x(itemp) = 0.0d0
         temp_y(itemp) = 0.0d0
         temp_z(itemp) = 0.0d0
      end do

      imolty = moltyp(i)
      itype = ntype(imolty,1)
!      if(i.eq.18) then
!         write(iou,*) 
!         write(iou,*) 'neigh_icnt',cnt,':',(ni(j),j=1,cnt)
!      end if
      do 10 j = 1,cnt
         atomj = ni(j)
         jmolty = moltyp(atomj)
         jtype = ntype(jmolty,1)

!     skip if i=H; go straight to l-j-i loop
         if (itype.ne.3) then

!     skip if Si-Si or O-O
            if(jtype.eq.itype) goto 10

!     skip if Si-H
            if(jtype.eq.3.and.itype.eq.1) goto 10

!     loop over other pairs with i as central atom
            do 20 k = j+1,cnt
               atomk = ni(k)
               kmolty = moltyp(atomk)
               ktype = ntype(kmolty,1)

               if(ktype.eq.itype) goto 20
               if(ktype.eq.3.and.itype.eq.1) goto 20
            
!     determine type
               if(jtype.eq.ktype) then 
                  ntang = jtype
!                  if(lwrite) write(iou,*) 'ntang determination a:',atomj,
!     &                 i,atomk,ntang,' (',i,')'
                  nta = 1
                  ntb = 1
                  if((nrij(j).gt.grij(ntang,1)).or. (nrij(k).gt.grij(ntang,1))) then
!                      if(lwrite) write(iou,*) 'skipped above'
                      goto 20
                   end if
               else
                  ntang = 4
!                  if(lwrite) write(iou,*) 'ntang determination a:',atomj,
!     &                 i,atomk,ntang,' (',i,')'
                  if(jtype.eq.1) then
                     nta = 1
                     ntb = 2
                     if((nrij(j).gt.grij(4,nta)).or. (nrij(k).gt.grij(4,ntb))) then
!                        if(lwrite) write(iou,*) 'skipped above'
!                        if(i.eq.18) 
!     &                       write(iou,*) 'catch distance Si-O-H',
!     &                       atomj,i,atomk,nrij(j),nrij(k)
                        goto 20
                     end if
                  else
                     nta = 2
                     ntb = 1
                     if((nrij(j).gt.grij(4,nta)).or. (nrij(k).gt.grij(4,ntb))) then
!                        if(i.eq.18)
!     &                       write(iou,*) 'catch distance H-O-Si',
!     &                       atomj,i,atomk,nrij(j),nrij(k)
                        goto 20
                     end if
                  end if
               end if

               thetac = (nxi(j)*nxi(k)+nyi(j)*nyi(k)+ nzi(j)*nzi(k))/(nrij(j)*nrij(k))
               if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
               if ( thetac .le. -1.0d0 ) thetac = -1.0d0
               
               p = (thetac-gtheta(ntang))**2
              g = dexp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta)) )+(ggamma(ntang,ntb)/(nrij(k)-grij(ntang,ntb))))
               
               vthreea = glambda(ntang)*p*g

               if(lwrite) then
!                  write(iou,*) 'details:',lupdate,' out of triad_en 1'
!$$$                  write(iou,*) 'thetac',thetac
                  write(iou,*) 'vthree jik',atomj,i,atomk,vthreea
!$$$                  write(iou,*) nxi(j),nyi(j),nzi(j)
!$$$                  write(iou,*) nxi(k),nyi(k),nzi(k)
!$$$                  write(iou,*) nrij(j),nrij(k)
!                  write(iou,*) g,thetac,p
!$$$                  write(iou,*) 
               end if
         
               vthree = vthree + vthreea
! actual neighbor counter
               if(.not.ltemp(j).and.lupdate) then
                  temp_cnt = temp_cnt+1
                  number(temp_cnt) = j
                  temp_nei(temp_cnt) = neighi(j)
                  temp_dist(temp_cnt) = ndiji(j)
                  temp_x(temp_cnt) = nxiji(j)
                  temp_y(temp_cnt) = nyiji(j)
                  temp_z(temp_cnt) = nziji(j)
                  if(lwrite) then
                     write(iou,*) 'temps a',temp_cnt,temp_nei(temp_cnt)
           write(iou,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(iou,*) 
                  end if
                  ltemp(j) = .true.
               end if
               if(.not.ltemp(k).and.lupdate) then
                  temp_cnt = temp_cnt+1
                  number(temp_cnt) = k
                  temp_nei(temp_cnt) = neighi(k)
                  temp_dist(temp_cnt) = ndiji(k)
                  temp_x(temp_cnt) = nxiji(k)
                  temp_y(temp_cnt) = nyiji(k)
                  temp_z(temp_cnt) = nziji(k)
                  if(lwrite) then
                     write(iou,*) 'temps b',temp_cnt,temp_nei(temp_cnt)
           write(iou,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(iou,*) 
                  end if
                  ltemp(k) = .true.
               end if
 20         continue
         end if

!     now starting loop to check i-j-m pairs
         if(jtype.ne.3.and..not.(itype.eq.3.and.jtype.eq.1)) then
            if(lwrite) then
               write(iou,*) 
            write(iou,*) 'j',atomj,' neigh_cnt(j)',neigh_cnt(atomj),':', (neighbor(m,atomj),m=1,neigh_cnt(atomj))
            end if
            do 30 m = 1,neigh_cnt(atomj)
               atomm = neighbor(m,atomj)
               mmolty = moltyp(atomm)
               mtype = ntype(mmolty,1)

               if(mtype.eq.jtype) goto 30
               if(mtype.eq.3.and.jtype.eq.1) goto 30
               if(atomm.eq.i) goto 30

!     determine type
               if(itype.eq.mtype) then 
                  ntang = mtype
!                  if(lwrite) write(iou,*) 'ntang determination b:',i,
!     &                 atomj,atomm,ntang,' (',i,')'
                  nta = 1
                  ntb = 1
                  if((nrij(j).gt.grij(ntang,1)) .or. (ndij(m,atomj).gt.grij(ntang,1))) goto 30
               else
                  ntang = 4
!                  if(lwrite) write(iou,*) 'ntang determination b:',i,
!     &                 atomj,atomm,ntang,' (',i,')'
                  if(itype.eq.1) then
                     nta = 1
                     ntb = 2
                     if((nrij(j).gt.grij(4,nta)).or. (ndij(m,atomj).gt.grij(4,ntb))) then
!                        if(lwrite) write(iou,*) 'skipped above'
!     & write(iou,*) 'catch distance Si-O-H',nrij(j),ndij(m,atomj)
                        goto 30
                     end if
                  else
                     nta = 2
                     ntb = 1
                     if((nrij(j).gt.grij(4,nta)).or. (ndij(m,atomj).gt.grij(4,ntb))) then
!                        if(lwrite) write(iou,*) 'skipped above'
!     & write(iou,*) 'catch distance H-O-Si',nrij(j),ndij(m,atomj)
                        goto 30
                     end if
                  end if
               end if
!     because the values saved are i-j not j-i must multiply all by -1 
               thetac = ((-nxi(j))*nxij(m,atomj)+(-nyi(j))* nyij(m,atomj)+(-nzi(j))*nzij(m,atomj))/ (nrij(j)*ndij(m,atomj))
               if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
               if ( thetac .le. -1.0d0 ) thetac = -1.0d0
               
               p = (thetac-gtheta(ntang))**2
              g = dexp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta)) )+(ggamma(ntang,ntb)/(ndij(m,atomj)- grij(ntang,ntb))))
               
               vthreea = glambda(ntang)*p*g

               if(lwrite) then
!                  write(iou,*) 'details:',lupdate,' triad_en 2'
!$$$                  write(iou,*) 'thetac',thetac,' ntang',ntang
                  write(iou,'(A13,I5,I5,I5,F15.7)') 'vthree ijm',i,atomj ,atomm,vthreea
!$$$                  write(iou,*) nxi(j),nyi(j),nzi(j)
!$$$                  write(iou,*) nxij(m,atomj),nyij(m,atomj),nzij(m,atomj)
!$$$                  write(iou,*) nrij(j),ndij(m,atomj)
!                  write(iou,*) g,p
!$$$                  write(iou,*)
               end if
              vthree = vthree + vthreea
! actual neighbor counter
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
                  write(iou,*) 'temps ijm',temp_cnt,temp_nei(temp_cnt)
           write(iou,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(iou,*) 
                end if
               end if

 30         continue
         end if
 10   continue

!     updating actual neighbors
!      if(lupdate) then
!         if(i.eq.18) then
!            write(iou,*) 'updating actual neighbors',temp_cnt
!            do itemp = 1,temp_cnt
!               write(iou,*) 'temp_cnt number',itemp,number(itemp),
!     &              neighi(number(itemp))
!            end do
!         end if
!         neigh_icnt=temp_cnt
!         do itemp=1,neigh_icnt
!            neighi(itemp) = temp_nei(itemp)
!            ndiji(itemp) = temp_dist(itemp)
!            nxiji(itemp) = temp_x(itemp)
!            nyiji(itemp) = temp_y(itemp)
!            nziji(itemp) = temp_z(itemp)
!         end do
!         do itemp = 1,neigh_icnt
!            write(iou,*) 'neigh_icnt number',itemp,neighi(itemp)
!            write(iou,*) 'dists',nxiji(itemp),nyiji(itemp),nziji(itemp)
!         end do
!      end if
      if(lwrite) write(iou,*) 'triad_en vthree',vthree

      return
      end
