MODULE energy_garofalini
  use var_type,only:dp
  use const_phys,only:qqfact
  use util_math,only:erfunc
  use sim_system
  implicit none
  private
  save
  public::garofalini,vthreebody,triad,triad_en

! GAROFALINI.INC
  integer::pair_max
  parameter(pair_max=50000)
  real,public::galpha(6),grho(6),gbeta(6),ga(6,3),gb(6,3),gc(6,3),glambda(4),grij(4,2),ggamma(4,2),gtheta(4),grijsq(4,2)
  real::dxij(pair_max),dyij(pair_max),dzij(pair_max),dij(pair_max),dik(pair_max),dxik(pair_max),dyik(pair_max),dzik(pair_max)
  integer::itr1(pair_max),itr2(pair_max),itr3(pair_max),ntr

contains
!     **************************************************************
! calculates the energy using the garofalini (SiO2/H2O) ***
! exp-6 potential modified using articles in suijtab    ***
! parameters are defined in suijtab.f  KE ANDERSON      ***
!     **************************************************************
  function garofalini(rijsq,ntij,qa,qb,aa,bb)

      real::rijsq,rij,hterm,coul,qa,qb
      real::hterma,garofalini
      integer::ntij,aa,bb,i

#ifdef __DEBUG__
      write(io_output,*) 'entering garofalini in ',myid,'. Input:',rijsq,ntij,qa,qb
#endif

      rij = sqrt(rijsq)
      hterm = 0.0E0_dp
      coul = 0.0E0_dp
      garofalini = 0.0E0_dp


! write(io_output,*) aa,bb,ntij,qa,qb,gbeta(ntij),galpha(ntij),grho(ntij)
! H term
      do i=1,3
         hterma =   ga(ntij,i)/(1+exp(gb(ntij,i)*(rij-gc(ntij,i))))
! write(io_output,*) i,hterma,' (',ntij,')'
         hterm = hterm + hterma
      end do

      coul = qa*qb*erfunc(rij/gbeta(ntij))/rij
! write(io_output,*) 'erfunc',coul*rij
      coul = coul * qqfact

      garofalini = galpha(ntij)*exp(-rij/grho(ntij)) + hterm  + coul

! write(io_output,*) 'i,j,v2',aa,bb,'H term:',hterm,' Coul term:'
!     &     ,coul,' the rest:',garofalini-hterm-coul,
!     &     ' Total:',garofalini
! write(io_output,*) '                     ',rij

#ifdef __DEBUG__
      write(io_output,*) 'leaving garofalini in ',myid
#endif
      return
  end function garofalini

  subroutine vthreebody(vthree)
    real(kind=dp)::vthree
      integer::ntang,nta,ntb,tri,i,j,k
      real::vthreea,thetac,g,p
      logical::lwrite

#ifdef __DEBUG__
      write(io_output,*) 'start Vthreebody in ',myid
#endif

      vthreea = 0.0E0_dp
      vthree = 0.0E0_dp
      nta = 0
      ntb = 0
      lwrite=.false.

      do 10 tri = 1,ntr
         i = itr1(tri)
         j = itr2(tri)
         k = itr3(tri)

! skip if i (central atom) is hydrogen
         if(ntype(moltyp(i),1).eq.3) goto 10
! determine type
         if(ntype(moltyp(j),1).eq.ntype(moltyp(k),1)) then
            ntang = ntype(moltyp(j),1)
            if(lwrite) then
! write(30,*) 'ntang determination:',j,i,k,ntype(moltyp(j),1)
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
! write(30,*) 'ntang determination:',j,i,k,ntype(moltyp(j),1)
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

! dij = sqrt(dijsq(tri))
! dik = sqrt(diksq(tri))

         thetac = (dxij(tri)*dxik(tri)+dyij(tri)*dyik(tri)+ dzij(tri)*dzik(tri))/(dij(tri)*dik(tri))
         if ( thetac .ge. 1.0E0_dp ) thetac = 1.0E0_dp
         if ( thetac .le. -1.0E0_dp ) thetac = -1.0E0_dp

         p = (thetac-gtheta(ntang))**2
         g = exp((ggamma(ntang,nta)/(dij(tri)-grij(ntang,nta))) + (ggamma(ntang,ntb)/(dik(tri)-grij(ntang,ntb))))

         vthreea = glambda(ntang)*p*g
         if(lwrite) then
! write(io_output,*) 'details:',ntang,' out of vthreebody',vthreea
! write(io_output,*) 'thetac',thetac,' gtheta',gtheta(ntang)
            write(69,'(A11,I10,I10,I10,F15.6)') 'vthree jik',j,i,k ,vthreea
            write(69,'(A11,3F10.6,F25.6)') 'ij',dxij(tri),dyij(tri) ,dzij(tri),dij(tri)
            write(69,'(A11,3F10.6,F25.6)') 'ik',dxik(tri),dyik(tri) ,dzik(tri),dik(tri)
! write(io_output,*) 'dij',dij(tri),' dik',dik(tri)
! write(io_output,*) 'g',g,' theta',p
! write(io_output,*)
         end if
! write(io_output,*) j,i,k,'out of vthreebody',vthreea

         vthree = vthree + vthreea
 10   continue

! clear accumlators for the next triad call
      do i=1,1000
         itr1(i) = 0
         itr2(i) = 0
         itr3(i) = 0
      end do
      ntr = 0

#ifdef __DEBUG__
      write(io_output,*) 'end Vthreebody in ',myid
#endif

      return
  end subroutine vthreebody

!**************************************************************************
  subroutine triad

      integer::i,j,k,ptr,ptr2

#ifdef __DEBUG__
      write(io_output,*) 'starting triad in ',myid
#endif

      ntr = 0

      do  10 i=1,nchain
! imolty = moltyp(i)
! do 20 ii=1,nunit(imolty)
! write(io_output,*) 'from triad, neigh_cnt(i)',i,neigh_cnt(i)
            do 30 ptr = 1,neigh_cnt(i)-1

               j = neighbor(ptr,i)

               do 40 ptr2 = ptr+1,neigh_cnt(i)
                  k = neighbor(ptr2,i)
                  ntr = ntr + 1
! write(65,*) 'triad',ntr,':',j,i,k
! write(io_output,*) 'ndists:',ndij(ptr,i),ndij(ptr2,i)

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

#ifdef __DEBUG__
      write(io_output,*) 'end triad in ',myid
#endif

      return
  end subroutine triad

!******************************************************************
  subroutine triad_en(i,vthree,cnt,ni,nrij,nxi,nyi,nzi,lupdate)
    real(kind=dp)::vthree
      integer::i,j,k,m,imolty,jmolty,kmolty,mmolty,atomj,atomk,atomm
      integer::nta,ntb,ntang,cnt,ni(maxneigh),temp_cnt,itemp,number(nmax)
      integer::temp_nei(nmax),itype,jtype,ktype ,mtype
      logical::ltemp(nmax),lupdate,lwrite
      real::vthreea,thetac,p,g,nrij(maxneigh),nxi(maxneigh),nyi(maxneigh),nzi(maxneigh)
      real::temp_dist(nmax),temp_x(nmax),temp_y(nmax),temp_z(nmax)

      lwrite=.false.


! this is to determine correct 3-body interactions for single particles
! particle moved: i
! number of pairs within cutoff: cnt
! identity of pair molecule: ni(cnt)   ---- set in energy

      vthree = 0.0E0_dp
      temp_cnt = 0
      do itemp = 1,nmax
         ltemp(itemp) = .false.
         temp_dist(itemp) = 0.0E0_dp
         temp_x(itemp) = 0.0E0_dp
         temp_y(itemp) = 0.0E0_dp
         temp_z(itemp) = 0.0E0_dp
      end do

      imolty = moltyp(i)
      itype = ntype(imolty,1)
! if(i.eq.18) then
! write(io_output,*)
! write(io_output,*) 'neigh_icnt',cnt,':',(ni(j),j=1,cnt)
! end if
      do 10 j = 1,cnt
         atomj = ni(j)
         jmolty = moltyp(atomj)
         jtype = ntype(jmolty,1)

! skip if i=H; go straight to l-j-i loop
         if (itype.ne.3) then

! skip if Si-Si or O-O
            if(jtype.eq.itype) goto 10

! skip if Si-H
            if(jtype.eq.3.and.itype.eq.1) goto 10

! loop over other pairs with i as central atom
            do 20 k = j+1,cnt
               atomk = ni(k)
               kmolty = moltyp(atomk)
               ktype = ntype(kmolty,1)

               if(ktype.eq.itype) goto 20
               if(ktype.eq.3.and.itype.eq.1) goto 20

! determine type
               if(jtype.eq.ktype) then
                  ntang = jtype
! if(lwrite) write(io_output,*) 'ntang determination a:',atomj,
!     &                 i,atomk,ntang,' (',i,')'
                  nta = 1
                  ntb = 1
                  if((nrij(j).gt.grij(ntang,1)).or. (nrij(k).gt.grij(ntang,1))) then
! if(lwrite) write(io_output,*) 'skipped above'
                      goto 20
                   end if
               else
                  ntang = 4
! if(lwrite) write(io_output,*) 'ntang determination a:',atomj,
!     &                 i,atomk,ntang,' (',i,')'
                  if(jtype.eq.1) then
                     nta = 1
                     ntb = 2
                     if((nrij(j).gt.grij(4,nta)).or. (nrij(k).gt.grij(4,ntb))) then
! if(lwrite) write(io_output,*) 'skipped above'
! if(i.eq.18)
!     &                       write(io_output,*) 'catch distance Si-O-H',
!     &                       atomj,i,atomk,nrij(j),nrij(k)
                        goto 20
                     end if
                  else
                     nta = 2
                     ntb = 1
                     if((nrij(j).gt.grij(4,nta)).or. (nrij(k).gt.grij(4,ntb))) then
! if(i.eq.18)
!     &                       write(io_output,*) 'catch distance H-O-Si',
!     &                       atomj,i,atomk,nrij(j),nrij(k)
                        goto 20
                     end if
                  end if
               end if

               thetac = (nxi(j)*nxi(k)+nyi(j)*nyi(k)+ nzi(j)*nzi(k))/(nrij(j)*nrij(k))
               if ( thetac .ge. 1.0E0_dp ) thetac = 1.0E0_dp
               if ( thetac .le. -1.0E0_dp ) thetac = -1.0E0_dp

               p = (thetac-gtheta(ntang))**2
              g = exp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta)) )+(ggamma(ntang,ntb)/(nrij(k)-grij(ntang,ntb))))

               vthreea = glambda(ntang)*p*g

               if(lwrite) then
! write(io_output,*) 'details:',lupdate,' out of triad_en 1'
!$$$                  write(io_output,*) 'thetac',thetac
                  write(io_output,*) 'vthree jik',atomj,i,atomk,vthreea
!$$$                  write(io_output,*) nxi(j),nyi(j),nzi(j)
!$$$                  write(io_output,*) nxi(k),nyi(k),nzi(k)
!$$$                  write(io_output,*) nrij(j),nrij(k)
! write(io_output,*) g,thetac,p
!$$$                  write(io_output,*)
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
                     write(io_output,*) 'temps a',temp_cnt,temp_nei(temp_cnt)
           write(io_output,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(io_output,*)
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
                     write(io_output,*) 'temps b',temp_cnt,temp_nei(temp_cnt)
           write(io_output,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(io_output,*)
                  end if
                  ltemp(k) = .true.
               end if
 20         continue
         end if

! now starting loop to check i-j-m pairs
         if(jtype.ne.3.and..not.(itype.eq.3.and.jtype.eq.1)) then
            if(lwrite) then
               write(io_output,*)
            write(io_output,*) 'j',atomj,' neigh_cnt(j)',neigh_cnt(atomj),':', (neighbor(m,atomj),m=1,neigh_cnt(atomj))
            end if
            do 30 m = 1,neigh_cnt(atomj)
               atomm = neighbor(m,atomj)
               mmolty = moltyp(atomm)
               mtype = ntype(mmolty,1)

               if(mtype.eq.jtype) goto 30
               if(mtype.eq.3.and.jtype.eq.1) goto 30
               if(atomm.eq.i) goto 30

! determine type
               if(itype.eq.mtype) then
                  ntang = mtype
! if(lwrite) write(io_output,*) 'ntang determination b:',i,
!     &                 atomj,atomm,ntang,' (',i,')'
                  nta = 1
                  ntb = 1
                  if((nrij(j).gt.grij(ntang,1)) .or. (ndij(m,atomj).gt.grij(ntang,1))) goto 30
               else
                  ntang = 4
! if(lwrite) write(io_output,*) 'ntang determination b:',i,
!     &                 atomj,atomm,ntang,' (',i,')'
                  if(itype.eq.1) then
                     nta = 1
                     ntb = 2
                     if((nrij(j).gt.grij(4,nta)).or. (ndij(m,atomj).gt.grij(4,ntb))) then
! if(lwrite) write(io_output,*) 'skipped above'
!     & write(io_output,*) 'catch distance Si-O-H',nrij(j),ndij(m,atomj)
                        goto 30
                     end if
                  else
                     nta = 2
                     ntb = 1
                     if((nrij(j).gt.grij(4,nta)).or. (ndij(m,atomj).gt.grij(4,ntb))) then
! if(lwrite) write(io_output,*) 'skipped above'
!     & write(io_output,*) 'catch distance H-O-Si',nrij(j),ndij(m,atomj)
                        goto 30
                     end if
                  end if
               end if
! because the values saved are i-j not j-i must multiply all by -1
               thetac = ((-nxi(j))*nxij(m,atomj)+(-nyi(j))* nyij(m,atomj)+(-nzi(j))*nzij(m,atomj))/ (nrij(j)*ndij(m,atomj))
               if ( thetac .ge. 1.0E0_dp ) thetac = 1.0E0_dp
               if ( thetac .le. -1.0E0_dp ) thetac = -1.0E0_dp

               p = (thetac-gtheta(ntang))**2
              g = exp((ggamma(ntang,nta)/(nrij(j)-grij(ntang,nta)) )+(ggamma(ntang,ntb)/(ndij(m,atomj)- grij(ntang,ntb))))

               vthreea = glambda(ntang)*p*g

               if(lwrite) then
! write(io_output,*) 'details:',lupdate,' triad_en 2'
!$$$                  write(io_output,*) 'thetac',thetac,' ntang',ntang
                  write(io_output,'(A13,I5,I5,I5,F15.7)') 'vthree ijm',i,atomj ,atomm,vthreea
!$$$                  write(io_output,*) nxi(j),nyi(j),nzi(j)
!$$$                  write(io_output,*) nxij(m,atomj),nyij(m,atomj),nzij(m,atomj)
!$$$                  write(io_output,*) nrij(j),ndij(m,atomj)
! write(io_output,*) g,p
!$$$                  write(io_output,*)
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
                  write(io_output,*) 'temps ijm',temp_cnt,temp_nei(temp_cnt)
           write(io_output,*) temp_x(temp_cnt),temp_y(temp_cnt),temp_z(temp_cnt)
                  write(io_output,*)
                end if
               end if

 30         continue
         end if
 10   continue

! updating actual neighbors
! if(lupdate) then
! if(i.eq.18) then
! write(io_output,*) 'updating actual neighbors',temp_cnt
! do itemp = 1,temp_cnt
! write(io_output,*) 'temp_cnt number',itemp,number(itemp),
!     &              neighi(number(itemp))
! end do
! end if
! neigh_icnt=temp_cnt
! do itemp=1,neigh_icnt
! neighi(itemp) = temp_nei(itemp)
! ndiji(itemp) = temp_dist(itemp)
! nxiji(itemp) = temp_x(itemp)
! nyiji(itemp) = temp_y(itemp)
! nziji(itemp) = temp_z(itemp)
! end do
! do itemp = 1,neigh_icnt
! write(io_output,*) 'neigh_icnt number',itemp,neighi(itemp)
! write(io_output,*) 'dists',nxiji(itemp),nyiji(itemp),nziji(itemp)
! end do
! end if
      if(lwrite) write(io_output,*) 'triad_en vthree',vthree

      return
  end subroutine triad_en
end MODULE energy_garofalini
