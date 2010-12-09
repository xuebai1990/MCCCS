      subroutine sumup( ovrlap, v, vinter,vtail,vintra,vvib,
     &                  vbend,vtg,vext,velect,vflucq,ibox,lvol)
                       

! sumup     
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
! John Stubbs, and Collin Wick and Ilja Siepmann  
!                     
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to 
!
! Free Software Foundation, Inc. 
! 59 Temple Place - Suite 330
! Boston, MA  02111-1307, USA.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!    *******************************************************************
!    ** calculates the total potential energy for a configuration.    **
!    **                                                               **
!    ** logical::ovrlap            true for substantial atom overlap   **
!    *******************************************************************
 
      implicit none

! *** common blocks ***
      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'poten.inc'
      include 'conver.inc' 
      include 'external.inc'
      include 'zeolite.inc'
      include 'connect.inc'
      include 'ewaldsum.inc'
      include 'fepsi.inc'
      include 'zeopoten.inc'
      include 'inputdata.inc'
      include 'qqlist.inc'
      include 'clusterbias.inc'
      include 'neigh.inc'
      include 'cell.inc'
      include 'nsix.inc'
      include 'peboco.inc'     
      include 'ipswpar.inc'
      include 'eepar.inc'
! kea include for garofalini potential
      include 'garofalini.inc'
      include 'tabulated.inc'
! RP added for MPI
      include 'mpif.h'
      include 'mpi.inc'
 
      logical::ovrlap, lvol 
      logical::lexplt,lqimol,lqjmol,lcoulo,lij2,liji,lqchgi,lexclude
      integer::i, imolty, ii, j, jmolty, jj, ntii, ntjj, ntij, iunit
     &     , ip1, ip2, ip3,ibox,nmcount,iii,jjj,ip
     &     ,iivib, jjvib, jjben, jjtor, it, ntj,itype,k, mmm
      real(8)::v, vinter, vintra, vtail, vvib, vbend, vtg, vext
     &     ,velect,vflucq,qqii
      real(8)::rcutsq,rminsq,rxui,ryui,rzui,rxuij,ryuij,rzuij
     &       ,rijsq,sr2, sr6, rho, thetac, theta 
     &     ,xaa1, yaa1, zaa1, xa1a2, ya1a2, za1a2, daa1, da1a2, dot
     &     ,vtorso, dzui, dz3, dz12, rhoz,xcc,ycc,zcc,tcc,spltor
     &     ,mmff,rij,vrecipsum,erfunc
     &     ,rbcut,ninesix,vwell, genlj
! tabulated potential variables
      real(8)::tabulated_vib, tabulated_bend, tabulated_vdW, 
     &     tabulated_elect, rbend, rbendsq

      real(8)::sx,sy,sz

!      real(8)::vtemp
      real(8)::vsc

      real(8)::coru,coruz,xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,qave
      real(8)::ljsami,ljpsur,ljmuir,exsami,exmuir,exzeo,exsix
      real(8)::exgrph,vintera,velecta,vol
      real(8)::xvec(numax,numax),yvec(numax,numax)
     &                ,zvec(numax,numax),distij(numax,numax),epsilon2
     &                ,sigma2, distij2(numax,numax)
      real(8)::slitpore, v_elect_field, field
      dimension lcoulo(numax,numax),lexclude(nmax)
! Neeraj & RP for MPI
      real(8)::sum_velect,sum_vinter,sum_vtail,sum_vintra
     &       ,sum_vflucq,sum_vvib,sum_vbend,sum_vtg,sum_vext,sum_vwell
     &      ,sum_sself,sum_correct,my_velect,sum_my_velect
      logical::all_ovrlap
! --------------------------------------------------------------------
      vintera = 0.0d0
      velecta = 0.0d0

!      write(iou,*) 'start SUMUP, box ', ibox
      ovrlap = .false.
! KM for MPI
      all_ovrlap = .false.
      sum_velect = 0.0d0
      sum_vinter = 0.0d0
      sum_vtail = 0.0d0
      sum_vintra = 0.0d0
      sum_vflucq = 0.0d0
      sum_vvib = 0.0d0
      sum_vbend = 0.0d0
      sum_vtg = 0.0d0
      sum_vext = 0.0d0
      sum_vwell = 0.0d0
      sum_sself = 0.0d0
      sum_correct = 0.0d0
      my_velect = 0.0d0
      sum_my_velect = 0.0d0

      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut = rcut(ibox)
      field = Elect_field(ibox)

      rminsq = rmin * rmin

      v = 0.0d0
      vinter = 0.0d0
      vintra = 0.0d0
      vtail = 0.0d0
      vtg = 0.0d0
      vbend = 0.0d0
      vvib = 0.0d0
      vext = 0.0d0
      velect = 0.0d0
      vflucq = 0.0d0
!kea - 3body garofalini term
      v3garo = 0.0d0
      vwell = 0.0d0
! *** check the molecule count ***
      nmcount = 0
      do i = 1, nchain
         if ( nboxi(i) .eq. ibox ) then
            nmcount=nmcount+1
         end if
         neigh_cnt(i) = 0
      end do
      if ( nmcount .ne. nchbox(ibox) ) then
         write(iou,*) 'SUMUP: nmcount ne nchbox', nmcount, nchbox
         call cleanup('')
      end if
 
! ###############################################################

! *******************************
! *** INTERCHAIN INTERACTIONS ***
! *******************************

! --- loop over all chains i 
! --- not if lgrand and ibox =2
! --- JLR 11-24-09 don't loop if box is ideal
!      if (.not.(lgrand.and.(ibox.eq.2))) then
      if (.not.(lgrand.and.ibox.eq.2) .and.
     &     .not.lideal(ibox) ) then
! --- END JLR 11-24-09
! RP added for MPI
!       do 100 i = 1, nchain - 1
         do 100 i = myid+1,nchain-1,numprocs
 
! ### check if i is in relevant box ###
            if ( nboxi(i) .eq. ibox ) then
               imolty = moltyp(i)
               lqimol = lelect(imolty)
               
               if ( lcutcm .and. lvol ) then
                  xcmi = xcm(i)
                  ycmi = ycm(i)
                  zcmi = zcm(i)
                  rcmi = rcmu(i)
               else
                  lij2 = .true.
               end if
               
               if ( nugrow(imolty) .eq. nunit(imolty) ) then
                  lexplt = .false.
               else
                  lexplt = .true.
               end if
               
! --- loop over all chains j with j>i 
               do 99 j = i + 1, nchain
! ### check for simulation box ###
                  if ( nboxi(j) .eq. ibox ) then


                     jmolty = moltyp(j)
                     lqjmol = lelect(jmolty)

                     if (lcutcm .and. lvol ) then
!                     --- check if ctrmas within rcmsq
                        rxuij = xcmi - xcm(j)
                        ryuij = ycmi - ycm(j)
                        rzuij = zcmi - zcm(j)
!                     --- minimum image the ctrmas pair separations
                        if ( lpbc ) call mimage (rxuij,ryuij,rzuij,ibox)

                        rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                        rij  = dsqrt(rijsq)
                        rcm = rbcut + rcmi + rcmu(j)
                        rcmsq = rcm*rcm

!                      if ( lneighbor .and. rcmsq .lt. rbsmax**2 
!     &                     .and. rcmsq .gt. rbsmin**2 ) then
!                         neigh_cnt(i,jmolty)=neigh_cnt(i,jmolty)+1
!                         neighbor(neigh_cnt(i,jmolty),i,jmolty)=j
!                         neigh_cnt(j,imolty)=neigh_cnt(j,imolty)+1
!                         neighbor(neigh_cnt(j,imolty),j,imolty)=i
!                      end if

                        if ( rijsq .gt. rcmsq ) then
                           if ( lqimol .and. lqjmol .and. lchgall ) then
                              lij2 = .false.
                              goto 98
                           else
                              goto 99
                           end if
                        else
                           lij2 = .true.
                        end if
                     end if

 98                  do ii = 1,nunit(imolty)
                        ntii = ntype(imolty,ii)
                        liji = lij(ntii)
                        lqchgi = lqchg(ntii)
                        rxui = rxu(i,ii)
                        ryui = ryu(i,ii)
                        rzui = rzu(i,ii)
                        
! --- loop over all beads jj of chain j 
                        do 97 jj = 1, nunit(jmolty) 
! --- check exclusion table
                           if ( lexclu(imolty,ii,jmolty,jj) ) goto 97
                           
                           ntjj = ntype(jmolty,jj)
                           if ( lij2 ) then
                              if ( (.not. (liji .and. lij(ntjj))) 
     &                             .and. (.not. (lqchgi .and. 
     &                             lqchg(ntjj)))) goto 97
                           else
                              if (.not. (lqchgi .and. lqchg(ntjj)))
     &                             goto 97
                           end if
                           if (lexpsix .or. lmmff) then
                              ntij = (ntii+ntjj)/2
                           elseif (lninesix) then
                              ntij = (ntii-1)*nxatom + ntjj
! Generalized Lennard Jones
                           elseif (lgenlj) then
                              ntij = (ntii-1)*nntype + ntjj
! KEA garofalini
                           elseif (lgaro) then
                              if (ntii.eq.ntjj) then
                                 ntij = ntii
                              else
                                 ntij = ntii+ntjj+1
                              end if
                           else
                              ntij = (ntii-1)*nntype + ntjj
                           end if
                           if (lexpee) rminsq =rminee(ntij)*rminee(ntij)
                         
                           rxuij = rxui - rxu(j,jj)
                           ryuij = ryui - ryu(j,jj)
                           rzuij = rzui - rzu(j,jj)
! *** minimum image the pair separations ***
                           if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)

                           rijsq = (rxuij*rxuij)+(ryuij*ryuij)
     &                          + (rzuij*rzuij)
                           rij   = dsqrt(rijsq)
                         
!                           if ( i .eq. 12 .and. ii .eq. 6 .and.
!     &                        j .eq. 95 .and. jj .eq. 1 ) then
!                           write(iou,*) 'CONTROL CONTROL CONTROL'
!                           write(iou,*) 'box',ibox,nboxi(i),nboxi(j)
!                           write(iou,*) 'i xyz',rxui,ryui,rzui
!                           write(iou,*) 'j xyz',rxu(j,jj),ryu(j,jj),
!     &                          rzu(j,jj)
!                           write(iou,*) 'r*uij',rxuij,ryuij,rzuij
!                           write(iou,*) 'dist2',rijsq
!                           write(iou,*) 'distance', dsqrt(rijsq)
!                        end if

                           if ( rijsq .lt. rminsq .and. .not.
     &                          (lexpand(imolty) .or. 
     &                          lexpand(jmolty))) then
                              if ( .not. lvol .and.myid.eq.0) then
                                 write(iou,*) 'overlap inter'
                                 write(iou,*)'rijsq rminsq',rijsq,rminsq
                                 write(iou,*) 'i ii', i, ii
                                 write(iou,*) 'i-pos', rxui,ryui,rzui
                                 write(iou,*) 'j jj', j, jj
                                 write(iou,*) 'j-pos', 
     &                                rxu(j,jj),ryu(j,jj),rzu(j,jj)
                              end if
                              ovrlap = .true.

!-------- RP added for MPI to compensate ovrlap
                              if(ovrlap .eq. .true.)then
                                 goto 199
                              end if
! -------------------------------
!                            return
                           elseif ( rijsq .lt. rcutsq .or. lijall) then

                              if (L_vdW_table.and.(.not.(
     &                             lexpand(imolty).or.lexpand(jmolty)))
     &                             )then
                                 call lininter_vdW(rij, 
     &                                tabulated_vdW, ntii, ntjj)
                                 
                                 vinter = vinter + tabulated_vdW
                                 
                              elseif (llj.and.(.not.(lexpand(imolty).or.
     &                                lexpand(jmolty)))) then
                                 if ( lij(ntii) .and. lij(ntjj) ) then
                                    sr2 = sig2ij(ntij) / rijsq
                                    epsilon2=epsij(ntij)        
                                    sr6 = sr2 * sr2 * sr2
                                    vinter = vinter 
     &                                   + sr6*(sr6-1.0d0)*epsilon2
                                 end if
                              elseif ( lsami ) then
                                 vinter = vinter + ljsami(rijsq,ntij)
                              elseif (lexpsix) then
                                 vinter = vinter + exsix(rijsq,ntij)
                              elseif (lmmff) then
                                 vinter = vinter + mmff(rijsq,ntij)
                              elseif (lninesix) then
                                 vinter = vinter + ninesix(rijsq,ntij)
!     Generalized Lennard Jones potential
                              elseif (lgenlj) then
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2=epsij(ntij)
                                 vinter = vinter + 
     &                                genlj(rijsq,sr2,epsilon2)
                              elseif ( lmuir ) then
                                 vinter = vinter + ljmuir(rijsq,ntij)
                              elseif ( lpsurf ) then
                                 vinter = vinter + ljpsur(rijsq,ntij)
! KEA garofalini potential
                              elseif ( lgaro) then
                                 vinter = vinter + garofalini(rijsq,ntij
     &                                ,qqu(i,ii),qqu(j,jj),i,j)
                                 if(lshift) then
                                    vinter = vinter-ecut(ntij)
                                 end if
                              else if (lshift) then
                                 sr2 = sig2ij(ntij) / rijsq
                                 sr6 = sr2 * sr2 * sr2
                                 vinter = vinter + 
     &                                sr6*(sr6-1.0d0)*epsij(ntij)-
     &                                ecut(ntij) 
                                 
                              elseif ( lfepsi ) then
                                 if ( lij(ntii) .and. lij(ntjj) ) then
                                    sr6 = rijsq*rijsq*rijsq
                                    
                                    if ( (.not. lqchg(ntii)) .and. 
     &                                   (.not. lqchg(ntjj)) ) then
                                       if ( nunit(imolty) .eq. 4 ) then
! *** TIP-4P structure (temperary use ???)
                                          qave=(qqu(i,4)+qqu(j,4))/2.0d0
                                       else
                                          qave=(qqu(i,4)+qqu(i,5)+
     &                                         qqu(j,4)+qqu(j,5))*0.85d0
                                       end if
                                    else
                                       qave =(qqu(i,ii)+qqu(j,jj))/2.0d0
                                    end if
                                    if ( lexpand(imolty) 
     &                                   .and. lexpand(jmolty)) then
                                       epsilon2=dsqrt(epsilon(imolty,ii)
     &                                      *epsilon(jmolty,jj))
                                    elseif (lexpand(imolty)) then
                                       epsilon2=dsqrt(epsilon(imolty,ii)
     &                                      *epsi(ntjj))
                                    elseif ( lexpand(jmolty) ) then
                                       
                                       epsilon2=dsqrt(epsi(ntii)*
     &                                      epsilon(jmolty,jj))
                                    else
                                       epsilon2=epsij(ntij)
                                    end if
                                    vinter = vinter + 
     &                                   ((aslope*(qave-a0)*(qave-a0)
     &                                   +ashift)/sr6 - (bslope*(qave-
     &                                   b0)*(qave-b0)+bshift))/
     &                                   sr6*epsilon2  
                                    
                                 end if      
                              else
                                 if ( lij(ntii) .and. lij(ntjj) ) then
                                    if ( lexpand(imolty) 
     &                                   .and. lexpand(jmolty)) then
                                       sigma2=(sigma(imolty,ii)+
     &                                      sigma(jmolty,jj))/2.0d0
                                       sr2 = sigma2*sigma2/rijsq
                                       epsilon2=dsqrt(epsilon(imolty,ii)
     &                                      *epsilon(jmolty,jj))
                                    elseif ( lexpand(imolty) ) then
                                       sigma2=(sigma(imolty,ii)+
     &                                      sigi(ntjj))/2.0d0
                                       sr2 = sigma2*sigma2/rijsq
                                       epsilon2=dsqrt(epsilon(imolty,ii)
     &                                      *epsi(ntjj))
                                    elseif ( lexpand(jmolty) ) then
                                       sigma2=(sigma(jmolty,jj)+
     &                                      sigi(ntii))/2.0d0
                                       sr2 = sigma2*sigma2/rijsq
                                       epsilon2=dsqrt(epsi(ntii)*
     &                                      epsilon(jmolty,jj))
                                    else
                                       sr2 = sig2ij(ntij) / rijsq
                                       epsilon2=epsij(ntij)
                                    end if
                                    sr6 = sr2 * sr2 * sr2
                                    vinter = vinter 
     &                                   + sr6*(sr6-1.0d0)*epsilon2
                                 end if
                              end if
                           end if
                           
! charge interactions
!kea - skip for garofalini; included in vinter
                           if(lgaro) then
                           elseif ( lchgall.and. lqchg(ntii) 
     &                             .and. lqchg(ntjj)) then
                              if ( lewald ) then
                                 velect = velect + qqu(i,ii)*qqu(j,jj)*
     &                                erfunc(calp(ibox)*rij)/
     &                              rij
                              else
                                 velect = velect + qqu(i,ii)*qqu(j,jj)/
     &                                rij
                              end if
                           elseif ( lqimol .and. lqjmol .and. 
     &                             lqchgi .and. lqchg(ntjj) ) then
                              
                              if (lewald) then               
                                 if (rijsq.lt.rcutsq) then
                                    velect = velect + qqu(i,ii)*
     &                                   qqu(j,jj)*erfunc(calp(ibox)*
     &                                   rij)/rij
                                 end if
                              else
! --- All-Atom charges (charge-group look-up table)
                                 iii = leaderq(imolty,ii)
                                 jjj = leaderq(jmolty,jj)
                                 if ( iii .eq. ii .and. jjj .eq.jj )then
!     --- set up the table
                                    if ( rijsq .lt. rcutsq ) then
                                       lcoulo(iii,jjj) = .true.
                                    else
                                       lcoulo(iii,jjj) = .false.
                                    end if
                                 end if
                                 if ( lcoulo(iii,jjj) ) then
                                    if (L_elect_table) then 
                                       call lininter_elect(rij, 
     &                                      tabulated_elect, ntii, ntjj)
                                       velect = velect + qqu(i,ii)*
     &                                      qqu(j,jj)*tabulated_elect
                                    else
                                       velect = velect + qqu(i,ii)
     &                                      *qqu(j,jj)/rij
                                    end if
                                 end if
                              end if
                           end if
                           
! End inter-charge loop

!cc  KM for MPI
!cc  all processors need to know neighbor information
!cc  lneighbor and lgaro will not work in parallel
!cc  calculation of neighbors assumes everything is sequential
                           if ( lneighbor .and. ii .eq. 1 .and. 
     &                          jj .eq. 1 .and. rijsq .lt. rbsmax**2 
     &                          .and. rijsq .gt. rbsmin**2 ) then
!                            neigh_cnt(i,jmolty)=neigh_cnt(i,jmolty)+1
!                            neighbor(neigh_cnt(i,jmolty),i,jmolty)=j
!                            neigh_cnt(j,imolty)=neigh_cnt(j,imolty)+1
!                            neighbor(neigh_cnt(j,imolty),j,imolty)=i
                            
                              neigh_cnt(i)=neigh_cnt(i)+1
                              neighbor(neigh_cnt(i),i)=j
                              neigh_cnt(j)=neigh_cnt(j)+1
                              neighbor(neigh_cnt(j),j)=i
                           elseif(lgaro) then
                              if((ntij.eq.4.and.rijsq.lt.grijsq(2,1))
     &                             .or.(ntij.eq.6.and.rijsq.lt.
     &                             grijsq(3,1)))
     &                             then
                                 write(64,*) 'neighbor',i,' (',
     &                                neigh_cnt(i)+1,')',j,' (',
     &                                neigh_cnt(j)+1,')'
                                 neigh_cnt(i)=neigh_cnt(i)+1
                                 neighbor(neigh_cnt(i),i)=j
                                 neigh_cnt(j)=neigh_cnt(j)+1
                                 neighbor(neigh_cnt(j),j)=i
                                 ndij(neigh_cnt(i),i) = rij
                                 ndij(neigh_cnt(j),j) = 
     &                                ndij(neigh_cnt(i),i)
                                 nxij(neigh_cnt(i),i) = rxuij
                                 nyij(neigh_cnt(i),i) = ryuij
                                 nzij(neigh_cnt(i),i) = rzuij
                                 nxij(neigh_cnt(j),j) = -rxuij
                                 nyij(neigh_cnt(j),j) = -ryuij
                                 nzij(neigh_cnt(j),j) = -rzuij
                              end if
                           end if
 97                     continue
                     end do
                  end if
 99            continue
            end if
            
 100     continue
! ----- Returning from ovrlap--------------

 199     continue

! KM don't check overlap until allreduce is finished
!      if(ovrlap .eq. .true.)then
!        write(iou,*)'521: in sumup ovrlap=',ovrlap,'myid=',myid
!      end if
! -----------------------------------------

         CALL MPI_ALLREDUCE(ovrlap,all_ovrlap,1,MPI_LOGICAL,MPI_LOR,
     &        MPI_COMM_WORLD,ierr)

         ovrlap = all_ovrlap
         if(ovrlap)then
!            write(iou,*)'530: in sumup ovrlap=',ovrlap,'myid=',myid
            return
         end if

         CALL MPI_ALLREDUCE(vinter, sum_vinter,1,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)

         CALL MPI_ALLREDUCE(velect, sum_velect,1,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)

         velect = sum_velect
         vinter = sum_vinter

         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &        .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro
     &        .and..not.L_vdW_table) then 
            vinter = 4.0d0 * vinter
         end if
         
! KEA garofalini 3 body potential
         if (lgaro) then
            call triad
            call vthreebody(v3garo)
         end if
         
         
         if (ltailc) then
!--     add tail corrections for the Lennard-Jones energy
            if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
               vol = cell_vol(ibox)
            else
               vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
            end if
            do jmolty = 1, nmolty
               rho = ncmt(ibox,jmolty) / vol
               do imolty = 1, nmolty
                  vtail = vtail + 
     &                 ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
               end do
            end do
!$$$            if (ibox .eq. 1 .and. lexzeo) then
!$$$               do jmolty = 1,zntype
!$$$                  rho = znum(jmolty)/vol
!$$$                  do imolty = 1, nmolty
!$$$                     vtail = vtail + ncmt(ibox,imolty)*
!$$$     &                    coruz(imolty,jmolty,rho,ibox)
!$$$                  end do
!$$$               end do
!$$$            end if
!-----
  
            vinter = vinter + vtail
!----
         end if
      end if

!$$$c      write(iou,*)
!$$$c      write(iou,*) '+++++++'
!$$$c      vtemp = velect
!$$$c      write(iou,*) 'direct space part:',velect*qqfact

      if ( ldielect ) then
         call dipole(ibox,0)
      end if
      
      if ( lewald ) then
         call recipsum(ibox,vrecipsum)
!---- update self terms and correction terms
         sself = 0.0d0
         correct = 0.0d0
! * combine to reduce numerical error
!         vsc = 0.0d0


! RP added for MPI
!         do i = 1,nchain
         do i = myid+1, nchain,numprocs          

            if (nboxi(i) .eq. ibox) then
               imolty = moltyp(i)
               do ii = 1,nunit(imolty)
                  sself = sself + qqu(i,ii)*qqu(i,ii)
! * 1.772.. is the square root of pi
!                  vsc = vsc - qqu(i,ii)*qqu(i,ii)*
!     &                 calp(ibox)/1.772453851d0
                  do jj = ii+1,nunit(imolty)
! * correct should only be calculated if ii and jj should NOT interact,
! * so only calculating it if lqinclu is false
!                    * this part is 1,2 and 1,3
                     if (.not. lqinclu(imolty,ii,jj)) then
                        rxuij = rxu(i,ii) - rxu(i,jj)
                        ryuij = ryu(i,ii) - ryu(i,jj)
                        rzuij = rzu(i,ii) - rzu(i,jj)
! --- JLR 11-17-09  need call to mimage for intrachain
                        if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
! --- END JLR 11-17-09
                        rij = dsqrt(rxuij*rxuij + ryuij*ryuij + 
     &                       rzuij*rzuij)
                        correct = correct + qqu(i,ii)*qqu(i,jj)*
     &                            (erfunc(calp(ibox)*rij)-1.0d0)/rij
!                        vsc = vsc + qqu(i,ii)*qqu(i,jj)*
!     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
                     elseif (lqinclu(imolty,ii,jj)) then
                        rxuij = rxu(i,ii) - rxu(i,jj)
                        ryuij = ryu(i,ii) - ryu(i,jj)
                        rzuij = rzu(i,ii) - rzu(i,jj)
! --- JLR 11-17-09  need call to mimage for intrachain
                        if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)
! --- END JLR 11-17-09
                        rij = dsqrt(rxuij*rxuij + ryuij*ryuij + 
     &                       rzuij*rzuij)
                        correct=correct+(1.0d0 - qscale2(imolty,ii,jj))
     &                           *qqu(i,ii)*qqu(i,jj)*
     &                           (erfunc(calp(ibox)*rij)-1.0d0)/rij
!                        vsc = vsc +
!     &                       (1.0d0 - qscale2(imolty,ii,jj))*qqu(i,ii)
!     *                             *qqu(i,jj)*
!     &                       (erfunc(calp(ibox)*rij)-1.0d0)/rij
                     end if
                  end do
               end do
            end if
         end do

! RP added for MPI

         CALL MPI_ALLREDUCE(correct, sum_correct,1,MPI_DOUBLE_PRECISION,
     &          MPI_SUM,MPI_COMM_WORLD,ierr)
         CALL MPI_ALLREDUCE(sself, sum_sself,1,MPI_DOUBLE_PRECISION,
     &          MPI_SUM,MPI_COMM_WORLD,ierr)

         sself = sum_sself
         correct = sum_correct

!            vdipole = (dipolex*dipolex+dipoley*dipoley+
!     &           dipolez*dipolez)*(2.0d0*onepi)/(3.0d0*
!     &           boxlx(ibox)**3.0d0)
!            write(iou,*) dipolex,dipoley,dipolez

         sself = -sself*calp(ibox)/dsqrt(onepi)
         velect = velect + sself + correct + vrecipsum/qqfact

!         velect = velect + vsc + vrecipsum/qqfact

      end if

!$$$c at this point velect contains all intermolecular charge interactions,
!$$$c plus the ewald self term and intramolecular corrections 
!$$$
!       write(iou,*)
!       write(iou,*) '== After Inter === velect is:',velect*qqfact
!$$$
!$$$       vtemp = velect

! ################################################################

! * have to recalculate ewald terms if volume changes
      if ( .not. lvol .or. (lvol .and. lewald) ) then

! *******************************
! *** INTRACHAIN INTERACTIONS ***
! *******************************
!         write(iou,*) 'starting intrachain'
! --- loop over all chains i 

! RP added for MPI
!         do i = 1, nchain
         do i = myid+1,nchain,numprocs

!            lcoulo(1,5) = .false.

! ### check if i is in relevant box ###
!          write(iou,*) 'nboxi(i),i,ibox',nboxi(i),i,ibox
            if ( nboxi(i) .eq. ibox ) then

               imolty = moltyp(i)

               do ii = 1, nunit(imolty)-1

!                write(iou,*) 'ntype(imolty,ii),ii',ntype(imolty,ii),ii

                  ntii = ntype(imolty,ii)
                  
                  rxui = rxu(i,ii)
                  ryui = ryu(i,ii)
                  rzui = rzu(i,ii)
                  
                  do jj = ii+1, nunit(imolty)
                     
                     if ( linclu(imolty,ii,jj) .or. 
     &                    lqinclu(imolty,ii,jj)) then
                        
                        ntjj = ntype(imolty,jj)

                        if (lexpsix .or. lmmff) then
                           ntij = (ntii+ntjj)/2
                        elseif (lninesix) then
                           ntij = (ntii-1)*nxatom + ntjj
                        elseif (lgenlj) then
                           ntij = (ntii-1)*nntype + ntjj
                        else
                           ntij = (ntii-1)*nntype + ntjj
                        end if
                        if (lexpee) rminsq = rminee(ntij)*rminee(ntij)
                        
                        rxuij = rxu(i,ii) - rxu(i,jj)
                        ryuij = ryu(i,ii) - ryu(i,jj)
                        rzuij = rzu(i,ii) - rzu(i,jj)
! --- JLR 11-17-09  need call to mimage for intrachain
                        if (lpbc) call mimage ( rxuij,ryuij,rzuij,ibox )
! --- END JLR 11-17-09
                        rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                        rij  = dsqrt(rijsq) 
                        if (linclu(imolty,ii,jj)) then

                           if ( rijsq .lt. rminsq .and. 
     &                          .not. lexpand(imolty)) then
                              if ( .not. lvol ) then
                                 write(iou,*) 'overlap intra'
                                 write(iou,*) 'rijsq rminsq', rijsq, 
     &                                rminsq
                                 write(iou,*) 'i ii', i, ii
                                 write(iou,*) 'i-pos', rxui,ryui,rzui
                                 write(iou,*) 'jj', jj
                                 write(iou,*) 'j-pos'
     &                                ,rxu(i,jj),ryu(i,jj),rzu(i,jj)
                              end if
                              ovrlap = .true.
!-------- RP added for MPI to compensate ovrlap
                              if(ovrlap .eq. .true.)then
                                 goto 299
                              end if
! -------------------------------
!                          return
                           elseif ( rijsq .lt. rcutsq .or. lijall) then
                          
                              if (L_vdW_table.or.L_bend_table.and.(.not.
     &                             (lexpand(imolty)))) then
                                 
                                 do mmm=1,inben(imolty,ii)
                                    if (ijben3(imolty,ii,mmm).eq.jj)then
                                       
                                       call lininter_bend(rij,
     &                                      tabulated_bend, 
     &                                      itben(imolty,ii,mmm))
                                       vintra = vintra + tabulated_bend
                                       
                                       goto 94
                                    end if
                                 end do
                                 
                                 call lininter_vdW(rij, tabulated_vdW, 
     &                                ntii, ntjj)
                                 vintra = vintra + tabulated_vdW
                                 
                              elseif (llj.and.(.not.(lexpand(imolty)
     &                                ))) then
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2=epsij(ntij)
                                 sr6 = sr2 * sr2 * sr2
                                 vintra = vintra 
     &                                + sr6*(sr6-1.0d0)*epsilon2
     &                                *ljscale(imolty,ii,jj)
                                 
                                 if (lainclu(imolty,ii,jj)) then
                                    vintra = vintra + 0.25d0 * 
     &                                   a15(a15type(imolty,ii,jj)) /
     &                                   ((rijsq**2)*(rijsq**2)*
     &                                   (rijsq**2))
                                 end if
                              elseif ( lsami ) then
                                 vintra = vintra + ljsami(rijsq,ntij)
                              elseif (lexpsix) then
                                 vintra = vintra + exsix(rijsq,ntij)
                              elseif ( lmuir ) then
                                 vintra = vintra + ljmuir(rijsq,ntij)
                              elseif (lmmff) then
                                 vintra = vintra + mmff(rijsq,ntij)
                              elseif (lninesix) then
                                 vintra = vintra + ninesix(rijsq,ntij)
                              elseif (lgenlj) then
                                 sr2 = sig2ij(ntij) / rijsq
                                 epsilon2=epsij(ntij)
                                 vintra = vintra + 
     &                                genlj(rijsq,sr2,epsilon2)
                              elseif ( lpsurf ) then
                                 vintra = vintra + ljpsur(rijsq,ntij)
                              else if (lshift) then
                                 sr2 = sig2ij(ntij) / rijsq
                                 sr6 = sr2 * sr2 * sr2
                                 vintra = vintra + 
     &                                (sr6*(sr6-1.0d0)*epsij(ntij)-
     &                                ecut(ntij))*ljscale(imolty,ii,jj) 
                              else
                                 if ( lexpand(imolty) ) then
                                    sigma2=(sigma(imolty,ii)+
     &                                   sigma(imolty,jj))/2.0d0
                                    sr2 = sigma2*sigma2/rijsq
                                    epsilon2 = dsqrt(epsilon(imolty,ii)
     &                                   *epsilon(imolty,jj))
                                 else
                                    sr2 = sig2ij(ntij) / rijsq
                                    epsilon2 = epsij(ntij)
                                 end if
                                 sr6 = sr2 * sr2 * sr2
                                 vintra = vintra + sr6*(sr6-1.0d0)
     &                                *epsilon2*ljscale(imolty,ii,jj)
                                 
! * OH 1-5 interaction
                                 if (lainclu(imolty,ii,jj)) then
                                    vintra = vintra + 0.25d0 * 
     &                                   a15(a15type(imolty,ii,jj)) /
     &                                   ((rijsq**2)*(rijsq**2)*
     &                                   (rijsq**2))
                                 end if
                              end if
                           end if
                        end if

! * calculate intramolecular charge interaction

 94                     if ( lchgall .and. lqchg(ntii)
     &                       .and. lqchg(ntjj)) then
                           if (lqinclu(imolty,ii,jj)) then
                              if ( lewald ) then
                                 my_velect = my_velect + 
     &                                qscale2(imolty,ii,jj)
     &                            * qqu(i,ii)*qqu(i,jj)*
     &                                erfunc(calp(ibox)*rij)/
     &                                rij
                              else
                                 my_velect = my_velect + 
     &                                qscale2(imolty,ii,jj)*qqu(i,ii)*
     &                                qqu(i,jj)/rij
                              end if
                           end if
                       
                        elseif ( lelect(imolty) .and.
     &                          (lqchg(ntii) .and. lqchg(ntjj)) ) then
                           if (lewald) then
                              if (lqinclu(imolty,ii,jj).and.
     &                             (rijsq.lt.rcutsq)) then
                                 my_velect = my_velect +
     &                                qscale2(imolty,ii,jj)*qqu(i,ii)*
     &                                qqu(i,jj)*erfunc(calp(ibox)*
     &                                rij)/rij
                              end if 
                           else
! --- All-Atom charges (charge-group look-up table)
                              iii = leaderq(imolty,ii)
                              jjj = leaderq(imolty,jj)
                              if ( iii .eq. ii .and. jjj .eq. jj ) then
! --- set up the table for neutral groups
                                 if ( rijsq .lt. rcutsq ) then
                                    lcoulo(iii,jjj) = .true.
                                 else
                                    lcoulo(iii,jjj) = .false.
                                 end if
                              end if
! *** set up table for neighboring groups- make sure they interact when 
! *** leaderqs are only 2 bonds apart.
                              if (.not. lqinclu(imolty,iii,jjj)) then
                                 lcoulo(iii,jjj)  = .true.
                              end if
                              if ( lcoulo(iii,jjj) ) then
                                 if (lqinclu(imolty,ii,jj)) then
                                    if (L_elect_table) then
!     check this                   if (rij.gt.rmin) then 
                                       call lininter_elect(rij, 
     &                                      tabulated_elect, ntii,ntjj)
                                       my_velect = my_velect + 
     &                                      qscale2(imolty,ii,jj)*
     &                                      qqu(i,ii)*qqu(i,jj)*
     &                                      tabulated_elect
!                                     end if
                                    else
                                       my_velect = my_velect + 
     &                                      qscale2(imolty,ii,jj)
     &                                      *qqu(i,ii)*qqu(i,jj)/rij
                                    end if
                                 end if
                              end if
                           end if
                        end if
                     end if
! end intramolecular charge
                  end do
               end do
            end if
         end do

! -----RP added for MPI----- Returning from ovrlap--------------

 299     continue
       
! KM don't check overlap until allreduce is finished
!       if(ovrlap .eq. .true.)then
!          write(iou,*)'924: in sumup ovrlap=',ovrlap,'myid=',myid
!       end if
       
         CALL MPI_ALLREDUCE(ovrlap,all_ovrlap,1,MPI_LOGICAL,MPI_LOR,
     &        MPI_COMM_WORLD,ierr)

         ovrlap = all_ovrlap
         if(ovrlap)then
!            write(iou,*)'941: in sumup ovrlap=',ovrlap,'myid=',myid
            return
         end if
! -----------------------------------------
         CALL MPI_ALLREDUCE(vintra, sum_vintra,1,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)

         CALL MPI_ALLREDUCE(my_velect, sum_my_velect,1
     &        ,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)
         
         vintra = sum_vintra
         velect = velect + sum_my_velect

!$$$       vtemp = velect - vtemp
!$$$  
!$$$c       write(iou,*) '== Intra Velect ===',vtemp*qqfact
!       write(iou,*) '== After Intra  === velect is:',velect*qqfact
!       write(iou,*) 'vintra ', vintra
!       write(iou,*) 'vinter ', vinter
!       write(6,*) 'test', vintra

         if ( .not. lsami .and. .not. lexpsix .and. .not. lmmff 
     &        .and. .not. lgenlj .and. .not. lninesix .and..not.lgaro
     &        .and..not.L_vdW_table) then
            vintra = 4.0d0 * vintra
         end if

! ################################################################

! *************************************
! *** INTRACHAIN FLUCQ INTERACTIONS ***
! *************************************
!c RP added for MPI
! removed at some point?
         do i = 1,nchain
!       do i = myid+1,nchain,numprocs

! --- calculate intramolecular flucq energy for chain i 
            if( nboxi(i) .eq. ibox ) then
               imolty = moltyp(i)
               if ( lelect(imolty) ) then

                  if ( lflucq(imolty) ) then
                     iunit = nunit(imolty)
                     do ii = 1, iunit
                        
                        ntii = ntype(imolty,ii)
                        qqii = qqu(i,ii)
                        do jj = ii, iunit
                           
                           if ( ii .eq. jj) then
                              vflucq = vflucq + xiq(ntii)*qqii 
     &                             + jayself(ntii)*qqii*qqii
                           else
                              ntjj = ntype(imolty,jj)
                              ntij = (ntii-1)*nntype + ntjj
                              
                              vflucq = vflucq 
     &                             + jayq(ntij)*qqii*qqu(i,jj)
                           end if
                        end do
                     end do
                     vflucq = vflucq - fqegp(imolty)
                  else
                     vflucq = 0.0d0
                  end if
               end if
            end if
         end do
         
!c RP added for MPI
!       CALL MPI_ALLREDUCE(vflucq, sum_vflucq,1,MPI_DOUBLE_PRECISION,
!     &        MPI_SUM,MPI_COMM_WORLD,ierr)

!       vflucq = sum_vflucq
! ------------------------------------------
! **************************************************
! *** CALCULATION OF VIB. + BEND. + TORS. ENERGY ***
! **************************************************

! NOTE here virtual coordinates can be used!!!
!
! RP added for MPI
!       do i = 1, nchain
         do i = myid+1,nchain,numprocs
          
            imolty = moltyp(i)

! ### check if i is in relevant box ###
            if ( nboxi(i) .eq. ibox ) then

! - branched and linear molecules with connectivity table -
! - go through entire chain -
! - calculate all bonds vectors and lengths
! - calculate all stretching, bending, and torsional potentials
! - that have an end-bead with an index smaller than the current bead
               do ii = 1, nunit(imolty)
                  rxui=rxu(i,ii)
                  ryui=ryu(i,ii)
                  rzui=rzu(i,ii)
                  do iivib = 1, invib(imolty,ii)
                     jj = ijvib(imolty,ii,iivib)
                     xvec(ii,jj) = rxu(i,jj) - rxui
                     yvec(ii,jj) = ryu(i,jj) - ryui
                     zvec(ii,jj) = rzu(i,jj) - rzui
                     distij2(ii,jj) = ( xvec(ii,jj)**2
     &                    + yvec(ii,jj)**2 + zvec(ii,jj)**2 )
                     distij(ii,jj) = dsqrt(distij2(ii,jj))
                     
                     if ( nunit(imolty) .ne. nugrow(imolty) )then
!                  --- account for explct atoms in opposite direction
                        xvec(jj,ii)   = -xvec(ii,jj)
                        yvec(jj,ii)   = -yvec(ii,jj)
                        zvec(jj,ii)   = -zvec(ii,jj)
                        distij(jj,ii) = distij(ii,jj)
                     end if
                  end do
               end do
               
! - stretching -
!             if ( brvibk(1) .gt. 0.01d0 .or. lninesix) then
               do j = 2, nunit(imolty)
                  do jjvib = 1, invib(imolty,j)
                     ip1 = ijvib(imolty,j,jjvib)
                     it  = itvib(imolty,j,jjvib)
                     if ( ip1. lt. j .and. L_vib_table) then
                        call lininter_vib(distij(ip1,j), 
     &                       tabulated_vib, it)
                        vvib = vvib + tabulated_vib
!                         write(2,*) 'TABULATED VVIB: ', tabulated_vib, 
!     &                        distij(ip1,j), ip1, j
                     end if
                     if ( ip1 .lt. j .and..not.L_vib_table) vvib = vvib
     &                    + brvibk(it) * (distij(ip1,j) - brvib(it))**2
                  end do
               end do
!     end if

! - bending -
! ### molecule with bond bending 
               do j = 2, nunit(imolty)
                  do jjben = 1, inben(imolty,j)
                     ip2 = ijben3(imolty,j,jjben)
                     if ( ip2 .lt. j ) then
                        ip1 = ijben2(imolty,j,jjben)
                        it  = itben(imolty,j,jjben)
                        thetac = ( xvec(ip1,j)*xvec(ip1,ip2) +
     &                       yvec(ip1,j)*yvec(ip1,ip2) +
     &                       zvec(ip1,j)*zvec(ip1,ip2) ) /
     &                       ( distij(ip1,j)*distij(ip1,ip2) )
                        if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
                        if ( thetac .le. -1.0d0 ) thetac = -1.0d0
                        
                        theta = dacos(thetac)

!$$$                      if (L_bend_table) then
!$$$                         rbendsq=distij2(ip1,j)+distij2(ip1,ip2)
!$$$     &                           -2.0d0*distij(ip1,j)*distij(ip1,ip2)
!$$$     &                           *thetac
!$$$                         rbend = dsqrt(rbendsq)
!$$$                         call lininter_bend(rbend, tabulated_bend, it)
!$$$                         vbend = vbend + tabulated_bend
!$$$                      else
                        vbend = vbend + 
     &                       brbenk(it) * (theta-brben(it))**2
!                      end if

!                        write(iou,*) 'j,ip1,ip2, it',j,ip1,ip2, it
!                        write(iou,*) 'bend energy, theta '
!     &                       ,brbenk(it) * (theta-brben(it))**2,theta

                     end if
                  end do
               end do
               
! - torsions -
! ### molecule with dihedral potenials ###
               do j = 2, nunit(imolty)
                  do jjtor = 1, intor(imolty,j)
                     ip3 = ijtor4(imolty,j,jjtor)
                     if ( ip3 .lt. j ) then
                        ip1 = ijtor2(imolty,j,jjtor)
                        ip2 = ijtor3(imolty,j,jjtor)
                        it  = ittor(imolty,j,jjtor)
!*** calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                        xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     &                       zvec(ip1,j) * yvec(ip1,ip2)
                        yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     &                       xvec(ip1,j) * zvec(ip1,ip2)
                        zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     &                       yvec(ip1,j) * xvec(ip1,ip2)
                        xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     &                       zvec(ip1,ip2) * yvec(ip3,ip2)
                        ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     &                       xvec(ip1,ip2) * zvec(ip3,ip2)
                        za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     &                       yvec(ip1,ip2) * xvec(ip3,ip2)
! *** calculate lengths of cross products ***
                        daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                        da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
! *** calculate dot product of cross products ***
                        dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                        thetac = - dot / ( daa1 * da1a2 )

                        if (thetac.gt.1.0d0) thetac=1.0d0
                        if (thetac.lt.-1.0d0) thetac=-1.0d0
!     KEA -- added for extending range to +/- 180
!     and for asymmetric potentials
                        if (L_tor_table) then
!     *** calculate cross product of cross products ***
                           xcc = yaa1*za1a2 - zaa1*ya1a2
                           ycc = zaa1*xa1a2 - xaa1*za1a2
                           zcc = xaa1*ya1a2 - yaa1*xa1a2
!     *** calculate scalar triple product ***
                           tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2)
     &                          + zcc*zvec(ip1,ip2)
                           theta = dacos (thetac)

                           if (tcc .lt. 0.0d0) theta = -theta
                           if (L_spline) then
                              call splint(theta,spltor,it)
                           elseif(L_linear) then
                              call lininter(theta,spltor,it)
                           end if
                           
                           vtg = vtg+spltor
                        else
                           vtg = vtg + vtorso( thetac, it )
                        end if
                     end if
                  end do
               end do
            end if
         end do
! RP added for MPI

         CALL MPI_ALLREDUCE(vbend, sum_vbend,1,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)

         CALL MPI_ALLREDUCE(vtg, sum_vtg,1,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)

         CALL MPI_ALLREDUCE(vvib, sum_vvib,1,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)

       
         vbend = sum_vbend
         vtg = sum_vtg
         vvib = sum_vvib
         
      end if
 
! ################################################################

! ***************************************************************
! *** CALCULATION OF INTERACTION ENERGY WITH EXTERNAL SURFACE ***
! ***************************************************************
 
! --- not if lgrand and ibox =2
!      if (.not.(lgrand.and.ibox.eq.2)) then
! --- for adsorption isotherms, don't calculate energy w/surface
! --- in box 2
!      if ( .not. (lslit.and.ibox.eq.2)) then
      if (ibox .eq. 1) then
         if ( ljoe .or. lsami .or. lmuir .or. 
     &        lexzeo .or. lgraphite .or. lslit ) then

! RP added for MPI
!         do i = 1, nchain
            do i = myid+1, nchain,numprocs

! ### check if i is in relevant box ###
               if ( nboxi(i) .eq. ibox ) then
                  imolty = moltyp(i)
                  do j = 1, nunit(imolty)
                     ntj = ntype(imolty,j)
                     if ( ljoe ) then
                        if ( extc12(ntj) .gt. 0.1d0 ) then
                           dzui = rzu(i,j) - extz0(ntj)
                           dz3  = dzui * dzui * dzui
                           dz12 = dz3**4
                           vext = vext + 
     &                          (extc12(ntj)/dz12) - (extc3(ntj)/dz3)  
                        end if
                     end if
		  
                     if (lslit) then
                        ntij = (ntj-1)*nntype + ntsubst
! -- calculate interaction with surface at the bottom of the box		  		  
                        vext = vext + slitpore(rzu(i,j),ntij)
! -- calculate interaction with the surface at the top of the box
                        dzui = boxlz(ibox)-rzu(i,j)
                        vext = vext +slitpore(dzui,ntij)
                     end if  
		  
                     if( lgraphite ) then
			ntij = (ntj-1)*nntype + ntsubst
			vext = vext + exgrph(rxu(i,j),ryu(i,j),
     &				      rzu(i,j),ntij)
                     end if
                     
                     if ( lsami ) 
     &                    vext = vext + exsami(rzu(i,j),ntj)
                     if ( lmuir ) 
     &                    vext = vext + exmuir(rzu(i,j),ntj)
                     if ( lexzeo ) vext = vext + 
     &                    exzeo(rxu(i,j),ryu(i,j),rzu(i,j),ntj)
                     
                  end do
               end if
            end do

! RP added for MPI
            CALL MPI_ALLREDUCE(vext, sum_vext,1,MPI_DOUBLE_PRECISION,
     &           MPI_SUM,MPI_COMM_WORLD,ierr)

            vext = sum_vext

         end if
!         if (ltailc .and. lexzeo) then
!     ---    add tailcorrections for the zeolite
!            rhoz=nzeo/(zeorx*zeory*zeorz)
!            vext=vext+nchbox(ibox)*coruz(iunit,rhoz,ibox)
!         end if
      end if

! **********************************************************************
! *** calculation of interaction energy with external electric field ***
! *** added 06/24/07 by KM
! **********************************************************************

      if (lelect_field) then
! RP added for MPI
!         do i = 1, nchain
         do i = myid+1, nchain,numprocs

            if (nboxi(i).eq.ibox) then
               if(lelect(moltyp(i))) then
                  do j = 1,nunit(moltyp(i)) 
                     vext = vext + v_elect_field(i,j,rzu(i,j),field)
                  end do
               end if 
            end if 
         end do
! RP added for MPI
         CALL MPI_ALLREDUCE(vext, sum_vext,1,MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)

       vext = sum_vext
       vext = vext * eXV_to_K 

      end if  

! --------------------------------------------------------------------
! calculation of additional gaussian potential needed in thermodynamic
! integration in stages b and c
! --------------------------------------------------------------------
                                                                                
      if (lmipsw) then
! RP added for MPI
!         do i = 1, nchain
         do i = myid+1,nchain,numprocs

            imolty = moltyp(i)
            if (lwell(imolty)) then
               rxui = xcm(i)
               ryui = ycm(i)
               rzui = zcm(i)
            do j = 1, nwell(imolty)*nunit(imolty)
               k = j-int(j/nunit(imolty))*nunit(imolty)
               if (k.eq.0) k = nunit(imolty)
               rxuij = rxui-rxwell(j,imolty)
               ryuij = ryui-rywell(j,imolty)
               rzuij = rzui-rzwell(j,imolty)
               call mimage(rxuij,ryuij,rzuij,ibox)
               rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
               rcm = rcut(ibox)+rcmu(i)
               rcmsq = rcm*rcm
               if (rijsq.lt.rcmsq) then
                  do ii = 1, nunit(imolty)
                     if (awell(ii,k,imolty).lt.1.0d-6) goto 666
                     rxui = rxu(i,ii)
                     ryui = ryu(i,ii)
                     rzui = rzu(i,ii)
                     rxuij = rxui-rxwell(j,imolty)
                     ryuij = ryui-rywell(j,imolty)
                     rzuij = rzui-rzwell(j,imolty)
                     call mimage(rxuij,ryuij,rzuij,ibox)
                     rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                     vwell = vwell-awell(ii,k,imolty)*dexp(-bwell*rijsq)
 666              end do
               end if
            end do
         end if
      end do
      
! RP added for MPI
      CALL MPI_ALLREDUCE(vwell, sum_vwell,1,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,MPI_COMM_WORLD,ierr)
      vwell = sum_vwell
      
      end if

! ----------------------------------------------------------------------------

!      write(6,*) 'self,corr:',
!     &     (velect-vrecipsum/qqfact)*qqfact
!      write(6,*) 'vsc, new self cor:',vsc*qqfact
!      write(6,*) 'recip space part :',vrecipsum
!      write(6,*) 'sc and recip:',vsc*qqfact + vrecipsum
      
      if (.not.L_elect_table) then
         velect = velect*qqfact
      end if

 
!      velect = qqfact*velect 

      v = vinter + vintra + vext + velect + vflucq + v3garo
     
!      write(iou,*) 'v in sumup',v
                                                                                
      vipsw = v
      vwellipsw = vwell
                                                                                
      if (lstagea) then
         v = (1.0d0-lambdais*(1.0d0-etais))*v
      elseif (lstageb) then
         v = etais*v+lambdais*vwell
      elseif (lstagec) then
         v = (etais+(1.0d0-etais)*lambdais)*v+(1.0d0-lambdais)*vwell
      end if
      
      v = v + vvib + vbend + vtg

      if ( .not. lvol.and.myid.eq.0 ) then
         write(iou,*)
         write(iou,*) 'sumup control'
         write(iou,*) 'number of chains', nmcount
         do i = 1, nmolty
            write(iou,*) 'number of chains of type',i,ncmt(ibox,i)
         end do
         write(iou,*) 'inter lj energy ', vinter
         write(iou,*) 'intra lj energy ', vintra
         if (ltailc) write(iou,*)
     &        'Tail correction ', vtail
         
!         if (ltailc.and.lexzeo) write(iou,*)
!     +        'Tail corr. zeol', nchbox(ibox)*coruz(iunit,rhoz,
!     +        ibox) 
         write(iou,*) 'bond vibration  ', vvib
         write(iou,*) 'bond bending    ', vbend
         write(iou,*) 'torsional       ', vtg
         write(iou,*) 'external        ', vext
         write(iou,*) 'coulombic energy', velect
!         write(iou,*) 'exact energy    ',
!     +        1.74756*1.67*831.441/3.292796
         write(iou,*) 'fluc Q energy   ', vflucq
         write(iou,*) 'well energy     ', vwellipsw
         if(lgaro) write(iou,*) '3-body garo     ', v3garo
         write(iou,*) 'total energy    ', v
      end if

!      write(iou,*) 'end SUMUP'

      return
      end
