      subroutine monper (acv,acpres,acsurf,acvolume,molfra,mnbox,asetel
     &       ,acdens,acmove,acnp,pres,nibox,nnn,nblock,lratio,lratv
     &       ,lprint,lmv,lrsave,lblock,lratfix,lsolute,acsolpar,
     &        acEnthalpy,acEnthalpy1)

! monper
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

! -----------------------------------------------------------------
! subroutine monper
! performs periodic updates and block averages
! -----------------------------------------------------------------
 
      implicit none

      include 'control.inc'
      include 'coord.inc'      
      include 'system.inc'
      include 'ensemble.inc'
      include 'cbmc.inc'
      include 'inputdata.inc'
      include 'blkavg.inc'
      include 'bnbsma.inc'
      include 'connect.inc'
      include 'coord2.inc'
      include 'cell.inc'
      include 'poten.inc'
      include 'mpi.inc'
 
      integer::nummol,ntii
      integer::nibox,im,nnn,ntot,nblock,imolty,m,mm,i,j,jjtor
     &     ,ibox,itype,itel,mnbox,zz,steps,igrow,jmolty,jbox
      integer::imend,itemp,bin,k,jjben,ip2,ip1,ip3,it,ii,jj,ivib
      logical::lratio,lratv,lprint,lmv,lrsave,lblock,lfq,lratfix
     &     ,lsolute,ovrlap
      double precision, dimension(nprop1,nbxmax,nbxmax):: acsolpar
      double precision, dimension(nbxmax):: acEnthalpy,acEnthalpy1
      real(8)::press1,press2,dp,dpp,debroglie,histrat
      real(8)::acv, molfra,acpres,acsurf,acvolume,asetel,acdens
     & ,histtot,acmove,acnp,dvalue,dnchoi,dnchoi1,dnchoih,dnunit,ratflcq
     & ,v,vintra,vinter,vext,velect,vewald,vtors,vtail,rho,coru
     & ,thetac,vbend,xvec,yvec,zvec,distij,theta,vtorso
     & ,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2,dot,daa1,da1a2
     & ,velect_intra,velect_inter

      real(8)::xcc,ycc,zcc,tcc,spltor,rcutmin      
 
      dimension acv(nener,nbxmax),lratfix(ntmax)
      dimension mnbox(nbxmax,ntmax)
      dimension acpres(nbxmax),acsurf(nbxmax),acvolume(nbxmax)
      dimension asetel(nbxmax,ntmax),acdens(nbxmax,ntmax)
     & ,molfra(nbxmax,ntmax)
      dimension lsolute(ntmax)
      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)
     &     ,distij(numax,numax)

      real(8)::ratrax,ratray,ratraz,rttrax,rttray,rttraz
     & ,rarotx,raroty,rarotz,rtrotx,rtroty,rtrotz,vol
     & ,ratvol,temmass,dn,pres(nbxmax)
! -------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
      DO im=1,nbxmax
         pres(im)=0.0d0
      end do
!      write(iou,*) 'begin MONPER'
! *** perform periodic operations  ***
 
! *** Deciding the minimum rcut for atom displacement. As we
! *** would have strecthing potential, this ought to be much
! *** smaller than the rcutmin

      rcutmin = 1.0d10

      do  im = 1,nbox
          if (rcutmin.gt.rcut(im)) then
             rcutmin = rcut(im)
          end if
      end do


      if ( lratio ) then
! *** adjust the atomic displacements
         if ( Abntrax .gt. 0.5d0 ) then
            ratrax = Abstrax / Abntrax
            rttrax = Armtrax * ratrax / tatra
            if ( rttrax .gt. 2.0d0 * rcutmin ) then
! --- maximum translational displacement
               Armtrax = 2.0d0*rcutmin

            elseif (rttrax .lt. 1.0d-10 ) then
! --- ratio must have been zero, so divide by 10
               Armtrax = Armtrax/10.0d0
               
            else
! --- accept new ratio
               Armtrax = rttrax
            end if
            
            if ( lneigh .and. Armtrax .ge. upnn) then
               Armtrax = upnn
               write(iou,*) '### problem : for target accept ',
     &              'ratio Armtrax should be smaller than upnn'
            end if
         end if


         if ( Abntray .gt. 0.5d0 ) then
            ratray = Abstray / Abntray
            rttray = Armtray * ratray / tatra
            if ( rttray .gt. 2.0d0 * rcutmin ) then
! --- maximum translational displacement
               Armtray = 2.0d0*rcutmin

            elseif (rttray .lt. 1.0d-10 ) then
! --- ratio must have been zero, so divide by 10
               Armtray = Armtray/10.0d0

            else
! --- accept new ratio
               Armtray = rttray
            end if

            if ( lneigh .and. Armtray .ge. upnn) then
               Armtray = upnn
               write(iou,*) '### problem : for target accept ',
     &              'ratio Armtray should be smaller than upnn'
            end if
         end if
        
         if ( Abntraz .gt. 0.5d0 ) then
            ratraz = Abstraz / Abntraz
            rttraz = Armtraz * ratraz / tatra
            if ( rttraz .gt. 2.0d0 * rcutin ) then
! --- maximum translational displacement
               Armtraz = 2.0d0*rcutin

            elseif (rttraz .lt. 1.0d-10 ) then
! --- ratio must have been zero, so divide by 10
               Armtraz = Armtraz/10.0d0

            else
! --- accept new ratio
               Armtraz = rttraz
            end if

            if ( lneigh .and. Armtraz .ge. upnn) then
               Armtraz = upnn
               write(iou,*) '### problem : for target accept ',
     &              'ratio Armtraz should be smaller than upnn'
            end if
         end if
 
! *** adjust maximum translational displacement for COM***
         imend = nbox
         do im=1,imend
            if (myid.eq.0) then
               write(iou,*) 'Box ',im
            end if
            do imolty = 1,nmolty
!              --- rmtrax
!              --- check whether any x translations have been done for 
!              --- molecule type imolty in box im
               if ( bntrax(imolty,im) .gt. 0.5d0 ) then

!                 --- compute possible new ratio
                  ratrax = bstrax(imolty,im) / bntrax(imolty,im)
                  rttrax = rmtrax(imolty,im) * ratrax / tatra

                  if ( rttrax .gt. 2.0d0 * rcut(im)) then 
!                    --- maximum translational displacement
                     rmtrax(imolty,im) = 2.0d0*rcut(im)
                     
                  elseif (rttrax .lt. 1.0d-10 ) then  
!                    --- ratio must have been zero, so divide by 10
                     rmtrax(imolty,im) = rmtrax(imolty,im)/10.0d0
                     
                  else 
!                    --- accept new ratio
                     rmtrax(imolty,im) = rttrax
                  end if

                  if ( lneigh .and. rmtrax(imolty,im) .ge. upnn) then
                     rmtrax(imolty,im) = upnn
                     write(iou,*) '### problem : for target accept '
     &                    ,'ratio rmtrax should be smaller than upnn'
                  end if
               end if

!              --- rmtray
!              --- check whether any y translations have been done for 
!              --- molecule type imolty in box im

               if ( bntray(imolty,im) .gt. 0.5d0 ) then

!                 --- compute possiblt new ratio
                  ratray = bstray(imolty,im) / bntray(imolty,im)
                  rttray = rmtray(imolty,im) * ratray / tatra

                  if ( rttray .gt. 2.0d0*rcut(im) ) then 
!                    --- maximum translational displacement
                     rmtray(imolty,im) = 2.0d0*rcut(im)
                     
                  elseif (rttray .eq. 0.0d0) then  
!                    --- ratio must have been zero, divide old by 10
                     rmtray(imolty,im) = rmtray(imolty,im)/10.0d0
                  else 
!                    --- accept new ratio
                     rmtray(imolty,im) = rttray
                  end if
                  if ( lneigh .and. rmtray(imolty,im) .ge. upnn) then
                     rmtray(imolty,im) = upnn
                     write(iou,*) '### problem : for target accept '
     &                    ,'ratio rmtray should be smaller than upnn'
                  end if
               end if

!              --- rmtraz
!              --- check whether any z translations have been done for 
!              --- molecule type imolty in box im

               if ( bntraz(imolty,im) .gt. 0.5d0 ) then

!                 --- compute possible new ratio
                  ratraz = bstraz(imolty,im) / bntraz(imolty,im)
                  rttraz = rmtraz(imolty,im) * ratraz / tatra

                  if ( rttraz .gt. 2.0d0*rcut(im) ) then 
!                    --- maximum translational displacement
                     rmtraz(imolty,im) = 2.0d0*rcut(im)

                  elseif ( rttraz .lt. 1.0d-10 ) then  
!                    --- ratio must have been zero, divide old by 10
                     rmtraz(imolty,im) = rmtraz(imolty,im)/10.0d0

                  else 
!                    --- accept new ratio
                     rmtraz(imolty,im) = rttraz
                  end if
!                 --- check neighbor list
                  if ( lneigh .and. rmtraz(imolty,im) .ge. upnn) then
                     rmtraz(imolty,im) = upnn
                     write(iou,*) '### problem : for target accept '
     &                    ,'ratio rmtraz should be smaller than upnn'
                  end if
               end if
 
!              --- rmrotx
!              --- check whether any x-axis rotations have been done for 
!              --- molecule type imolty in box im

               if ( bnrotx(imolty,im) .gt. 0.5d0 ) then

!                 --- compute possible new ratio
                  rarotx = bsrotx(imolty,im) / bnrotx(imolty,im)
                  rtrotx = rmrotx(imolty,im) * rarotx / tarot

                  if ( rtrotx .lt. 1.0d-10 )  then
!                    --- ratio was zero, divide old by 10
                     rmrotx(imolty,im) = rmrotx(imolty,im)/10.d0

                  elseif (rtrotx .gt. 3.1415d0) then
!                    --- maximum rotational displacement
                     rmrotx(imolty,im) = 3.1415d0

                  else
!                    --- accept trial ratio
                     rmrotx(imolty,im) = rtrotx
                  end if

!                 --- check neighbour list
                  if ( lneigh .and. rmrotx(imolty,im) .ge. upnndg) then
                     rmrotx(imolty,im) = upnndg
                     write(iou,*) '### problem : for target accept '
     &                    ,'ratio rmrotx should be smaller than upnndg'
                  end if
               end if

!              --- rmroty
!              --- check whether any y-axis rotations have been done for 
!              --- molecule type imolty in box im
               if ( bnroty(imolty,im) .gt. 0.5d0 ) then
                  
!                 --- compute possible new ratio
                  raroty = bsroty(imolty,im) / bnroty(imolty,im)
                  rtroty = rmroty(imolty,im) * raroty / tarot

                  if ( rtroty .lt. 1.0d-10)  then
!                    --- ratio was zero, divide old by 10
                     rmroty(imolty,im) = rmroty(imolty,im)/10.0d0
                  elseif (rtroty .gt. 3.1415d0) then
!                    --- maximum rotational displacement
                     rmroty(imolty,im) = 3.1415d0
                  else
!                    --- accept trial ratio
                     rmroty(imolty,im) = rtroty
                  end if

!                 --- check neighbour list
                  if ( lneigh .and. rmroty(imolty,im) .ge. upnndg) then
                     rmroty(imolty,im) = upnndg
                     write(iou,*) '### problem : for target accept'
     &                    ,' ratio rmroty should be smaller than upnndg'
                  end if
               end if
               
!              --- rmrotz
!              --- check whether any y-axis rotations have been done for 
!              --- molecule type imolty in box im

               if ( bnrotz(imolty,im) .gt. 0.5d0 ) then

!                 --- compute possible new ratio
                  rarotz = bsrotz(imolty,im) / bnrotz(imolty,im)
                  rtrotz = rmrotz(imolty,im) * rarotz / tarot
                  
                  if (rtrotz .eq. 0.0d0)  then
!                    --- ratio was zero, divide old by 10
                     rmrotz(imolty,im) = rmrotz(imolty,im)/10.0d0
                  elseif (rtrotz .gt. 3.1415d0) then
!                    --- maximum rotational displacement
                     rmrotz(imolty,im) = 3.1415d0
                  else
!                    --- accept trial ratio
                     rmrotz(imolty,im) = rtrotz
                  end if
!                 --- check neighbour list
                  if ( lneigh .and. rmrotz(imolty,im) .ge. upnndg) then
                     rmrotz(imolty,im) = upnndg
                     write(iou,*) '### problem : for target accept '
     &                    ,'ratio rmrotz should be smaller than upnndg'
                  end if
               end if

! KM for MPI
! only processor 0 writes to output files (except for error messages)
               if (myid.eq.0) then
! *** write some ratio update information ***
                  write(iou,57) imolty,bntrax(imolty,im),
     &                 bntray(imolty,im), bntraz(imolty,im), 
     &                 bnrotx(imolty,im), bnroty(imolty,im), 
     &                 bnrotz(imolty,im)
!              write(iou,58) ratrax, ratray, ratraz,rarotx, raroty,rarotz
                  write(iou,59) rmtrax(imolty,im), rmtray(imolty,im)
     &                 , rmtraz(imolty,im), rmrotx(imolty,im)
     &                 , rmroty(imolty,im), rmrotz(imolty,im)
               end if

 57            format(' Type ',i1,' bn ',6(f7.0,1x))
!     58            format(' ratio       ',6(f7.4,1x))
 59            format(' max.displ.  ',6(f7.4,1x))
               
               acntrax(imolty,im) = acntrax(imolty,im)+bntrax(imolty,im)
               acntray(imolty,im) = acntray(imolty,im)+bntray(imolty,im)
               acntraz(imolty,im) = acntraz(imolty,im)+bntraz(imolty,im)
               acstrax(imolty,im) = acstrax(imolty,im)+bstrax(imolty,im)
               acstray(imolty,im) = acstray(imolty,im)+bstray(imolty,im)
               acstraz(imolty,im) = acstraz(imolty,im)+bstraz(imolty,im)

               acnrotx(imolty,im) = acnrotx(imolty,im)+bnrotx(imolty,im)
               acnroty(imolty,im) = acnroty(imolty,im)+bnroty(imolty,im)
               acnrotz(imolty,im) = acnrotz(imolty,im)+bnrotz(imolty,im)
               acsrotx(imolty,im) = acsrotx(imolty,im)+bsrotx(imolty,im)
               acsroty(imolty,im) = acsroty(imolty,im)+bsroty(imolty,im)
               acsrotz(imolty,im) = acsrotz(imolty,im)+bsrotz(imolty,im)

               bstrax(imolty,im) = 0.0d0
               bstray(imolty,im) = 0.0d0
               bstraz(imolty,im) = 0.0d0
               bntrax(imolty,im) = 0.0d0
               bntray(imolty,im) = 0.0d0
               bntraz(imolty,im) = 0.0d0
               bsrotx(imolty,im) = 0.0d0
               bsroty(imolty,im) = 0.0d0
               bsrotz(imolty,im) = 0.0d0
               bnrotx(imolty,im) = 0.0d0
               bnroty(imolty,im) = 0.0d0
               bnrotz(imolty,im) = 0.0d0
            end do
         end do

!        ---  adjust maximum charge displacement for fluc Q
         lfq = .false.
         do i = 1,nmolty
            do im = 1,imend
               if ( bnflcq(i,im) .gt. 0.5d0 ) then
                  lfq = .true.
                  ratflcq = bsflcq(i,im)/(bnflcq(i,im)*taflcq)
                  if ( ratflcq .lt. 0.1d0 ) then
                     rmflcq(i,im) = rmflcq(i,im) * 0.1d0
                  else
                     rmflcq(i,im) = rmflcq(i,im) * ratflcq
                  end if
               end if
!              --- accumulate flcq info for final output
               bsflcq2(i,im) = bsflcq2(i,im) + bsflcq(i,im)
               bnflcq2(i,im) = bnflcq2(i,im) + bnflcq(i,im)
!              --- rezero flcq
               bsflcq(i,im) = 0.0d0
               bnflcq(i,im) = 0.0d0
            end do
         end do
         if ( lfq.and.myid.eq.0 ) then
!           --- write out information about fluctuating charge success
            write(iou,*) 'Box:   rmflcq for moltyps'
            do im =1,imend
               write(iou,*) im,(rmflcq(i,im),i=1,nmolty)
            end do
         end if
         
      end if  ! end if (lratio)

      do imolty = 1, nmolty
         if (lratfix(imolty)) then
!     *** readjust fixed end point data to assure optimum efficiency ***

            counttot = counttot + counthist
            histrat = dble(counthist) / dble(counttot) 

!     --- reset counthist
            counthist = 0

            do j = 1, iring(imolty) 
               do k = 1, iring(imolty)
                  if (j.eq.k) goto 150
                  histtot = 0
                  do bin = 1, maxbin
                     histtot = histtot + hist(j,k,bin)
                  end do

                  if (histtot.eq.0) goto 150
                  
!     *** normalize and multiply hist by its weighting using above condition
                  do bin = 1, maxbin
                     
                     hist(j,k,bin) = hist(j,k,bin)*dble(histrat)
     &                    / histtot
                                         
!     *** add weighed hist to hist

                     probf(j,k,bin) = probf(j,k,bin) + hist(j,k,bin)
                     
!     *** reset hist to zero for next iteration
                     hist(j,k,bin) = 0
                  end do
!     *** renormalize new distribution
                  histtot = 0
                  do bin = 1, maxbin
                     histtot = histtot + probf(j,k,bin)
                  end do
                  do bin = 1, maxbin
                     probf(j,k,bin) = probf(j,k,bin) / histtot
                  end do
 150              continue
               end do
            end do
         end if

! KM for MPI
! only do anything for lsolute if monper is not called from readdat (nnn.ne.0)
!     *** calculate energy and write out movie for lsolute
         if (lsolute(imolty).and.nnn.ne.0) then
            do ibox = 1, nbox
               do k = 1, ncmt(ibox,imolty)
                  i = parbox(k,ibox,imolty)

!     --- set coords for energy and write out conformations

                  if (myid.eq.0) write(11,*) imolty,ibox,nunit(imolty)

                     do j = 1, nunit(imolty)
                        rxuion(j,1) = rxu(i,j)                     
                        ryuion(j,1) = ryu(i,j)
                        rzuion(j,1) = rzu(i,j)
                        
                        if (myid.eq.0) write(11,*) ntype(imolty,j)
     &                       ,rxuion(j,1),ryuion(j,1),rzuion(j,1)
     &                       ,qqu(j,1)

                     end do

                  call energy(i,imolty,v,vintra,vinter,vext,velect
     &                 ,vewald,1,ibox,1,nunit(imolty),.true.,ovrlap
     &                 ,.false.,vtors,.false.,.false.)
                  if (ovrlap) write(iou,*) 
     &                 '*** DISASTER, OVERLAP IN MONPER'

                  if (lsolid(ibox).and..not.lrect(ibox)) then
                     vol = cell_vol(ibox) 

                  end if

!     --- tail corrections
                  vtail = 0.0d0
                  do jmolty = 1, nmolty
                     
                     rho = ncmt(ibox,jmolty) / 
     &                    ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
                     
                     vtail = vtail + ncmt(ibox,imolty)
     &                    * coru(imolty,jmolty,rho,ibox)
                  end do

!     --- bending energy
                  do ii = 1, nunit(imolty)
                     do ivib = 1, invib(imolty,ii)
                        jj = ijvib(imolty,ii,ivib)
                        xvec(ii,jj) = rxu(i,jj) - rxu(i,ii)
                        yvec(ii,jj) = ryu(i,jj) - ryu(i,ii)
                        zvec(ii,jj) = rzu(i,jj) - rzu(i,ii)
                        distij(ii,jj) = dsqrt(xvec(ii,jj)**2
     &                       + yvec(ii,jj)**2 + zvec(ii,jj)**2)
                        
                        xvec(jj,ii)   = -xvec(ii,jj)
                        yvec(jj,ii)   = -yvec(ii,jj)
                        zvec(jj,ii)   = -zvec(ii,jj)
                        distij(jj,ii) = distij(ii,jj)
                     end do
                  end do
                     
                  vbend = 0.0d0

                  do j = 3, nunit(imolty)
                     do jjben = 1, inben(imolty,j)
                        ip2 = ijben3(imolty,j,jjben)
                        if ( ip2 .lt. j ) then
                           ip1 = ijben2(imolty,j,jjben)
                           it  = itben(imolty,j,jjben)
                           if (brbenk(it).gt.0.001d0) then

                              thetac = ( xvec(ip1,j)*xvec(ip1,ip2) +
     &                             yvec(ip1,j)*yvec(ip1,ip2) +
     &                             zvec(ip1,j)*zvec(ip1,ip2) ) /
     &                             ( distij(ip1,j)*distij(ip1,ip2) )
                              if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
                              if ( thetac .le. -1.0d0 ) thetac = -1.0d0
                              
                              theta = dacos(thetac)
                              vbend = vbend + 
     &                             brbenk(it) * (theta-brben(it))**2
                           end if

                        end if
                     end do
                  end do

!     --- torsional energy
                  vtors = 0.0d0
                  do j = 4, nunit(imolty)
                     do jjtor = 1, intor(imolty,j)
                        ip3 = ijtor4(imolty,j,jjtor)
                        if ( ip3 .lt. j ) then
                           ip1 = ijtor2(imolty,j,jjtor)
                           ip2 = ijtor3(imolty,j,jjtor)
                           it  = ittor(imolty,j,jjtor)
!***  calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                           xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     &                          zvec(ip1,j) * yvec(ip1,ip2)
                           yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     &                          xvec(ip1,j) * zvec(ip1,ip2)
                           zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     &                          yvec(ip1,j) * xvec(ip1,ip2)
                           xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     &                          zvec(ip1,ip2) * yvec(ip3,ip2)
                           ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     &                          xvec(ip1,ip2) * zvec(ip3,ip2)
                           za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     &                          yvec(ip1,ip2) * xvec(ip3,ip2)
!     *** calculate lengths of cross products ***
                           daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                           da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
!     *** calculate dot product of cross products ***
                           dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                           thetac = - dot / ( daa1 * da1a2 )

                           if (thetac.gt.1.0d0) thetac=1.0d0
                           if (thetac.lt.-1.0d0) thetac=-1.0d0
!     KEA -- added for extending range to +/- 180
                           if (L_tor_table) then
!     *** calculate cross product of cross products ***
                              xcc = yaa1*za1a2 - zaa1*ya1a2
                              ycc = zaa1*xa1a2 - xaa1*za1a2
                              zcc = xaa1*ya1a2 - yaa1*xa1a2
!     *** calculate scalar triple product ***
                              tcc = xcc*xvec(ip1,ip2)+ycc*yvec(ip1,ip2)
     &                             + zcc*zvec(ip1,ip2)
                              theta = dacos(thetac)
                              if (tcc .lt. 0.0d0) theta = -theta
                              if (L_spline) then
                                 call splint(theta,spltor,it)
                              elseif(L_linear) then
                                 call lininter(theta,spltor,it)
                              end if

                              vtors = vtors + spltor
                           else
                              vtors = vtors + vtorso( thetac, it )
                           end if
                        end if
                     end do
                  end do
                 
                  solcount(ibox,imolty) = solcount(ibox,imolty) + 1
                  avsolinter(ibox,imolty) = 
     &                 avsolinter(ibox,imolty) + vinter / 2.0d0 + vtail
                  avsolintra(ibox,imolty) = avsolintra(ibox,imolty) 
     &                 + vintra 
                  avsoltor(ibox,imolty) = avsoltor(ibox,imolty) + vtors
                  
                  avsolbend(ibox,imolty) = avsolbend(ibox,imolty) 
     &                 + vbend

                  avsolelc(ibox,imolty) = avsolelc(ibox,imolty) 
     &                 + velect + vewald

               end do

            end do
            
         end if   
      end do

      if ( lgibbs .or. lnpt ) then
         if ( lratv ) then
!     *** adjust maximum volume displacement ***
            do ibox = 1, nbox

               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  do j = 1,9
                     if ( bnhmat(ibox,j) .gt. 0.5d0 ) then
                        ratvol = bshmat(ibox,j) / bnhmat(ibox,j)
                        if (ratvol .eq. 0.0d0) then
                           rmhmat(ibox,j) = rmhmat(ibox,j) * 0.1d0
                        else
                           rmhmat(ibox,j) = rmhmat(ibox,j)*ratvol/tavol
                        end if
                     end if
                  end do
               else
                  if ( bnvol(ibox) .gt. 0.5d0 ) then
                     ratvol = bsvol(ibox) / bnvol(ibox)
                     if ( ratvol .eq. 0.0d0 ) then
                        rmvol(ibox) = rmvol(ibox) * 0.1d0
                     else
                        rmvol(ibox) = rmvol(ibox) * ratvol / tavol
                     end if
                  end if
               end if

            end do

! KM for MPI            
            if (myid.eq.0) then
               do ibox = 1, nbox
                  if (lsolid(ibox) .and. .not. lrect(ibox)) then
                     do j = 1,9
                        write(iou,56) bnhmat(ibox,j),bshmat(ibox,j),
     &                       rmhmat(ibox,j)
                     end do
                  else
                     write(iou,60) bnvol(ibox),bsvol(ibox),rmvol(ibox)
                  end if
               end do
            end if	

            do ibox = 1, nbox

               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  do j = 1,9
                     acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
                     acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
                     bshmat(ibox,j) = 0.0d0
                     bnhmat(ibox,j) = 0.0d0
                  end do
               else
                  acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
                  acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
                  bnvol(ibox) = 0.0d0
                  bsvol(ibox) = 0.0d0
               end if

            end do
         end if

 56      format(' h-matrix change:  bn =',f8.1,
     &        '   bs =',f8.1,'   max.displ. =',e12.5)
 60      format(' volume change:  bn =',f8.1,
     &        '   bs =',f8.1,'   max.displ. =',e12.5)
      end if

! KM for MPI

      if ( lprint.and.myid.eq.0 ) then
!     *** write out runtime information ***
         ntot = nnn + nnstep
         write(iou,'(i6,i8,e12.4,f10.3,f12.1,15i4)') nnn,ntot,
     &        vbox(1),boxlx(1),pres(1)
     &        ,(ncmt(1,imolty),imolty=1,nmolty)
         if ( lgibbs ) then
            do ibox = 2, nbox
               write(iou,'(14x,e12.4,f10.3,f12.1,15i4)') 
     &              vbox(ibox),boxlx(ibox),pres(ibox)
     &              ,(ncmt(ibox,imolty),imolty=1,nmolty)
            end do
         end if
      end if

! KM for MPI      
      if ( lmv .and.myid.eq.0) then
!     *** write out the movie configurations ***
         write(10,*) nnn
         do ibox = 1, nbox
            write(10,*) (ncmt(ibox,zz),zz=1,nmolty)

            if (lsolid(ibox) .and. .not. lrect(ibox)) then
               write(10,*) (hmat(ibox,zz),zz=1,9)
            else
               write(10,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
            end if
         end do
         do m = 1, nchain
            imolty = moltyp(m)
            write(10,1015) m, imolty, nunit(imolty), nboxi(m)
     &           ,xcm(m),ycm(m),zcm(m)
            do mm = 1, nunit(imolty)
               write(10,1012) rxu(m,mm), ryu(m,mm), rzu(m,mm)
     &              , qqu(m,mm),ntype( imolty, mm )
            end do
         end do
 1012    format(4(1x,f14.6),i5)
 1015    format(4(1x,i5),3(1x,f16.6))
      end if

! KM for MPI
      if ( lmv .and. L_movie_xyz.and.myid.eq.0) then
         do ibox = 1,nbox
           nummol = 0
           do i = 1,nchain
              if (nboxi(i).eq.ibox) then
                  nummol = nummol + nunit(moltyp(i))
              end if
           end do
           write(210+ibox,*) nummol
           write(210+ibox,*)
           do i = 1,nchain
              if(nboxi(i).eq.ibox) then
                 imolty = moltyp(i)
                 do ii = 1,nunit(imolty)
                    ntii = ntype(imolty,ii)
                    write(210+ibox,'(a4,5x,3f15.4)') chemid(ntii),
     &               rxu(i,ii), ryu(i,ii), rzu(i,ii)
                 end do
              end if
           end do
         end do
      end if

      if ( lrsave .and.myid.eq.0) then
!     *** write out the restart configurations to SAFETY-file ***
         open(unit=88,file="save-config")
         write (88,*) nnn + nnstep
         write (88,*) Armtrax, Armtray, Armtraz 
         do im=1,nbox
            do imolty=1,nmolty
               write (88,*) rmtrax(imolty,im), rmtray(imolty,im)
     &              , rmtraz(imolty,im)
               write (88,*) rmrotx(imolty,im), rmroty(imolty,im)
     &              , rmrotz(imolty,im)
            end do
         end do
         do im=1,nbox
            write (88,*) (rmflcq(i,im),i=1,nmolty)
         end do
         if ( lgibbs .or. lgrand .or. lnpt) then
! KM 01/10 
            write (88,*) (rmvol(ibox),ibox=1,nbox)
            do im = 1,nbox

               if (lsolid(im) .and. .not. lrect(im)) then
                  write(88,*) (rmhmat(im,zz),zz=1,9)
                  write(88,*) (hmat(im,zz),zz=1,9)
               else
                  write (88,*) boxlx(im),boxly(im),boxlz(im)
               end if
            end do
         end if
         write (88,*) nchain
         write (88,*) nmolty
         write (88,*) (nunit(i),i=1,nmolty)
         write (88,*) (moltyp(i),i=1,nchain)
         write (88,*) (nboxi(i),i=1,nchain)
         do i = 1, nmolty
            if ( lexpand(i) ) write(88,*) eetype(i)
         end do
         do i = 1, nmolty
            if ( lexpand(i) ) write(88,*) rmexpc(i)
         end do
         do i = 1, nchain
            imolty = moltyp(i)
            do j = 1, nunit(imolty)
               write (88,*) rxu(i,j), ryu(i,j), rzu(i,j),qqu(i,j)
            end do
         end do
! * rewind doesn't always finish writing
!         rewind(88)
         close(88)
      end if
         
!     ***    calculation of block averages ***
      if ( lblock )  then
         nblock = nblock + 1
         do ibox=1,nibox

! 1                                           = specific density
! 2                                           = pressure
! 3                    to (2+nener)           = energies
! 1+(2+nener)          to (2+nener) +  nmolty = chemical potential
! 1+(2+nener)+  nmolty to (2+nener) +2*nmolty = square-end-to-end-length
! 1+(2+nener)+2*nmolty to (2+nener) +3*nmolty = number density
! 1+(2+nener)+3*nmolty to (2+nener) +4*nmolty = mol fraction
! 4+nener+4*nmolty+1                      = Enthalpy
! ---------------------------------------------------------

! - specific density
            if ( lpbcz ) then
               temmass = 0.0d0
               do itype = 1, nmolty
                  temmass = temmass + masst(itype)*acdens(ibox,itype)
               end do
               dp = temmass / 0.6022045d0
            else
               temmass = 0.0d0
               do itype = 1, nmolty
                  temmass = temmass + acdens(ibox,itype)
               end do
               dp = temmass
            end if
            call update(nblock,1,ibox,dp,acmove)

! - pressure
            call update(nblock,2,ibox,acpres(ibox),acnp)

! - surface tension
            itel = 2+nener+4*nmolty+1
            call update(nblock,itel,ibox,acsurf(ibox),acnp)

! - box volume
            itel = 4+nener+4*nmolty
            call update(nblock,itel,ibox,acvolume(ibox),acmove)

! - energies
            do j=1,nener
               itel=j+2
               call update(nblock,itel,ibox,acv(j,ibox),acmove)
            end do

! - chemical potential
            do itype = 1, nmolty
               itel = (2+nener) + itype
               dnchoi = dble(nchoi(itype))
               dnchoi1 = dble(nchoi1(itype))
!                 --- determine how many steps it takes to grow the 
!                 --- molecule not counting the first inserted bead
               igrow = nugrow(itype)
! --- need the first schedule call for rigid molecules, else it will seg fault
               if (lrigid(itype))then
                  call schedule(igrow,itype,steps,1,0,4)
               else
                  call schedule(igrow,itype,steps,1,0,2)
               end if
               dnunit = dble(steps)
               dnchoih = dble(nchoih(itype))
               debroglie = 17.458d0/( dsqrt(masst(itype)/beta ))
               dpp = acchem(ibox,itype)/( dnchoih*
     &              dnchoi1*(dnchoi**(dnunit))
     &              *debroglie*debroglie*debroglie )
               dvalue = bnchem(ibox,itype)
               call update(nblock,itel,ibox,dpp,dvalue)
            end do

! - square end-to-end length
            do itype = 1, nmolty
               itel = (2 + nener) + nmolty + itype
               if ( mnbox(ibox,itype) .eq. 0 ) then
!                 avoid division by zero in update
                  dn = 1.0d0
               else
                  dn = dble(mnbox(ibox,itype))
               end if
               call update(nblock,itel,ibox,asetel(ibox,itype),dn)
            end do

! - number density
            do itype = 1, nmolty
               itel = (2 + nener) + (2*nmolty) + itype
               dpp = acdens(ibox,itype)
               call update(nblock,itel,ibox,dpp,acmove)
            end do
! - mol fraction
            do itype = 1, nmolty
               itel = (2 + nener) + (3*nmolty) + itype
               dpp = molfra(ibox,itype)
               call update(nblock,itel,ibox,dpp,acmove)
            end do

! - Enthalpy
            itel = 4+nener+4*nmolty+1
!            write(iou,*) acEnthalpy(ibox)
            call update(nblock,itel,ibox,acEnthalpy(ibox),acnp)
            itel = 4+nener+4*nmolty+2
            call update(nblock,itel,ibox,acEnthalpy1(ibox),acnp)
         end do
         do ibox = 1,nibox-1
            do jbox = ibox+1,nibox 
               do i = 1,nprop1
                  call update1(nblock,i,acsolpar(i,ibox,jbox),acnp,
     &                         ibox,jbox)  
               end do
            end do
         end do
      end if
      
!      write(iou,*) 'end MONPER'
      return
      end
