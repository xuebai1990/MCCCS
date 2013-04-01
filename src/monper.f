! -----------------------------------------------------------------
! subroutine monper
! performs periodic updates and block averages
! -----------------------------------------------------------------
subroutine monper (acv,acpres,acsurf,acvolume,molfra,mnbox,asetel ,acdens,acmove,acnp,pres,nibox,nnn,nblock,lratio,lratv ,lprint,lmv,lrsave,lblock,lratfix,lsolute,acsolpar, acEnthalpy,acEnthalpy1)
  use const_phys,only:debroglie_factor
  use sim_particle,only:check_neighbor_list
  use sim_system
  use sim_cell
  use energy_intramolecular,only:U_bonded
  use energy_pairwise,only:energy,coru
  use moves_simple
  use moves_volume
  use moves_cbmc,only:schedule
  use transfer_swap,only:acchem,bnchem
  implicit none

  integer::nummol,ntii
  integer::nibox,im,nnn,ntot,nblock,imolty,m,mm,i,j,ibox,itype,itel,mnbox,zzz,steps,igrow,jmolty,jbox
  integer::imend,bin,k,ii
  logical::lratio,lratv,lprint,lmv,lrsave,lblock,lfq,lratfix,lsolute,ovrlap
  real, dimension(nprop1,nbxmax,nbxmax)::acsolpar
  real, dimension(nbxmax)::acEnthalpy,acEnthalpy1
  real::dpr,dpp,debroglie,histrat
  real::acv,molfra,acpres,acsurf,acvolume,asetel,acdens,histtot,acmove,acnp,dvalue,ratflcq,v,vintra,vinter,vext,velect,vewald,vtors,vtail,rho,vbend,rcutmin,vvib

  dimension acv(nener,nbxmax),lratfix(ntmax)
  dimension mnbox(nbxmax,ntmax)
  dimension acpres(nbxmax),acsurf(nbxmax),acvolume(nbxmax)
  dimension asetel(nbxmax,ntmax),acdens(nbxmax,ntmax) ,molfra(nbxmax,ntmax)
  dimension lsolute(ntmax)

  real::ratrax,ratray,ratraz,rttrax,rttray,rttraz,rarotx,raroty,rarotz,rtrotx,rtroty,rtrotz,vol,ratvol,temmass,dn,pres(nbxmax)
! -------------------------------------------------------------------
#ifdef __DEBUG__
  write(io_output,*) 'begin MONPER in ',myid
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
  pres=0.0E0_dp
! perform periodic operations  ***
 
! Deciding the minimum rcut for atom displacement. As we
! would have strecthing potential, this ought to be much
! smaller than the rcutmin

  rcutmin = 1.0E10_dp

  do  im = 1,nbox
     if (rcutmin.gt.rcut(im)) then
        rcutmin = rcut(im)
     end if
  end do

  if ( lratio ) then
! adjust the atomic displacements
     if ( Abntrax .gt. 0.5E0_dp ) then
        ratrax = Abstrax / Abntrax
        rttrax = Armtrax * ratrax / tatra
        if ( rttrax .gt. 2.0E0_dp * rcutmin ) then
! maximum translational displacement
           Armtrax = 2.0E0_dp*rcutmin
        else if (rttrax .lt. 1.0E-10_dp ) then
! ratio must have been zero, so divide by 10
           Armtrax = Armtrax/10.0E0_dp
        else
! accept new ratio
           Armtrax = rttrax
        end if
     end if

     if ( Abntray .gt. 0.5E0_dp ) then
        ratray = Abstray / Abntray
        rttray = Armtray * ratray / tatra
        if ( rttray .gt. 2.0E0_dp * rcutmin ) then
! maximum translational displacement
           Armtray = 2.0E0_dp*rcutmin
        else if (rttray .lt. 1.0E-10_dp ) then
! ratio must have been zero, so divide by 10
           Armtray = Armtray/10.0E0_dp
        else
! accept new ratio
           Armtray = rttray
        end if
     end if
        
     if ( Abntraz .gt. 0.5E0_dp ) then
        ratraz = Abstraz / Abntraz
        rttraz = Armtraz * ratraz / tatra
        if ( rttraz .gt. 2.0E0_dp * rcutin ) then
! maximum translational displacement
           Armtraz = 2.0E0_dp*rcutin
        else if (rttraz .lt. 1.0E-10_dp ) then
! ratio must have been zero, so divide by 10
           Armtraz = Armtraz/10.0E0_dp
        else
! accept new ratio
           Armtraz = rttraz
        end if
     end if
 
! adjust maximum translational displacement for COM***
     imend = nbox
     do im=1,imend
        if (myid.eq.0) then
           write(io_output,*) 'Box ',im
        end if
        do imolty = 1,nmolty
! rmtrax
! check whether any x translations have been done for 
! molecule type imolty in box im
           if ( bntrax(imolty,im) .gt. 0.5E0_dp ) then
! compute possible new ratio
              ratrax = bstrax(imolty,im) / bntrax(imolty,im)
              rttrax = rmtrax(imolty,im) * ratrax / tatra

              if ( rttrax .gt. 2.0E0_dp * rcut(im)) then 
! maximum translational displacement
                 rmtrax(imolty,im) = 2.0E0_dp*rcut(im)
              else if (rttrax .lt. 1.0E-10_dp ) then  
! ratio must have been zero, so divide by 10
                 rmtrax(imolty,im) = rmtrax(imolty,im)/10.0E0_dp
              else 
! accept new ratio
                 rmtrax(imolty,im) = rttrax
              end if
           end if

! rmtray
! check whether any y translations have been done for 
! molecule type imolty in box im
           if ( bntray(imolty,im) .gt. 0.5E0_dp ) then
! compute possiblt new ratio
              ratray = bstray(imolty,im) / bntray(imolty,im)
              rttray = rmtray(imolty,im) * ratray / tatra

              if ( rttray .gt. 2.0E0_dp*rcut(im) ) then 
! maximum translational displacement
                 rmtray(imolty,im) = 2.0E0_dp*rcut(im)
              else if (rttray .eq. 0.0E0_dp) then  
! ratio must have been zero, divide old by 10
                 rmtray(imolty,im) = rmtray(imolty,im)/10.0E0_dp
              else 
! accept new ratio
                 rmtray(imolty,im) = rttray
              end if
           end if

! rmtraz
! check whether any z translations have been done for 
! molecule type imolty in box im
           if ( bntraz(imolty,im) .gt. 0.5E0_dp ) then
! compute possible new ratio
              ratraz = bstraz(imolty,im) / bntraz(imolty,im)
              rttraz = rmtraz(imolty,im) * ratraz / tatra

              if ( rttraz .gt. 2.0E0_dp*rcut(im) ) then 
! maximum translational displacement
                 rmtraz(imolty,im) = 2.0E0_dp*rcut(im)
              else if ( rttraz .lt. 1.0E-10_dp ) then  
! ratio must have been zero, divide old by 10
                 rmtraz(imolty,im) = rmtraz(imolty,im)/10.0E0_dp
              else 
! accept new ratio
                 rmtraz(imolty,im) = rttraz
              end if
           end if
 
! rmrotx
! check whether any x-axis rotations have been done for 
! molecule type imolty in box im
           if ( bnrotx(imolty,im) .gt. 0.5E0_dp ) then
! compute possible new ratio
              rarotx = bsrotx(imolty,im) / bnrotx(imolty,im)
              rtrotx = rmrotx(imolty,im) * rarotx / tarot

              if ( rtrotx .lt. 1.0E-10_dp )  then
! ratio was zero, divide old by 10
                 rmrotx(imolty,im) = rmrotx(imolty,im)/10.E0_dp
              else if (rtrotx .gt. 3.1415E0_dp) then
! maximum rotational displacement
                 rmrotx(imolty,im) = 3.1415E0_dp
              else
! accept trial ratio
                 rmrotx(imolty,im) = rtrotx
              end if
           end if

! rmroty
! check whether any y-axis rotations have been done for 
! molecule type imolty in box im
           if ( bnroty(imolty,im) .gt. 0.5E0_dp ) then
! compute possible new ratio
              raroty = bsroty(imolty,im) / bnroty(imolty,im)
              rtroty = rmroty(imolty,im) * raroty / tarot

              if ( rtroty .lt. 1.0E-10_dp)  then
! ratio was zero, divide old by 10
                 rmroty(imolty,im) = rmroty(imolty,im)/10.0E0_dp
              else if (rtroty .gt. 3.1415E0_dp) then
! maximum rotational displacement
                 rmroty(imolty,im) = 3.1415E0_dp
              else
! accept trial ratio
                 rmroty(imolty,im) = rtroty
              end if
           end if
               
! rmrotz
! check whether any y-axis rotations have been done for 
! molecule type imolty in box im
           if ( bnrotz(imolty,im) .gt. 0.5E0_dp ) then
! compute possible new ratio
              rarotz = bsrotz(imolty,im) / bnrotz(imolty,im)
              rtrotz = rmrotz(imolty,im) * rarotz / tarot
                  
              if (rtrotz .eq. 0.0E0_dp)  then
! ratio was zero, divide old by 10
                 rmrotz(imolty,im) = rmrotz(imolty,im)/10.0E0_dp
              else if (rtrotz .gt. 3.1415E0_dp) then
! maximum rotational displacement
                 rmrotz(imolty,im) = 3.1415E0_dp
              else
! accept trial ratio
                 rmrotz(imolty,im) = rtrotz
              end if
           end if

           if (lneigh) call check_neighbor_list(im,imolty)
           call averageMaximumDisplacement(im,imolty)

! KM for MPI
! only processor 0 writes to output files (except for error messages)
           if (myid.eq.0) then
! write some ratio update information ***
              write(io_output,"(' Type ',i2,' bn ',6(f7.0,1x))") imolty ,bntrax(imolty,im),bntray(imolty,im), bntraz(imolty ,im),bnrotx(imolty,im), bnroty(imolty,im) ,bnrotz(imolty,im)
! write(io_output,"(' ratio       ',6(f7.4,1x))") ratrax, ratray, ratraz,rarotx, raroty,rarotz
              write(io_output,"(' max.displ.  ',6(f9.4,1x))") rmtrax(imolty,im), rmtray(imolty,im), rmtraz(imolty ,im), rmrotx(imolty,im), rmroty(imolty,im), rmrotz(imolty,im)
           end if

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

           bstrax(imolty,im) = 0.0E0_dp
           bstray(imolty,im) = 0.0E0_dp
           bstraz(imolty,im) = 0.0E0_dp
           bntrax(imolty,im) = 0.0E0_dp
           bntray(imolty,im) = 0.0E0_dp
           bntraz(imolty,im) = 0.0E0_dp
           bsrotx(imolty,im) = 0.0E0_dp
           bsroty(imolty,im) = 0.0E0_dp
           bsrotz(imolty,im) = 0.0E0_dp
           bnrotx(imolty,im) = 0.0E0_dp
           bnroty(imolty,im) = 0.0E0_dp
           bnrotz(imolty,im) = 0.0E0_dp
        end do
     end do

! adjust maximum charge displacement for fluc Q
     lfq = .false.
     do i = 1,nmolty
        do im = 1,imend
           if ( bnflcq(i,im) .gt. 0.5E0_dp ) then
              lfq = .true.
              ratflcq = bsflcq(i,im)/(bnflcq(i,im)*taflcq)
              if ( ratflcq .lt. 0.1E0_dp ) then
                 rmflcq(i,im) = rmflcq(i,im) * 0.1E0_dp
              else
                 rmflcq(i,im) = rmflcq(i,im) * ratflcq
              end if
           end if
! accumulate flcq info for final output
           bsflcq2(i,im) = bsflcq2(i,im) + bsflcq(i,im)
           bnflcq2(i,im) = bnflcq2(i,im) + bnflcq(i,im)
! rezero flcq
           bsflcq(i,im) = 0.0E0_dp
           bnflcq(i,im) = 0.0E0_dp
        end do
     end do
     if ( lfq.and.myid.eq.0 ) then
! write out information about fluctuating charge success
        write(io_output,*) 'Box:   rmflcq for moltyps'
        do im =1,imend
           write(io_output,*) im,(rmflcq(i,im),i=1,nmolty)
        end do
     end if

  end if  ! end if (lratio)

  do imolty = 1, nmolty
     if (lratfix(imolty)) then
! readjust fixed end point data to assure optimum efficiency ***
        counttot = counttot + counthist
        histrat = dble(counthist) / dble(counttot) 

! reset counthist
        counthist = 0
        do j = 1, iring(imolty) 
           do k = 1, iring(imolty)
              if (j.eq.k) goto 150
              histtot = 0
              do bin = 1, maxbin
                 histtot = histtot + hist(j,k,bin)
              end do

              if (histtot.eq.0) goto 150

! normalize and multiply hist by its weighting using above condition
              do bin = 1, maxbin
                 hist(j,k,bin) = hist(j,k,bin)*dble(histrat) / histtot
! add weighed hist to hist
                 probf(j,k,bin) = probf(j,k,bin) + hist(j,k,bin)
! reset hist to zero for next iteration
                 hist(j,k,bin) = 0
              end do
! renormalize new distribution
              histtot = 0
              do bin = 1, maxbin
                 histtot = histtot + probf(j,k,bin)
              end do
              do bin = 1, maxbin
                 probf(j,k,bin) = probf(j,k,bin) / histtot
              end do
150           continue
           end do
        end do
     end if

! KM for MPI
! only do anything for lsolute if monper is not called from readdat (nnn.ne.0)
! calculate energy and write out movie for lsolute
     if (lsolute(imolty).and.nnn.ne.0) then
        do ibox = 1, nbox
           do k = 1, ncmt(ibox,imolty)
              i = parbox(k,ibox,imolty)
! set coords for energy and write out conformations
              if (myid.eq.0) write(11,*) imolty,ibox,nunit(imolty)
              do j = 1, nunit(imolty)
                 rxuion(j,1) = rxu(i,j)                     
                 ryuion(j,1) = ryu(i,j)
                 rzuion(j,1) = rzu(i,j)
                 if (myid.eq.0) write(11,*) ntype(imolty,j) ,rxuion(j,1),ryuion(j,1),rzuion(j,1) ,qqu(j,1)
              end do

              call energy(i,imolty,v,vintra,vinter,vext,velect ,vewald,1,ibox,1,nunit(imolty),.true.,ovrlap ,.false.,vtors,.false.,.false.,.false.)
              if (ovrlap) write(io_output,*)  '*** DISASTER, OVERLAP IN MONPER'

              if (ltailc) then
! tail corrections
                 if (lsolid(ibox).and..not.lrect(ibox)) then
                    vol = cell_vol(ibox) 
                 else
                    vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
                 end if
                 vtail = 0.0E0_dp
                 do jmolty = 1, nmolty
                    rho = ncmt(ibox,jmolty) / vol
                    vtail = vtail + ncmt(ibox,imolty) * coru(imolty,jmolty,rho,ibox)
                 end do
              end if

              call U_bonded(i,imolty,vvib,vbend,vtors)

              solcount(ibox,imolty) = solcount(ibox,imolty) + 1
              avsolinter(ibox,imolty) =  avsolinter(ibox,imolty) + vinter / 2.0E0_dp + vtail
              avsolintra(ibox,imolty) = avsolintra(ibox,imolty)  + vintra 
              avsoltor(ibox,imolty) = avsoltor(ibox,imolty) + vtors
              avsolbend(ibox,imolty) = avsolbend(ibox,imolty)  + vbend
              avsolelc(ibox,imolty) = avsolelc(ibox,imolty)  + velect + vewald
           end do
        end do
     end if
  end do

  if ( lgibbs .or. lnpt ) then
     if ( lratv ) then
! adjust maximum volume displacement ***
        do ibox = 1, nbox
           if (lsolid(ibox) .and. .not. lrect(ibox)) then
              do j = 1,9
                 if ( bnhmat(ibox,j) .gt. 0.5E0_dp ) then
                    ratvol = bshmat(ibox,j) / bnhmat(ibox,j)
                    if (ratvol .eq. 0.0E0_dp) then
                       rmhmat(ibox,j) = rmhmat(ibox,j) * 0.1E0_dp
                    else
                       rmhmat(ibox,j) = rmhmat(ibox,j)*ratvol/tavol
                    end if
                 end if
              end do
           else
              if ( bnvol(ibox) .gt. 0.5E0_dp ) then
                 ratvol = bsvol(ibox) / bnvol(ibox)
                 if ( ratvol .eq. 0.0E0_dp ) then
                    rmvol(ibox) = rmvol(ibox) * 0.1E0_dp
                 else
                    rmvol(ibox) = rmvol(ibox) * ratvol / tavol
                    if (rmvol(ibox).gt.(0.10E0_dp*boxlx(ibox)*boxly(ibox)*boxlz(ibox))) then
                       rmvol(ibox)=0.1E0_dp*(boxlx(ibox)*boxly(ibox)*boxlz(ibox))
                    end if
                 end if
              end if
           end if
        end do

! KM for MPI            
        if (myid.eq.0) then
           do ibox = 1, nbox
              if (lsolid(ibox) .and. .not. lrect(ibox)) then
                 do j = 1,9
                    write(io_output,"(' h-matrix change:  bn =',f8.1, '   bs =',f8.1,'   max.displ. =',e12.5)") bnhmat(ibox,j),bshmat(ibox,j), rmhmat(ibox,j)
                 end do
              else
                 write(io_output,"(' volume change:  bn =',f8.1, '   bs =',f8.1,'   max.displ. =',e12.5)") bnvol(ibox),bsvol(ibox),rmvol(ibox)
              end if
           end do
        end if

        do ibox = 1, nbox
           if (lsolid(ibox) .and. .not. lrect(ibox)) then
              do j = 1,9
                 acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
                 acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
                 bshmat(ibox,j) = 0.0E0_dp
                 bnhmat(ibox,j) = 0.0E0_dp
              end do
           else
              acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
              acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
              bnvol(ibox) = 0.0E0_dp
              bsvol(ibox) = 0.0E0_dp
           end if
        end do
     end if
  end if

! KM for MPI
  if ( lprint.and.myid.eq.0 ) then
! write out runtime information ***
     ntot = nnn + nnstep
     write(io_output,'(i6,i8,e12.4,f10.3,f12.1,15i4)') nnn,ntot, vbox(1),boxlx(1),pres(1) ,(ncmt(1,imolty),imolty=1,nmolty)
     if ( lgibbs ) then
        do ibox = 2, nbox
           write(io_output,'(14x,e12.4,f10.3,f12.1,15i4)')  vbox(ibox),boxlx(ibox),pres(ibox) ,(ncmt(ibox,imolty),imolty=1,nmolty)
        end do
     end if
  end if

! KM for MPI      
  if ( lmv .and.myid.eq.0) then
! write out the movie configurations ***
     write(10,*) nnn
     do ibox = 1, nbox
        write(10,*) (ncmt(ibox,zzz),zzz=1,nmolty)

        if (lsolid(ibox) .and. .not. lrect(ibox)) then
           write(10,*) (hmat(ibox,zzz),zzz=1,9)
        else
           write(10,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
        end if
     end do
     do m = 1, nchain
        imolty = moltyp(m)
        write(10,'(4(1x,i5),3(1x,f16.6))') m, imolty, nunit(imolty), nboxi(m),xcm(m),ycm(m),zcm(m)
        do mm = 1, nunit(imolty)
           write(10,'(4(1x,f14.6),i5)') rxu(m,mm), ryu(m,mm), rzu(m ,mm), qqu(m,mm),ntype( imolty, mm )
        end do
     end do
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
                 write(210+ibox,'(a4,5x,3f15.4)') chemid(ntii), rxu(i,ii), ryu(i,ii), rzu(i,ii)
              end do
           end if
        end do
     end do
  end if

  if ( lrsave .and.myid.eq.0) then
! write out the restart configurations to SAFETY-file ***
     call dump('save-config')
  end if
         
! calculation of block averages ***
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

! specific density
        if ( lpbcz ) then
           temmass = 0.0E0_dp
           do itype = 1, nmolty
              temmass = temmass + masst(itype)*acdens(ibox,itype)
           end do
           dpr = temmass / 0.6022045E0_dp
        else
           temmass = 0.0E0_dp
           do itype = 1, nmolty
              temmass = temmass + acdens(ibox,itype)
           end do
           dpr = temmass
        end if
        call update(nblock,1,ibox,dpr,acmove)

! pressure
        call update(nblock,2,ibox,acpres(ibox),acnp)

! surface tension
        itel = 2+nener+4*nmolty+1
        call update(nblock,itel,ibox,acsurf(ibox),acnp)

! box volume
        itel = 4+nener+4*nmolty
        call update(nblock,itel,ibox,acvolume(ibox),acmove)

! energies
        do j=1,nener
           itel=j+2
           call update(nblock,itel,ibox,acv(j,ibox),acmove)
        end do

! chemical potential
        do itype = 1, nmolty
           itel = (2+nener) + itype
           debroglie = debroglie_factor*sqrt(beta/masst(itype))
           ! determine how many steps it takes to grow the 
           ! molecule not counting the first inserted bead
           igrow = nugrow(itype)
           if (lrigid(itype))then
              call schedule(igrow,itype,steps,1,0,4)
              dpp = acchem(ibox,itype)/( real(nchoi1(itype)*nchoir(itype)*(nchoi(itype)**steps)*nchoih(itype),dp) * (debroglie**3) )
           else
              call schedule(igrow,itype,steps,1,0,2)
              dpp = acchem(ibox,itype)/( real(nchoi1(itype)*(nchoi(itype)**steps)*nchoih(itype),dp) * (debroglie**3) )
           end if
           dvalue = bnchem(ibox,itype)
           call update(nblock,itel,ibox,dpp,dvalue)
        end do

! square end-to-end length
        do itype = 1, nmolty
           itel = (2 + nener) + nmolty + itype
           if ( mnbox(ibox,itype) .eq. 0 ) then
! avoid division by zero in update
              dn = 1.0E0_dp
           else
              dn = dble(mnbox(ibox,itype))
           end if
           call update(nblock,itel,ibox,asetel(ibox,itype),dn)
        end do

! number density
        do itype = 1, nmolty
           itel = (2 + nener) + (2*nmolty) + itype
           dpp = acdens(ibox,itype)
           call update(nblock,itel,ibox,dpp,acmove)
        end do
! mol fraction
        do itype = 1, nmolty
           itel = (2 + nener) + (3*nmolty) + itype
           dpp = molfra(ibox,itype)
           call update(nblock,itel,ibox,dpp,acmove)
        end do

! Enthalpy
        itel = 4+nener+4*nmolty+1
! write(io_output,*) acEnthalpy(ibox)
        call update(nblock,itel,ibox,acEnthalpy(ibox),acnp)
        itel = 4+nener+4*nmolty+2
        call update(nblock,itel,ibox,acEnthalpy1(ibox),acnp)
     end do
     do ibox = 1,nibox-1
        do jbox = ibox+1,nibox 
           do i = 1,nprop1
              call update1(nblock,i,acsolpar(i,ibox,jbox),acnp, ibox,jbox)  
           end do
        end do
     end do
  end if

#ifdef __DEBUG__      
  write(io_output,*) 'end MONPER in ',myid
#endif
  return
end subroutine monper

subroutine averageMaximumDisplacement(ibox,imol)
  use sim_system
  implicit none
  integer,intent(in)::ibox,imol

  if (numberDimensionIsIsotropic(ibox).eq.3) then
     rmtrax(imol,ibox)=(rmtrax(imol,ibox)+rmtray(imol,ibox)+rmtraz(imol,ibox))/3
     rmtray(imol,ibox)=rmtrax(imol,ibox)
     rmtraz(imol,ibox)=rmtrax(imol,ibox)
     rmrotx(imol,ibox)=(rmrotx(imol,ibox)+rmroty(imol,ibox)+rmrotz(imol,ibox))/3
     rmroty(imol,ibox)=rmrotx(imol,ibox)
     rmrotz(imol,ibox)=rmrotx(imol,ibox)
  else if (numberDimensionIsIsotropic(ibox).eq.2) then
     rmtrax(imol,ibox)=(rmtrax(imol,ibox)+rmtray(imol,ibox))/2
     rmtray(imol,ibox)=rmtrax(imol,ibox)
     rmrotx(imol,ibox)=(rmrotx(imol,ibox)+rmroty(imol,ibox))/2
     rmroty(imol,ibox)=rmrotx(imol,ibox)
  end if
end subroutine averageMaximumDisplacement
