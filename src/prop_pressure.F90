MODULE prop_pressure
  use var_type,only:dp
  use const_phys,only:qqfact
  use util_math,only:erfunc
  use util_mp,only:mp_sum
  use sim_system
  use sim_cell
  use energy_kspace,only:recippress,calp
  use energy_pairwise
  implicit none
  private
  save
  public::pressure

contains
!> \brief Calculates the pressure for a configuration, unit [kPa]
  subroutine pressure(press,surf,ibox)
    use const_math,only:sqrtpi
    use const_phys,only:k_B
    real,intent(out)::press,surf
    integer,intent(in)::ibox

    real::rbcut,rcutsq,calpi,calpisq,pxx,pyy,pzz,xcmi,ycmi,zcmi,rcmi,fxcmi,fycmi,fzcmi,rxuij,ryuij,rzuij,rijsq,rcm,rcmsq,rxui,ryui,rzui,rij,fij,rs2,rs1,rs6,rs7,sr1,sr7,sr2,srij,sigma2,epsilon2,sr6,qave,repress,rpxx,rpyy,rpzz,rpxy,rpyx,rpxz,rpzx,rpyz,rpzy,vol,volsq,rhosq,pwell
    integer::i,imolty,j,jmolty,ii,jj,ntii,ntjj,ntij,iii,jjj,k
    logical::lqimol,lexplt,lij2,lqjmol,lcoulo(numax,numax)
! --------------------------------------------------------------------
    if (lsolid(ibox).and.(.not.lrect(ibox))) then
       vol = cell_vol(ibox)
    else
       vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
    end if

    if (lideal(ibox)) then
       press = k_B*1.0E27_dp * (nchbox(ibox)+ghost_particles(ibox))/beta/vol
       surf = 0.0_dp
       return
    end if

    if (lpbc) call setpbc(ibox)

    rbcut  = rcut(ibox)
    rcutsq = rbcut*rbcut
    calpi=calp(ibox)
    calpisq=calpi*calpi

    press = 0.0E0_dp
    pxx = 0.0E0_dp
    pyy = 0.0E0_dp
    pzz = 0.0E0_dp
    pips = 0.0E0_dp

! *******************************
! INTERCHAIN INTERACTIONS ***
! *******************************
    ! if(LSOLPAR.and.(ibox.eq.2)) then
    !    press= 1.380662E4_dp * ( ( nchbox(ibox) / beta) - ( press/3.0E0_dp ) ) / ( boxlx(ibox)*boxly(ibox)*boxlz(ibox) )
    !    surf = 0.0E0_dp
    !    return
    ! end if

    ! loop over all chains i
    ! RP added for MPI
    do i = myid+1, nchain - 1,numprocs
       ! check if i is in relevant box ###
       if ( nboxi(i) .eq. ibox ) then
          imolty = moltyp(i)
          lqimol = lelect(imolty)
          if (nugrow(imolty).eq.nunit(imolty)) then
             lexplt = .false.
          else
             lexplt = .true.
          end if
          xcmi = xcm(i)
          ycmi = ycm(i)
          zcmi = zcm(i)
          if (lcutcm) then
             rcmi = rcmu(i)
          else
             lij2 = .true.
          end if

          ! loop over all chains j with j>i
          do j = i + 1, nchain
             ! check for simulation box ###
             if ( nboxi(j) .eq. ibox ) then
                jmolty = moltyp(j)
                lqjmol = lelect(jmolty)
                fxcmi = 0.0E0_dp
                fycmi = 0.0E0_dp
                fzcmi = 0.0E0_dp
                if ( lcutcm ) then
                   ! check if ctrmas within rcmsq
                   rxuij = xcmi - xcm(j)
                   ryuij = ycmi - ycm(j)
                   rzuij = zcmi - zcm(j)
                   ! minimum image the ctrmas pair separations
                   if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                   rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                   rcm = rbcut + rcmi + rcmu(j)
                   rcmsq = rcm*rcm
                   if (rijsq .le. rcmsq) then
                      lij2 = .true.
                   else if (lqimol .and. lqjmol .and. lchgall) then
                      lij2 = .false.
                   else
                      cycle
                   end if
                end if

                ! loop over all beads ii of chain i
                do ii = 1, nunit(imolty)
                   ntii = ntype(imolty,ii)
                   rxui = rxu(i,ii)
                   ryui = ryu(i,ii)
                   rzui = rzu(i,ii)

                   ! loop over all beads jj of chain j
                   bead2: do jj = 1, nunit(jmolty)
                      ! check exclusion table
                      if (lexclu(imolty,ii,jmolty,jj)) cycle bead2

                      ntjj = ntype(jmolty,jj)
                      if ( lij2 ) then
                         if ((.not.(lij(ntii).and.lij(ntjj))) .and.(.not.(lqchg(ntii).and.lqchg(ntjj)))) cycle bead2
                      else
                         if (.not.(lqchg(ntii).and.lqchg(ntjj))) cycle bead2
                      end if

                      ntij=type_2body(ntii,ntjj)

                      rxuij = rxui - rxu(j,jj)
                      ryuij = ryui - ryu(j,jj)
                      rzuij = rzui - rzu(j,jj)
                      ! minimum image the pair separations ***
                      if ( lpbc ) call mimage(rxuij,ryuij,rzuij,ibox)
                      rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                      rij=sqrt(rijsq)

                      ! compute whether the charged groups will interact & fij
                      fij = 0.0E0_dp
                      !> \todo when lgaro=.true.
                      if (lqimol.and.lqjmol.and.lqchg(ntii).and.lqchg(ntjj) ) then
                         if (.not.lewald) then
                            if (.not.lchgall) then
                               iii = leaderq(imolty,ii)
                               jjj = leaderq(jmolty,jj)
                               if (iii.eq.ii.and.jjj.eq.jj) then
                                  ! set up the charge-interaction table
                                  if (rijsq.lt.rcutsq) then
                                     lcoulo(iii,jjj) = .true.
                                  else
                                     lcoulo(iii,jjj) = .false.
                                  end if
                               end if
                            end if
                            !> \todo tabulated potential
                            if (lchgall.or.lcoulo(iii,jjj) ) then
                               fij = -qqfact*qqu(i,ii)*qqu(j,jj)/rijsq/rij
                            end if
                         else if (lchgall.or.rijsq.lt.rcutsq) then
                               fij = -qqfact*qqu(i,ii)*qqu(j,jj)/rijsq*(2.0_dp*calpi*exp(-calpisq*rijsq)/sqrtpi+erfunc(calpi*rij)/rij)
                         end if
                      end if

                      if (rijsq.lt.rcutsq.or.lijall) then
                         !> \todo when lsami, lmuir, lpsur, lgaro is .true.
                         if ( lexpsix ) then
                            fij = fij + (-6.0_dp*aexsix(ntij)/(rijsq*rijsq*rijsq)+bexsix(ntij)*cexsix(ntij)*rij*exp(cexsix(ntij)*rij))/rijsq
                         else if (lmmff) then
                            rs2 = rijsq/sigisq(ntij)
                            rs1 = sqrt(rs2)
                            rs6 = rs2*rs2*rs2
                            rs7 = rs1*rs6
                            sr1 = 1.07_dp/(rs1+0.07_dp)
                            sr7 = sr1**7.0_dp
                            fij = fij - 7.0_dp*epsimmff(ntij)*rs1*sr7*(sr1*(1.12_dp/(rs7+0.12_dp)-2.0_dp)/1.07_dp+1.12_dp*rs6/((rs7+0.12_dp)*(rs7+0.12_dp)))/rijsq
                         else if (lninesix) then
                            fij = fij + 72.0_dp*epsnx(ntij)/(rij*rzero(ntij))*(rzero(ntij)/rij)**7*(1.0_dp-(rzero(ntij)/rij)**3)
                         else if (lgenlj) then
                            sr2 = sig2ij(ntij)/rijsq
                            srij = sqrt(sr2)
                            if (rij.le.rij*srij*2.0_dp**(2.0_dp/n0)) then
                               fij = fij - 4.0E0_dp*epsij(ntij) * (n0*(srij)**n0-n0/2.0_dp*srij**(n0/2.0E0_dp))/rijsq
                            else
                               fij = fij - epsij(ntij)*(2.0_dp*n1*srij**(2.0E0_dp*n1)*2.0_dp**(4.0_dp*n1/n0)-n1*srij**n1*2.0_dp**(2.0_dp*n1/n0+1.0E0_dp))/rijsq
                            end if
                         else
                            if (lexpand(imolty).and.lexpand(jmolty)) then
                               sigma2=(sigma_f(imolty,ii)+sigma_f(jmolty,jj))/2.0_dp
                               sr2 = sigma2*sigma2/rijsq
                               epsilon2=sqrt(epsilon_f(imolty,ii)*epsilon_f(jmolty,jj))
                            else if (lexpand(imolty)) then
                               sigma2=(sigma_f(imolty,ii)+sigi(ntjj))/2.0_dp
                               sr2 = sigma2*sigma2/rijsq
                               epsilon2=sqrt(epsilon_f(imolty,ii)*epsi(ntjj))
                            else if (lexpand(jmolty)) then
                               sigma2=(sigma_f(jmolty,jj)+sigi(ntii))/2.0_dp
                               sr2 = sigma2*sigma2/rijsq
                               epsilon2=sqrt(epsilon_f(jmolty,jj)*epsi(ntii))
                            else
                               sr2 = sig2ij(ntij) / rijsq
                               epsilon2 = epsij(ntij)
                            end if

                            if (lfepsi) then
                               sr6 = 1.0_dp/rijsq
                               sr6 = sr6**3
                               if ((.not.lqchg(ntii)).and.(.not.lqchg(ntjj))) then
                                  if (nunit(imolty).eq.4) then
                                     !> \bug TIP-4P structure (temperary use?)
                                     qave = (qqu(i,4)+qqu(j,4))/2.0_dp
                                  else
                                     qave = (qqu(i,4)+qqu(i,5)+qqu(j,4)+qqu(j,5))*0.85_dp
                                  end if
                               else
                                  qave=(qqu(i,ii)+qqu(j,jj))/2.0_dp
                               end if

                               fij = fij + 48.0_dp*epsilon2*sr6*(-sr6*(aslope*(qave-a0)*(qave-a0)+ashift)+0.5_dp*(bslope*(qave-b0)*(qave-b0)+ bshift))/rijsq
                            else
                               sr6 = sr2**3
                               fij = fij + 48.0_dp*sr6*(-sr6+0.5_dp)*epsilon2/rijsq
                            end if
                         end if
                      end if
                      fxcmi = fxcmi + fij * rxuij
                      fycmi = fycmi + fij * ryuij
                      fzcmi = fzcmi + fij * rzuij
                   end do bead2
                end do

                ! calculate distance between c-o-m ---
                rxuij = xcmi - xcm(j)
                ryuij = ycmi - ycm(j)
                rzuij = zcmi - zcm(j)
                ! minimum image the pair separations ***
                if (lpbc) call mimage(rxuij,ryuij,rzuij,ibox)

                press = press +  fxcmi*rxuij + fycmi*ryuij + fzcmi*rzuij

                ! for surface tension
                ! this is correct for the coulombic part and for LJ.  Note sign difference!
                pxx = pxx - fxcmi*rxuij
                pyy = pyy - fycmi*ryuij
                pzz = pzz - fzcmi*rzuij
                pips(1,2) = pips(1,2) - rxuij*fycmi
                pips(1,3) = pips(1,3) - rxuij*fzcmi
                pips(2,1) = pips(2,1) - ryuij*fxcmi
                pips(2,3) = pips(2,3) - ryuij*fzcmi
                pips(3,1) = pips(3,1) - rzuij*fxcmi
                pips(3,2) = pips(3,2) - rzuij*fycmi
             end if
          end do
       end if
    end do

    call mp_sum(press,1,groupid)
    call mp_sum(pxx,1,groupid)
    call mp_sum(pyy,1,groupid)
    call mp_sum(pzz,1,groupid)
    call mp_sum(pips,size(pips),groupid)
! ################################################################

    if (lewald) then
       ! Compute the reciprocal space contribution
       ! by using the thermodynamic definition
       call recippress(ibox,repress,rpxx,rpyy,rpzz,rpxy,rpyx,rpxz,rpzx,rpyz,rpzy)
       press = press - repress
       pxx = pxx + rpxx
       pyy = pyy + rpyy
       pzz = pzz + rpzz
       pips(1,2) = pips(1,2) + qqfact*rpxy
       pips(1,3) = pips(1,3) + qqfact*rpxz
       pips(2,1) = pips(2,1) + qqfact*rpyx
       pips(2,3) = pips(2,3) + qqfact*rpyz
       pips(3,1) = pips(3,1) + qqfact*rpzx
       pips(3,2) = pips(3,2) + qqfact*rpzy
    end if

    pips(1,1) = pxx
    pips(2,2) = pyy
    pips(3,3) = pzz

    ! Compute the Gaussian well contribution to the pressure
    pwell = 0.0E0_dp
    pwellips = 0.0E0_dp

    ! KM for MPI - comment:
    ! this could likely be parallelized in the future
    do i = 1,nchain
       imolty = moltyp(i)
       if (lwell(imolty)) then
          rxui = xcm(i)
          ryui = ycm(i)
          rzui = zcm(i)
          do j = 1, nwell(imolty)*nunit(imolty)
             k = j - int(j/nunit(imolty))*nunit(imolty)
             if (k.eq.0) k = nunit(imolty)
             rxuij = rxui-rxwell(j,imolty)
             ryuij = ryui-rywell(j,imolty)
             rzuij = rzui-rzwell(j,imolty)
             call mimage (rxuij,ryuij,rzuij,ibox)
             rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
             rcm = rcut(ibox)+rcmu(i)
             rcmsq = rcm*rcm
             if (rijsq.lt.rcmsq) then
                do ii = 1, nunit(imolty)
                   if (awell(ii,k,imolty).lt.1.0E-6_dp) cycle
                   rxui = rxu(i,ii)
                   ryui = ryu(i,ii)
                   rzui = rzu(i,ii)
                   rxuij = rxui-rxwell(j,imolty)
                   ryuij = ryui-rywell(j,imolty)
                   rzuij = rzui-rzwell(j,imolty)
                   call mimage (rxuij,ryuij,rzuij,ibox)
                   rijsq = rxuij*rxuij+ryuij*ryuij+rzuij*rzuij
                   fij = 2.0E0_dp*awell(ii,k,imolty)*bwell* exp(-bwell*rijsq)
                   pwell = pwell+fij*rijsq
                   pwellips(1,1) = pwellips(1,1)+fij*rxuij*rxuij
                   pwellips(2,2) = pwellips(2,2)+fij*ryuij*ryuij
                   pwellips(3,3) = pwellips(3,3)+fij*rzuij*rzuij
                   pwellips(1,2) = pwellips(1,2)+fij*rxuij*ryuij
                   pwellips(1,3) = pwellips(1,3)+fij*rxuij*rzuij
                   pwellips(2,3) = pwellips(2,3)+fij*ryuij*rzuij
                end do
             end if
          end do
       end if
    end do
    pwellips(2,1) = pwellips(1,2)
    pwellips(3,1) = pwellips(1,3)
    pwellips(3,2) = pwellips(2,3)

    pipsw = -(press/3.0E0_dp)/vol
    pwellipsw = -(pwell/3.0E0_dp)/vol

    do i = 1, 3
       do j = 1, 3
          pips(i,j) = pips(i,j)/vol
          pwellips(i,j) = -pwellips(i,j)/vol
       end do
    end do

    if (lstagea) then
       press = (1.0E0_dp-lambdais*(1.0E0_dp-etais))*press
    else if (lstageb) then
       press = etais*press+lambdais*pwell
    else if (lstagec) then
       press = (etais+(1.0E0_dp-etais)*lambdais)*press +(1.0E0_dp-lambdais)*pwell
    end if

    press = ((nchbox(ibox)+ghost_particles(ibox))/beta - press/3.0_dp)/vol
    surf = pzz - 0.5_dp*(pxx + pyy)
    ! divide by surface area and convert from K to put surf in mN/m
    surf = k_B*1.0E23_dp*surf/(2.0_dp*boxlx(ibox)*boxly(ibox))

    ! tail correction and impulsive force contribution to pressure
    if (ltailc.or.(.not.lshift.and.numberDimensionIsIsotropic(ibox).eq.3)) then
       ! add tail corrections for the Lennard-Jones energy
       ! Not adding tail correction for the ghost particles
       ! as they are ideal (no interaction) Neeraj.
       volsq = vol*vol
       do imolty=1, nmolty
          do jmolty=1, nmolty
             rhosq = ncmt(ibox,imolty)*ncmt(ibox,jmolty)/volsq
             press=press + corp(imolty,jmolty,rhosq,ibox)
          end do
       end do
    end if

    press = k_B*1.0E27_dp * press

    return
  end subroutine pressure

!> \brief Impulsive force and tail corrections to pressure
!>
!> @todo Impulsive force correction only done for LJ 12-6 potential
  function corp(imolty,jmolty,rhosq,ibox)
    use const_math,only:onepi
    real::corp,rci3,rhosq, epsilon2, sigma2
    real::rci1
    integer::imolty,jmolty,ii,jj, ntii, ntjj, ntij ,ibox

    corp = 0.0E0_dp
    do ii = 1, nunit(imolty)
       ntii = ntype(imolty,ii)
       do jj = 1, nunit(jmolty)
          ntjj = ntype(jmolty,jj)
          ntij=type_2body(ntii,ntjj)
          if (lexpsix .and. ltailc) then
             corp = corp + consp(ntij)
          else if (lmmff .and. ltailc) then
             corp = corp+((-2.0E0_dp)/3.0E0_dp)*onepi*epsimmff(ntij)*sigimmff(ntij)**3.0E0_dp*corp_cons(ntij)
          else if (lninesix .and. ltailc) then
             corp = corp + 16.0E0_dp * onepi * epsnx(ntij) * rzero(ntij)**3 * (0.5E0_dp*(rzero(ntij)/rcut(ibox))**6 - (rzero(ntij)/rcut(ibox))**3)
          else if (lgenlj .and. ltailc) then
             rci3 = sig2ij(ntij)**(3.0E0_dp/2.0E0_dp) / rcut(ibox)**3
             rci1 = rci3 **(1.0E0_dp/3.0E0_dp)
             if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                sigma2 = (sigma_f(imolty,ii)+sigma_f(jmolty,jj))**2/4.0E0_dp
                epsilon2 = sqrt(epsilon_f(imolty,ii) *epsilon_f(jmolty,jj))
             else if ( lexpand(imolty) ) then
                sigma2 = (sigma_f(imolty,ii)+sigi(ntjj))**2/4.0E0_dp
                epsilon2 = sqrt(epsilon_f(imolty,ii)*epsi(ntjj))
             else if ( lexpand(jmolty) ) then
                sigma2 = (sigma_f(jmolty,jj)+sigi(ntii))**2/4.0E0_dp
                epsilon2 = sqrt(epsilon_f(jmolty,jj)*epsi(ntii))
             else
                sigma2 = sig2ij(ntij)
                epsilon2 = epsij(ntij)
             end if
             corp = corp + (2.0E0_dp/3.0E0_dp) * onepi * epsilon2 * sigma2 ** (1.50E0_dp) * n1&
              * ( (( (2.0E0_dp**((4.0E0_dp*n1/n0)+1.0E0_dp))/(2.0E0_dp*n1-3.0E0_dp)) * rci1 **(2.0E0_dp*n1-3.0E0_dp) )&
              - ( ((2.0E0_dp**((2.0E0_dp*n1/n0)+1.0E0_dp))/(n1-3.0E0_dp))* rci1 **(n1-3.0E0_dp) )  )
          else if (.not.(lexpsix.or.lmmff.or.lninesix.or.lgenlj)) then
             if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                sigma2 = (sigma_f(imolty,ii)+sigma_f(jmolty,jj))**2/4.0E0_dp
                epsilon2 = sqrt(epsilon_f(imolty,ii)* epsilon_f(jmolty,jj))
             else if ( lexpand(imolty) ) then
                sigma2 = (sigma_f(imolty,ii)+sigi(ntjj))**2/4.0E0_dp
                epsilon2 = sqrt(epsilon_f(imolty,ii)*epsi(ntjj))
             else if ( lexpand(jmolty) ) then
                sigma2 = (sigma_f(jmolty,jj)+sigi(ntii))**2/4.0E0_dp
                epsilon2 = sqrt(epsilon_f(jmolty,jj)*epsi(ntii))
             else
                sigma2 = sig2ij(ntij)
                epsilon2 = epsij(ntij)
             end if
             rci3 = (sqrt(sigma2)/rcut(ibox))**3
             if (ltailc) then ! both impulsive force and tail corrections
                corp = corp+8.0E0_dp*onepi*epsilon2*sigma2**(1.5E0_dp)*(rci3*rci3*rci3*7.0E0_dp/9.0E0_dp-rci3)
             else ! only impulsive force corrections
                corp=corp+(8.0E0_dp/3.0E0_dp)*onepi*epsilon2*sigma2**(1.5E0_dp)*(rci3*rci3*rci3 - rci3)
             end if
          end if
       end do
    end do

    corp=corp*rhosq

    return
  end function corp
end MODULE prop_pressure
