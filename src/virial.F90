!> \brief Computes the 2nd virial coefficient
!>
!> B2(T) = -2Pi Int 0toInf [ Exp[-beta*u(r)] -1] r^2 dr \n
!> Using the trapazoid method of numerical integration give \n
!> B2(T) = -2*Pi*stepvir* Sum(i=2,n-1)[ <Exp[-beta*u(r)]-1> ri^2 \n
!> 1/2( <Exp[-beta*u(r1)]-1> r1 + <Exp[-beta*u(rn)]-1> rn
!> \author Marcus Martin 1-15-97
!> revised by Robert Hembree 2-10-2016
!> \brief Computes the 2nd virial coefficient
!>
!> Uses a trapezoidal integration to integrate the virial coefficient of a
!> specific conformation of two molecules. The contribution to B2 is then weighted
!> according to the configurational intramolecular energy.
!> I.e. we compute the boltzmann weighted (via intramolecular energy) B2 virail
!> coefficient
!> subroutine parameters:
!>      binvirnumerator: The numerator in the average B2 calculation
!>      binvirdenominator: the denominator of the average B2 calculation.
subroutine virial(binvirnumerator,binvirdenominator)
  use const_math,only:onepi,twopi
  use util_runtime,only:err_exit
  use sim_system
  use energy_pairwise,only:U2,type_2body,energy
  use energy_intramolecular,only:U_bonded
  implicit none
      integer::i, imolty, ii, j, jmolty, jj, ntii, ntjj , ntij,nnn,ip,itemp,iii
      real::xdiff,ydiff ,zdiff,dvircm
      real::binvir
      real::mass_t, factor,corr,vprev,deri_u
      real::binvirnumerator,binvirdenominator
      real::Uintramol1,Uintramol2,vvbend,vvib,vvtors,uIntermolecular
      real::xcmsep,ycmsep,zcmsep,deviation,their_distance
      real::mayerterm,boltzfact,smallexpfact,fullexpfact,integralvalue
      real::storefirstval
      real::vEnergy(nEnergy),vmol1(nEnergy),vmol2(nEnergy)
      logical::olp=.false., firstval=.true.
      integer::iunit
#ifdef __DEBUG__
      write(io_output,*) 'start VIRIAL in ',myid
#endif

      firstval=.true.
      ! start by comptuing the intramolecular energy.
      ! this is used to properly weight the integral
      call U_bonded(1,moltyp(1),vvib,vvbend,vvtors)
      Uintramol1=vvib+vvbend+vvtors
      call U_bonded(2,moltyp(2),vvib,vvbend,vvtors)
      Uintramol2=vvib+vvbend+vvtors
      ! We now need to include the nonbonded portion of the intramolecular
      ! energies. It is very important that we don't have the molecules in the
      ! same box.
      ! This is done to distinguish the inter- and intramolecular portions of
      ! the electric and ewald interactions. Being in individual boxes we won't
      ! see them interact with each other.
      ! also note I will later want to subtract the intramolecular portion of
      ! these energies from the energy calculated when the molecules are in the
      ! same box.

      call energy(1,moltyp(1),vmol1,1,nboxi(1),1,nunit(1),.true.,olp,.false.,.false.,.false.,.false.)
      call energy(2,moltyp(2),vmol2,1,nboxi(2),1,nunit(2),.true.,olp,.false.,.false.,.false.,.false.)
      Uintramol1 = Uintramol1+vmol1(ivElect)+vmol1(ivIntraLJ)
      Uintramol2 = Uintramol2+vmol2(ivElect)+vmol2(ivIntraLJ)

      if ( nboxi(1) .eq. nboxi(2) ) then
         !write(io_output,*) 'particles found in same box'
         call err_exit(__FILE__,__LINE__,'',myid+1)
      end if
      olp=.false.
      ! calculate the differences in their COM in each direction
      xdiff = xcm(2) - xcm(1)
      ydiff = ycm(2) - ycm(1)
      zdiff = zcm(2) - zcm(1)
      imolty = moltyp(1)
      jmolty = moltyp(2)
      iunit = nunit(imolty)
      ! old code retained for future development
      mass_t = 0.0E0_dp
      do ii = 1, iunit
         mass_t = mass_t + mass(ntype(imolty,ii))
      end do
      mass_t = mass_t/1000E0_dp
      factor = -(6.6260755E-34_dp)**2*6.0221367E23_dp*1E20_dp /  (24.0E0_dp*onepi*mass_t*1.380658E-23_dp*twopi)
      deviation = starvir
      integralvalue=0.0E0_dp
      do while (deviation > rmin)
        !set the "trial" location of the chain2
        do jj=1,nunit(jmolty)
            ! center them ontop of each other in the y- and z- directions.
            ! Slowly translate molecule 2 down the x axis until we reach begin
            ! to overlap. If we overlap then the energy will be infinity.
            rxuion(jj,2) = rxu(2,jj)-xdiff+deviation
            ryuion(jj,2) = ryu(2,jj)-ydiff
            rzuion(jj,2) = rzu(2,jj)-zdiff
        end do
        ! compute the energy of the system as if the second molecule was in hte
        ! same box as the first molecule at a distance of deviation away.
        call energy(2,jmolty,vEnergy,2,nboxi(1),1,nunit(jmolty),.false.,olp,.false.,.false.,.false.,.false.)
        uIntermolecular = vEnergy(ivInterLJ)+vEnergy(ivTail)+vEnergy(ivElect)+vEnergy(iv3body)+vEnergy(ivFlucq)
        ! subtract out the intramolecular contrubutions to the electric
        ! energies that were calculated before leaving only the
        ! intermolecular contributions. Only the intramolecular components
        ! for mol2 were calculated here so only those need to be removed.
        uIntermolecular = uIntermolecular-vmol2(ivElect)
        ! lets do the soft cut stuff
        smallexpfact = -uIntermolecular/virtemp
        fullexpfact = -(Uintramol1+Uintramol2)/virtemp
        if(smallexpfact < -1.0*softcut) then
            mayerterm=0.0
        else if (smallexpfact > softcut.or.olp) then
            ! essentially these are states where the energy appears to be
            ! infinite.
            mayerterm = 0.0E0_dp
        else
            mayerterm = exp(-uIntermolecular/virtemp)
        end if
        ! the 2* accounts for the double counting of terms in the trapezoidal
        ! method. These are later removed by dividing by 2.
        integralvalue=integralvalue+2*(1-mayerterm)*deviation**2
        if(firstval) then
          !stores the first term because the first and last terms have half
          !value of the rest of the terms in the trapezoidal integration scheme.
          firstval = .false.
          storefirstval = integralvalue/2.0
        end if
        deviation = deviation-stepvir
      end do
      ! handles the end points of the trapezoidal rule which have half of the
      ! value of all of the other points.
      integralvalue = integralvalue-storefirstval
      integralvalue = integralvalue-(1.0-mayerterm)*(deviation+stepvir)**2
      ! now multiply by the constant terms in the integral... 2pi*dr and divide
      ! by two for the trapezoidal rule
      integralvalue = integralvalue*twopi*stepvir/2.0

      if(fullexpfact < -1.0*softcut) then
          boltzfact=0.0
      else if (fullexpfact > softcut) then
          call err_exit(__FILE__,__LINE__,'',myid+1)
      else
          boltzfact = exp(-(fullexpfact)/virtemp)
      end if
      fullexpfact = -(Uintramol1+Uintramol2)/virtemp
      binvirdenominator=binvirdenominator+boltzfact
      binvirnumerator = binvirnumerator+boltzfact*integralvalue
#ifdef __DEBUG__
      write(io_output,*) "B2: ", binvirnumerator/binvirdenominator
      write(io_output,*) 'end VIRIAL in ',myid
#endif

end subroutine virial
