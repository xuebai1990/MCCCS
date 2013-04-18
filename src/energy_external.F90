MODULE energy_external
  use var_type,only:dp
  use const_math,only:twopi,fourpi
  use util_math,only:mbessel
  use sim_system,only:sig2ij,epsij,qqu,lcorreg
  implicit none
  private
  save
  public::U_ext,init_energy_external

  real,allocatable::extc12(:),extc3(:),extz0(:)

! Slitpore
  real,parameter::a1=2.460E0_dp,delta=3.40E0_dp& !< a1, delta are in Angstroms
   ,rsol=0.114E0_dp !< rsol [=] 1/A^3
  integer,parameter::ntsubst=190
! Joe Hautman's parameters are read from standard input, to be used with ljoe = .true.
! AT PRESENT no parameters for polymeric surfactants (see LJPSUR)

contains
!> \brief calculates the surface energy of a bend with a featureless
!> graphite surface
  function slitpore(z,ntij)
	integer::ntij
	real::vgs,z
	real::slitpore
	real::sig
!	real::coef1,coef2,coef3,coef4

	sig = sqrt(sig2ij(ntij))

!	coef1 = twopi*rsol*epsij(ntij)*sig2ij(ntij)*delta
!	coef2 = 2.0/5.0*(sig/z)**10
!	coef3 = (sig/z)**4
!	coef4 = sig2ij(ntij)**2/(3*delta*(z+0.61*delta)**3)

	vgs = twopi*rsol*epsij(ntij)*sig2ij(ntij)*delta*(2.0/5.0*(sig/z)**10 		-(sig/z)**4	-(sig2ij(ntij)**2/(3*delta*(z+0.61*delta)**3)))

    slitpore = vgs
  end function slitpore

!> \brief calculates the energy of a bead with a graphite surface
  function exgrph(x,y,z,ntij)
        real::aa,aa2
        real::a1sq
        real::e0,e1
        real::exgrph
        real::fxy
        real::x,y,z
        real::bb,cc,dd
        real::k2,k5
        real::sz2
        real::zzz
        integer::ntij

        exgrph = 0.0E0_dp
        e0 = 0.0E0_dp
        e1 = 0.0E0_dp
        fxy = 0.0E0_dp

! write(81,*) ntij,sqrt(sig2ij(ntij))

        sz2 = sig2ij(ntij)/(z**2)

        aa = twopi * rsol * delta * sig2ij(ntij)

        e0 = aa*epsij(ntij)*((2.0E0_dp/5.0E0_dp)*(sz2**5) - (sz2**2) - (sig2ij(ntij)**2/(3.0E0_dp*delta*(0.61*delta+z)**3)))

! write(82,*) e0,aa,delta,z
        if ( lcorreg ) then
                a1sq = a1**2

                aa2 = (sig2ij(ntij)/a1sq)**3

                bb = aa2*fourpi*epsij(ntij)/sqrt(3.0E0_dp)

! bb = fourpi*epsij(ntij)*sig2ij(ntij)**3/
!     +                 (sqrt(3.0E0_dp)*a1**6)

                cc = aa2/(30.0E0_dp*(twopi/sqrt(3.0E0_dp))**5)

! cc = sig2ij(ntij)**6/
!     +         (30.0E0_dp*a1**6*(twopi/sqrt(3.0E0_dp)**5))

                dd = 2.0E0_dp*(twopi/sqrt(3.0E0_dp))**2

                zzz = fourpi*z/(sqrt(3.0E0_dp)*a1)

                k2 = mbessel(zzz,2.0E0_dp)

                k5 = mbessel(zzz,5.0E0_dp)
! write(84,*) zzz,k2,k5
                e1 = bb*(cc * k5 * (a1/z)**5 - dd * k2 * (a1/z)**2)
! write(82,*) bb,cc,dd,e1,k2,k5
                fxy = -2.0E0_dp*(cos(twopi*(x/a1 + y/sqrt(3.0E0_dp)/a1)) + cos(twopi*(x/a1 - y/sqrt(3.0E0_dp)/a1)) + cos(fourpi*y/sqrt(3.0E0_dp)/a1))

! write(82,'(6g12.5)') x,y,z,fxy,e1,e0
                exgrph = e0 + e1*fxy
        else
! write(83,'(6g12.5)') x,y,z,e0
                exgrph = e0
        end if
  end function exgrph

!> \brief calculates interaction of molecule i with an external field E
!> ********************************************
!> \par Units
!> E in V/A, q in e, rz in A \n
!> E*q*rz = V*e \n
!> 1 V*e = 11600 K
!> ********************************************
!> \author added 06/24/07 by KM
  function v_elect_field(i, j, rzfield,E)
      real::v_elect_field,rzfield, E
      integer::i,j
      v_elect_field = -E*rzfield*qqu(i,j)
! write(io_output,*) 'E ', E, ' exfield ', exfield
      return
  end function v_elect_field

  function U_ext(ibox,i,j,ntj)
    use const_phys,only:eXV_to_K
    use sim_system,only:rxu,ryu,rzu,nntype,lelect_field,Elect_field,ljoe,lslit,lgraphite,lsami,lmuir,lexzeo,io_output,nchain,boxlz
    use energy_sami,only:exsami,exmuir
    use zeolite
    real::U_ext
    integer,intent(in)::ibox,i,j,ntj
    integer::ntij
    real::dzui,dz3,dz12,vtmp

    U_ext=0.0E0_dp

! **********************************************************************
! calculation of interaction energy with external electric field ***
! added 06/24/07 by KM
! **********************************************************************
    if (lelect_field) then
       U_ext = U_ext + v_elect_field(i,j,rzu(i,j),Elect_field(ibox)) * eXV_to_K
    end if

    if (ibox.ne.1) return

    if ( ljoe ) then
       if ( extc12(ntj) .gt. 0.1E0_dp ) then
          dzui = rzu(i,j) - extz0(ntj)
          dz3  = dzui * dzui * dzui
          dz12 = dz3**4
          U_ext = U_ext +  (extc12(ntj)/dz12) - (extc3(ntj)/dz3)
       end if
    end if

    if (lslit) then
       ! Carbon slitpore
       ntij = (ntj-1)*nntype + ntsubst
       ! calculate interaction with surface at the bottom of the box
       U_ext = U_ext + slitpore(rzu(i,j),ntij)
       ! calculate interaction with the surface at the top of the box
       dzui = boxlz(ibox)-rzu(i,j)
       U_ext = U_ext +slitpore(dzui,ntij)
    end if

    if( lgraphite ) then
       ntij = (ntj-1)*nntype + ntsubst
       U_ext = U_ext + exgrph(rxu(i,j),ryu(i,j),rzu(i,j),ntij)
    end if

    if ( lsami )  U_ext = U_ext + exsami(rzu(i,j),ntj)
    if ( lmuir )  U_ext = U_ext + exmuir(rzu(i,j),ntj)
    if ( lexzeo ) then
       vtmp=exzeo(rxu(i,j),ryu(i,j),rzu(i,j),ntj,ignoreTable=.false.)
       if (i.le.nchain.and.abs(vtmp).gt.1E5_dp) write(io_output,*) i,j,rxu(i,j),ryu(i,j),rzu(i,j),vtmp
       U_ext = U_ext + vtmp
    end if

  end function U_ext

  subroutine init_energy_external(io_input,lprint)
    use util_runtime,only:err_exit
    use sim_system,only:nntype,nbox,ljoe,lmuir,io_output,lelect_field,Elect_field
    use energy_sami
    INTEGER,INTENT(IN)::io_input
    LOGICAL,INTENT(IN)::lprint
    integer::jerr
    namelist /E_field/ Elect_field

    if (lelect_field) then
       !> read namelist E_field
       Elect_field=0.0_dp

       rewind(io_input)
       read(UNIT=io_input,NML=E_field,iostat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'reading namelist: E_field',jerr)

       if (lprint) then
          write(io_output,*) 'Electric field in z direction: ',Elect_field(1:nbox),' [V/A]'
       end if
    end if

    if ( ljoe ) then
       allocate(extc12(1:nntype),extc3(1:nntype),extz0(1:nntype),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_pairwise: nonbond allocation failed',jerr)

       ! STANDARD METHYL GROUP
       extc12(4) = 3.41E7_dp
       extc3(4)  = 20800.0E0_dp
       extz0(4)  = 0.86E0_dp

       ! STANDARD METHYLENE GROUP
       extc12(5) = 2.80E7_dp
       extc3(5)  = 17100.0E0_dp
       extz0(5)  = 0.86E0_dp

       ! Methane
       extc12(3) = 3.41E7_dp
       extc3(3)  = 20800.0E0_dp
       extz0(3)  = 0.86E0_dp

       ! Martin's methyl (CH3)
       extc12(18) = 3.41E7_dp
       extc3(18)  = 20800.0E0_dp
       extz0(18)  = 0.86E0_dp
    end if

    ! calculate constants for lmuir external potential ***
    if ( lmuir ) then
       sigpri = 0.715E0_dp * sqrt( 3.8E0_dp * 3.93E0_dp )
       c9ch2 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*47.0E0_dp) ) * sigpri**9
       c3ch2 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*47.0E0_dp) ) * sigpri**3
       c9ch3 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*114.0E0_dp) ) * sigpri**9
       c3ch3 = 4.0E0_dp * ( 1.43E0_dp * sqrt(80.0E0_dp*114.0E0_dp) ) * sigpri**3
       zprmin = ( 3.0E0_dp**(1/6.0E0_dp) ) * sigpri
       v2prmin = c9ch2 / zprmin**9 - c3ch2 / zprmin**3
       v3prmin = c9ch3 / zprmin**9 - c3ch3 / zprmin**3
       betac2 = beta1 - v2prmin
       betac3 = beta1 - v3prmin
       if (lprint) then
          write(io_output,*) 'external potential for Langmuir monolayers used'
          write(io_output,*) 'zprmin',zprmin
          write(io_output,*) 'v2prmin',v2prmin,'v3prmin',v3prmin
       end if
    end if
  end subroutine init_energy_external
end MODULE energy_external
