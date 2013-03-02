MODULE energy_intramolecular
  use var_type,only:default_string_length
  use const_math,only:onepi,twopi,raddeg
  use util_runtime,only:err_exit
  use util_memory,only:reallocate
  use util_string,only:uppercase,integer_to_string
  use util_files,only:readLine
  use util_search,only:LookupTable,initiateTable,addToTable
  use sim_system,only:io_output,brvib,brvibk,brben,brbenk,myid,L_spline,L_linear
  implicit none
  private
  save
  public::init_energy_bonded,vtorso,lininter_vib,lininter_bend,U_torsion,U_bonded,allocate_energy_bonded,bonds,angles,dihedrals

  integer,parameter::torsion_nParameter(8)=(/4,10,3,2,5,10,4,5/)
  integer,allocatable::vib_type(:),ben_type(:),torsion_type(:) !< type 0: dummy torsion type for setting up interaction table
                                       !< type 1: OPLS potential (three terms), angle in protein convention (trans is 180 deg)
                                       !< type 2: Ryckaert-Bellemans potential, angle in polymer convention (trans is 0 deg)
                                       !< type 3: periodic type, angle in protein convention (trans is 180 deg)
                                       !< type 4: harmonic type, angle in polymer convention (trans is 0 deg)
                                       !< type 5: OPLS potential (four terms), angle in protein convention (trans is 180 deg)
                                       !< type 6: nine-term Fourier cosine series, angle in protein convention (trans is 180 deg)
                                       !< type 7: Ryckaert-Bellemans potential (three terms), angle in polymer convention (trans is 0 deg)
                                       !< type 8: Ryckaert-Bellemans potential (four terms), angle in polymer convention (trans is 0 deg)
  real,allocatable::vtt(:,:)

  integer,allocatable::vibsplits(:),bendsplits(:),splpnts(:)
  real,allocatable::vib(:,:),tabvib(:,:),bend(:,:),tabbend(:,:),deg(:,:),tabtorso(:,:),torderiv2(:,:)
  integer::ntabvib,ntabbend,nttor

  type(LookupTable)::bonds,angles,dihedrals

  real,allocatable::rxvec(:,:),ryvec(:,:),rzvec(:,:),distij2(:,:),distanceij(:,:)
contains
  subroutine init_energy_bonded(io_ff)
    integer,intent(in)::io_ff
    integer,parameter::initial_size=20
    integer::i,n,jerr
    character(LEN=default_string_length)::line_in

    call initiateTable(bonds,initial_size)
    call initiateTable(angles,initial_size)
    call initiateTable(dihedrals,initial_size)

    call init_tabulated_potential_bonded()

    ! Looking for section BONDS
    REWIND(io_ff)
    CYCLE_READ_BONDS:DO
       call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) then
          exit cycle_read_bonds
       end if

       if (UPPERCASE(line_in(1:5)).eq.'BONDS') then
          allocate(vib_type(1:initial_size),brvib(1:initial_size),brvibk(1:initial_size),stat=jerr)
          if (jerr.ne.0) then
             call err_exit(__FILE__,__LINE__,'init_intramolecular: bonds allocation failed',jerr)
          end if
          brvib=0.0d0
          brvibk=0.0d0
          n=0
          do
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) then
                call err_exit(__FILE__,__LINE__,'Reading section BONDS',jerr)
             end if
             if (UPPERCASE(line_in(1:9)).eq.'END BONDS') exit
             n=n+1
             read(line_in,*) i
             i=addToTable(bonds,i,expand=.true.)
             if (i.gt.ubound(vib_type,1)) then
                call reallocate(vib_type,1,2*ubound(vib_type,1))
                call reallocate(brvib,1,2*ubound(brvib,1))
                call reallocate(brvibk,1,2*ubound(brvibk,1))
             end if
             read(line_in,*) jerr,vib_type(i),brvib(i),brvibk(i)
          end do
          exit cycle_read_bonds
       end if
    END DO CYCLE_READ_BONDS

    ! do j=1,nvmax
    ! if (brvib(j).gt.0) then
    ! write(101,'(I3,1X,I1,1X,F7.5,1X,G13.6)') j,1,brvib(j),brvibk(j)
    ! end if
    ! end do

    ! Looking for section ANGLES
    REWIND(io_ff)
    CYCLE_READ_ANGLES:DO
       call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) then
          exit cycle_read_angles
       end if

       if (UPPERCASE(line_in(1:6)).eq.'ANGLES') then
          allocate(ben_type(1:initial_size),brben(1:initial_size),brbenk(1:initial_size),stat=jerr)
          if (jerr.ne.0) then
             call err_exit(__FILE__,__LINE__,'init_intramolecular: angles allocation failed',jerr)
          end if
          brben=0.0d0
          brbenk=0.0d0
          n=0
          do
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) then
                call err_exit(__FILE__,__LINE__,'Reading section ANGLES',jerr)
             end if
             if (UPPERCASE(line_in(1:10)).eq.'END ANGLES') exit
             n=n+1
             read(line_in,*) i
             i=addToTable(angles,i,expand=.true.)
             if (i.gt.ubound(ben_type,1)) then
                call reallocate(ben_type,1,2*ubound(ben_type,1))
                call reallocate(brben,1,2*ubound(brben,1))
                call reallocate(brbenk,1,2*ubound(brbenk,1))
             end if
             read(line_in,*) jerr,ben_type(i),brben(i),brbenk(i)
             brben(i) = brben(i) * onepi / 180.0d0
          end do
          exit cycle_read_angles
       end if
    END DO CYCLE_READ_ANGLES

    ! do j=1,nvmax
    ! if (brben(j).ne.0) then
    ! write(102,'(I3,1X,I1,1X,F8.4,1X,G13.6)') j,1,brben(j),brbenk(j)
    ! end if
    ! end do

    ! Looking for section DIHEDRALS
    REWIND(io_ff)
    CYCLE_READ_DIHEDRALS:DO
       call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
       if (jerr.ne.0) then
          exit cycle_read_dihedrals
       end if

       if (UPPERCASE(line_in(1:9)).eq.'DIHEDRALS') then
          allocate(vtt(0:9,1:initial_size),torsion_type(1:initial_size),stat=jerr)
          if (jerr.ne.0) then
             call err_exit(__FILE__,__LINE__,'init_intramolecular: dihedrals allocation failed',jerr)
          end if
          vtt=0.0d0
          n=0
          do
             call readLine(io_ff,line_in,skipComment=.true.,iostat=jerr)
             if (jerr.ne.0) then
                call err_exit(__FILE__,__LINE__,'Reading section DIHEDRALS',jerr)
             end if
             if (UPPERCASE(line_in(1:13)).eq.'END DIHEDRALS') exit
             n=n+1
             read(line_in,*) i
             i=addToTable(dihedrals,i,expand=.true.)
             if (i.gt.UBOUND(torsion_type,1)) then
                call reallocate(torsion_type,1,2*UBOUND(torsion_type,1))
                call reallocate(vtt,0,9,1,2*UBOUND(vtt,2))
             end if
             read(line_in,*) jerr,torsion_type(i),vtt(0:torsion_nParameter(torsion_type(i))-1,i)
             if (torsion_type(i).eq.3) then
                vtt(2,i)=vtt(2,i)*onepi/180.0d0
             else if (torsion_type(i).eq.1.or.torsion_type(i).eq.5) then
                ! convert OPLS to Ryckaert-Bellemans because the latter is more efficient
                vtt(0,i)=vtt(0,i)+vtt(1,i)+2.0d0*vtt(2,i)+vtt(3,i)
                vtt(1,i)=-vtt(1,i)+3.0d0*vtt(3,i)
                vtt(2,i)=-2.0d0*vtt(2,i)+8.0d0*vtt(4,i)
                vtt(3,i)=-4.0d0*vtt(3,i)
                vtt(4,i)=-8.0d0*vtt(4,i)
                if (torsion_type(i).eq.1) then
                   torsion_type(i)=7
                else if (torsion_type(i).eq.5) then
                   torsion_type(i)=8
                end if
             end if
          end do
          exit cycle_read_dihedrals
       end if
    END DO CYCLE_READ_DIHEDRALS

    ! do j=1,ntormax
    ! if (vtt0(j).ne.0.or.vtt1(j).ne.0.or.vtt2(j).ne.0.or.vtt3(j).ne.0.or.vtt4(j).ne.0.or.vtt5(j).ne.0.or.vtt6(j).ne.0.or.vtt7(j).ne.0.or.vtt8(j).ne.0.or.vtt9(j).ne.0) then
    ! write(103,'(I3,1X,I1,10(1X,F13.7))') j,1,vtt0(j),vtt1(j),vtt2(j),vtt3(j),vtt4(j),vtt5(j),vtt6(j),vtt7(j),vtt8(j),vtt9(j)
    ! end if
    ! end do

    return
  end subroutine init_energy_bonded

!DEC$ ATTRIBUTES FORCEINLINE :: vtorso
  function vtorso(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,itype)
    use sim_system,only:L_tor_table
    real::vtorso
    real,intent(in)::xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3
    integer,intent(in)::itype
    real::thetac,theta,tac2,tac3,tac4,tac5,tac6,tac7,tac8,tac9

    if (torsion_type(itype).eq.0) then
       ! type 0: dummy torsion type for setting up interaction table
       vtorso=0.0d0
    else if (L_tor_table) then
       call dihedral_angle(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,thetac,theta,.true.)
       vtorso=inter_tor(theta,itype)
    else
       call dihedral_angle(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,thetac,theta,.false.)

        if (torsion_type(itype).eq.7) then
           !< type 7: Ryckaert-Bellemans potential (three terms), angle in polymer convention (trans is 0 deg)
           tac2 = thetac*thetac
           tac3 = tac2*thetac
           vtorso = vtt(0,itype)+vtt(1,itype)*thetac+vtt(2,itype)*tac2+vtt(3,itype)*tac3
        else if (torsion_type(itype).eq.8) then
           !< type 8: Ryckaert-Bellemans potential (four terms), angle in polymer convention (trans is 0 deg)
           tac2 = thetac*thetac
           tac3 = tac2*thetac
           tac4 = tac3*thetac
           vtorso = vtt(0,itype)+vtt(1,itype)*thetac+vtt(2,itype)*tac2+vtt(3,itype)*tac3+vtt(4,itype)*tac4
        else if (torsion_type(itype).eq.2) then
           ! type 2: Ryckaert-Bellemans potential, angle in polymer convention (trans is 0 deg)
           tac2 = thetac*thetac
           tac3 = tac2*thetac
           tac4 = tac3*thetac
           tac5 = tac4*thetac
           tac6 = tac5*thetac
           tac7 = tac6*thetac
           tac8 = tac7*thetac
           tac9 = tac8*thetac
           vtorso = vtt(0,itype)+vtt(1,itype)*thetac+vtt(2,itype)*tac2+vtt(3,itype)*tac3+vtt(4,itype)*tac4+vtt(5,itype)*tac5+vtt(6,itype)*tac6+vtt(7,itype)*tac7+vtt(8,itype)*tac8+vtt(9,itype)*tac9
        else if (torsion_type(itype).eq.3) then
           ! type 3: periodic type, angle in protein convention (trans is 180 deg)
           theta=theta+onepi
           vtorso=vtt(0,itype)*(1+dcos(vtt(1,itype)*theta-vtt(2,itype)))
        else if (torsion_type(itype).eq.4) then
           ! type 4: harmonic type, angle in polymer convention (trans is 0 deg)
           tac2=theta-vtt(1,itype)
           vtorso=vtt(0,itype)*tac2*tac2
        else if (torsion_type(itype).eq.1) then
           ! type 1: OPLS potential (three terms), angle in protein convention (trans is 180 deg)
           theta=theta+onepi
    ! remember: 1 + cos( theta+onepi ) = 1 - cos( theta )
           vtorso = vtt(0,itype) + vtt(1,itype)*(1.0d0-thetac) + vtt(2,itype)*(1.d0-dcos(2.d0*theta)) + vtt(3,itype)*(1.d0+dcos(3.d0*theta))
        else if (torsion_type(itype).eq.5) then
           ! type 5: OPLS potential (four terms), angle in protein convention (trans is 180 deg)
           theta=theta+onepi
           vtorso = vtt(0,itype) + vtt(1,itype)*(1.0d0-thetac) + vtt(2,itype)*(1.d0-dcos(2.d0*theta)) + vtt(3,itype)*(1.d0+dcos(3.d0*theta)) + vtt(4,itype)*(1.d0-dcos(4.d0*theta))
        else if (torsion_type(itype).eq.6) then
           ! type 6: nine-term Fourier cosine series, angle in protein convention (trans is 180 deg)
           theta=theta+onepi
           vtorso=vtt(0,itype)-vtt(1,itype)*thetac+vtt(2,itype)*dcos(2.0d0*theta)+vtt(3,itype)*dcos(3.0d0*theta)+vtt(4,itype)*dcos(4.0d0*theta)+vtt(5,itype)*dcos(5.0d0*theta)+vtt(6,itype)*dcos(6.0d0*theta)+vtt(7,itype)*dcos(7.0d0*theta)+vtt(8,itype)*dcos(8.0d0*theta)+vtt(9,itype)*dcos(9.0d0*theta)
        else
           call err_exit(__FILE__,__LINE__,'you picked a non-defined torsional type',myid+1)
        end if
     end if

    return
  end function vtorso

!> \brief Calculate the dihedral angle and its cosine
!>
!> The dihedral is formed between vectors (1,2), (2,3), and (3,4) using polymer convention (trans is 0 degree)
!DEC$ ATTRIBUTES FORCEINLINE :: dihedral_angle
  subroutine dihedral_angle(xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3,thetac,theta,extended)
    real,intent(in)::xvec1,yvec1,zvec1,xvec2,yvec2,zvec2,xvec3,yvec3,zvec3
    real,intent(out)::thetac,theta
    logical,intent(in)::extended !< whether to extend the dihedral angle to the range of -180 -- +180 degree

    real::x12,y12,z12,x23,y23,z23,d12,d23,dot,tcc,xcc,ycc,zcc

! calculate cross products d_a x d_a-1
    x12 = yvec1 * zvec2 - zvec1 * yvec2
    y12 = zvec1 * xvec2 - xvec1 * zvec2
    z12 = xvec1 * yvec2 - yvec1 * xvec2

! calculate cross products d_a-1 x d_a-2
    x23 = yvec2 * zvec3 - zvec2 * yvec3
    y23 = zvec2 * xvec3 - xvec2 * zvec3
    z23 = xvec2 * yvec3 - yvec2 * xvec3

! calculate lengths of cross products ***
    d12 = dsqrt ( x12*x12 + y12*y12 + z12*z12 )
    d23 = dsqrt ( x23*x23 + y23*y23 + z23*z23 )

! Addition for table look up for Torsion potential
! calculate dot product of cross products ***
    dot = x12*x23 + y12*y23 + z12*z23
    thetac = - (dot / ( d12 * d23 ))

    if (thetac.gt.1.0d0) thetac=1.0d0
    if (thetac.lt.-1.0d0) thetac=-1.0d0
    theta = dacos(thetac)

    if (extended) then
       ! calculate cross product of cross products ***
       xcc = y12*z23 - z12*y23
       ycc = z12*x23 - x12*z23
       zcc = x12*y23 - y12*x23
       ! calculate scalar triple product ***
       tcc = xcc*xvec2 + ycc*yvec2 + zcc*zvec2
       ! determine angle between -180 and 180, not 0 to 180
       if (tcc .lt. 0.0d0) theta = -theta
    end if

    return
  end subroutine dihedral_angle

  function lininter_vib(r,typ) result(tabulated_vib)
    use util_math,only:polint
    use util_search,only:indexOf,LOCATE
    real::tabulated_vib
    real,intent(in)::r
    integer,intent(in)::typ
    integer::low,high

    low=locate(vib(:,typ),vibsplits(typ),r,2)
    high=low+1
    if (vib(low,typ).gt.r.or.vib(high,typ).lt.r) then
       write(io_output,*) 'problem in lininter_vib!'
       write(io_output,*) 'len', r, ' vibtyp', typ
       write(io_output,*) 'low ', low, vib(low, typ)
       write(io_output,*) 'high ',high, vib(high, typ)
       write(io_output,*)
    end if
    call polint(vib(low:high,typ),tabvib(low:high,typ),2,r,tabulated_vib)
    return
  end function lininter_vib

  function lininter_bend(r,typ) result(tabulated_bend)
    use util_math,only:polint
    use util_search,only:indexOf,LOCATE
    real::tabulated_bend
    real,intent(in)::r
    integer,intent(in)::typ
    integer::low,high

    low=locate(bend(:,typ),bendsplits(typ),r,2)
    high=low+1
    if (bend(low,typ).gt.r.or.bend(high,typ).lt.r) then
       write(io_output,*) 'problem in lininter_bend!'
       write(io_output,*) 'r', r, ' bendtyp', typ
       write(io_output,*) 'low ', low, bend(low, typ)
       write(io_output,*) 'high ',high, bend(high, typ)
       write(io_output,*)
    end if
    call polint(bend(low:high,typ),tabbend(low:high,typ),2,r,tabulated_bend)
    return
  end function lininter_bend

  function inter_tor(thetarad,typ) result(tabulated_tor)
    use util_math,only:polint,splint
    use util_search,only:indexOf,LOCATE
    real::tabulated_tor
    real,intent(in)::thetarad
    integer,intent(in)::typ
    real::theta
    integer::low,high

    theta=raddeg*thetarad
    if (L_spline) then
       call splint(deg(:,typ),tabtorso(:,typ),torderiv2(:,typ),splpnts(typ),theta,tabulated_tor)
    else if (L_linear) then
       low=locate(deg(:,typ),splpnts(typ),theta,2)
       high=low+1
       if (deg(low,typ).gt.theta.or.deg(high,typ).lt.theta) then
          write(io_output,*) 'problem in inter_tor_linear!'
          write(io_output,*) 'theta ',thetarad, ' [rad], ',theta, ' [deg]. tortyp', typ
          write(io_output,*) 'low ', low, deg(low, typ)
          write(io_output,*) 'high ',high, deg(high, typ)
          write(io_output,*)
       end if
       call polint(deg(low:high,typ),tabtorso(low:high,typ),2,theta,tabulated_tor)
    end if
    return
  end function inter_tor

! branched and linear molecules with connectivity table -
! go through entire chain -
! calculate all bonds vectors and lengths
!DEC$ ATTRIBUTES FORCEINLINE :: calc_connectivity
  subroutine calc_connectivity(i,imolty)
    use sim_system,only:nunit,nugrow,rxu,ryu,rzu,ijvib,invib
    integer,intent(in)::i,imolty

    real::rxui,ryui,rzui
    integer::ii,iivib,jj

    do ii = 1, nunit(imolty)
       rxui=rxu(i,ii)
       ryui=ryu(i,ii)
       rzui=rzu(i,ii)
       do iivib = 1, invib(imolty,ii)
          jj = ijvib(imolty,ii,iivib)
          rxvec(ii,jj) = rxu(i,jj) - rxui
          ryvec(ii,jj) = ryu(i,jj) - ryui
          rzvec(ii,jj) = rzu(i,jj) - rzui
          distij2(ii,jj) = ( rxvec(ii,jj)**2 + ryvec(ii,jj)**2 + rzvec(ii,jj)**2 )
          distanceij(ii,jj) = dsqrt(distij2(ii,jj))

          if ( nunit(imolty) .ne. nugrow(imolty) )then
! account for explct atoms in opposite direction
             rxvec(jj,ii)   = -rxvec(ii,jj)
             ryvec(jj,ii)   = -ryvec(ii,jj)
             rzvec(jj,ii)   = -rzvec(ii,jj)
             distanceij(jj,ii) = distanceij(ii,jj)
          end if
       end do
    end do
  end subroutine calc_connectivity

!DEC$ ATTRIBUTES FORCEINLINE :: U_torsion
  function U_torsion(i,imolty,ist,lupdate_connectivity) result(vtg)
    use sim_system,only:nunit,intor,ittor,ijtor2,ijtor3,ijtor4,L_tor_table
    real::vtg
    integer,intent(in)::i,imolty,ist
    logical,intent(in)::lupdate_connectivity

    integer::j,jjtor,ip1,ip2,ip3

    if (lupdate_connectivity) call calc_connectivity(i,imolty)

    vtg=0.0d0
    do j = ist, nunit(imolty)
       do jjtor = 1, intor(imolty,j)
          ip3 = ijtor4(imolty,j,jjtor)
          if ( ip3 .lt. j ) then
             ip1 = ijtor2(imolty,j,jjtor)
             ip2 = ijtor3(imolty,j,jjtor)
             vtg = vtg + vtorso(rxvec(j,ip1),ryvec(j,ip1),rzvec(j,ip1),rxvec(ip1,ip2),ryvec(ip1,ip2),rzvec(ip1,ip2),rxvec(ip2,ip3),ryvec(ip2,ip3),rzvec(ip2,ip3),ittor(imolty,j,jjtor))
          end if
       end do
    end do
  end function U_torsion

! calculate all stretching, bending, and torsional potentials
! that have an end-bead with an index smaller than the current bead
!DEC$ ATTRIBUTES FORCEINLINE :: U_bonded
  subroutine U_bonded(i,imolty,vvib,vbend,vtg)
    use sim_system,only:nunit,invib,itvib,ijvib,inben,itben,ijben2,ijben3,L_vib_table
    real,intent(out)::vvib,vbend,vtg
    integer,intent(in)::i,imolty

    real::theta,thetac
    integer::j,jjvib,ip1,ip2,it,jjben

    call calc_connectivity(i,imolty)

! stretching -
    vvib=0.0d0
    do j = 2, nunit(imolty)
       do jjvib = 1, invib(imolty,j)
          ip1 = ijvib(imolty,j,jjvib)
          it  = itvib(imolty,j,jjvib)
          if ( ip1.lt. j .and. L_vib_table) then
             vvib = vvib + lininter_vib(distanceij(ip1,j),it)
! write(io_output,*) 'TABULATED VVIB: ', tabulated_vib,
!   &         distanceij(ip1,j), ip1, j
          end if
          if ( ip1 .lt. j .and..not.L_vib_table) vvib = vvib + brvibk(it) * (distanceij(ip1,j) - brvib(it))**2
       end do
    end do

! bending -
! molecule with bond bending
    vbend=0.0d0
    do j = 2, nunit(imolty)
       do jjben = 1, inben(imolty,j)
          ip2 = ijben3(imolty,j,jjben)
          if ( ip2 .lt. j ) then
             ip1 = ijben2(imolty,j,jjben)
             it  = itben(imolty,j,jjben)
             thetac = ( rxvec(ip1,j)*rxvec(ip1,ip2) + ryvec(ip1,j)*ryvec(ip1,ip2) + rzvec(ip1,j)*rzvec(ip1,ip2) ) / ( distanceij(ip1,j)*distanceij(ip1,ip2) )
             if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
             if ( thetac .le. -1.0d0 ) thetac = -1.0d0

             theta = dacos(thetac)

             ! if (L_bend_table) then
             ! rbendsq=distij2(ip1,j)+distij2(ip1,ip2)-2.0d0*distanceij(ip1,j)*distanceij(ip1,ip2)*thetac
             ! rbend = dsqrt(rbendsq)
             ! vbend = vbend + lininter_bend(rbend,it)
             ! else
             vbend = vbend +  brbenk(it) * (theta-brben(it))**2
             ! end if

! write(io_output,*) 'j,ip1,ip2, it',j,ip1,ip2, it
! write(io_output,*) 'bend energy, theta ',brbenk(it) * (theta-brben(it))**2,theta
          end if
       end do
    end do

! torsions -
! molecule with dihedral potenials ###
    vtg=U_torsion(i,imolty,2,.false.)

  end subroutine U_bonded

!> \brief Read in tabulated potential for bonded interactions (stretching, bending, and torsion) and set up linear interpolation
  subroutine read_tabulated_potential_bonded(file_tab,ntab,r,tab,splits,lists)
    use util_files,only:get_iounit
    character(LEN=*),intent(in)::file_tab
    integer,intent(out)::ntab
    integer,allocatable,intent(inout)::splits(:)
    real,allocatable,intent(inout)::r(:,:),tab(:,:)
    type(LookupTable),intent(inout)::lists

    integer::io_tab,mmm,t,i,jerr

    io_tab=get_iounit()
    open(unit=io_tab,access='sequential',action='read',file=file_tab,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open tabulated potential file: '//file_tab,myid+1)
    end if

    read(io_tab,*) ntab
    do mmm=1,ntab
       read(io_tab,*) t
       t=addToTable(lists,t,expand=.true.)
       if (t.gt.size(splits)) then
          call reallocate(splits,1,2*size(splits))
          call reallocate(r,1,size(r,1),1,2*size(r,2))
          call reallocate(tab,1,size(tab,1),1,2*size(tab,2))
       end if
       i=1
       do
          if (i.gt.size(r,1)) then
             call reallocate(r,1,2*size(r,1),1,size(r,2))
             call reallocate(tab,1,2*size(tab,1),1,size(tab,2))
          end if
          read(io_tab,*,end=100) r(i,t), tab(i,t)
          if (r(i,t).eq.1000) exit
          ! write(io_tab+10,*) i,v(i,t),tab(i,t)
          i=i+1
       end do
100    splits(t)=i-1
    end do
    close(io_tab)
  end subroutine read_tabulated_potential_bonded

! L_spline: Requires file (fort.40) running from -195 to 195 in degree steps
! (Extra 15 degrees on each side required so that second derivatives are
! reasonable for the degrees of interest)
! L_linear: Requires a file (fort.40) running from -180 to 180 in 1/4 degree intervals
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc  Calculates vibrational and 1-3 nonbonded 'bending' potential
!ccc  using linear interpolation between two points
!ccc    must specify equilibrium bond length in suvibe, force
!ccc  constant must be zero
!ccc  requires a file (fort.41) that starts with 0.0 (not 0.5)
!ccc  fort.41: number of tabulated potentials, potential number from
!ccc  suvibe, number of points per angstrom, tabulated potential
!ccc  (repeat last three parts for each additional potential)
!ccc  KM 12/02/08
!ccc    must include 1-3 interactions
!ccc  must specify equilibrium angle in suvibe, force constant
!ccc  must be very small but non-zero
!ccc  requires a file (fort.42) with distances in A
!ccc  fort.42: number of tabulated potentials, potential number from
!ccc  suvibe, number of points per degree, tabulated potential
!ccc  (repeat last three parts for each additional potential,
!ccc   separated by 1000 1000)
!ccc  make sure potential does not go up to infinity!
!ccc  KM 12/03/08
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine init_tabulated_potential_bonded()
    use util_math,only:spline
    use sim_system,only:L_tor_table,L_vib_table,L_bend_table
    integer,parameter::initial_size=10,grid_size=1500
    integer::jerr,ttor

    if (L_tor_table) then
       allocate(splpnts(1:initial_size),deg(1:grid_size,1:initial_size),tabtorso(1:grid_size,1:initial_size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_bonded: allocation failed for vib_table 1',myid+1)
       call read_tabulated_potential_bonded('fort.40',nttor,deg,tabtorso,splpnts,dihedrals)
       if (L_spline) then
          if (myid.eq.0) write(io_output,*) 'using spline interpolation'
          allocate(torderiv2(1:grid_size,1:initial_size),stat=jerr)
          if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_bonded: allocation failed for tor_table 2',myid+1)
          do ttor=1,nttor
             call spline(deg(:,ttor),tabtorso(:,ttor),splpnts(ttor),1.0d31,1.0d31,torderiv2(:,ttor))
          end do
       else if (L_linear) then
          if (myid.eq.0) write(io_output,*) 'using linear interpolation'
       end if
    end if

    if (L_vib_table) then
       allocate(vibsplits(1:initial_size),vib(1:grid_size,1:initial_size),tabvib(1:grid_size,1:initial_size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_bonded: allocation failed for vib_table',myid+1)
       call read_tabulated_potential_bonded('fort.41',ntabvib,vib,tabvib,vibsplits,bonds)
       if (myid.eq.0)  write(io_output,'(/,A)') 'using linear interpolation for vibrations'
    end if

    if (L_bend_table) then
       allocate(bendsplits(1:initial_size),bend(1:grid_size,1:initial_size),tabbend(1:grid_size,1:initial_size),stat=jerr)
       if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'init_tabulated_potential_bonded: allocation failed for bend_table',myid+1)
       call read_tabulated_potential_bonded('fort.42',ntabbend,bend,tabbend,bendsplits,angles)
       if (myid.eq.0) write(io_output,*) 'using linear interpolation for 1-3 nonbonded bending'
    end if
  end subroutine init_tabulated_potential_bonded

  subroutine allocate_energy_bonded()
    use sim_system,only:numax
    integer::jerr
    allocate(rxvec(numax,numax),ryvec(numax,numax),rzvec(numax,numax),distij2(numax,numax),distanceij(numax,numax),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'allocate_energy_bonded: allocation failed',jerr)
    end if
  end subroutine allocate_energy_bonded
end MODULE energy_intramolecular
