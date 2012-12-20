MODULE sim_cell
  use var_type,only:RealPtr
  use sim_system,only:nbxmax
  implicit none
  private
  save
  public::CellType,CellMaskType,matops,setpbc,mimage

  type CellType
     logical::ortho,solid ! ortho=.true. if the simulation box is orthorhombic
                          ! solid=.true. if the simulation box can change shape
     logical,dimension(3)::pbc
     real::cut,vol,calp,boxl(3),ang(3),height(3),hmat(3,3),hmati(3,3)
  end type CellType

  type CellMaskType
     logical,pointer::ortho,solid ! ortho=.true. if the simulation box is orthorhombic
                          ! solid=.true. if the simulation box can change shape
     logical,dimension(3)::pbc
     real,pointer::cut,vol,calp
     type(RealPtr)::boxl(3),ang(3),height(3),hmat(3,3),hmati(3,3)
  end type CellMaskType

  real,target,public::hmat(nbxmax,9),hmati(nbxmax,9),cell_length(nbxmax,3),min_width(nbxmax,3),cell_vol(nbxmax),cell_ang(nbxmax,3)

! PEBOCO.INC
  real::bx,by,bz,hbx,hby,hbz,bxi,byi,bzi

contains
!> \brief Calculates for non cubic simulation cell:
!>       boxlengths: cell_length,
!>       minimum boxwidths: min_length,
!>       boxvolume: cell_vol,
!>       inverse H matrix
!> \author Neeraj Rai (in Merck Apr 2005)
  subroutine matops(ibox)
    use util_runtime,only:err_exit
    use sim_system,only:boxlx,boxly,boxlz
    integer,intent(in)::ibox
    real::abx,aby,abz,bcx,bcy,bcz,cax,cay,caz
    real::elem(9),inv_vol,adj(9),cosa,cosb,cosg
    integer::i

    do i=1,9
       elem(i)=hmat(ibox,i)
    end do

    !     -- calculating the length of cell vectors

    cell_length(ibox,1)=sqrt(elem(1)*elem(1)+elem(2)*elem(2)+ elem(3)*elem(3))
    cell_length(ibox,2)=sqrt(elem(4)*elem(4)+elem(5)*elem(5)+ elem(6)*elem(6))
    cell_length(ibox,3)=sqrt(elem(7)*elem(7)+elem(8)*elem(8)+ elem(9)*elem(9))

    boxlx(ibox) = cell_length(ibox,1)
    boxly(ibox) = cell_length(ibox,2)
    boxlz(ibox) = cell_length(ibox,3)

    !    -- calculating cross product of cell vectors

    abx=elem(2)*elem(6)-elem(3)*elem(5)
    aby=elem(3)*elem(4)-elem(1)*elem(6)
    abz=elem(1)*elem(5)-elem(2)*elem(4)
    bcx=elem(5)*elem(9)-elem(6)*elem(8)
    bcy=elem(6)*elem(7)-elem(4)*elem(9)
    bcz=elem(4)*elem(8)-elem(5)*elem(7)
    cax=elem(8)*elem(3)-elem(2)*elem(9)
    cay=elem(1)*elem(9)-elem(3)*elem(7)
    caz=elem(2)*elem(7)-elem(1)*elem(8)

    !    -- calculating cell volume

    cell_vol(ibox) = elem(1)*bcx+elem(2)*bcy+elem(3)*bcz

    if(abs(cell_vol(ibox)).lt.1D-16) then
       call err_exit('Volume of cell negligible, check input H matrix')
    end if

    inv_vol = 1.0d0/cell_vol(ibox)

    !    -- calculating minimum cell widths

    min_width(ibox,1) = cell_vol(ibox)/sqrt(bcx*bcx+bcy*bcy+ bcz*bcz)
    min_width(ibox,2) = cell_vol(ibox)/sqrt(cax*cax+cay*cay+ caz*caz)
    min_width(ibox,3) = cell_vol(ibox)/sqrt(abx*abx+aby*aby+ abz*abz)

    !    -- calculating adjoint for inverting the h-matrix

    adj(1)=elem(5)*elem(9)-elem(6)*elem(8)
    adj(2)=elem(3)*elem(8)-elem(2)*elem(9)
    adj(3)=elem(2)*elem(6)-elem(3)*elem(5)
    adj(4)=elem(6)*elem(7)-elem(4)*elem(9)
    adj(5)=elem(1)*elem(9)-elem(3)*elem(7)
    adj(6)=elem(3)*elem(4)-elem(1)*elem(6)
    adj(7)=elem(4)*elem(8)-elem(5)*elem(7)
    adj(8)=elem(2)*elem(7)-elem(1)*elem(8)
    adj(9)=elem(1)*elem(5)-elem(2)*elem(4)

    !     --  inverting the matrix


    hmati(ibox,1) = inv_vol * adj(1)
    hmati(ibox,2) = inv_vol * adj(2)
    hmati(ibox,3) = inv_vol * adj(3)
    hmati(ibox,4) = inv_vol * adj(4)
    hmati(ibox,5) = inv_vol * adj(5)
    hmati(ibox,6) = inv_vol * adj(6)
    hmati(ibox,7) = inv_vol * adj(7)
    hmati(ibox,8) = inv_vol * adj(8)
    hmati(ibox,9) = inv_vol * adj(9)

    !   ---  calculating alpha, beta and gamma using the dot product rule

    !   ---  cell_ang(ibox,1)=alpha;2= beta; 3= gamma
    !   ---  In crystallography Literature angle between b & c is alpha (1),
    !        angle between c & a is Beta (2) and angle between a & b is gamma (3)


    cosa = (elem(4)*elem(7)+elem(5)*elem(8)+elem(6)*elem(9))/ (cell_length(ibox,2)*cell_length(ibox,3))

    cosb = (elem(1)*elem(7)+elem(2)*elem(8)+elem(3)*elem(9))/ (cell_length(ibox,1)*cell_length(ibox,3))

    cosg = (elem(4)*elem(1)+elem(5)*elem(2)+elem(6)*elem(3))/ (cell_length(ibox,2)*cell_length(ibox,1))

    cell_ang(ibox,1) = dacos(cosa)
    cell_ang(ibox,2) = dacos(cosb)
    cell_ang(ibox,3) = dacos(cosg)

    return
  end subroutine matops

  subroutine setpbc(ibox)
    use sim_system,only:boxlx,boxly,boxlz,lpbcx,lpbcy,lpbcz,lfold
    integer,intent(in)::ibox

    if ( lpbcx ) then
       bx = boxlx(ibox)
       if ( lfold ) then
          hbx = 0.5d0 * bx
       else
          bxi = 1.0d0 / bx
       end if
    end if

    if ( lpbcy ) then
       by = boxly(ibox)
       if ( lfold ) then
          hby = 0.5d0 * by
       else
          byi = 1.0d0 / by
       end if
    end if

    if ( lpbcz ) then
       bz = boxlz(ibox)
       if ( lfold ) then
          hbz = 0.5d0 * bz
       else
          bzi = 1.0d0 / bz
       end if
    end if

    return
  end subroutine setpbc

  pure subroutine mimage(rxuij,ryuij,rzuij,ibox)
    use sim_system,only:lsolid,lrect,lpbcx,lpbcy,lpbcz,lfold,ltwice
    !$$$      include 'peboco.inc'
    !$$$      include 'control.inc'
    !$$$      include 'system.inc'
    !$$$      include 'cell.inc'
    integer,intent(in)::ibox
    real,intent(inout)::rxuij,ryuij,rzuij

    real::hsx,hsy,hsz,sx,sy,sz

    if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
       ! *** non-hexagonal box
       hsx = 0.5d0*hmat(ibox,1)
       hsy = 0.5d0*hmat(ibox,5)
       hsz = 0.5d0*hmat(ibox,9)
       ! sx, sy, sz are the coordinates of vector (rxuij,ryuij,rzuij) in the
       ! basis (n1,n2,n3), which is the transpose of the H matrix, and is the
       ! transforming matrix from basis (n1,n2,n3) to the canonical basis
       ! (e1,e2,e3), where e1, e2, e3 are the three perpendicular unit vector.
       sx = rxuij*hmati(ibox,1)+ryuij*hmati(ibox,4) +rzuij*hmati(ibox,7)
       sy = rxuij*hmati(ibox,2)+ryuij*hmati(ibox,5) +rzuij*hmati(ibox,8)
       sz = rxuij*hmati(ibox,3)+ryuij*hmati(ibox,6) +rzuij*hmati(ibox,9)

       !         if ( sx .gt. 0.5d0 ) then
       !            sx = sx-1d0
       !         else if ( sx .lt. -0.5d0 ) then
       !            sx = sx+1d0
       !         end if
       !         if ( sy .gt. 0.5d0 ) then
       !            sy = sy-1d0
       !         else if ( sy .lt. -0.5d0 ) then
       !            sy = sy+1d0
       !         end if
       !         if ( sz .gt. 0.5d0 ) then
       !            sz = sz-1d0
       !         else if ( sz .lt. -0.5d0 ) then
       !            sz = sz+1d0
       !         end if
       !         sx = sx-nint(sx)
       !         sy = sy-nint(sy)
       !         sz = sz-nint(sz)
       !         rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
       !         ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)+sz*hmat(ibox,8)
       !         rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)+sz*hmat(ibox,9)

       !         print*,rxuij,ryuij,rzuij

       ! Here it implies that in the H matrix, the first vector must only have
       ! the x component, the second vector the x,y components, and only the
       ! third can have all the x,y,z components
       if ( rzuij .gt. hsz ) then
          rzuij=rzuij-hmat(ibox,9)
          sz=sz-1d0
          if ( rzuij .gt. hsz ) then
             rzuij=rzuij-hmat(ibox,9)
             sz=sz-1d0
          end if
       else if ( rzuij .lt. -hsz ) then
          rzuij=rzuij+hmat(ibox,9)
          sz=sz+1d0
          if ( rzuij .lt. -hsz ) then
             rzuij=rzuij+hmat(ibox,9)
             sz=sz+1d0
          end if
       end if

       ryuij=sy*hmat(ibox,5)+sz*hmat(ibox,8)
       if ( ryuij .gt. hsy ) then
          ryuij=ryuij-hmat(ibox,5)
          sy=sy-1d0
          if ( ryuij .gt. hsy ) then
             ryuij=ryuij-hmat(ibox,5)
             sy=sy-1d0
          end if
       else if ( ryuij .lt. -hsy ) then
          ryuij=ryuij+hmat(ibox,5)
          sy=sy+1d0
          if ( ryuij .lt. -hsy ) then
             ryuij=ryuij+hmat(ibox,5)
             sy=sy+1d0
          end if
       end if

       rxuij=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
       if ( rxuij .gt. hsx ) then
          rxuij=rxuij-hmat(ibox,1)
          if ( rxuij .gt. hsx ) then
             rxuij=rxuij-hmat(ibox,1)
          end if
       else if ( rxuij .lt. -hsx ) then
          rxuij=rxuij+hmat(ibox,1)
          if ( rxuij .lt. -hsx ) then
             rxuij=rxuij+hmat(ibox,1)
          end if
       end if
       !         print*, sqrt(rxuij**2+ryuij**2+rzuij**2)
       !         print*
    else
       ! *** orthorhombic box
       if ( lpbcx ) then
          if ( lfold ) then
             if ( rxuij .gt. hbx ) then
                rxuij=rxuij-bx
             else
                if (rxuij.lt.-hbx) rxuij=rxuij+bx
             end if
          else
             !            rxuij = rxuij - bx*anint(rxuij*bxi)
             rxuij = rxuij - bx*dint(rxuij*bxi+dsign(0.5d0,rxuij))
          end if
       end if

       if ( lpbcy ) then
          if ( lfold ) then
             if ( ryuij .gt. hby ) then
                ryuij=ryuij-by
             else
                if (ryuij.lt.-hby) ryuij=ryuij+by
             end if
          else
             !           ryuij  = ryuij - by*anint(ryuij*byi)
             ryuij = ryuij - by*dint(ryuij*byi+dsign(0.5d0,ryuij))
          end if
       end if

       if ( lpbcz ) then
          if ( lfold ) then
             if (rzuij.gt.hbz) then
                rzuij=rzuij-bz
             else
                if (rzuij.lt.-hbz) rzuij=rzuij+bz
             end if
          else
             !            rzuij  = rzuij - bz*anint(rzuij*bzi)
             rzuij = rzuij - bz*dint(rzuij*bzi+dsign(0.5d0,rzuij))
          end if
       end if

       ! --- JLR 11-17-09
       ! --- for orginal RPLC setup, mimage must be applied twice or energy errors can result
       ! --- this may be removed in future if larger system size is used!!!
       if (ltwice(ibox)) then !everthing is folded again
          if ( rxuij .gt. hbx ) then
             rxuij=rxuij-bx
          else
             if (rxuij.lt.-hbx) rxuij=rxuij+bx
          end if
          if ( ryuij .gt. hby ) then
             ryuij=ryuij-by
          else
             if (ryuij.lt.-hby) ryuij=ryuij+by
          end if
       end if
       ! --- END JLR 11-17-09
    end if
    return
  end subroutine mimage
end MODULE sim_cell
