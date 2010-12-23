      subroutine align_planes(isplist,imol_b,imol_a,itype_b,itype_a,
     &     xb,yb,zb)

!--- align_planes.f
!
!  designed with swatching rigid planar PAH molecules in mind
!  it will work for other rigid stuff, but it's not tailored for that
!
!  1) Translates so swatch bead 1 of the swathced molecule is in 
!     same position as swatch bead 1 of the other molecule
!  2) Rotate the swathched molecule so its swatch bead1--bead2 vector
!     is aligned with the swatch bead1--bead2 of the other molecule
!  3) Rotate the swatched molecule again so its swatch bead2--bead3 
!     vector is parallel with the the bead2--bead3 of the other molecule
!
!  RESULTS:
!  A) The planes defined by the three swatched beads in each molecule 
!     will become coplanar.  Thus, if both molecules are planar, then they 
!     will definitely become coplanar.
!  B) If both molecules have the same 1--2--3 "bond" angles and
!     1--2 "bond" lengths, the 1--2--3 bonds will be exactly aligned
!     (where 1,2,3 refer to the three beads you swatched). If not, 
!     best of luck; I hope you know what you are trying to do :)
!
!  Written by Jake L. Rafferty on the fine day of 2.28.07
!
! --  a subscripts refer to molecule that exists, the "other"
! --  b subscripts refer to molecule that being swatched in, the "self"
!


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
!$$$      include 'control.inc'
!$$$      include 'conver.inc'
!$$$      include 'coord.inc'
!$$$      include 'swtcmove.inc'

!     --- INPUT VARIABLES ---
      integer(KIND=normal_int)::isplist, imol_b, imol_a, itype_b,
     & itype_a
   
!     --- OUTPUT VARIABLES ---
      real(KIND=double_precision)::xb(numax), yb(numax), zb(numax)

!     --- LOCAL VARIABLES ---
      integer(KIND=normal_int)::i, imoltype_b, nunit_b, nunit_a
      integer(KIND=normal_int)::ia_bead1, ia_bead2, ia_bead3
      integer(KIND=normal_int)::ib_bead1, ib_bead2, ib_bead3

      real(KIND=double_precision)::xa(3),   ya(3),   za(3)
      real(KIND=double_precision)::xorigin, yorigin, zorigin
      real(KIND=double_precision)::xcross,  ycross,  zcross
      real(KIND=double_precision)::xtemp,   ytemp,   ztemp
      real(KIND=double_precision)::dxa,     dya,     dza
      real(KIND=double_precision)::dxb,     dyb,     dzb

      real(KIND=double_precision)::rnorm, d
      real(KIND=double_precision)::gamma_a, gamma_b
      real(KIND=double_precision)::theta, cos_theta, sin_theta

      

!      write(iou,*) 'BEGIN align_planes'

      imoltype_b = moltyp(imol_b)
      nunit_b = nunit(imoltype_b)

!     Find the three beads on each that are being swatched
      ia_bead1=splist(isplist,1,itype_a)
      ia_bead2=splist(isplist,2,itype_a)
      ia_bead3=splist(isplist,3,itype_a)
      ib_bead1=splist(isplist,1,itype_b)
      ib_bead2=splist(isplist,2,itype_b)
      ib_bead3=splist(isplist,3,itype_b)

!CC -- FIRST ROTATION

!CC -- STEP 1, translate bead 1 of both molecules to origin

!     save position of bead_1 on molecule_a
      xorigin = rxu(imol_a,ia_bead1)
      yorigin = ryu(imol_a,ia_bead1)
      zorigin = rzu(imol_a,ia_bead1)
!     use this bead as the origin
      xa(1) = 0.0d0
      ya(1) = 0.0d0
      za(1) = 0.0d0
      xa(2) = rxu(imol_a,ia_bead2) - xorigin
      ya(2) = ryu(imol_a,ia_bead2) - yorigin
      za(2) = rzu(imol_a,ia_bead2) - zorigin
      xa(3) = rxu(imol_a,ia_bead3) - xorigin
      ya(3) = ryu(imol_a,ia_bead3) - yorigin
      za(3) = rzu(imol_a,ia_bead3) - zorigin

!     translate molecule be to the origin
      do i=1,nunit_b
         xb(i) = rxu(imol_b,i) - rxu(imol_b,ib_bead1)
         yb(i) = ryu(imol_b,i) - ryu(imol_b,ib_bead1)
         zb(i) = rzu(imol_b,i) - rzu(imol_b,ib_bead1)
      end do

               
!     Get first rotation vector -- the vector orthogonal to 
!     both of the 1--2 vectors, i.e. the cross product

!     take cross product of 1--2 vectors
      xcross = ya(2)*zb(ib_bead2) - za(2)*yb(ib_bead2)
      ycross = za(2)*xb(ib_bead2) - xa(2)*zb(ib_bead2)
      zcross = xa(2)*yb(ib_bead2) - ya(2)*xb(ib_bead2)
!     normalize cross product to unit length
      rnorm = sqrt(xcross*xcross + ycross*ycross + zcross*zcross)
      xcross = xcross/rnorm
      ycross = ycross/rnorm
      zcross = zcross/rnorm
!     find projection of this unit vector onto yz plane
      d = sqrt(ycross*ycross + zcross*zcross)

!CC -- STEP 2, rotate space about x-axis so that rotation lies in yz plane

      do i=2,3
         ytemp = (ya(i)*zcross - za(i)*ycross)/d
         ztemp = (ya(i)*ycross + za(i)*zcross)/d
         ya(i) = ytemp
         za(i) = ztemp
      end do

      do i=1,nunit_b
         ytemp = (yb(i)*zcross - zb(i)*ycross)/d
         ztemp = (yb(i)*ycross + zb(i)*zcross)/d
         yb(i) = ytemp
         zb(i) = ztemp
      end do

!CC -- STEP 3, rotate space so that rotation axis is along positive z-axis

      do i=2,3
         xtemp = xa(i)*d - za(i)*xcross
         ztemp = xa(i)*xcross + za(i)*d
         xa(i) = xtemp
         za(i) = ztemp
      end do

      do i=1,nunit_b
         xtemp = xb(i)*d - zb(i)*xcross
         ztemp = xb(i)*xcross + zb(i)*d
         xb(i) = xtemp
         zb(i) = ztemp
      end do

!CC -- STEP 4, Rotate by theta that about z-axis
      
!     first we need the angle between the two 1--2 vectors

!     angle of vector for molecule_a with x axis
      gamma_a = datan(abs(ya(2)/xa(2)))
      if (ya(2).lt.0.0) then
         if (xa(2).lt.0.0) then
            gamma_a = gamma_a + onepi
         else
            gamma_a = twopi - gamma_a
         end if
      else
         if (xa(2).lt.0.0) gamma_a = onepi - gamma_a
      end if

!     angle of vector for molecule_b with x axis
      gamma_b = datan(abs(yb(ib_bead2)/xb(ib_bead2)))
      if (yb(ib_bead2).lt.0.0) then
         if (xb(ib_bead2).lt.0.0) then
            gamma_b = gamma_b + onepi
         else
            gamma_b = twopi - gamma_b
         end if
      else
         if (xb(ib_bead2).lt.0.0) gamma_b = onepi - gamma_b
      end if

      theta = gamma_a - gamma_b
      cos_theta = dcos(theta)
      sin_theta = dsin(theta)

!     now rotate molecule_b by theta
      do i=1,nunit_b
         xtemp = xb(i)*cos_theta - yb(i)*sin_theta
         ytemp = xb(i)*sin_theta + yb(i)*cos_theta
         xb(i) = xtemp
         yb(i) = ytemp
      end do

!CC -- INVERT STEP 3

      do i=2,3
         xtemp = xa(i)*d + za(i)*xcross
         ztemp = -xa(i)*xcross + za(i)*d
         xa(i) = xtemp
         za(i) = ztemp
      end do

      do i=1,nunit_b
         xtemp = xb(i)*d + zb(i)*xcross
         ztemp = -xb(i)*xcross + zb(i)*d
         xb(i) = xtemp
         zb(i) = ztemp
      end do

!CC --- INVERT STEP 2

      do i=2,3
         ytemp = (ya(i)*zcross + za(i)*ycross)/d
         ztemp = (-ya(i)*ycross + za(i)*zcross)/d
         ya(i) = ytemp
         za(i) = ztemp
      end do

      do i=1,nunit_b
         ytemp = (yb(i)*zcross + zb(i)*ycross)/d
         ztemp = (-yb(i)*ycross + zb(i)*zcross)/d
         yb(i) = ytemp
         zb(i) = ztemp
      end do


!CC --- NOW PERFORM THE SECOND ROTATION
!CC --- THIS ONE IS ABOUT THE 1--2 VECTOR OF MOLECULE_a

!CC --- STEP 1, already done in last rotation

!CC --- STEP 2, rotate space about x-axis so that rotation lies in yz plane

!     The unit rotation axis
      rnorm = dsqrt( xa(2)*xa(2) + ya(2)*ya(2) + za(2)*za(2) )
      xcross = xa(2)/rnorm
      ycross = ya(2)/rnorm
      zcross = za(2)/rnorm
!     find projection of this unit vector onto yz plane
      d = sqrt(ycross*ycross + zcross*zcross)

      do i=2,3
         ytemp = (ya(i)*zcross - za(i)*ycross)/d
         ztemp = (ya(i)*ycross + za(i)*zcross)/d
         ya(i) = ytemp
         za(i) = ztemp
      end do

      do i=1,nunit_b
         ytemp = (yb(i)*zcross - zb(i)*ycross)/d
         ztemp = (yb(i)*ycross + zb(i)*zcross)/d
         yb(i) = ytemp
         zb(i) = ztemp
      end do

!CC -- STEP 3, rotate space so that rotation axis is along positive z-axis

      do i=2,3
         xtemp = xa(i)*d - za(i)*xcross
         ztemp = xa(i)*xcross + za(i)*d
         xa(i) = xtemp
         za(i) = ztemp
      end do

      do i=1,nunit_b
         xtemp = xb(i)*d - zb(i)*xcross
         ztemp = xb(i)*xcross + zb(i)*d
         xb(i) = xtemp
         zb(i) = ztemp
      end do

!CC -- STEP 4, Rotate by theta that about z-axis

!     This time theta is angle between the two planes
!     which correponds to the 2--3 vectors
      
      dxa = xa(3) - xa(2)
      dya = ya(3) - ya(2)

      dxb = xb(ib_bead3) - xb(ib_bead2)
      dyb = yb(ib_bead3) - yb(ib_bead2)
!     angle of vector for molecule_a with x axis
      gamma_a = datan(abs(dya/dxa))
      if (dya.lt.0.0) then
         if (dxa.lt.0.0) then
            gamma_a = gamma_a + onepi
         else
            gamma_a = twopi - gamma_a
         end if
      else
         if (dxa.lt.0.0) gamma_a = onepi - gamma_a
      end if

!     angle of vector for molecule_b with x axis
      gamma_b = datan(abs(dyb/dxb))
      if (dyb.lt.0.0) then
         if (dxb.lt.0.0) then
            gamma_b = gamma_b + onepi
         else
            gamma_b = twopi - gamma_b
         end if
      else
         if (dxb.lt.0.0) gamma_b = onepi - gamma_b
      end if

      theta = gamma_a - gamma_b
      cos_theta = dcos(theta)
      sin_theta = dsin(theta)

!     now rotate molecule_b by theta
      do i=1,nunit_b
         xtemp = xb(i)*cos_theta - yb(i)*sin_theta
         ytemp = xb(i)*sin_theta + yb(i)*cos_theta
         xb(i) = xtemp
         yb(i) = ytemp
      end do

!CC -- INVERT STEP 3

      do i=1,nunit_b
         xtemp = xb(i)*d + zb(i)*xcross
         ztemp = -xb(i)*xcross + zb(i)*d
         xb(i) = xtemp
         zb(i) = ztemp
      end do

!CC --- INVERT STEP 2

      do i=1,nunit_b
         ytemp = (yb(i)*zcross + zb(i)*ycross)/d
         ztemp = (-yb(i)*ycross + zb(i)*zcross)/d
         yb(i) = ytemp
         zb(i) = ztemp
      end do

!CC --- INVERT STEP 1

      do i=1,nunit_b
         xb(i) = xb(i) + xorigin
         yb(i) = yb(i) + yorigin
         zb(i) = zb(i) + zorigin
      end do

!CC --- FINALLY WE ARE DONE

!     Write a check
!      open(unit=90, file='align.xyz',status='unknown')
!      nunit_a = nunit(moltyp(imol_a))
!      write(90,*) nunit_a+nunit_b
!      write(90,*) 'Plane Alignment Test'
!      do i=1,nunit_a
!         write(90,*) 'C ',rxu(imol_a,i),ryu(imol_a,i),rzu(imol_a,i)
!      end do
!      do i=1,nunit_b
!         write(90,*) 'O ',xb(i),yb(i),zb(i)
!      end do
!      close(90)
!      write(iou,*) 'END align_planes'

      return

      end
