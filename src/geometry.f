      subroutine geometry(lnew,iw,i,imolty,angstart,iuprev,glist
     &     ,bondlen,bendang,phi,vvibtr,vbbtr, maxlen, wei_bend )

!     ***********************************************************************
!     ** determines the new geometry of the bond lengths and angles to be  **
!     ** rotated on the cone for rosenb.f                                  **
!     ** for old computes the rosenbluth weight for bending and determines **
!     ** the old bond lenghts, angles, and phi for growth                  **
!     ** bondlen(count) is the bondlengths from the grow bead to count     **
!     ** bendang(count) is the bond angle between iuprev,iufrom, and count **
!     ** phi(count) is the angle around the cone between count and 1       **
!     ** written by M.G. Martin 7-10-98 from geomnew and geomold           **
!     ** last modified by Neeraj Rai on 12/23/2008 for CG models           **
!     ***********************************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'connect.inc'
!$$$      include 'conver.inc'
!$$$      include 'rosen.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'tabulated.inc'

!     --- variables passed to/from the subroutine

      logical::lnew
      integer(KIND=int)::iw,imolty,angstart,iuprev,glist,i
      real(KIND=double_precision)::bondlen,bendang,phi,vvibtr,vbbtr,maxlen
      real(KIND=double_precision)::wei_bend

      dimension glist(numax)
      dimension bondlen(numax),bendang(numax),phi(numax)

!     --- local variables
      
      integer(KIND=int)::count,ntogrow,iugrow,iufrom,iv,juvib
     &     ,jtvib,iu2back,ib,iulast,type,aaa,iuone

      real(KIND=double_precision)::equil,kforce,thetaone,thetatwo
     &     ,vvib,length,angle,phione,phitwo,vangle,random,vphi

!     --- new variables
      integer(KIND=int)::ibend,nchben_a,nchben_b,start,start_ang
      real(KIND=double_precision)::bsum_try,rsint,ang_trial,bfactor,rbf,bs

      dimension ang_trial(nchbn_max),bfactor(nchbn_max)

!     --- variables from geomold
      real(KIND=double_precision)::rxui,ryui,rzui,rxuij,ryuij,rzuij
     &     ,xvecprev,yvecprev,zvecprev,distprev,xvecgrow
     &     ,yvecgrow,zvecgrow,distgrow,anglec
      real(KIND=double_precision)::xub,yub,zub,dum,ux,uy,uz,alpha,gamma
      real(KIND=double_precision)::tabulated_vib, tabulated_bend
      real(KIND=double_precision)::rbend, rbendsq

! Neeraj: Adding for the lookup table for CG model

      real(KIND=double_precision)::distprev2,distgrow2
      real(KIND=double_precision)::lengtha,lengthb,lengtha2,lengthb2,
     &                 lengthc,lengthc2,lengthFP,lengthFP2

!     --- assign nchben_a and ncben_b
      nchben_a = nchbna(imolty)
      nchben_b = nchbnb(imolty)

!      write(iou,*) 'START GEOMETRY'

!     --- initialize trial energies
      vvibtr = 0.0d0
      vbbtr = 0.0d0

!     --- assign grownum and growfrom to local variables
      ntogrow = grownum(iw)
      iufrom = growfrom(iw)

      if ( .not. lnew ) then
!        --- OLD store r*ui positions for unit iufrom
         rxui = rxu(i,iufrom)
         ryui = ryu(i,iufrom)
         rzui = rzu(i,iufrom)
      end if

!     --- Begin Bond length selection based on Boltzmann rejection

!     --- determine the bond lengths of the beads to be grown
      maxlen = 0.0d0

      do count = 1,ntogrow
         iugrow = growlist(iw,count)
         glist(count) = iugrow
            
!        --- determine the vibration (bond) type as jtvib
         do iv = 1, invib(imolty,iugrow)
            juvib = ijvib(imolty,iugrow,iv)
            if ( juvib .eq. iufrom ) then
               jtvib = itvib(imolty,iugrow,iv)
            end if
         end do

         if ( lnew ) then
!           --- compute bond length
            equil = brvib(jtvib)
            kforce = brvibk(jtvib)
            call bondlength( jtvib,equil,kforce,beta,length,vvib )

         else
!           --- compute bond length
            rxuij = rxu(i,iugrow) - rxui
            ryuij = ryu(i,iugrow) - ryui
            rzuij = rzu(i,iugrow) - rzui

!           --- do not need mimage for intramolecular
            length = dsqrt(rxuij*rxuij + ryuij*ryuij + rzuij*rzuij)
         
!           --- compute vibration energy

            if (L_vib_table) then
               call lininter_vib(length,tabulated_vib, jtvib)
               vvib = tabulated_vib
            else
               equil = brvib(jtvib)
               kforce = brvibk(jtvib)
               vvib = kforce * (length-equil )**2
            end if
         end if
            
!        --- adjust maximum bond length of those being grown
         if ( length .gt. maxlen ) maxlen = length

!        --- assign bondlength and add up vibrational energy
         bondlen(count) = length
         vvibtr = vvibtr + vvib
      end do

      if ( growprev(iw) .eq. 0 ) then
!        --- need to choose one bead to grow on unit sphere
         iuprev = growlist(iw,1)
         angstart = 2
      else
         iuprev = growprev(iw)
         angstart = 1
      end if

      if (L_bend_table) then
         if ( growprev(iw) .eq. 0 ) then
             iuprev = growlist(iw,1)
             lengthFP = bondlen(1)
             lengthFP2 = lengthFP*lengthFP
         else
             iuprev = growprev(iw)
             if(.not.lnew) then
                rxuij = rxu(i,iuprev) - rxu(i,iufrom)
                ryuij = ryu(i,iuprev) - ryu(i,iufrom)
                rzuij = rzu(i,iuprev) - rzu(i,iufrom)
                lengthFP = dsqrt(rxuij*rxuij+ryuij*ryuij+rzuij*rzuij)
                lengthFP2 = lengthFP*lengthFP
              else
                rxuij = rxnew(iuprev) - rxnew(iufrom)
                ryuij = rynew(iuprev) - rynew(iufrom)
                rzuij = rznew(iuprev) - rznew(iufrom)
                lengthFP = dsqrt(rxuij*rxuij+ryuij*ryuij+rzuij*rzuij)
                lengthFP2 = lengthFP*lengthFP
              end if
         end if
      end if     


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     --- Begin Bond angle biased selection ---                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if ( .not. lnew ) then
!        --- compute the vector from iufrom to iuprev
         xvecprev = rxu(i,iuprev) - rxui
         yvecprev = ryu(i,iuprev) - ryui
         zvecprev = rzu(i,iuprev) - rzui
         distprev = dsqrt( xvecprev*xvecprev + yvecprev*yvecprev 
     &        + zvecprev*zvecprev )
      end if

!     --- initialize wei_bend
      wei_bend = 1.0d0

!     --- determine the iugrow-iufrom-iuprev angles 
      do count = angstart,ntogrow
         iugrow = growlist(iw,count)
         do ib = 1, inben(imolty,iugrow)
            iulast = ijben2(imolty,iugrow,ib)
            if ( iulast .eq. iufrom ) then
               iu2back = ijben3(imolty,iugrow,ib)
               if ( iu2back .eq. iuprev ) then
                  type = itben(imolty,iugrow,ib)
                  equil = brben(type)
                  kforce = brbenk(type)
               end if
            end if
         end do


         if ( kforce .gt. 0.1d0 ) then
!           --- flexible bond angle
!           --- initialize bsum_try
            bsum_try = 0.0d0

            if ( .not. lnew ) then
!              --- first ibend is the old conformation
               xvecgrow = rxu(i,iugrow) - rxui
               yvecgrow = ryu(i,iugrow) - ryui
               zvecgrow = rzu(i,iugrow) - rzui
               distgrow = bondlen(count)
               distgrow2 = distgrow*distgrow 
!              --- dot product divided by lengths gives cos(angle)
               anglec = ( xvecprev*xvecgrow + yvecprev*yvecgrow 
     &              + zvecprev*zvecgrow ) / (distprev*distgrow)
               angle = dacos(anglec)
               if (L_bend_table) then
                  lengthc2 = lengthFP2 + distgrow2 - 
     &                       2.0d0*lengthFP*distgrow*anglec
                  lengthc = dsqrt(lengthc2)
                  call lininter_bend(lengthc, tabulated_bend, type)
                  vangle = tabulated_bend
                else
                  vangle = kforce * (angle - equil)**2
               end if

               ang_trial(1) = angle
               bfactor(1) = dexp( -beta*vangle )
               bsum_try = bsum_try + bfactor(1)
!              --- skip first ibend in next loop
               start = 2

            else
!              --- new conformation start at 1
               start = 1
            end if
         
            if (L_bend_table) then
                  distgrow = bondlen(count)
                  distgrow2 = distgrow*distgrow
            end if 
!           --- compute trial angles and energies
            do ibend = start,nchben_a
!               --- choose the angle uniformly on sin(theta)
               rsint = 2.0d0*random() - 1.0d0
               angle = dacos(rsint)
               ang_trial(ibend) = angle

!              --- calculate the bond angle energy
               if (L_bend_table) then
                  lengthc2 = lengthFP2 + distgrow2 -
     &                       2.0d0*lengthFP*distgrow*dcos(angle)
                  lengthc = dsqrt(lengthc2)   
                  call lininter_bend(lengthc,tabulated_bend, type)
                  vangle = tabulated_bend
               else
                  vangle = kforce * (angle - equil)**2
               end if
               bfactor(ibend) = dexp(-beta*vangle)
               bsum_try = bsum_try + bfactor(ibend)
            end do

            if ( lnew ) then
!              --- select one of the trial sites via bias
               rbf = random()*bsum_try
               bs = 0.0d0
               do ibend = 1,nchben_a
                  bs = bs + bfactor(ibend)
                  if ( rbf .lt. bs ) then
                     angle = ang_trial(ibend)
                     vangle = dlog(bfactor(ibend))/(-beta)
                     goto 10
                  end if
               end do
 10            continue
            else
!              --- select the old conformation
               angle = ang_trial(1)
               vangle = dlog(bfactor(1))/(-beta)
            end if

!           --- propagate the rosenbluth weight
            wei_bend = wei_bend * bsum_try/dble(nchben_a)

         else
!           --- fixed bond angle
            angle = equil
            vangle = 0.0d0
         end if
         
         bendang(count) = angle

         vbbtr = vbbtr + vangle
         
      end do

! Neeraj  iugrow-iufrom-iuprev bend angle has been selected!

      if ( lnew ) then
!        --- assign phi(angstart) to 0.0
         phi(angstart) = 0.0d0
!        --- skip angstart in the loop below
         start_ang = angstart+1
      else
!        --- set up the cone using iuprev
         xub = -xvecprev/distprev
         yub = -yvecprev/distprev
         zub = -zvecprev/distprev
         
         call cone(1,xub,yub,zub,dum,dum,dum,dum,dum)
         
!        --- need to determine angstart
         start_ang = angstart
      end if
         
!     --- determine the angles of the grown beads relative to anglestart

      do count = start_ang,ntogrow
         iugrow = growlist(iw,count)

         lengtha = bondlen(count)
         lengtha2 = lengtha * lengtha

!        --- initialize bsum_try
         bsum_try = 0.0d0

         if ( .not. lnew ) then
!           --- compute vector from iufrom to iugrow
            xvecgrow = rxu(i,iugrow) - rxui
            yvecgrow = ryu(i,iugrow) - ryui
            zvecgrow = rzu(i,iugrow) - rzui
            distgrow = bondlen(count)

!           --- turn this into a unit vector
            ux = xvecgrow/distgrow
            uy = yvecgrow/distgrow
            uz = zvecgrow/distgrow
            alpha = bendang(count)

            call cone(3,dum,dum,dum,alpha,gamma,ux,uy,uz)

            phitwo = gamma

!     ******************************************************************
!     * if iinit = 2 then it creates a unit vector that has an angle   *
!     * of alpha from the +z direction (previous vector) and an angle  *
!     * of gamma (0,2Pi) around the cone circle and returns this as    *
!     * ux,uy,uz                                                       *
!     ******************************************************************

            call cone(2,dum,dum,dum,alpha,gamma,ux,uy,uz)

            if ( angstart .eq. count ) then
!              --- this is the first phi, no need to compute bias
               vphi = 0.0d0
               goto 20
            else
               vphi = 0.0d0
               do aaa = angstart,count-1
                  iuone = growlist(iw,aaa)
                  lengthb = bondlen(aaa)
                  lengthb2 = lengthb * lengthb
                  phione = phi(aaa)
                  thetaone = bendang(aaa)
                  thetatwo = bendang(count)
                  call coneangle(thetaone,phione,thetatwo,phitwo
     &                 ,angle)
            
                  do ib = 1, inben(imolty,iugrow)
                     iulast = ijben2(imolty,iugrow,ib)
                     if ( iulast .eq. iufrom ) then
                        iu2back = ijben3(imolty,iugrow,ib)
                        if ( iu2back .eq. iuone ) then
                           type = itben(imolty,iugrow,ib)
!                          --- calculate the bond angle energy
                           if (L_bend_table) then
                              lengthc2 = lengtha2 + lengthb2 - 
     &                                 2.0d0*lengtha*lengthb*dcos(angle)
                              lengthc = dsqrt(lengthc2)
                              call lininter_bend(lengthc, 
     &                             tabulated_bend, type)
                              vphi = vphi + tabulated_bend
                           else
                              vphi = vphi + brbenk(type) 
     &                             * (angle - brben(type))**2
                           end if
                        end if
                     end if
                  end do
               end do

               ang_trial(1) = phitwo
               bfactor(1) = dexp( -beta * vphi )
               bsum_try = bsum_try + bfactor(1)

            end if
!           --- skip first ibend for OLD
            start = 2
         else
!           --- perform all ibend for NEW
            start = 1
         end if

!        --- compute trial energies and weights
         do ibend = start,nchben_b
!           --- determine a value of phitwo
            phitwo = random()*twopi
            vphi = 0.0d0
            do aaa = angstart,count-1
               lengthb = bondlen(aaa)
               lengthb2 = lengthb*lengthb
               iuone = growlist(iw,aaa)
               phione = phi(aaa)
               thetaone = bendang(aaa)
               thetatwo = bendang(count)
               call coneangle(thetaone,phione,thetatwo,phitwo
     &              ,angle)
            
               do ib = 1, inben(imolty,iugrow)
                  iulast = ijben2(imolty,iugrow,ib)
                  if ( iulast .eq. iufrom ) then
                     iu2back = ijben3(imolty,iugrow,ib)
                     if ( iu2back .eq. iuone ) then
                        type = itben(imolty,iugrow,ib)
!                       --- calculate the bond angle energy
                        if (L_bend_table) then
                           lengthc2 = lengtha2 + lengthb2 - 
     &                                2.0d0*lengtha*lengthb*dcos(angle)
                           lengthc = dsqrt(lengthc2) 
                           call lininter_bend (lengthc, 
     &                          tabulated_bend, type)
                           vphi = vphi + tabulated_bend
                        else
                           vphi = vphi + brbenk(type) 
     &                          * (angle - brben(type))**2
                        end if
                     end if
                  end if
               end do
            end do

!           --- store the boltzmann factors and phi
            bfactor(ibend) = dexp(-beta*vphi)
            ang_trial(ibend) = phitwo
            bsum_try = bsum_try + bfactor(ibend)
         end do
   
         if ( lnew ) then
!           --- select a value of phitwo in a biased fashion
            rbf = random()*bsum_try
            bs = 0.0d0
            do ibend = 1,nchben_b
               bs = bs + bfactor(ibend)
               if ( rbf .lt. bs ) then
                  phitwo = ang_trial(ibend)
                  vphi = dlog( bfactor(ibend) )/(-beta)
                  goto 15
               end if
            end do
 15         continue
         else
!           --- select the OLD value of phitwo
            phitwo = ang_trial(1)
            vphi = dlog( bfactor(1) )/(-beta)
         end if

!        --- propagate angle weight
         wei_bend = wei_bend * (bsum_try/dble(nchben_b))

 20      continue
!        --- store the angle for thetatwo
         phi(count) = phitwo
         vbbtr = vbbtr + vphi
      end do

       

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     --- End Bond angle biased selection ---           c 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      write(iou,*) 'FINISH GEOMETRY'

      return
      end

      subroutine bendangle(equil,kforce,beta,angle,vangle )


!     *************************************************************
!     *** choose a bond angle                                   ***
!     *** equil is equilibrium bond angle                       ***
!     *** kforce is the force constant                          ***
!     *** beta is 1/kT                                          ***
!     *** angle is the angle returned by this subroutine        ***
!     *** vangle is the energy of the angle                     ***
!     *** M.G. Martin                                           ***
!     *************************************************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none

      real(KIND=double_precision)::equil, kforce, beta, angle, vangle, rr, v1, v2
     &     ,random, tabulated_bend

!      write(2,*) 'start BENDANGLE'

      if ( kforce .gt. 0.1d0 ) then
!        --- find a vector inside the unit sphere
 80      v1 = 2.0d0*random()-1.0d0
         v2 = 2.0d0*random()-1.0d0
         rr = v1*v1 + v2*v2
         if (rr .ge. 1.0d0 ) goto 80
         
!        --- select angle from a gaussian distribution
         angle = equil + v1*dsqrt( (-dlog(rr))
     &        /( kforce*beta*rr) )
         
         if (angle .le. 0.0d0 .or. angle .ge. 3.1415 ) then
            write(2,*) 'chose angle outside of 0,Pi in bendangle'
            goto 80
         end if

!        --- correct for the phase space of the angle
         if ( random() .gt. dsin(angle) ) goto 80
!         if (L_bend_table) then
!            call lininter_bend(angle, tabulated_bend, type)
!         type isn't specified in this subroutine, but
!         I don't think its called anywhere anymore
         vangle = kforce * (angle - equil)**2
      else
!        --- fixed bond angle
         angle = equil
         vangle = 0.0d0
      end if

!      write(2,*) 'end BENDANGLE'

      return
      end
