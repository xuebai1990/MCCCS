      subroutine inclus( inclnum,inclmol,inclbead,inclsign,ncarbon, ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2)

      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_runtime,only:err_exit
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'common.inc'
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'connect.inc'
!$$$      include 'poten.inc'

      integer(KIND=normal_int)::m,n,nb,mb,imolty,ioffset
      integer(KIND=normal_int)::inclnum,inclmol,inclbead,inclsign ,ncarbon
      integer(KIND=normal_int)::ainclnum,ainclmol,ainclbead,a15t
      dimension inclmol(ntmax*numax*numax),inclsign(ntmax*numax*numax)
      dimension inclbead(ntmax*numax*numax,2),ncarbon(ntmax)
      dimension ainclmol(ntmax*numax*numax)
      dimension ainclbead(ntmax*numax*numax,2)
      dimension a15t(ntmax*numax*numax)
      
! -- variables added (3/24/05) for variable 1-4 interactions 	
      real(KIND=double_precision)::ofscale,ofscale2
      dimension ofscale(ntmax*numax*numax),ofscale2(ntmax*numax*numax)

! ----------------------------------------------------------------

!c !!!!This is modified to work only for TATB NR-2007!!!


! - triple loop over all types of molecules -
      do imolty = 1, nmolty

         do m = 1, nunit(imolty)
            do n = 1, nunit(imolty)
               linclu(imolty,m,n) = .true.
               lqinclu(imolty,m,n) = .true.	      
! * by default, dont want any 1-5 r^12 interactions
               lainclu(imolty,m,n) = .false.
	       ljscale(imolty,m,n) = 1.0
	       qscale2(imolty,m,n) = 1.0
            end do
         end do
         
! - double loop over all units -
         do m = 1, nunit(imolty)

! - exclude all self interactions -
            linclu(imolty,m,m) = .false.
            lqinclu(imolty,m,m) = .false.
! - exclude all directly bonded beads (vibrations) -
            do n = 1, invib(imolty,m)
               nb = ijvib(imolty,m,n)
               linclu(imolty,m,nb) = .false.
               lqinclu(imolty,m,nb) = .false.
            end do

! - exclude carbons around a quaternary center for explct
            if (invib(imolty,m) .eq. 4) then
               do n = 1,4
                  do nb = 1,4
                     linclu(imolty,ijvib(imolty,m,n) ,ijvib(imolty,m,nb))=.false.
                     lqinclu(imolty,ijvib(imolty,m,n) ,ijvib(imolty,m,nb))=.false.
                  end do
               end do
            end if

! - exclude all next-nearest neighbor bonded beads (bending) -
            do n = 1, inben(imolty,m)
               nb = ijben3(imolty,m,n)
               linclu(imolty,m,nb) = .false.
               lqinclu(imolty,m,nb) = .false.
            end do
! - exclude all third-nearest neighbor bonded beads (torsions) -
            do n = 1, intor(imolty,m)
               nb = ijtor4(imolty,m,n)
               linclu(imolty,m,nb) = .false.
! * dont set lqinclu since we want 1-4 interactions, unless 1q14scale is F
               if (.not.lq14scale(imolty)) then
                  lqinclu(imolty,m,nb) = .false.
               else
                  qscale2(imolty,m,nb) = qscale(imolty) 
                  qscale2(imolty,nb,m) = qscale(imolty)
               end if
            end do

         end do

!     - include or exclude additional beads accoring to incl
         do n = 1,inclnum
            if ( inclmol(n) .eq. imolty ) then
               m = inclbead(n,1)
               nb = inclbead(n,2)
               if ( inclsign(n) .eq. 1 ) then
                  linclu(imolty,m,nb) = .true.
                  linclu(imolty,nb,m) = .true.
                  lqinclu(imolty,m,nb) = .true.
                  lqinclu(imolty,nb,m) = .true.
		  ljscale(imolty,m,nb) = ofscale(n)
		  ljscale(imolty,nb,m) = ofscale(n)
		  qscale2(imolty,m,nb) = ofscale2(n)
		  qscale2(imolty,nb,m) = ofscale2(n)
               elseif (inclsign(n) .eq. -1 ) then
                  linclu(imolty,m,nb) = .false.
                  linclu(imolty,nb,m) = .false.
                  lqinclu(imolty,m,nb) = .false.
                  lqinclu(imolty,nb,m) = .false.
               else
                  write(iou,*) 'INCLUS: n,inclsign(n)',n,inclsign(n)
                  call err_exit('inclusign must be 1 or -1')
               end if
            end if

         end do

! * add in 1-5 interactions according to aincl
         do m = 1,ainclnum
            if ( ainclmol(m) .eq. imolty ) then
               mb = ainclbead(m,1)
               nb = ainclbead(m,2)
                  lainclu(imolty,mb,nb) = .true.
                  lainclu(imolty,nb,mb) = .true.
                  a15type(imolty,mb,nb) = a15t(m)
                  a15type(imolty,nb,mb) = a15t(m)
            end if
         end do

! - exclude all hydrogens that have their carbons excluded 
         if ( ncarbon(imolty) .lt. nunit(imolty) ) then
            if (ncarbon(imolty) .eq. 3 .and. nunit(imolty) .eq. 8) then
! - ethane with bead 3 being hydrogen
               ioffset = 0
            else
               ioffset = 1
            end if
            do m = ncarbon(imolty)+ioffset,nunit(imolty)
! - hydrogens only have one vibration and that is to the C atom               
               mb = ijvib(imolty,m,1)
               do nb = 1, ncarbon(imolty)
                  if ( .not. linclu(imolty,mb,nb) ) then
                     linclu(imolty,m,nb) = .false.
                     linclu(imolty,nb,m) = .false.
                     lqinclu(imolty,m,nb) = .false.
                     lqinclu(imolty,nb,m) = .false.
                  end if
               end do
               do n=m+1,nunit(imolty)
                  nb = ijvib(imolty,n,1)
                  linclu(imolty,n,nb) = .false.
                  linclu(imolty,nb,n) = .false.
                  lqinclu(imolty,n,nb) = .false.
                  lqinclu(imolty,nb,n) = .false.
                  if ( .not. linclu(imolty,mb,nb) ) then
                     linclu(imolty,m,n) = .false.
                     linclu(imolty,n,m) = .false.
                     linclu(imolty,m,nb) = .false.
                     linclu(imolty,nb,m) = .false.
                     linclu(imolty,n,mb) = .false.
                     linclu(imolty,mb,n) = .false.
                     lqinclu(imolty,m,n) = .false.
                     lqinclu(imolty,n,m) = .false.
                     lqinclu(imolty,m,nb) = .false.
                     lqinclu(imolty,nb,m) = .false.
                     lqinclu(imolty,n,mb) = .false.
                     lqinclu(imolty,mb,n) = .false.
                  end if
               end do
            end do
         end if

! * exclude charge interactions if lqchg is false
         do m = 1,nunit(imolty)
            do n = m+1,nunit(imolty)
               if (.not. lqchg(ntype(imolty,m)) .or. .not. lqchg(ntype(imolty,n))) then
                  lqinclu(imolty,m,n) = .false.
                  lqinclu(imolty,n,m) = .false.
               end if
            end do
         end do

! * self consistency check *

         do m = 1,nunit(imolty)
            do n = m+1,nunit(imolty)

               if (linclu(imolty,m,n) .neqv. linclu(imolty,n,m)) then
                  linclu(imolty,m,n) = .false.
                  linclu(imolty,n,m) = .false.
               end if

               if (lqinclu(imolty,m,n) .neqv. lqinclu(imolty,n,m)) then
                  lqinclu(imolty,m,n) = .false.
                  lqinclu(imolty,n,m) = .false.
               end if

            end do
         end do

!! Removing this part for TATB (NR-2007)

         if (lrigid(imolty)) then
! - dont include rigid beads

! - there will be no intramolecular forces between rigid beads
! - or beads connected one away from a rigid bead          
            do m = 1, invib(imolty,riutry(imolty,1))
               mb = ijvib(imolty,riutry(imolty,1),m)
               do n = riutry(imolty,1), nunit(imolty)
                  linclu(imolty,n,mb) = .false.
                  linclu(imolty,mb,n) = .false.
                  lqinclu(imolty,n,mb) = .false.
                  lqinclu(imolty,mb,n) = .false.
               end do
            end do

 
            do m = riutry(imolty,1), nunit(imolty)
               do n = riutry(imolty,1), nunit(imolty)
                  linclu(imolty,m,n) = .false.
                  lqinclu(imolty,m,n) = .false.
               end do
            end do
         end if

         if (myid.eq.0) then
            write(iou,*) 
            write(iou,*) 'INCLUSION TABLE'

            do m = 1, nunit(imolty)
               write(iou,*) m, (linclu(imolty,m,n),n=1,nunit(imolty))
            end do
            write(iou,*) 

            write(iou,*) 
            write(iou,*) 'CHARGE INCLUSION TABLE'

            do m = 1, nunit(imolty)
               write(iou,*) m,  (lqinclu(imolty,m,n),n=1,nunit(imolty))
            end do

!  400  format (<nunit(imolty)> F5.2)
            write(iou,*) 
	 
            write(iou,*) '1-4 LJ SCALING FACTORS'
            do m = 1, nunit(imolty)
               write(iou,*) m,  (ljscale(imolty,m,n),n=1,nunit(imolty))
            end do
! 500  format (i5,<nunit(imolty)> F5.2)
	 
            write(iou,*)
            write(iou,*) '1-4 CHARGE SCALING FACTORS'
            do m = 1, nunit(imolty)
               write(iou,*) m,  (qscale2(imolty,m,n),n=1,nunit(imolty))
            end do
         end if
! 600  format (i5,<nunit(imolty)> F5.2)

! * not really that important to write out
!         write(iou,*) 
!         write(iou,*) '1-5 OH INTERACTION TABLE'
!         do m = 1, nunit(imolty)
!            write(iou,*) m, (lainclu(imolty,m,n),n=1,nunit(imolty))
!         end do
!         write(iou,*) 

      end do
!      write(iou,*) 'finished inclus'
      return
      end






