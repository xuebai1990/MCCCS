      subroutine recipsum(ibox,vrecip)
!    *********************************************************************
!    ** calculates the total reciprocal space ewald-sum term for volume **
!    ** moves, written in 1998 by Bin Chen.                             **
!    ** rewritten in 2001 by Bin Chen.                                  **
!    ** rewritten again, probably by Bin.                               **
!    *********************************************************************
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
!$$$      include 'system.inc'
!$$$      include 'coord.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'poten.inc'
!$$$      include 'conver.inc'
!$$$      include 'cell.inc'
!$$$! RP added for MPI
!$$$      include 'mpif.h'
!$$$      include 'mpi.inc'

      integer(KIND=normal_int)::ibox,i,ii,imolty,ncount

! * from h-matrix formulation
      integer(KIND=normal_int)::l,m,n,m_min,n_min,kmaxl
     & ,kmaxm,kmaxn

      real(KIND=double_precision)::alpsqr4,vol,ksqr,sumr,sumi,arg
     & ,vrecip,bx1,by1,bz1
     & ,hmatik(9),kx1,ky1,kz1,hmaxsq,calpi
!      real(KIND=double_precision)::sum_sumr,sum_sumi

! RP added for calculating time for communication step
      integer(KIND=normal_int)::mystart,myend,blocksize
      integer(KIND=normal_int)::ncount_arr(numprocmax)
     & ,ncount_displs(numprocmax)
      real(KIND=double_precision)::sum_vrecip,kx_arr(vectormax)
     & ,ky_arr(vectormax),kz_arr(vectormax),kx_one(vectormax)
     & ,ky_one(vectormax),kz_one(vectormax),ssumi_arr(vectormax)
     & ,prefact_arr(vectormax),ssumr_one(vectormax),ssumi_one(vectormax)
     & ,prefact_one(vectormax),ssumr_arr(vectormax)
 
! KM for MPI
      do i=1,numprocmax
         ncount_displs(i) = 0
         ncount_arr(i) = 0
      end do
      

! *** Set up the reciprocal space vectors ***

      ncount = 0
      vrecip = 0.0d0

      calpi = calp(ibox)

      if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
         bx1 = boxlx(ibox)
         by1 = boxly(ibox)
         bz1 = boxlz(ibox)
         hmat(ibox,1) = bx1
         hmat(ibox,5) = by1
         hmat(ibox,9) = bz1
         do i = 1,9
            hmatik(i) = 0.0d0
         end do
         hmatik(1) = twopi/bx1
         hmatik(5) = twopi/by1
         hmatik(9) = twopi/bz1
         kmaxl = dint(bx1*calpi)+1
         kmaxm = dint(by1*calpi)+1
         kmaxn = dint(bz1*calpi)+1
      else
         do i = 1,9
            hmatik(i) = twopi*hmati(ibox,i)
         end do
         kmaxl = dint(hmat(ibox,1)*calpi)+2
         kmaxm = dint(hmat(ibox,5)*calpi)+2
         kmaxn = dint(hmat(ibox,9)*calpi)+2
      end if
    
      alpsqr4 = 4.0d0*calpi*calpi
         
      vol = hmat(ibox,1)*
     &     (hmat(ibox,5)*hmat(ibox,9) - hmat(ibox,8)*hmat(ibox,6))
     &     + hmat(ibox,4)*
     &     (hmat(ibox,8)*hmat(ibox,3) - hmat(ibox,2)*hmat(ibox,9))
     &     + hmat(ibox,7)*
     &     (hmat(ibox,2)*hmat(ibox,6) - hmat(ibox,5)*hmat(ibox,3))

      vol = vol/(4.0d0*onepi)

      hmaxsq = alpsqr4*onepi*onepi
! RP added for MPI        
      blocksize = kmaxl/numprocs
       
      mystart = myid * blocksize
      if(myid .eq. (numprocs-1))then
         myend = kmaxl
      else 
         myend = (myid + 1) * blocksize - 1
      end if
! *** generate the reciprocal-space
! here -kmaxl,-kmaxl+1,...,-1 are skipped, so no need to divide by 2 for the prefactor
      do l = mystart,myend
!      do l = 0,kmaxl 
        if ( l .eq. 0 ) then
            m_min = 0
         else
            m_min = -kmaxm
         end if
         do m = m_min, kmaxm
            if (l .eq. 0 .and. m .eq. 0) then
               n_min = 1
            else
               n_min = -kmaxn
            end if
            do n = n_min, kmaxn
               kx1 = dble(l)*hmatik(1)+dble(m)*hmatik(2)+
     &              dble(n)*hmatik(3)
               ky1 = dble(l)*hmatik(4)+dble(m)*hmatik(5)+
     &              dble(n)*hmatik(6)
               kz1 = dble(l)*hmatik(7)+dble(m)*hmatik(8)+
     &              dble(n)*hmatik(9)
               ksqr = kx1*kx1+ky1*ky1+kz1*kz1
!               if ( ksqr .lt. hmaxsq ) then
! --- sometimes these are about equal, which can cause different
! --- behavior on 32 and 64 bit machines without this .and. statement
               if ( ksqr .lt. hmaxsq .and.
     &              abs(ksqr-hmaxsq) .gt. 1d-9 ) then
                  ncount = ncount + 1
                  kx_arr(ncount) = kx1
                  ky_arr(ncount) = ky1
                  kz_arr(ncount) = kz1
                  prefact_arr(ncount) =
     &                 exp(-ksqr/alpsqr4)/(ksqr*vol)
!     *** sum up q*cos and q*sin ***
                  sumr = 0.0d0
                  sumi = 0.0d0
!                  do i = myid+1,nchain,numprocs
                  do i = 1,nchain
                     imolty = moltyp(i)
                     if (.not.lelect(imolty).or.nboxi(i).ne.ibox ) cycle
                     do ii = 1,nunit(imolty)
                        if ( lqchg(ntype(imolty,ii)) ) then
                           arg=kx1*rxu(i,ii)+ky1*ryu(i,ii)+kz1*rzu(i
     &                          ,ii)
                           sumr = sumr + cos(arg)*qqu(i,ii)
                           sumi = sumi + sin(arg)*qqu(i,ii)
                        end if
                     end do
                  end do

                  ssumr_arr(ncount) = sumr
                  ssumi_arr(ncount) = sumi
! *** Potential energy ***
                  vrecip = vrecip + (sumr*sumr + sumi*sumi)
     &                 * prefact_arr(ncount)
               end if
            end do
         end do
      end do
    
       CALL MPI_ALLREDUCE(vrecip,sum_vrecip,1,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,
     &     ierr)
 
      vrecip = sum_vrecip

       CALL MPI_ALLGATHER(ncount,1,MPI_INTEGER,ncount_arr,
     &       1,MPI_INTEGER,MPI_COMM_WORLD,ierr)

       ncount_displs(1) = 0
       do i = 2,numprocs
           ncount_displs(i) = ncount_displs(i-1) + ncount_arr(i-1)
       end do
 
      CALL MPI_ALLGATHERV(kx_arr,ncount,MPI_DOUBLE_PRECISION,kx_one,
     &         ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &         MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHERV(ky_arr,ncount,MPI_DOUBLE_PRECISION,ky_one,
     &        ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHERV(kz_arr,ncount,MPI_DOUBLE_PRECISION,
     &       kz_one,
     &       ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &       MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHERV(ssumr_arr,ncount,MPI_DOUBLE_PRECISION,
     &          ssumr_one,
     &      ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &      MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHERV(ssumi_arr,ncount,MPI_DOUBLE_PRECISION,
     &          ssumi_one,
     &      ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &      MPI_COMM_WORLD,ierr)
       CALL MPI_ALLGATHERV(prefact_arr,ncount,MPI_DOUBLE_PRECISION,
     &           prefact_one,
     &      ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &      MPI_COMM_WORLD,ierr)       

      ncount = 0
      do i = 1,numprocs
        ncount = ncount + ncount_arr(i)
       end do
      do i = 1, ncount
        kx(i,ibox) = kx_one(i)
        ky(i,ibox) = ky_one(i)
        kz(i,ibox) = kz_one(i)
        ssumr(i,ibox) = ssumr_one(i)
        ssumi(i,ibox) = ssumi_one(i)
        prefact(i,ibox) = prefact_one(i)
      end do

      vrecip = vrecip*qqfact
!      write(iou,*) 'in recipsum:',ssumr(100,ibox),ibox
! *** safety check ***
!      write(iou,*) 'A total of ',ncount,' vectors are used'
      if ( ncount .gt. vectormax ) call cleanup
     &     ('choose a larger vectormax')

      numvect(ibox) = ncount
      return

      end
