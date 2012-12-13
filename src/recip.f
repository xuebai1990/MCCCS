      subroutine recip(ibox,vrecipnew,vrecipold,type)

!    *********************************************************************
!    ** calculates the reciprocal ewald-sum term for trans, rot, flucq, **
!    ** swatch and swap moves, and update the reciprocal ewald-sum.     **
!    ** rewritten on June 25/99 by Bin Chen.                            **
!    *********************************************************************

      use sim_system
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
!$$$      include 'coord.inc'
!$$$      include 'coord2.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'poten.inc'
!$$$! RP added for MPI     
!$$$      include 'mpif.h'
!$$$      include 'mpi.inc'

      integer(KIND=normal_int)::ic,izz,ii,imolty,ibox,ncount,type
      real(KIND=double_precision)::vrecipnew,vrecipold,sumr(2),sumi(2) ,arg       
! RP added for MPI
      integer(KIND=normal_int)::mystart,myend,blocksize,i
      real(KIND=double_precision)::ssumrn_arr(vectormax) ,ssumin_arr(vectormax),ssumrn_one(vectormax) ,ssumin_one(vectormax)
      integer(KIND=normal_int)::countn,ncount_displs(numprocmax) ,ncount_arr(numprocmax)   
    
!      if (LSOLPAR.and.(ibox.eq.2))then
!        return
!      end if

! KM for MPI
      do i=1,numprocmax
         ncount_displs(i) = 0
         ncount_arr(i) = 0
      end do

      ncount = numvect(ibox)
      if ( type .eq. 1 ) then
         
! *** recalculate the reciprocal space part for one-particle move, translation,
! *** rotation, swap, flucq, and swatch.
! *** old conformation izz = 1 (which is 0 for swap inserted molecule)
! *** new conformation izz = 2 (which is 0 for swap removed molecule)

!         write(iou,*) 'in recip:',moltion(1),moltion(2)
!         do izz = 1,2
!            imolty = moltion(izz)
!            do ii = 1, nunit(imolty)
!               write(iou,*) rxuion(ii,izz),ryuion(ii,izz),rzuion(ii,izz),
!     &              qquion(ii,izz)
!            end do
!         end do
      
! RP added for MPI
         blocksize = ncount/numprocs
         mystart = myid * blocksize + 1
         if(myid .eq. (numprocs-1))then
            myend = ncount
         else 
            myend = (myid + 1) * blocksize
         end if
         countn = myend - mystart + 1
         do 30 ic = mystart,myend
!         do 30 ic = 1,ncount
            do 20 izz = 1,2
! --- izz = 1: old configuration 
! --- izz = 2: new configuration

               sumr(izz) = 0.0d0
               sumi(izz) = 0.0d0
               imolty = moltion(izz)
               do ii = 1, nunit(imolty)
                  if ( lqchg(ntype(imolty,ii)) ) then
                     arg = kx(ic,ibox)*rxuion(ii,izz) + ky(ic,ibox)*ryuion(ii,izz) + kz(ic,ibox)*rzuion(ii,izz)
                     sumr(izz) = sumr(izz) +  qquion(ii,izz)*dcos(arg)
                     sumi(izz) = sumi(izz) +  qquion(ii,izz)*dsin(arg)
                  end if
               end do
 20         continue
!          ssumrn(ic,ibox) = ssumr(ic,ibox) - sumr(1)
!     &           + sumr(2)
!            ssumin(ic,ibox) = ssumi(ic,ibox) - sumi(1)
!     &           + sumi(2)


! RP added for MPI
            ssumrn_arr(ic-mystart + 1) = ssumr(ic,ibox) - sumr(1) + sumr(2)
            ssumin_arr(ic-mystart + 1) = ssumi(ic,ibox) - sumi(1) + sumi(2)

 30      continue
        
          CALL MPI_ALLGATHER(countn,1,MPI_INTEGER,ncount_arr, 1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         
         ncount_displs(1) = 0
         do i = 2,numprocs
           ncount_displs(i) = ncount_displs(i-1) + ncount_arr(i-1)    
         end do

         CALL MPI_ALLGATHERV(ssumrn_arr,ncount_arr(myid + 1), MPI_DOUBLE_PRECISION, ssumrn_one,ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
     
        CALL MPI_ALLGATHERV(ssumin_arr,ncount_arr(myid+1), MPI_DOUBLE_PRECISION, ssumin_one,ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION, MPI_COMM_WORLD,ierr)
     
       do i = 1,ncount
        ssumrn(i,ibox) = ssumrn_one(i)
        ssumin(i,ibox) = ssumin_one(i)
      end do
!----------------------------------------------------------------------    
         vrecipnew = 0.0d0
         vrecipold = 0.0d0
         do ic = 1,ncount
            vrecipnew = vrecipnew + (ssumrn(ic,ibox)* ssumrn(ic,ibox) + ssumin(ic,ibox)* ssumin(ic,ibox))*prefact(ic,ibox)
            vrecipold = vrecipold + (ssumr(ic,ibox)* ssumr(ic,ibox) + ssumi(ic,ibox)* ssumi(ic,ibox))*prefact(ic,ibox)
         end do
       
         vrecipnew = vrecipnew*qqfact
         vrecipold = vrecipold*qqfact

      elseif (type .eq. 2) then

! *** update the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumrn(ic,ibox)
            ssumi(ic,ibox) = ssumin(ic,ibox)
         end do

      elseif (type .eq. 3) then

! *** store the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumro(ic,ibox) = ssumr(ic,ibox)
            ssumio(ic,ibox) = ssumi(ic,ibox)
         end do

      elseif (type .eq. 4) then

! *** restore the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumr(ic,ibox) = ssumro(ic,ibox)
            ssumi(ic,ibox) = ssumio(ic,ibox)
         end do

      end if

!      write(iou,*) 'in recip:',ssumr(100,ibox),ibox,ssumrn(100,ibox)

      return 
      end

