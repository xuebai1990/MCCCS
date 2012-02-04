      subroutine recip(ibox,vrecipnew,vrecipold,type)

c recip
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
c John Stubbs, and Collin Wick and Ilja Siepmann  
c                     
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to 
c
c Free Software Foundation, Inc. 
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c    *********************************************************************
c    ** calculates the reciprocal ewald-sum term for trans, rot, flucq, **
c    ** swatch and swap moves, and update the reciprocal ewald-sum.     **
c    ** rewritten on June 25/99 by Bin Chen.                            **
c    *********************************************************************

      implicit none
      integer ic,zz,ii,imolty,ibox,ncount,type
      double precision vrecipnew,vrecipold,sumr(2),sumi(2),arg
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ewaldsum.inc'
      include 'poten.inc'
c RP added for MPI     
      include 'mpif.h'
      include 'mpi.inc'
       
c RP added for MPI
      integer mystart,myend,blocksize,i
      double precision ssumrn_arr(vectormax),ssumin_arr(vectormax),
     &      ssumrn_one(vectormax),ssumin_one(vectormax)
      integer countn,ncount_displs(numprocmax),ncount_arr(numprocmax)   
    
!      if (LSOLPAR.and.(ibox.eq.2))then
!        return
!      endif

c KM for MPI
      do i=1,numprocmax
         ncount_displs(i) = 0
         ncount_arr(i) = 0
      enddo

      ncount = numvect(ibox)
      if ( type .eq. 1 ) then
         
c *** recalculate the reciprocal space part for one-particle move, translation,
c *** rotation, swap, flucq, and swatch.
c *** old conformation zz = 1 (which is 0 for swap inserted molecule)
c *** new conformation zz = 2 (which is 0 for swap removed molecule)

c         write(iou,*) 'in recip:',moltion(1),moltion(2)
c         do zz = 1,2
c            imolty = moltion(zz)
c            do ii = 1, nunit(imolty)
c               write(iou,*) rxuion(ii,zz),ryuion(ii,zz),rzuion(ii,zz),
c     &              qquion(ii,zz)
c            enddo
c         enddo
      
c RP added for MPI
         blocksize = ncount/numprocs
         mystart = myid * blocksize + 1
         if(myid .eq. (numprocs-1))then
            myend = ncount
         else 
            myend = (myid + 1) * blocksize
         endif
         countn = myend - mystart + 1
         do 30 ic = mystart,myend
!         do 30 ic = 1,ncount
            do 20 zz = 1,2
c --- zz = 1: old configuration 
c --- zz = 2: new configuration

               sumr(zz) = 0.0d0
               sumi(zz) = 0.0d0
               imolty = moltion(zz)
               do ii = 1, nunit(imolty)
                  if ( lqchg(ntype(imolty,ii)) ) then
                     arg = kx(ic,ibox)*rxuion(ii,zz) +
     &                    ky(ic,ibox)*ryuion(ii,zz) +
     &                    kz(ic,ibox)*rzuion(ii,zz)
                     sumr(zz) = sumr(zz) + 
     &                    qquion(ii,zz)*dcos(arg)
                     sumi(zz) = sumi(zz) + 
     &                    qquion(ii,zz)*dsin(arg)
                  endif
               enddo
 20         continue
!          ssumrn(ic,ibox) = ssumr(ic,ibox) - sumr(1)
!     &           + sumr(2)
!            ssumin(ic,ibox) = ssumi(ic,ibox) - sumi(1)
!     &           + sumi(2)


c RP added for MPI
            ssumrn_arr(ic-mystart + 1) = ssumr(ic,ibox) - sumr(1)
     &           + sumr(2)
            ssumin_arr(ic-mystart + 1) = ssumi(ic,ibox) - sumi(1)
     &           + sumi(2)

 30      continue
        
          CALL MPI_ALLGATHER(countn,1,MPI_INTEGER,ncount_arr,
     &       1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
         
         ncount_displs(1) = 0
         do i = 2,numprocs
           ncount_displs(i) = ncount_displs(i-1) + ncount_arr(i-1)    
         enddo

         CALL MPI_ALLGATHERV(ssumrn_arr,ncount_arr(myid + 1),
     &     MPI_DOUBLE_PRECISION,
     &     ssumrn_one,ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &     MPI_COMM_WORLD,ierr)
     
        CALL MPI_ALLGATHERV(ssumin_arr,ncount_arr(myid+1),
     &     MPI_DOUBLE_PRECISION,
     &     ssumin_one,ncount_arr,ncount_displs,MPI_DOUBLE_PRECISION,
     &     MPI_COMM_WORLD,ierr)
     
       do i = 1,ncount
        ssumrn(i,ibox) = ssumrn_one(i)
        ssumin(i,ibox) = ssumin_one(i)
      enddo
c----------------------------------------------------------------------    
         vrecipnew = 0.0d0
         vrecipold = 0.0d0
         do ic = 1,ncount
            vrecipnew = vrecipnew + (ssumrn(ic,ibox)*
     &           ssumrn(ic,ibox) + ssumin(ic,ibox)*
     &           ssumin(ic,ibox))*prefact(ic,ibox)
            vrecipold = vrecipold + (ssumr(ic,ibox)*
     &           ssumr(ic,ibox) + ssumi(ic,ibox)*
     &           ssumi(ic,ibox))*prefact(ic,ibox)
         enddo
       
         vrecipnew = vrecipnew*qqfact
         vrecipold = vrecipold*qqfact

      elseif (type .eq. 2) then

c *** update the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumrn(ic,ibox)
            ssumi(ic,ibox) = ssumin(ic,ibox)
         enddo

      elseif (type .eq. 3) then

c *** store the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumro(ic,ibox) = ssumr(ic,ibox)
            ssumio(ic,ibox) = ssumi(ic,ibox)
         enddo

      elseif (type .eq. 4) then

c *** restore the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumr(ic,ibox) = ssumro(ic,ibox)
            ssumi(ic,ibox) = ssumio(ic,ibox)
         enddo

      endif

c      write(iou,*) 'in recip:',ssumr(100,ibox),ibox,ssumrn(100,ibox)

      return 
      end

