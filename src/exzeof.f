      function exzeof(i,j,k,idi)

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
!$$$      include 'grid.inc'
!$$$      include 'zeopoten.inc'
!$$$      include 'zeolite.inc'
!$$$      include 'control.inc'
!$$$      include 'external.inc'
!$$$      include 'system.inc'
!$$$      include 'cell.inc'
!$$$      include 'poten.inc'      
!$$$      include 'ewaldsum.inc'
!$$$      include 'mpif.h'
!$$$      include 'mpi.inc'
      real(KIND=double_precision)::exzeof,xi,yi,zi,r2,rcutsq,xr,yr,zr,r2i,r6,vljnew,vqnew
     &     ,overflow=1.0d8,rminsq,r,recipzeo,erfunc,i,j,k,sx,sy,sz
      integer(KIND=normal_int)::izeo,idi,idj,ntij,layer,ii,jj,kk,ibox=1

!     calculate the Lennard-Jones interactions, include as many layers
!     of neighboring unit cells as needed for the specified precision

      exzeof=0.
      rminsq=rmin*rmin
      
      if (ltailc.and.lij(idi)) then
         vljnew=eps+1.
         layer=0
         do while (abs(vljnew).gt.eps)
            if (layer.gt.nlayermax) nlayermax=layer
            vljnew=0.
            do izeo=1,nzeo
               if (lunitcell(izeo)) then
                  idj=idzeo(izeo)
                  ntij = (idi - 1) * nntype + idj
                  do ii=-layer,layer
                     do jj=-layer,layer
                        do kk=-layer,layer
                           if (abs(ii).eq.layer .or. abs(jj).eq.layer
     &                          .or.abs(kk).eq.layer) then

                              sx = zeox(izeo)*hmati(ibox,1)+zeoy(izeo)
     &                             *hmati(ibox,4)+zeoz(izeo)*hmati(ibox
     &                             ,7)+dble(ii-i)/nx
                              sy = zeox(izeo)*hmati(ibox,2)+zeoy(izeo)
     &                             *hmati(ibox,5)+zeoz(izeo)*hmati(ibox
     &                             ,8)+dble(jj-j)/ny
                              sz = zeox(izeo)*hmati(ibox,3)+zeoy(izeo)
     &                             *hmati(ibox,6)+zeoz(izeo)*hmati(ibox
     &                             ,9)+dble(kk-k)/nz

                              xr=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz
     &                             *hmat(ibox,7)
                              yr=sy*hmat(ibox,5)+sz*hmat(ibox,8)
                              zr=sz*hmat(ibox,9)

                              r2=xr*xr+yr*yr+zr*zr
                              if (r2.le.rminsq) then
                                 exzeof=overflow
                                 return
                              end if
                              r2i=sig2ij(ntij)/r2
                              r6=r2i*r2i*r2i
                              if (lshift) then
                                 vljnew=vljnew+4.*(epsij(ntij)*(r6-1.0)
     &                                *r6-ecut(ntij))
                              else
                                 vljnew=vljnew+4.*epsij(ntij)*(r6-1.0)
     &                                *r6
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
            exzeof=exzeof+vljnew
            layer=layer+1
         end do
      end if

      if ((.not.ltailc.and.lij(idi)).or.(lqchg(idi))) then
         vljnew=0.
         vqnew=0.
         rcutsq = rcut(ibox)*rcut(ibox)
         sx=dble(i)/nx
         sy=dble(j)/ny
         sz=dble(k)/nz
         xi=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
         yi=sy*hmat(ibox,5)+sz*hmat(ibox,8)
         zi=sz*hmat(ibox,9)
         do izeo=1,nzeo
            idj=idzeo(izeo)
            ntij = (idi - 1) * nntype + idj
            xr=xi-zeox(izeo)
            yr=yi-zeoy(izeo)
            zr=zi-zeoz(izeo)
            call mimage(xr,yr,zr,ibox)
            r2=xr*xr+yr*yr+zr*zr
            if (r2.le.rminsq) then
               exzeof=overflow
               return
            end if
            if (r2 .lt. rcutsq) then
               if (.not.ltailc.and.lij(idi)) then
                  r2i=sig2ij(ntij)/r2
                  r6=r2i*r2i*r2i
                  if (lshift) then     
                     vljnew=vljnew+4.*(epsij(ntij)*(r6-1.0)*r6
     &                    -ecut(ntij))
                  else
                     vljnew=vljnew+4.*epsij(ntij)*(r6-1.0)*r6
                  end if               
               end if
               if (lqchg(idi)) then
                  r=dsqrt(r2)
                  if (lewald) then
                     vqnew=vqnew+qelect(idi)*qelect(idj)
     &                    *erfunc(calp(ibox)*r)/r
                  else
                     vqnew=vqnew+qelect(idi)*qelect(idj)/r
                  end if
               end if
            end if
         end do

         if (lqchg(idi).and.lewald) vqnew=vqnew+recipzeo(xi,yi,zi
     &        ,qelect(idi))

         exzeof=exzeof+vljnew+vqnew*qqfact
      end if

      return
      end
 
      function recipzeo(xi,yi,zi,qi)

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
!$$$      include 'zeolite.inc'
!$$$      include 'system.inc'
!$$$      include 'coord.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'poten.inc'
!$$$      include 'conver.inc'
!$$$      include 'cell.inc'
!$$$      include 'mpif.h'
!$$$      include 'mpi.inc'

      integer(KIND=normal_int)::ibox=1,i,ii

! * from h-matrix formulation
      integer(KIND=normal_int)::l,m,n,m_min,m_max,n_min,n_max,kmaxl,kmaxm,kmaxn
      integer(KIND=normal_int)::mystart,myend,blocksize
      real(KIND=double_precision)::alpsqr4,vol,ksqr,sum,arg,recipzeo,xi,yi,zi,qi
     &     ,hmatik(9),kx1,ky1,kz1,hmaxsq,calpi

! *** Set up the reciprocal space vectors ***

      recipzeo = 0.0d0

      calpi = calp(ibox)

      do i = 1,9
         hmatik(i) = twopi*hmati(ibox,i)
      end do
      if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
         kmaxl = dint(hmat(ibox,1)*calpi)+1
         kmaxm = dint(hmat(ibox,5)*calpi)+1
         kmaxn = dint(hmat(ibox,9)*calpi)+1
      else
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

      vol = vol/(8.0d0*onepi)

      hmaxsq = alpsqr4*onepi*onepi

! *** generate the reciprocal-space
! here -kmaxl,-kmaxl+1,...,-1 are skipped, so no need to divide by 2 for the prefactor
!      do l = myid,kmaxl,numprocs       
      do l = 0,kmaxl 
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
!     *** sum up q*cos and q*sin ***
                  sum = 0.0d0
                  do i = 1,nzeo
                     arg=kx1*(xi-zeox(i))+ky1*(yi-zeoy(i))+kz1*(zi
     &                    -zeoz(i))
                     sum=sum+qelect(idzeo(i))*cos(arg)
                  end do
                  recipzeo=recipzeo+sum*exp(-ksqr/alpsqr4)/ksqr
! *** Potential energy ***
               end if
            end do
         end do
      end do

!      CALL MPI_ALLREDUCE(recipzeo,sumrecipzeo,1,MPI_DOUBLE_PRECISION
!     &     ,MPI_SUM,MPI_COMM_WORLD,ierr)
!      recipzeo = sumrecipzeo

      recipzeo=recipzeo*qi/vol

      return

      end
