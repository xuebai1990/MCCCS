      function exzeof(xi,yi,zi,idi)

      use grid
      implicit none 
      real(8)::exzeof,xi,yi,zi,r2,rcutsq,xr,yr,zr,r2i,r6,vljnew,vqnew
     &     ,overflow=10.0_8**15,rminsq,r,recipzeo,erfunc
      integer::j,idi,idj,ntij,layer,ii,jj,kk
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'control.inc'
      include 'external.inc'
      include 'system.inc'
      include 'poten.inc'      
      include 'ewaldsum.inc'
      include 'mpif.h'
      include 'mpi.inc'

!     calculate the Lennard-Jones and Coulombic interactions, include as
!     many layers of neighboring unit cells as needed for the specified
!     precision

      exzeof=0.
      rminsq=rmin*rmin
      
      if (ltailc.and.lij(idi)) then
         vljnew=100
         layer=0
         do while (vljnew.gt.eps)
            if (layer.gt.nlayermax) nlayermax=layer
            vljnew=0.
            do j=1,nzeo
               if (zeox(j).lt.zunitx .and. zeoy(j).lt.zunity .and.
     &              zeoz(j).lt.zunitz) then
                  idj=idzeo(j)
                  ntij = (idi - 1) * nntype + idj
                  do ii=-layer,layer
                     do jj=-layer,layer
                        do kk=-layer,layer
                           if (abs(ii).eq.layer .or. abs(jj).eq.layer
     &                          .or.abs(kk).eq.layer) then
                              xr=zeox(j)+ii*zunitx-xi
                              yr=zeoy(j)+jj*zunity-yi
                              zr=zeoz(j)+kk*zunitz-zi
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
         rcutsq = rcut(1)*rcut(1)
         do j=1,nzeo
            idj=idzeo(j)
            ntij = (idi - 1) * nntype + idj
            xr=xi-zeox(j)
            xr=xr-zeorx*anint(xr*zeorxi)
            yr=yi-zeoy(j)
            yr=yr-zeory*anint(yr*zeoryi)
            zr=zi-zeoz(j)
            zr=zr-zeorz*anint(zr*zeorzi)
            r2=xr*xr+yr*yr+zr*zr
            if (r2.le.rminsq) then
               exzeof=overflow
               return
            else if (r2 .lt. rcutsq) then
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
                     vqnew=vqnew+qelect(idi)*qelect(idj)*erfunc(calp(1)
     &                    *r)/r
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

      implicit none

      include 'control.inc'
      include 'zeolite.inc'
      include 'system.inc'
      include 'coord.inc'
      include 'ewaldsum.inc'
      include 'poten.inc'
      include 'conver.inc'
      include 'cell.inc'
      include 'mpif.h'
      include 'mpi.inc'

      integer::ibox=1,i,ii

! * from h-matrix formulation
      integer::l,m,n,m_min,m_max,n_min,n_max,kmaxl,kmaxm,kmaxn
      integer::mystart,myend,blocksize
      real(8)::alpsqr4,vol,ksqr,sum,arg,recipzeo,xi,yi,zi,qi
     &     ,hmatik(9),kx1,ky1,kz1,hmaxsq,calpi

! *** Set up the reciprocal space vectors ***

      recipzeo = 0.0d0
      calpi = calp(ibox)

      if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
         hmat(ibox,1) = boxlx(ibox)
         hmat(ibox,5) = boxly(ibox)
         hmat(ibox,9) = boxlz(ibox)
         do i = 1,9
            hmatik(i) = 0.0d0
         end do
         hmatik(1) = twopi/hmat(ibox,1)
         hmatik(5) = twopi/hmat(ibox,5)
         hmatik(9) = twopi/hmat(ibox,9)
         kmaxl = dint(hmat(ibox,1)*calpi)+1
         kmaxm = dint(hmat(ibox,5)*calpi)+1
         kmaxn = dint(hmat(ibox,9)*calpi)+1
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

      vol = vol/(8.0d0*onepi)

      hmaxsq = alpsqr4*onepi*onepi

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
      recipzeo=recipzeo*qi/vol

      return

      end
