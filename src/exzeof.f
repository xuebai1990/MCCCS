      function exzeof(xi,yi,zi,idi)

      use grid
      implicit none 
      real(8)::exzeof,xi,yi,zi,r2,rcutsq,xr,yr,zr,r2i,r6,vljnew,vqnew
     &     ,overflow=10.0_8**15,rminsq
      integer::j,idi,idj,ntij,layer,ii,jj,kk
      include 'zeopoten.inc'
      include 'zeolite.inc'
      include 'control.inc'
      include 'external.inc'
      include 'system.inc'
      include 'poten.inc'      
      include 'mpif.h'
      include 'mpi.inc'

!     calculate the Lennard-Jones and Coulombic interactions, include as
!     many layers of neighboring unit cells as needed for the specified
!     precision

      exzeof=0.
      rminsq=rmin*rmin
      
      if ((ltailc.and.lij(idi)).or.(lewald.and.lqchg(idi))) then
         if (ltailc.and.lij(idi)) then
            vljnew=100.
         else
            vljnew=0.
         end if
         if (lewald.and.lqchg(idi)) then
            vqnew=100.
         else
            vqnew=0.
         end if
         layer=0
         do while (vljnew.gt.eps .or. vqnew.gt.eps)
            if (layer.gt.nlayermax) nlayermax=layer
            vljnew=0.
            vqnew=0.
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
                                 if (myid.eq.0) write(iou,*) 'OF!'
                                 exzeof=overflow
                                 return
                              end if
                              r2i=sig2ij(ntij)/r2
                              r6=r2i*r2i*r2i
                              if (ltailc.and.lij(idi)) then
                                 if (lshift) then
                                    vljnew=vljnew+4.*(epsij(ntij)*(r6
     &                                   -1.0)*r6-ecut(ntij))
                                 else
                                    vljnew=vljnew+4.*epsij(ntij)*(r6
     &                                   -1.0)*r6
                                 end if
                              end if
                              if (lewald.and.lqchg(idi)) then
                                 vqnew=vqnew+qelect(idi)*qelect(idj)
     &                                /dsqrt(r2)
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
            vqnew=vqnew*qqfact
            exzeof=exzeof+vljnew+vqnew
            layer=layer+1
            if (myid.eq.0) write(iou,*) 'layer ',layer,'...'
         end do
      end if

      vljnew=0.
      vqnew=0.
      if ((.not.ltailc.and.lij(idi)).or.(.not.lewald.and.lqchg(idi)))
     &     then
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
               if (myid.eq.0) write(iou,*) '2:overflow!'
               exzeof=overflow
               return
            else if (r2 .lt. rcutsq) then
               r2i=sig2ij(ntij)/r2
               r6=r2i*r2i*r2i
               if (.not.ltailc.and.lij(idi)) then
                  if (lshift) then     
                     vljnew=vljnew+4.*(epsij(ntij)*(r6-1.0)*r6
     &                    -ecut(ntij))
                  else
                     vljnew=vljnew+4.*epsij(ntij)*(r6-1.0)*r6
                  end if               
               end if
               if (.not.lewald.and.lqchg(idi)) then
                  vqnew=vqnew+qelect(idi)*qelect(idj)/dsqrt(r2)
               end if
            end if
         end do
      end if
      exzeof=exzeof+vljnew+vqnew*qqfact

      return
      end
 
