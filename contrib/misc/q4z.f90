subroutine Q4z(bin_width)

  use fort10_vars !stores all variables from fort.10

  implicit none

  integer :: iframe,ichain,imolty,ibin,nbin,ibox,i
  integer :: jchain,jmolty,jbox,jbin,j
  integer, dimension(1:4) :: ineighbor
  integer, dimension(:), allocatable :: ncount

  real :: bin_width,halfbox, xl, yl, onethird, qsum, cosphi
  real :: dx, dy, dz, dr, z
  real :: dr4, dr3, dr2, dr1
  real :: xi,yi,zi,dxi,dyi,dzi,dri  
  real :: xj,yj,zj,dxj,dyj,dzj,drj
  real :: xO,yO,zO 
  real :: rsum, ravg
  real, dimension(1:4) :: rtetra
  real, dimension(:), allocatable :: hist_angle, hist_dist


  write(6,*) '   ...Begin Q4z.f90'

  nbin = int(boxlz(1,1)/bin_width)+1

  allocate(ncount(1:nbin))
  allocate(hist_angle(1:nbin))
  allocate(hist_dist(1:nbin))

  ncount = 0
  hist_angle = 0
  hist_dist = 0

  xl = boxlx(1,1)
  yl = boxly(1,1)

  onethird = 1.0/3.0

  loop_frames : do iframe = 1,nframe
 
     loop_ichain : do ichain = 1,nchain

        imolty = moltyp(ichain)
        ibox = molbox(iframe,ichain)

        limolty: if ( (imolty.eq.1) .and. (ibox.eq.1) ) then

           dr1 = 1000.0
           dr2 = 1000.0
           dr3 = 1000.0
           dr4 = 1000.0
           ineighbor = 0

           loop_jchain: do jchain = 1,nchain

              jmolty = moltyp(jchain)
              jbox = molbox(iframe,jchain)

              ljmolty: if ( (jmolty.eq.1) .and. (jbox.eq.1) .and. &
                   (ichain.ne.jchain) ) then

                 !compute O-O distance
                 dx = xbead(iframe,ichain,1) -  &
                      xbead(iframe,jchain,1)
                 dy = ybead(iframe,ichain,1) -  &
                      ybead(iframe,jchain,1)
                 dz = zbead(iframe,ichain,1) -  &
                      zbead(iframe,jchain,1)
                 dx = dx - anint(dx/xl)*xl 
                 dy = dy - anint(dy/yl)*yl 
                 !note no PBC in z-direction
                 dr = sqrt(dx*dx + dy*dy + dz*dz)

                 ldr4: if(dr.lt.dr4) then
                    dr4 = dr
                    ineighbor(4) = jchain
                    ldr3: if(dr.lt.dr3) then
                       dr4 = dr3
                       dr3 = dr
                       ineighbor(4) = ineighbor(3)
                       ineighbor(3) = jchain
                       ldr2: if (dr.lt.dr2) then
                          dr3 = dr2
                          dr2 = dr
                          ineighbor(3) = ineighbor(2)
                          ineighbor(2) = jchain
                          ldr1: if (dr.lt.dr1) then
                             dr2 = dr1
                             dr1 = dr
                             ineighbor(2) = ineighbor(1)
                             ineighbor(1) = jchain
                          endif ldr1
                       endif ldr2
                    endif ldr3
                 endif ldr4

              endif ljmolty
           enddo loop_jchain

           !compute order parameters
           qsum = 0.0
           loopi: do i = 1,3
              loopj: do j = i+1,4
                 
                 xi = xbead(iframe,ineighbor(i),1)
                 yi = ybead(iframe,ineighbor(i),1)
                 zi = zbead(iframe,ineighbor(i),1)

                 xj = xbead(iframe,ineighbor(j),1)
                 yj = ybead(iframe,ineighbor(j),1)
                 zj = zbead(iframe,ineighbor(j),1)
                 
                 xO = xbead(iframe,ichain,1)
                 yO = ybead(iframe,ichain,1)
                 zO = zbead(iframe,ichain,1)

                 dxi = xi - xO
                 dyi = yi - yO
                 dzi = zi - zO
                 dxi = dxi - anint(dxi/xl)*xl 
                 dyi = dyi - anint(dyi/yl)*yl 
                 dri = sqrt(dxi*dxi + dyi*dyi + dzi*dzi)

                 dxj = xj - xO
                 dyj = yj - yO
                 dzj = zj - zO
                 dxj = dxj - anint(dxj/xl)*xl 
                 dyj = dyj - anint(dyj/yl)*yl 
                 drj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj)

                 cosphi = (dxi*dxj + dyi*dyj + dzi*dzj)/(dri*drj)
                 qsum = qsum + (cosphi+onethird)**2

                 rtetra(i) = dri
                 rtetra(j) = drj

              enddo loopj
           enddo loopi

           qsum = 1.0 - 3.0*qsum/8.0

           ravg = 0.0
           rsum = 0.0
           do i = 1,4
              ravg = ravg + rtetra(i)
           enddo
           ravg = ravg/4.0
           do i = 1,4
              rsum = rsum + (rtetra(i)-ravg)**2/(4*ravg*ravg)
           enddo
           
           rsum = 1.0 - onethird*rsum
         

           z = zcom(iframe,ichain)
           ibin = int(z/bin_width) + 1
           if (ibin.lt.nbin) then
              ncount(ibin) = ncount(ibin) + 1
              hist_angle(ibin) = hist_angle(ibin) + qsum
              hist_dist(ibin) = hist_dist(ibin) + rsum
           else
              stop 'q4z.f: ibin.gt.nbin'
           endif

        endif limolty
     enddo loop_ichain
  enddo loop_frames

  open(unit=107,file='tetra-angle.dat',status='unknown')
  open(unit=108,file='tetra-dist.dat',status='unknown')

  halfbox = boxlz(1,1)/2.0
  do ibin = 1,nbin
     z = ibin*bin_width-halfbox-bin_width
     if (ncount(ibin).gt.0) then
        write(107,*) z,hist_angle(ibin)/real(ncount(ibin)), ' 0.0 ', ncount(ibin)
        write(108,*) z,hist_dist(ibin)/real(ncount(ibin)), ' 0.0 ', ncount(ibin)
     else
        write(107,*) z,' 0.0  0.0 ', ncount(ibin)
        write(108,*) z,' 0.0  0.0 ', ncount(ibin)
     endif
  enddo
  close(107)

end subroutine Q4z
