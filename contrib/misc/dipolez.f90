subroutine dipolez(nbin,target_box,target_molty)

  use fort10_vars !stores all variables from fort.10

  implicit none

  integer :: iframe,ichain,imolty,ibin,nbin,ibox,target_box,target_molty
  integer, dimension(:), allocatable :: ncount

  real :: bin_width,halfbox,z,dx,dy,dz,dr,cos,p2
  real, dimension(:), allocatable :: hist_cos
  real, dimension(:), allocatable :: hist_p2

  write(6,*) '   Calculating orientation and orientation parameter ...'

  bin_width = boxlz(1,target_box)/nbin

  allocate(ncount(1:nbin))
  allocate(hist_cos(1:nbin))
  allocate(hist_p2(1:nbin))
  ncount = 0
  hist_cos = 0
  hist_p2 = 0

  loop_chains : do ichain = 1,nchain
     if (moltyp(ichain).ne.target_molty) cycle loop_chains
     loop_frames : do iframe = 1,nframe
        lmolty: if (molbox(iframe,ichain).eq.target_box) then
           dx = xbead(iframe,ichain,2) - xbead(iframe,ichain,1)
           dy = ybead(iframe,ichain,2) - ybead(iframe,ichain,1)
           dz = zbead(iframe,ichain,2) - zbead(iframe,ichain,1)
           cos = abs(dz)/sqrt( dx*dx + dy*dy + dz*dz )
           p2 = 1.5*cos*cos - 0.5
           ibin = int(zcom(iframe,ichain)/bin_width) + 1
           if (ibin.eq.nbin+1) then
              ibin = nbin
           elseif (ibin .gt. nbin+1) then
              write(6,*) 'ibin =',ibin, '  nbin =',nbin
              write(6,*) 'COMz =',zcom(iframe,ichain)
              stop 'dipolez.f: ibin > nbin'
           endif
           ncount(ibin) = ncount(ibin) + 1
           hist_cos(ibin) = hist_cos(ibin) + cos
           hist_p2(ibin) = hist_p2(ibin) + p2
        end if lmolty
     end do loop_frames
  end do loop_chains

  !shift by half of boxlength
  halfbox = bin_width/2 !+ boxlz(1,target_box)/2.0

  open(unit=102,file='cos_dipole.dat',status='unknown')
  open(unit=103,file='p2_dipole.dat',status='unknown')
  do ibin = 1,nbin
     z = ibin*bin_width-halfbox
     if (ncount(ibin).gt.0) then
        write(102,*) z,hist_cos(ibin)/real(ncount(ibin)), ' 0.0 ', ncount(ibin)
        write(103,*) z,hist_p2(ibin)/real(ncount(ibin)), ' 0.0 ', ncount(ibin)
     else
        write(102,*) z,' 0.0  0.0  0'
        write(103,*) z,' 0.0  0.0  0'
     endif
  enddo
  close(102)
  close(103)

end subroutine dipolez
