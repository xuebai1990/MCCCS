subroutine rhoz(nbin,target_box,target_molty)

  use fort10_vars !stores all variables from fort.10

  implicit none

  integer :: target_box,target_molty,iframe,ichain,imolty,ibin,nbin,ibox

  real :: bin_width,halfbox,z,bin_vol
  real, dimension(:), allocatable :: hist_rhoz

  write(6,*) '   Calculating density profile (as a function of z) ...'

  bin_width = boxlz(1,target_box)/nbin

  allocate(hist_rhoz(1:nbin))
  hist_rhoz = 0

  loop_chains : do ichain = 1,nchain
     if (moltyp(ichain) .ne. target_molty) cycle loop_chains
     loop_frames : do iframe = 1,nframe
        lmolty: if (molbox(iframe,ichain) .eq. target_box) then
           ibin = int(zcom(iframe,ichain)/bin_width) + 1
           if (ibin.eq.nbin+1) then
              ibin = nbin;
           elseif (ibin .gt. nbin+1) then
              write(6,*) 'ibin =',ibin,'  nbin =',nbin
              write(6,*) 'COMz =',zcom(iframe,ichain)
              stop 'rhoz.f: ibin.gt.nbin'
           endif
           hist_rhoz(ibin) = hist_rhoz(ibin) + 1.0
        end if lmolty
     end do loop_frames
  end do loop_chains

  !divide by volume of bin and number of frames
  hist_rhoz = hist_rhoz/nframe/(bin_width*boxlx(1,target_box)*boxly(1,target_box))

  !convert from 1/A^3 to g/mL
  !1 ethane molecule = 30.0694 amu = 4.993258x10^-23 g
  !1 A^3 = 10^-24 cm^3
!  hist_rhoz = hist_rhoz*4.993258d-23/1.0e-24
  hist_rhoz = hist_rhoz*(3.207**3)

  !shift by half of boxlength
  halfbox = bin_width/2.0 + boxlz(1,target_box)/2.0

  open(unit=101,file='rhoz.dat',status='unknown')
  do ibin = 1,nbin
     write(101,*) ibin*bin_width-halfbox,hist_rhoz(ibin), ' 0.0 ', nframe 
  enddo
  close(101)

end subroutine rhoz
     

