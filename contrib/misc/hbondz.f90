subroutine hbondz(bin_width,rOO,rOH)

  use fort10_vars !stores all variables from fort.10

  implicit none

  integer :: iframe,ichain,imolty,ibin,nbin,ibox,i
  integer :: jchain,jmolty,jbox,jbin
  integer, dimension(:), allocatable :: ncount

  real :: bin_width,halfbox,rOO,rOH,z,xl,yl
  real :: dx, dy, dz, dr
  real, dimension(:), allocatable :: hist_hbond
  real, dimension(:), allocatable :: hist_donor
  real, dimension(:), allocatable :: hist_acceptor

  logical :: lhbond_Hi_Oj, lhbond_Hj_Oi

  write(6,*) '   ...Begin hbondz.f90'

  nbin = int(boxlz(1,1)/bin_width)+1

  allocate(ncount(1:nbin))
  allocate(hist_hbond(1:nbin))
  allocate(hist_donor(1:nbin))
  allocate(hist_acceptor(1:nbin))

  ncount = 0
  hist_hbond = 0
  hist_donor = 0
  hist_acceptor = 0


  xl = boxlx(1,1)
  yl = boxly(1,1)


  loop_frames : do iframe = 1,nframe

     !count the number of waters in each bin
     loop_count: do ichain = 1,nchain
        imolty = moltyp(ichain)
        ibox = molbox(iframe,ichain)
        if ( (imolty.eq.1) .and. (ibox.eq.1) ) then
           z = zcom(iframe,ichain)
           ibin = int(z/bin_width) + 1
           if (ibin.lt.nbin) then
              ncount(ibin) = ncount(ibin) + 1
           endif
        endif
     enddo loop_count

     loop_ichain : do ichain = 1,nchain-1

        imolty = moltyp(ichain)
        ibox = molbox(iframe,ichain)

        limolty: if ( (imolty.eq.1) .and. (ibox.eq.1) ) then

           loop_jchain: do jchain = ichain+1,nchain

              jmolty = moltyp(jchain)
              jbox = molbox(iframe,jchain)

              ljmolty: if ( (jmolty.eq.1) .and. (jbox.eq.1) ) then

                 lhbond_Hi_Oj = .false.
                 lhbond_Hj_Oi = .false.

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

                 lOO: if (dr.lt.rOO) then

                    !compute OH distance ichain-H to jchain-O
                    do i = 2,3
                       dx = xbead(iframe,ichain,i) -  &
                            xbead(iframe,jchain,1)
                       dy = ybead(iframe,ichain,i) -  &
                            ybead(iframe,jchain,1)
                       dz = zbead(iframe,ichain,i) -  &
                            zbead(iframe,jchain,1)
                       dx = dx - anint(dx/xl)*xl 
                       dy = dy - anint(dy/yl)*yl 
                       dr = sqrt(dx*dx + dy*dy + dz*dz)
                       if (dr.lt.rOH) lhbond_Hi_Oj = .true.
                    enddo

                    !compute OH distance ichain-O to jchain-H
                    do i = 2,3
                       dx = xbead(iframe,ichain,1) -  &
                            xbead(iframe,jchain,i)
                       dy = ybead(iframe,ichain,1) -  &
                            ybead(iframe,jchain,i)
                       dz = zbead(iframe,ichain,1) -  &
                            zbead(iframe,jchain,i)
                       dx = dx - anint(dx/xl)*xl 
                       dy = dy - anint(dy/yl)*yl 
                       dr = sqrt(dx*dx + dy*dy + dz*dz)
                       if (dr.lt.rOH) lhbond_Hj_Oi = .true.
                    enddo

                 end if lOO

                 z = zcom(iframe,ichain)
                 ibin = int(z/bin_width) + 1

                 z = zcom(iframe,jchain)
                 jbin = int(z/bin_width) + 1

                 libin: if (ibin.le.nbin) then
                    if (lhbond_Hi_Oj) then
                       hist_hbond(ibin) = hist_hbond(ibin) + 1.0
                       hist_donor(ibin) = hist_donor(ibin) + 1.0
                    endif
                    if (lhbond_Hj_Oi) then
                       hist_hbond(ibin) = hist_hbond(ibin) + 1.0
                       hist_acceptor(ibin) = hist_acceptor(ibin) + 1.0
                    endif
                 else
                    stop 'hbondz: ibin.gt.nbin'
                 endif libin

                 ljbin: if (jbin.le.nbin) then
                    if (lhbond_Hi_Oj) then
                       hist_hbond(jbin) = hist_hbond(jbin) + 1.0
                       hist_acceptor(jbin) = hist_acceptor(jbin) + 1.0
                    endif
                    if (lhbond_Hj_Oi) then
                       hist_hbond(jbin) = hist_hbond(jbin) + 1.0
                       hist_donor(jbin) = hist_donor(jbin) + 1.0
                    endif
                 else
                    stop 'hbondz: jbin.gt.nbin'
                 endif ljbin

              end if ljmolty
           end do loop_jchain
        end if limolty
     end do loop_ichain
  end do loop_frames

  open(unit=104,file='hbondz.dat',status='unknown')
  open(unit=105,file='donorz.dat',status='unknown')
  open(unit=106,file='acceptz.dat',status='unknown')

  halfbox = boxlz(1,1)/2.0
  do ibin = 1,nbin
     z = ibin*bin_width-halfbox-bin_width
     if (ncount(ibin).gt.0) then
        write(104,*) z,hist_hbond(ibin)/real(ncount(ibin)), ' 0.0 ', ncount(ibin)
        write(105,*) z,hist_donor(ibin)/real(ncount(ibin)), ' 0.0 ', ncount(ibin)
        write(106,*) z,hist_acceptor(ibin)/real(ncount(ibin)), ' 0.0 ', ncount(ibin)
     else
        write(104,*) z,' 0.0  0.0 ', ncount(ibin)
        write(105,*) z,' 0.0  0.0 ', ncount(ibin)
        write(106,*) z,' 0.0  0.0 ', ncount(ibin)
     endif
  enddo
  close(104)
  close(105)
  close(106)

end subroutine hbondz

