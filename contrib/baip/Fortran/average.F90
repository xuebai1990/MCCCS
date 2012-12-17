!> \brief Reads x-y plots and averages the y-values weighted
!>       by the number of cycles used to generate them
!>
!> Files 1 through N must be named fort.101, fort.102 ... fort.(100+N)
!> \par EXAMPLE INPUT
!> 3 => number of plots
!> 200 => number of x-y pairs in plot 1
!> 220 => number of x-y pairs in plot 2
!> 180 => number of x-y pairs in plot 3
program average
  implicit none

  !INPUT VARIABLES
  integer::nplot,maxval
  integer,allocatable::nval(:)
  real,allocatable::weight(:)

  !VARIABLES READ FROM FORT.10X
  integer,dimension(:,:),allocatable::ncount
  real,dimension(:,:),allocatable::yval
  character(LEN=32),dimension(:),allocatable::xval

  !LOCAL VARIABLES
  integer::ierr,iplot,ival,ifile
  real::dev,rdum
  real,allocatable::sum_sqdev(:)
  character(len=128)::valStr

  !OUTPUT
  integer,allocatable::ntcount(:)
  real,dimension(:),allocatable::yavg,stdev,err_mean

  open(unit=100,access='sequential',action='read',file='average.cfg',form='formatted',iostat=ierr,status='old')
  read(100,*) nplot
  allocate(nval(nplot))
  maxval=0
  do iplot=1,nplot
     read(100,*) nval(iplot)
     maxval=max(maxval,nval(iplot))
  end do
  close(100)

  allocate(xval(1:maxval))
  allocate(yval(1:nplot,1:maxval))
  allocate(ncount(1:nplot,1:maxval))
  allocate(sum_sqdev(1:maxval))
  allocate(ntcount(1:maxval))
  allocate(yavg(1:maxval))
  allocate(stdev(1:maxval))
  allocate(err_mean(1:maxval))

  yval=0.0
  ncount=0
  yavg=0.0
  ntcount=0

  !READ IN FORT.10X FILES
  do iplot = 1,nplot
     ifile = 100 + iplot
     do ival = 1,nval(iplot)
        read(ifile,*) xval(ival), valStr, rdum, ncount(iplot,ival)
        if (valStr.eq.'NaN') cycle
        read(valStr,*) yval(iplot,ival)
        ntcount(ival) = ntcount(ival) + ncount(iplot,ival)
        yavg(ival) = yavg(ival) + yval(iplot,ival) * ncount(iplot,ival)
     enddo
  enddo

  where (ntcount.ne.0) yavg = yavg / ntcount

  !compute standard deviation in result
  sum_sqdev = 0.
  err_mean = 0.
  do ival = 1,maxval
     do iplot = 1,nplot
        dev = ncount(iplot,ival)*(yval(iplot,ival)-yavg(ival))**2
        sum_sqdev(ival) = sum_sqdev(ival) + dev
     enddo
     if (ntcount(ival).ne.0) then
        stdev(ival) = sqrt(sum_sqdev(ival)/ntcount(ival))
     end if
     err_mean(ival) = stdev(ival)/sqrt(real(nplot))
  enddo

  !write out final result
  open(unit=99,access='sequential',action='write',file='plotavg.xy',form='formatted',iostat=ierr,status='unknown')
  do ival = 1,maxval
        write(99,*) trim(xval(ival)), yavg(ival), err_mean(ival), ntcount(ival)
  enddo

end program average
