module util_timings
  use var_type,only:double_precision
  implicit none
  private
  public::time_date_str,time_init,time_now
  real(kind=double_precision),save::time_start
  include 'common.inc'

contains
! *****************************************************************************
!> \brief returns a datum in human readable format using a standard Fortran routine
!> \par History
!>      10.2009 created [Joost VandeVondele]
!> taken from CP2K
! *****************************************************************************
  function time_date_str()
    CHARACTER::time_date_str*23,date*8,time*10
    CALL DATE_AND_TIME(date=date, time=time)
    time_date_str=date(1:4)//"-"//date(5:6)//"-"//date(7:8)//" "//time(1:2)//":"//time(3:4)//":"//time(5:10)
  end function time_date_str

  subroutine time_init()
! subroutine cpu_time(real time) is available in Fortran 95
! and later, while double precision function MPI_WTIME() is
! available using MPI, whether or not running with multiple
! processes
!    call cpu_time(time_start)
    time_start=MPI_WTIME()
!!$   time_start=omp_get_wtime()
  end subroutine time_init

  function time_now()
    real(kind=double_precision)::time_now
!    call cpu_time(time_now)
    time_now=MPI_WTIME()
!!$   time_now=omp_get_wtime()
    time_now=time_now-time_start
  end function time_now

end module util_timings
