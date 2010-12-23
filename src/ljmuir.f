       function ljmuir ( rijsq, ntij )

!    *************************************************
!    ** calculates SAMI's 12+3 energy for headgroup **
!    **            and normal LJ energy for tail    **
!    *************************************************
 
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
!$$$      include 'coord.inc'
!$$$      include 'poten.inc'

      real(KIND=double_precision)::ljmuir, rijsq, sr, sr2, sr6, epshead, sighead
      integer(KIND=int)::ntij

! --- attention: eps_hh / 4 used, since later multiplied by 4 --- 
!      parameter (epshead=27.67204d0,sighead=4.22d0)
      parameter (epshead=27.7204d0,sighead=6.5d0)

! --------------------------------------------------------------------

!       write(iou,*) 'sig2ij',sig2ij
!       write(iou,*) 'epsij',epsij

      if ( ntij .eq. 1 ) then
         sr = sighead / dsqrt( rijsq )
!       write(iou,*) 'sr',sr,'v',4.0d0*epshead*sr**3*(sr**9+1.0d0)
         ljmuir = epshead * sr**3 * ( sr**9 + 1.0d0 )
      else
         sr2 = sig2ij(ntij) / rijsq
         sr6 = sr2 * sr2 * sr2
         ljmuir = epsij(ntij) * sr6 * ( sr6 - 1.0d0)
!         if (ljmuir .gt. 100.0d0)
!     &       write(18,*) sig2ij(ntij),rijsq,'sr',dsqrt(sr2),'v',
!     &            4.0d0 * epsij(ntij) * sr6 * ( sr6 - 1.0d0)
      end if

      return

! ----------------------------------------------------------------------------

      end
