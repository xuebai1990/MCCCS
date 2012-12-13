      subroutine init_vars
 
!    *******************************************************************
!    *** initializes some variables that otherwise cause errors      ***
!    *** Matt McGrath, October 1, 2009                               ***
!    *******************************************************************
 
      use sim_system
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
!$$$      include 'system.inc'
!$$$      include 'poten.inc'
!$$$      include 'conver.inc' 
!$$$      include 'external.inc'
!$$$      include 'zeolite.inc'
!$$$      include 'connect.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'qqlist.inc'
!$$$      include 'clusterbias.inc'
!$$$      include 'neigh.inc'
!$$$      include 'cell.inc'
!$$$      include 'nsix.inc'
!$$$      include 'peboco.inc'     
!$$$      include 'ipswpar.inc'
!$$$      include 'eepar.inc'
!$$$      include 'cbmc.inc'
!$$$! KM 01/10 remove analysis
!$$$!      include 'gor.inc'
!$$$      include 'blkavg.inc'
!$$$! kea include for garofalini potential
!$$$      include 'garofalini.inc'
!$$$      include 'tabulated.inc'
! **********************************************************************
! local variables

      integer(KIND=normal_int)::ibend,ibox,i,imolty,iunit,iprop,jblock ,itype


! **********************************************************************
      nmolty1=0
      leemove=.false.
      iratipsw=0
      acipsw=0.0d0
! KM remove analysis
!      nframe=0.0d0
      DO ibend=1,nvmax
         brben(ibend)=0.0d0
      end do
      DO itype=1,nntype
         lij(itype)=.FALSE.
      end do
      DO itype=1,nntype*nntype
         sig2ij(itype)=0.0d0
         epsij(itype)=0.0d0
      end do
      DO iprop=1,nprop
         DO ibox=1,nbxmax
            DO jblock=1,blockm
               baver(iprop,ibox,jblock)=0.0d0
            end do
         end do
      end do
      DO ibox=1,nbxmax
         rmvol(ibox)=0.0d0
         boxlx(ibox)=0.0d0
         boxly(ibox)=0.0d0
         boxlz(ibox)=0.0d0
         DO i=1,9
            hmat(ibox,i)=0.0d0
         end do
         lsolid(ibox)=.FALSE.
         lrect(ibox)=.FALSE.
      end do
      upnn=0.0d0

      DO imolty=1,ntmax
         llplace(imolty)=.FALSE.
         lwell(imolty)=.FALSE.
         nrig(imolty)=0
      end do

      DO iunit=1,numax
         lsave(iunit)=.false.
      end do

      return
      end 

