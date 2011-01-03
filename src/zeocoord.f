      subroutine zeocoord

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
!$$$      include 'grid.inc'
!$$$      include 'zeolite.inc'
!$$$      include 'zeopoten.inc'
!$$$      include 'system.inc'
!$$$      include 'cell.inc'
!$$$      include 'mpi.inc'
      integer(KIND=normal_int)::count,frac,izeo,bonding(8),atomtype,ibox
     & =1
      real(KIND=double_precision)::wzeo,charge,alpha,betaang,gamma,sx,sy
     & ,sz
      character(LEN=4)::atom

      open (unit = 47, file = 'zeolite.cssr')
      if (myid.eq.0) write(16,"(/
     & ,' READING ZEOLITE LATTICE FROM FILE zeolite.cssr:',/
     & ,' --------------------------------------------------')")

      read(47,*)    zeorx,zeory,zeorz,alpha,betaang,gamma
      if (abs(alpha-90)>eps.or.abs(betaang-90)>eps.or.abs(gamma-90)>eps)
     &     then
         lsolid(ibox)=.true.
         lrect(ibox)=.false.
      else
         lsolid(ibox)=.false.
      end if
      if (myid.eq.0) write(16,"(' box angles                        =
     & ',3f10.3,' degrees')") alpha,betaang,gamma
      alpha=alpha*onepi/180
      betaang=betaang*onepi/180
      gamma=gamma*onepi/180
      hmat(ibox,1)=zeorx
      hmat(ibox,2)=0.
      hmat(ibox,3)=0.
      hmat(ibox,4)=zeory*cos(gamma)
      hmat(ibox,5)=zeory*sin(gamma)
      hmat(ibox,6)=0.
      hmat(ibox,7)=zeorz*cos(betaang)
      hmat(ibox,8)=zeorz*sqrt(cos(betaang)**2*cos(gamma)**2+cos(alpha)
     & **2-2*cos(alpha)*cos(betaang)*cos(gamma))/sin(gamma)
      hmat(ibox,9)=zeorz*sqrt(1-cos(alpha)**2-cos(betaang)**2-cos(gamma)
     &     **2+2*cos(alpha)*cos(betaang)*cos(gamma))/sin(gamma)
      where(abs(hmat(ibox,:)).lt.eps) hmat(ibox,:)=0
      call matops(ibox)

      if (myid.eq.0) then
         write(16,"(  ' box dimensions                    = ',3f10.3,
     &    ' Angstrom',/, ' simulation volume                 =
     &    ', f10.1, 20x,' Angst**3')") zeorx,zeory,zeorz,cell_vol(ibox)
      end if

      read(47,*)   nzeo,frac,nx,ny,nz

      if (myid.eq.0) write(16,"(' number of zeolite atoms           =
     & ',i10)") nzeo
      if (nzeo.gt.nzeomax)
     &     call cleanup('** zeocoord: nzeo gt nzeomax **')

!     === Converting to absolute coordinates within [0,ri>
      ztype=(/177,178/)
      znum=0
      do izeo = 1,nzeo
         read(47,'(i5,1x,a4,2x,3(f9.5,1x))') count,atom,zeox(izeo)
     &    ,zeoy(izeo),zeoz(izeo)
         sx = zeox(izeo)*hmati(ibox,1)+zeoy(izeo)*hmati(ibox,4)
     &        +zeoz(izeo)*hmati(ibox,7)
         sy = zeox(izeo)*hmati(ibox,2)+zeoy(izeo)*hmati(ibox,5)
     &        +zeoz(izeo)*hmati(ibox,8)
         sz = zeox(izeo)*hmati(ibox,3)+zeoy(izeo)*hmati(ibox,6)
     &        +zeoz(izeo)*hmati(ibox,9)
         sx = sx - floor(sx)
         sy = sy - floor(sy)
         sz = sz - floor(sz)
         zeox(izeo)=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
         zeoy(izeo)=sy*hmat(ibox,5)+sz*hmat(ibox,8)
         zeoz(izeo)=sz*hmat(ibox,9)
         if (sx*nx.lt.1.and.sy*ny.lt.1.and.sz*nz.lt.1) then
            lunitcell(izeo)=.true.
         else
            lunitcell(izeo)=.false.
         end if
         idzeo(izeo)=atomtype(zntype,atom)
!         if (myid.eq.0) write(16,*) zeox(izeo),zeoy(izeo),
!     &        zeoz(izeo),idzeo(izeo)
      end do

!     === Calculate zeolite density from assumption no of Si = 0.5* no O
      wzeo = znum(1)*28.09 + znum(2)*16.00
      if (myid.eq.0) write(16,"(' mass zeolite                      =
     & ',e10.3,' grams',/ ,' one adsorbed molecule in sim box  =
     & ',f10.3 ,' mmol/gram')") wzeo/6.023e23,1000.0/wzeo

      if (myid.eq.0) then
         write(16,"(  ' number of zeolite cells           = ',3i5,/,
     &    ' --------------------------------------------------')") nx,ny
     &    ,nz
         write(16,*) 'number of Si:',znum(1),'number of O:',znum(2)
      end if

      return
      end


