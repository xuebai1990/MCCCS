      subroutine suzeo()

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none 
!$$$      include 'grid.inc'
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'zeopoten.inc'
!$$$      include 'zeolite.inc'
!$$$      include 'external.inc'
!$$$      include 'mpi.inc'
      integer(KIND=int)::imol,iunit,igtype,idi,jerr,i,j,k,
     &     ngrxt,ngryt,ngrzt,ibox=1
      character(LEN=default_path_length)::filename
      real(KIND=double_precision)::zunitxt,zunityt,zunitzt,exzeof,xi,yi,zi

! === load force field
!      call forcefield(rczeo)

! === tabulation of the zeolite potential
      if (lzgrid) then
         call setpbc(ibox)
! find all bead types and store them in an array
         gntype=0
         do imol=1,nmolty
            do iunit=1,nunit(imol)
               do igtype=1,gntype
                  if (gtable(igtype).eq.ntype(imol,iunit)) exit
               end do
               if (igtype.gt.gntype) then
                  gntype=igtype
                  gtable(gntype)=ntype(imol,iunit)
               end if
            end do
         end do

         zunitx = zeorx/nx
         zunity = zeory/ny
         zunitz = zeorz/nz
!     --- volume minimal box
!     ideal grid size: dgr
!         dgr = 0.2
!         dgr=(dgr)**(1.d00/3.d00)
!     --- find closest values for x and y
         ngrx=int(zunitx/dgr)
         ngry=int(zunity/dgr)
         ngrz=int(zunitz/dgr)
!     if (ngrid.gt.maxtab) call cleanup('ngrid.gt.maxtab')
         dgrx=zunitx/ngrx
         dgry=zunity/ngry
         dgrz=zunitz/ngrz
         if (myid.eq.0) write(iou,1000) zunitx,zunity,zunitz,
     &        ngrx,dgrx,ngry,dgry,ngrz,dgrz


         allocate(egrid(0:ngrx-1,0:ngry-1,0:ngrz-1,gntype),stat=jerr)
         if (jerr.ne.0) call cleanup('allocate failed')

         nlayermax=0

         do igtype=1,gntype
            idi=gtable(igtype)
            write(filename,'(I3.3,A)'),idi,'.ztb'
            open(91,file=filename,form='binary',iostat=jerr,
     &           status='old')
            if (jerr.eq.0) then
! ---    read zeolite table from disk
               if (myid.eq.0) write(iou,*) 'read in ',filename
               read(91) zunitxt,zunityt,zunitzt,ngrxt,ngryt,ngrzt
               if (abs(zunitxt-zunitx).gt.eps .or. abs(zunityt-zunity)
     &              .gt.eps .or. abs(zunitzt-zunitz).gt.eps .or.
     &              ngrxt.ne.ngrx .or. ngryt.ne.ngry .or. ngrzt.ne.ngrz)
     &              call cleanup('problem in zeolite potential table')
               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        read(91) egrid(i,j,k,igtype)
                     end do
                  end do
               end do
            else
! make a tabulated potential of the zeolite
               if (myid.eq.0) write(iou,*) 'make new ',filename
               if (myid.eq.0) open(91,file=filename,form='binary')
               if (myid.eq.0) write(91) zunitx,zunity,zunitz,ngrx, ngry
     &              ,ngrz
               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        egrid(i,j,k,igtype)=exzeof(dble(i)/ngrx
     &                       ,dble(j)/ngry,dble(k)/ngrz,idi)
                        if (myid.eq.0) write(91) egrid(i,j,k,igtype)
                     end do
                  end do
               end do
               if (myid.eq.0) write(iou,*) 'maxlayer = ',nlayermax
!     call ztest(idi)
            end if
            if (myid.eq.0) then
               close(91)
               close(16)
            end if
         end do        
      end if

      return

 1000 format(' Tabulated zeolite potential: ',/,
     &     ' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/
     &     '         x-dir           : ',i12,'  size: ',f7.4,/,
     &     '         y-dir           : ',i12,'  size: ',f7.4,/,
     &     '         z-dir           : ',i12,'  size: ',f7.4,/)

      end  
