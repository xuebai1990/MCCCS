module zeolite
    use global_data
    use var_type
    use const_math
    use const_phys
    use util_math 
    use util_files
    use util_string
    implicit none
    private
    save
    public zeocoord,suzeo,exzeo
      real(KIND=double_precision)::zeorx,zeory,zeorz,zunitx,zunity,zunitz
      integer(KIND=normal_int)::nx,ny,nz,nzeo,zntype
      logical,allocatable::lunitcell(:)
      integer(KIND=normal_int),allocatable::idzeo(:)
      real(KIND=double_precision),allocatable::zeox(:),zeoy(:),zeoz(:),zradiisq(:)
      integer(KIND=normal_int),allocatable::ztype(:),znum(:)
      character(LEN=default_string_length),allocatable::zname(:)
      integer(KIND=normal_int),parameter::gntypemax=32
      integer(KIND=normal_int)::gntype,gtable(0:gntypemax)
      integer(KIND=normal_int)::ngrx,ngry,ngrz,nlayermax
      real(KIND=double_precision),allocatable::egrid(:,:,:,:)
      real(KIND=double_precision),parameter::dgr=0.2,require_precision=1.0E-2,overflow=1.0E+8_double_precision,upperlimit=1.0E+7_double_precision
      integer(KIND=normal_int),parameter::ibox=1
  contains
    subroutine zeocoord(file_zeocoord,lhere)
      character(LEN=*),intent(in)::file_zeocoord
      logical,intent(out)::lhere(:)
      integer(KIND=normal_int)::io_zeocoord,jerr
      integer(KIND=normal_int)::pos,i
      real(KIND=double_precision)::alpha,betaang,gamma,sx,sy,sz
      character(LEN=4)::atom

      io_zeocoord=get_iounit()
      open(unit=io_zeocoord,access='sequential',action='read',file=file_zeocoord,form='formatted',iostat=jerr,status='old')
      if (jerr.ne.0) then
         call cleanup('cannot open zeolite coordinate file')
      end if

      read(io_zeocoord,*)    zeorx,zeory,zeorz,alpha,betaang,gamma,nx,ny,nz

      if (myid.eq.0) write(iou,"(/,' READING ZEOLITE LATTICE FROM FILE zeolite.cssr:',/&
                                  ,' --------------------------------------------------',/&
                                  ,' box dimensions                    = ',3f10.3,' Angstrom',/&
                                  ,' box angles                        = ',3f10.3,' degrees',/&
                                  ,' number of zeolite cells           = ',3i5)") zeorx,zeory,zeorz,alpha,betaang,gamma,nx,ny,nz

      if (abs(alpha-90)>eps.or.abs(betaang-90)>eps.or.abs(gamma-90)>eps) then
         lsolid(ibox)=.true.
         lrect(ibox)=.false.
      else
         lsolid(ibox)=.false.
      end if

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
      hmat(ibox,8)=zeorz*sqrt(cos(betaang)**2*cos(gamma)**2+cos(alpha)**2-2*cos(alpha)*cos(betaang)*cos(gamma))/sin(gamma)
      hmat(ibox,9)=zeorz*sqrt(1-cos(alpha)**2-cos(betaang)**2-cos(gamma)**2+2*cos(alpha)*cos(betaang)*cos(gamma))/sin(gamma)

      where(abs(hmat(ibox,:)).lt.eps) hmat(ibox,:)=0

      call matops(ibox)

      if (myid.eq.0) write(iou,"(' simulation volume                 = ', f10.1, 20x,' Angst**3')") cell_vol(ibox)

      read(io_zeocoord,*) nzeo,zntype
      allocate(zeox(nzeo),zeoy(nzeo),zeoz(nzeo),idzeo(nzeo),lunitcell(nzeo),zname(zntype),zradiisq(zntype),ztype(zntype),znum(zntype),stat=jerr)
      if (jerr.ne.0) call cleanup('zeocoord: allocation failed')

      if (myid.eq.0) write(iou,"(' number of zeolite atoms           = ',i10,/&
                                ,' number of atomtypes in the lattice= ',i10,/)") nzeo,zntype
      do i=1,zntype
         read(io_zeocoord,*) zname(i),ztype(i),zradiisq(i)
         lhere(ztype(i))=.true.
         zradiisq(i)=zradiisq(i)*zradiisq(i)
      end do

!     === Converting to absolute coordinates within [0,ri>
      znum=0
      do i = 1,nzeo
         read(io_zeocoord,'(i5,1x,a4,2x,3(f9.5,1x))') pos,atom,zeox(i),zeoy(i),zeoz(i)

         sx = zeox(i)*hmati(ibox,1)+zeoy(i)*hmati(ibox,4)+zeoz(i)*hmati(ibox,7)
         sy = zeox(i)*hmati(ibox,2)+zeoy(i)*hmati(ibox,5)+zeoz(i)*hmati(ibox,8)
         sz = zeox(i)*hmati(ibox,3)+zeoy(i)*hmati(ibox,6)+zeoz(i)*hmati(ibox,9)
         sx = sx - floor(sx)
         sy = sy - floor(sy)
         sz = sz - floor(sz)
         zeox(i)=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
         zeoy(i)=sy*hmat(ibox,5)+sz*hmat(ibox,8)
         zeoz(i)=sz*hmat(ibox,9)

         if (sx*nx.lt.1.and.sy*ny.lt.1.and.sz*nz.lt.1) then
            lunitcell(i)=.true.
         else
            lunitcell(i)=.false.
         end if

         pos=str_search(zname,zntype,atom)
         if (pos.eq.0) call cleanup('** atomtype: unknown atomtype **')

         idzeo(i)=pos
         znum(pos)=znum(pos)+1

!         if (myid.eq.0) write(iou,*) zeox(i),zeoy(i),zeoz(i),idzeo(i)
      end do

      close(io_zeocoord)

      return
    end subroutine zeocoord

    subroutine suzeo()
      integer(KIND=normal_int)::io_ztb,imol,iunit,igtype,idi,jerr,i,j,k,oldi,oldj,oldk,ngrxt,ngryt,ngrzt
      character(LEN=default_path_length)::file_ztb
      real(KIND=double_precision)::zunitxt,zunityt,zunitzt,wzeo

! === tabulation of the zeolite potential
      if (lzgrid) then
         call setpbc(ibox)
! find all bead types and store them in an array
         gntype=0
         gtable(0)=0
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
!     --- find closest values for x and y
         ngrx=int(zunitx/dgr)
         ngry=int(zunity/dgr)
         ngrz=int(zunitz/dgr)
         allocate(egrid(0:ngrx-1,0:ngry-1,0:ngrz-1,0:gntype),stat=jerr)
         if (jerr.ne.0) call cleanup('suzeo: allocation failed')

!     === Calculate zeolite density
         wzeo=dot_product(znum,mass(ztype))
         if (myid.eq.0) write(iou,"(' Tabulated zeolite potential: ',/&
                                   ,' --------------------------------------------------',/&
                                   ,<zntype>(' number of ',A,':',I5,3X),/&
                                   ,' mass zeolite                      = ',e12.5,' grams',/ &
                                   ,' one adsorbed molecule in sim box  = ',e12.5 ,' mmol/gram',/&
                                   ,' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/&
                                   ,'         x-dir           : ',i12,'  size: ',f7.4,/&
                                   ,'         y-dir           : ',i12,'  size: ',f7.4,/&
                                   ,'         z-dir           : ',i12,'  size: ',f7.4,/)") (trim(zname(i)),znum(i),i=1,zntype),wzeo/6.023e23,1000.0/wzeo,zunitx,zunity,zunitz,ngrx,zunitx/ngrx,ngry,zunity/ngry,ngrz,zunitz/ngrz

         nlayermax=0
         io_ztb=get_iounit()

         do igtype=0,gntype
            if (igtype.eq.0.and..not.lewald) cycle
            idi=gtable(igtype)
            write(file_ztb,'(I3.3,A)'),idi,'.ztb'
            open(unit=io_ztb,access='sequential',action='read',file=file_ztb,form='binary',iostat=jerr, status='old')
            if (jerr.eq.0) then
! ---    read zeolite table from disk
               if (myid.eq.0) write(iou,*) 'read in ',trim(file_ztb)
               read(io_ztb) zunitxt,zunityt,zunitzt,ngrxt,ngryt,ngrzt
               if (abs(zunitxt-zunitx).gt.eps .or. abs(zunityt-zunity) .gt.eps .or. abs(zunitzt-zunitz).gt.eps .or. ngrxt.ne.ngrx .or. ngryt.ne.ngry .or. ngrzt.ne.ngrz) call cleanup('problem in zeolite potential table')
               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        read(io_ztb,end=100) egrid(i,j,k,igtype)
                     end do
                  end do
               end do
            end if

100         if (jerr.eq.0.and.i.lt.ngrx) then
               oldi=i
               oldj=j
               oldk=k
               jerr=1
               close(io_ztb)
               write(iou,*) oldi,oldj,oldk
            else
               oldi=0
               oldj=0
               oldk=0
            end if

            if (jerr.ne.0) then
! make a tabulated potential of the zeolite
               if (myid.eq.0) write(iou,*) 'make new ',trim(file_ztb)
               if (myid.eq.0) open(unit=io_ztb,access='sequential',action='write',file=file_ztb,form='binary',iostat=jerr,status='unknown')
               if (jerr.ne.0) then
                  call cleanup('cannot create file for tabulated potential')
               end if
               if (myid.eq.0) write(io_ztb) zunitx,zunity,zunitz,ngrx,ngry,ngrz

               do i=0,ngrx-1
                  do j=0,ngry-1
                     do k=0,ngrz-1
                        if (i.gt.oldi.or.(i.eq.oldi.and.j.ge.oldj)) egrid(i,j,k,igtype)=exzeof(dble(i)/ngrx,dble(j)/ngry,dble(k)/ngrz,idi)
                        if (myid.eq.0) write(io_ztb) egrid(i,j,k,igtype)
                     end do
                  end do
               end do
               if (myid.eq.0.and.ltailc) write(iou,*) 'maxlayer = ',nlayermax
!     call ztest(idi)
            end if

            if (idi.ne.0.and.lewald.and.lqchg(idi)) then
               where (egrid(:,:,:,igtype).lt.overflow) egrid(:,:,:,igtype)=egrid(:,:,:,igtype)+qelect(idi)*egrid(:,:,:,0)
            end if
            if (myid.eq.0) then
               close(io_ztb)
            end if
         end do        
      end if

      return
    end subroutine suzeo

    function exzeof(i,j,k,idi)
      real(KIND=double_precision)::exzeof
      real(KIND=double_precision),intent(in)::i,j,k
      integer(KIND=normal_int),intent(in)::idi
      
      integer(KIND=normal_int)::izeo,idj,ntij,layer,ii,jj,kk
      real(KIND=double_precision)::sx,sy,sz,xi,yi,zi,xr,yr,zr,r2,rcutsq,r,r2i,r6,vljnew,vqnew

      exzeof=0.
      
      if (idi.eq.0.or.(.not.ltailc).or.(.not.lewald.and.lqchg(idi))) then
         vljnew=0.
         vqnew=0.
         rcutsq = rcut(ibox)*rcut(ibox)
         sx=dble(i)/nx
         sy=dble(j)/ny
         sz=dble(k)/nz
         xi=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
         yi=sy*hmat(ibox,5)+sz*hmat(ibox,8)
         zi=sz*hmat(ibox,9)
         do izeo=1,nzeo
            idj=ztype(idzeo(izeo))
            ntij = (idi - 1) * nntype + idj
            xr=xi-zeox(izeo)
            yr=yi-zeoy(izeo)
            zr=zi-zeoz(izeo)
            call mimage(xr,yr,zr,ibox)
            r2=xr*xr+yr*yr+zr*zr
            if (r2.le.zradiisq(idzeo(izeo))) then
               exzeof=overflow
               return
            end if
            if (r2 .lt. rcutsq) then
               if (idi.eq.0) then
                  if (lqchg(idj)) then
                     r=dsqrt(r2)
                     vqnew=vqnew+qelect(idj)*erfunc(calp(ibox)*r)/r
                  end if
               else
                  if (.not.lewald.and.lqchg(idi).and.lqchg(idj)) then
                     r=dsqrt(r2)
                     vqnew=vqnew+qelect(idi)*qelect(idj)/r
                  end if
                  if (.not.ltailc.and.lij(idi).and.lij(idj)) then
                     r2i=sig2ij(ntij)/r2
                     r6=r2i*r2i*r2i
                     if (lshift) then     
                        vljnew=vljnew+4.*(epsij(ntij)*(r6-1.0)*r6-ecut(ntij))
                     else
                        vljnew=vljnew+4.*epsij(ntij)*(r6-1.0)*r6
                     end if
                  end if
               end if
            end if
         end do

         if (idi.eq.0) vqnew=vqnew+recipzeo(xi,yi,zi,1._double_precision)

         exzeof=vljnew+vqnew*qqfact
      end if

! Calculate the Lennard-Jones interactions, include as many layers
! of neighboring unit cells as needed for the specified precision
      if (idi.ne.0.and.ltailc) then
         vljnew=require_precision+1.
         layer=0
         do while (abs(vljnew).gt.require_precision)
            if (layer.gt.nlayermax) nlayermax=layer
            vljnew=0.
            do izeo=1,nzeo
               idj=ztype(idzeo(izeo))
               if (lunitcell(izeo)) then
                  ntij = (idi - 1) * nntype + idj
                  do ii=-layer,layer
                     do jj=-layer,layer
                        do kk=-layer,layer
                           if (abs(ii).eq.layer .or. abs(jj).eq.layer .or.abs(kk).eq.layer) then

                              sx = zeox(izeo)*hmati(ibox,1)+zeoy(izeo)*hmati(ibox,4)+zeoz(izeo)*hmati(ibox,7)+dble(ii-i)/nx
                              sy = zeox(izeo)*hmati(ibox,2)+zeoy(izeo)*hmati(ibox,5)+zeoz(izeo)*hmati(ibox,8)+dble(jj-j)/ny
                              sz = zeox(izeo)*hmati(ibox,3)+zeoy(izeo)*hmati(ibox,6)+zeoz(izeo)*hmati(ibox,9)+dble(kk-k)/nz

                              xr=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
                              yr=sy*hmat(ibox,5)+sz*hmat(ibox,8)
                              zr=sz*hmat(ibox,9)

                              r2=xr*xr+yr*yr+zr*zr
                              if (r2.le.zradiisq(idzeo(izeo))) then
                                 exzeof=overflow
                                 return
                              end if
                              if (lij(idi).and.lij(idj)) then
                                 r2i=sig2ij(ntij)/r2
                                 r6=r2i*r2i*r2i
                                 if (lshift) then
                                    vljnew=vljnew+4.*(epsij(ntij)*(r6-1.0)*r6-ecut(ntij))
                                 else
                                    vljnew=vljnew+4.*epsij(ntij)*(r6-1.0)*r6
                                 end if
                              end if
                           end if
                        end do
                     end do
                  end do
               end if
            end do
            exzeof=exzeof+vljnew
            layer=layer+1
         end do
      end if

      return
    end function exzeof
 
    function recipzeo(xi,yi,zi,qi)
      real(KIND=double_precision)::recipzeo
      real(KIND=double_precision),intent(in)::xi,yi,zi,qi

      integer(KIND=normal_int)::i,l,m,n,m_min,n_min,kmaxl,kmaxm,kmaxn
      real(KIND=double_precision)::alpsqr4,vol,ksqr,sum,arg,hmatik(9),kx1,ky1,kz1,hmaxsq,calpi

! *** Set up the reciprocal space vectors ***
      recipzeo = 0.0d0

      calpi = calp(ibox)

      do i = 1,9
         hmatik(i) = twopi*hmati(ibox,i)
      end do
      if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
         kmaxl = dint(hmat(ibox,1)*calpi)+1
         kmaxm = dint(hmat(ibox,5)*calpi)+1
         kmaxn = dint(hmat(ibox,9)*calpi)+1
      else
         kmaxl = dint(hmat(ibox,1)*calpi)+2
         kmaxm = dint(hmat(ibox,5)*calpi)+2
         kmaxn = dint(hmat(ibox,9)*calpi)+2
      end if
    
      alpsqr4 = 4.0d0*calpi*calpi
         
      vol = hmat(ibox,1)* (hmat(ibox,5)*hmat(ibox,9) - hmat(ibox,8)*hmat(ibox,6)) + hmat(ibox,4)* (hmat(ibox,8)*hmat(ibox,3) - hmat(ibox,2)*hmat(ibox,9)) + hmat(ibox,7)* (hmat(ibox,2)*hmat(ibox,6) - hmat(ibox,5)*hmat(ibox,3))

      vol = vol/(8.0d0*onepi)

      hmaxsq = alpsqr4*onepi*onepi

! *** generate the reciprocal-space
! here -kmaxl,-kmaxl+1,...,-1 are skipped, so no need to divide by 2 for the prefactor
      do l = 0,kmaxl 
        if ( l .eq. 0 ) then
            m_min = 0
         else
            m_min = -kmaxm
         end if
         do m = m_min, kmaxm
            if (l .eq. 0 .and. m .eq. 0) then
               n_min = 1
            else
               n_min = -kmaxn
            end if
            do n = n_min, kmaxn
               kx1 = dble(l)*hmatik(1)+dble(m)*hmatik(2)+ dble(n)*hmatik(3)
               ky1 = dble(l)*hmatik(4)+dble(m)*hmatik(5)+ dble(n)*hmatik(6)
               kz1 = dble(l)*hmatik(7)+dble(m)*hmatik(8)+ dble(n)*hmatik(9)
               ksqr = kx1*kx1+ky1*ky1+kz1*kz1
!               if ( ksqr .lt. hmaxsq ) then
! --- sometimes these are about equal, which can cause different
! --- behavior on 32 and 64 bit machines without this .and. statement
               if ( ksqr .lt. hmaxsq .and. abs(ksqr-hmaxsq) .gt. 1d-9 ) then
                  sum = 0.0d0
                  do i = 1,nzeo
                     arg=kx1*(xi-zeox(i))+ky1*(yi-zeoy(i))+kz1*(zi -zeoz(i))
                     sum=sum+qelect(ztype(idzeo(i)))*cos(arg)
                  end do
                  recipzeo=recipzeo+sum*exp(-ksqr/alpsqr4)/ksqr
               end if
            end do
         end do
      end do

      recipzeo=recipzeo*qi/vol

      return

    end function recipzeo
    
    function exzeo(xi,yi,zi,idi)
      real(KIND=double_precision)::exzeo
      real(KIND=double_precision),intent(in)::xi,yi,zi
      integer(KIND=normal_int),intent(in)::idi

      integer(KIND=normal_int),parameter::m=2,mt=2*m+1,mst=-m
      integer(KIND=normal_int)::j,j0,jp,k,k0,kp,l,l0,lp,igtype
      real(KIND=double_precision)::yjtmp(mst:m),yktmp(mst:m),yltmp(mst:m)
      real(KIND=double_precision)::xt(mst:m),yt(mst:m),zt(mst:m),sx,sy,sz,xr,yr,zr
      
! --- determine cell parameters
      sx=(xi*hmati(ibox,1)+yi*hmati(ibox,4)+zi*hmati(ibox,7))*nx
      sy=(xi*hmati(ibox,2)+yi*hmati(ibox,5)+zi*hmati(ibox,8))*ny
      sz=(xi*hmati(ibox,3)+yi*hmati(ibox,6)+zi*hmati(ibox,9))*nz
      sx=sx-floor(sx)
      sy=sy-floor(sy)
      sz=sz-floor(sz)

      if (.not.lzgrid) then
         exzeo=exzeof(sx,sy,sz,idi)
         if (lewald.and.lqchg(idi)) exzeo=exzeo+exzeof(sx,sy,sz,0)*qelect(idi)
      else
!     calculation using a grid
!         write(iou,*) 'entering exzeo. xi yi zi idi',xi,yi,zi,idi

         do igtype=1,gntype
            if (gtable(igtype).eq.idi) exit
         end do
         if (igtype.gt.gntype) then
            call cleanup('exzeo: no such bead type')
         end if
         
         xr=sx*hmat(ibox,1)/nx+sy*hmat(ibox,4)/ny+sz*hmat(ibox,7)/nz
         yr=sy*hmat(ibox,5)/ny+sz*hmat(ibox,8)/nz
         zr=sz*hmat(ibox,9)/nz
         j = sx*ngrx
         k = sy*ngry
         l = sz*ngrz

! ---    test if in the reasonable regime
         exzeo=upperlimit
         if (egrid(j,k,l,igtype).ge.exzeo) return
! --     block m*m*m centered around: j,k,l
! ---  set up hulp array: (allow for going beyond unit cell
!      for polynom fitting)
         do l0=mst,m
            lp=l+l0
            sz=dble(lp)/ngrz/nz
! ---    store x,y,z values around xi,yi,zi in arrays
            zt(l0)=sz*hmat(ibox,9)
            if (lp.lt.0)    lp=lp+ngrz
            if (lp.ge.ngrz) lp=lp-ngrz
            do k0=mst,m
	       kp=k+k0
               sy=dble(kp)/ngry/ny
               yt(k0)=sy*hmat(ibox,5)+sz*hmat(ibox,8)
	       if (kp.lt.0)    kp=kp+ngry
	       if (kp.ge.ngry) kp=kp-ngry
               do j0=mst,m
                  jp=j+j0
                  sx=dble(jp)/ngrx/nx
                  xt(j0)=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox ,7)
                  if (jp.lt.0)    jp=jp+ngrx
                  if (jp.ge.ngrx) jp=jp-ngrx
                  yjtmp(j0)=egrid(jp,kp,lp,igtype)
                  if (yjtmp(j0).ge.exzeo) return
               end do
               call polint(xt,yjtmp,mt,xr,yktmp(k0))
	    end do
            call polint(yt,yktmp,mt,yr,yltmp(l0))
         end do
         call polint(zt,yltmp,mt,zr,exzeo)
      end if
      return
    end function exzeo

end module zeolite
