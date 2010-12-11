      subroutine recippress(ibox,repress,pxx,pyy,pzz,pxy,pyx,pxz,pzx,
     &                      pyz,pzy)

      implicit none
!    ********************************************************************
!    ** calculates the reciprocal space contribution to pressure using **
!    ** thermodynamic definition. See J. Chem. Phys. Vol. 109 P2791.   **
!    ** written in 1998 by Bin Chen.                                   **
!    ** modified to calculate surface tension, 11/24/03 JMS            **
!    ********************************************************************

      include 'control.inc'
      include 'ewaldsum.inc'
      include 'coord.inc'
      include 'poten.inc'
! ---RP added for MPI
      include 'mpi.inc'
      include 'mpif.h'

      integer::ncount,ibox,i,ii,imolty
      real(8)::factor,repress,repressx,repressy,repressz
     &     ,recipintra,piix,piiy,piiz,xcmi,ycmi,zcmi,arg

      real(8)::pxx,pyy,pzz,intraxx,intrayy,intrazz,intraxy
     &     ,intraxz,intrazy,intrayz,intrayx,intrazx,pxy,pyx,pyz,pzy
     &     ,pxz,pzx

!----RP added for MPI
      integer::blocksize,mystart,myend
      real(8)::sum_repressx,sum_repressy,sum_repressz,
     &   sum_pxy,sum_pxz,sum_pyz
      real(8)::sum_recipintra,sum_intraxx,sum_intrayy,
     & sum_intrazz,sum_intraxy,sum_intrazy,sum_intraxz,sum_intrayz,
     & sum_intrazx,sum_intrayx

      repress  = 0.0d0
      repressx = 0.0d0
      repressy = 0.0d0
      repressz = 0.0d0
      recipintra = 0.0d0
      pxy = 0.0d0
      pxz = 0.0d0
      pyx = 0.0d0
      pyz = 0.0d0
      pzx = 0.0d0
      pzy = 0.0d0

      intraxx = 0.0d0
      intrayy = 0.0d0
      intrazz = 0.0d0
      intraxy = 0.0d0
      intrazy = 0.0d0
      intraxz = 0.0d0
      intrazx = 0.0d0
      intrayz = 0.0d0
      intrayx = 0.0d0

! KM for MPI
      sum_repressx = 0.0d0
      sum_repressy = 0.0d0
      sum_repressz = 0.0d0
      sum_pxy = 0.0d0
      sum_pxz = 0.0d0
      sum_pyz = 0.0d0
      sum_recipintra = 0.0d0
      sum_intraxx = 0.0d0
      sum_intrayy = 0.0d0
      sum_intrazz = 0.0d0
      sum_intraxy = 0.0d0
      sum_intrazy = 0.0d0
      sum_intraxz = 0.0d0
      sum_intrayz = 0.0d0
      sum_intrazx = 0.0d0
      sum_intrayx = 0.0d0

! RP for MPI
      do ncount = myid+1,numvect(ibox),numprocs
!      do ncount = 1, numvect(ibox)
         factor = prefact(ncount,ibox)*(ssumr(ncount,ibox)*
     &        ssumr(ncount,ibox) + ssumi(ncount,ibox)*
     &        ssumi(ncount,ibox))
         repressx = repressx + factor*(1.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kx(ncount,ibox)*kx(ncount,ibox))
         repressy = repressy + factor*(1.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*ky(ncount,ibox)*ky(ncount,ibox))
         repressz = repressz + factor*(1.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kz(ncount,ibox)*kz(ncount,ibox))
         pxy = pxy + factor*(0.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kx(ncount,ibox)*ky(ncount,ibox))
         pxz = pxz + factor*(0.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*kx(ncount,ibox)*kz(ncount,ibox))
         pyz = pyz + factor*(0.0d0 - (1.0d0/(4.0d0*calp(ibox)
     &        *calp(ibox)) + 1.0d0/(kx(ncount,ibox)*kx(ncount,ibox)+
     &        ky(ncount,ibox)*ky(ncount,ibox)+kz(ncount,ibox)*
     &        kz(ncount,ibox)))*2.0d0*ky(ncount,ibox)*kz(ncount,ibox))
      
       end do

! -- RP added for MPI
       CALL MPI_ALLREDUCE(repressx,sum_repressx,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(repressy,sum_repressy,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(repressz,sum_repressz,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pxy,sum_pxy,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pxz,sum_pxz,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(pyz,sum_pyz,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

      repressx = sum_repressx
      repressy = sum_repressy
      repressz = sum_repressz
      pxy = sum_pxy
      pxz = sum_pxz
      pyz = sum_pyz

      repress = repressx + repressy + repressz
! * keep x,y,z separate for surface tension calculation
      pxx = repressx
      pyy = repressy
      pzz = repressz
      pyx = pxy
      pzx = pxz
      pzy = pyz

! --- the intramolecular part should be substracted
! RP for MPI
      do 10 i = myid+1, nchain, numprocs
! ### check if i is in relevant box ###
         if ( nboxi(i) .eq. ibox ) then

            imolty = moltyp(i)
            if ( .not. lelect(imolty) ) goto 10
            xcmi = xcm(i)
            ycmi = ycm(i)
            zcmi = zcm(i)  

! --- loop over all beads ii of chain i 
            do ii = 1, nunit(imolty)
               
! --- compute the vector of the bead to the COM (p)
               
               piix = rxu(i,ii) - xcmi
               piiy = ryu(i,ii) - ycmi
               piiz = rzu(i,ii) - zcmi
             
               do ncount = 1,numvect(ibox)

! --- compute the dot product of k and r
                  
                  arg = kx(ncount,ibox)*rxu(i,ii) +
     &                 ky(ncount,ibox)*ryu(i,ii) +
     &                 kz(ncount,ibox)*rzu(i,ii)
                  factor = prefact(ncount,ibox)*2.0d0*
     &                 (-ssumr(ncount,ibox)*dsin(arg)
     &                 +ssumi(ncount,ibox)*dcos(arg))*qqu(i,ii)
                  recipintra = recipintra + factor*
     &                 (kx(ncount,ibox)*piix+ky(ncount,ibox)*piiy
     &                 +kz(ncount,ibox)*piiz)
! * keep x,y and z separate for surface tension calculation
                  intraxx = intraxx + factor*(kx(ncount,ibox)*piix)
                  intrayy = intrayy + factor*(ky(ncount,ibox)*piiy)
                  intrazz = intrazz + factor*(kz(ncount,ibox)*piiz)
                  intraxy = intraxy + factor*(kx(ncount,ibox)*piiy)
                  intraxz = intraxz + factor*(kx(ncount,ibox)*piiz)
                  intrayx = intrayx + factor*(ky(ncount,ibox)*piix)
                  intrayz = intrayz + factor*(ky(ncount,ibox)*piiz)
                  intrazx = intrazx + factor*(kz(ncount,ibox)*piix)
                  intrazy = intrazy + factor*(kz(ncount,ibox)*piiy)

               end do
            end do
         end if
 10   continue
 
       CALL MPI_ALLREDUCE(recipintra,sum_recipintra,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intraxx,sum_intraxx,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intrayy,sum_intrayy,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intrazz,sum_intrazz,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intraxy,sum_intraxy,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intraxz,sum_intraxz,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intrayx,sum_intrayx,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intrayz,sum_intrayz,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intrazx,sum_intrazx,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
       CALL MPI_ALLREDUCE(intrazy,sum_intrazy,1,
     &    MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 
      recipintra = sum_recipintra
      intraxx = sum_intraxx
      intrayy = sum_intrayy
      intrazz = sum_intrazz
      intraxy = sum_intraxy
      intraxz = sum_intraxz
      intrayx = sum_intrayx
      intrayz = sum_intrayz
      intrazx = sum_intrazx
      intrazy = sum_intrazy

      repress = (repress + recipintra)*qqfact
      
      pxx = (pxx + intraxx)*qqfact
      pyy = (pyy + intrayy)*qqfact
      pzz = (pzz + intrazz)*qqfact

      pxy = pxy + intraxy
      pyx = pyx + intrayx
      pxz = pxz + intraxz
      pzx = pzx + intrazx
      pyz = pyz + intrayz
      pzy = pzy + intrazy

!      write(iou,*) 'internal part:',intraxx,intrayy,intrazz
      
      return
      end
