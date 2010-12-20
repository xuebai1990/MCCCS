      subroutine dipole(ibox,mtype)

      implicit none
      integer::ibox,mtype,i,imolty,zz,ii
      real(8)::dipox(2),dipoy(2),dipoz(2)
      include 'control.inc'
      include 'coord.inc'
      include 'coord2.inc'
      include 'ewaldsum.inc'
      
      if ( mtype .eq. 0 ) then
! in sumup, initiliaze dipole to be 0 and then sum up all the dipoles

         dipolex(ibox) = 0.0d0
         dipoley(ibox) = 0.0d0
         dipolez(ibox) = 0.0d0
         
         do i = 1, nchain
            if ( nboxi(i) .eq. ibox ) then
               imolty = moltyp(i)
               do ii = 1,nunit(imolty)
                  dipolex(ibox) = dipolex(ibox) 
     &                 + qqu(i,ii)*rxu(i,ii)
                  dipoley(ibox) = dipoley(ibox) 
     &                 + qqu(i,ii)*ryu(i,ii)
                  dipolez(ibox) = dipolez(ibox) 
     &                 + qqu(i,ii)*rzu(i,ii)
               end do
            end if
         end do

      elseif(mtype .eq. 1 ) then
         
! *** calculate the dipole moment after the traslation, rotation and 
! *** charge move

         do zz = 1,2
            dipox(zz) = 0.0d0
            dipoy(zz) = 0.0d0
            dipoz(zz) = 0.0d0
            imolty = moltion(zz)
            do i = 1,nunit(imolty)
               dipox(zz) = dipox(zz) + 
     &              qquion(i,zz)*rxuion(i,zz)
               dipoy(zz) = dipoy(zz) +
     &              qquion(i,zz)*ryuion(i,zz)
               dipoz(zz) = dipoz(zz) +
     &              qquion(i,zz)*rzuion(i,zz)
            end do
         end do
         dipolex(ibox) = dipolex(ibox) - dipox(1) + dipox(2) 
         dipoley(ibox) = dipoley(ibox) - dipoy(1) + dipoy(2)
         dipolez(ibox) = dipolez(ibox) - dipoz(1) + dipoz(2)

      elseif(mtype .eq. 2) then

! *** store the old dipole moment

         dipolexo = dipolex(ibox) 
         dipoleyo = dipoley(ibox)
         dipolezo = dipolez(ibox)
         
      elseif(mtype .eq. 3) then

! *** restore the old dipole moment

         dipolex(ibox) = dipolexo 
         dipoley(ibox) = dipoleyo
         dipolez(ibox) = dipolezo

      end if
      
       
      return
      end
