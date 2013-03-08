      subroutine dipole(ibox,mtype)

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
!$$$      include 'coord2.inc'
!$$$      include 'ewaldsum.inc'
      integer(KIND=normal_int)::ibox,mtype,i,imolty,zzz,ii
      real(KIND=double_precision)::dipox(2),dipoy(2),dipoz(2)
      
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

         do zzz = 1,2
            dipox(zzz) = 0.0d0
            dipoy(zzz) = 0.0d0
            dipoz(zzz) = 0.0d0
            imolty = moltion(zzz)
            do i = 1,nunit(imolty)
               dipox(zzz) = dipox(zzz) + 
     &              qquion(i,zzz)*rxuion(i,zzz)
               dipoy(zzz) = dipoy(zzz) +
     &              qquion(i,zzz)*ryuion(i,zzz)
               dipoz(zzz) = dipoz(zzz) +
     &              qquion(i,zzz)*rzuion(i,zzz)
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
