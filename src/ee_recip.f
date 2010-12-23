      subroutine ee_recip(ibox,vrecipnew,vrecipold,type)

!    *********************************************************************
!    ** calculates the reciprocal ewald-sum term for trans, rot, flucq, **
!    ** swatch and swap moves, and update the reciprocal ewald-sum.     **
!    ** rewritten on June 25/99 by Bin Chen.                            **
!    *********************************************************************

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
!$$$      include 'poten.inc'
      integer(KIND=normal_int)::ic,zzz,ii,imolty,ibox,ncount,type
      real(KIND=double_precision)::vrecipnew,vrecipold,sumr(2),sumi(2)
     & ,arg

      ncount = numvect(ibox)

      if ( type .eq. 1 ) then
         
! *** recalculate the reciprocal space part for one-particle move, translation,
! *** rotation, swap, flucq, and swatch.
! *** old conformation zzz = 1 (which is 0 for swap inserted molecule)
! *** new conformation zzz = 2 (which is 0 for swap removed molecule)

!         write(iou,*) 'in recip:',moltion(1),moltion(2)
!         do zzz = 1,2
!            imolty = moltion(zzz)
!            do ii = 1, nunit(imolty)
!               write(iou,*) rxuion(ii,zzz),ryuion(ii,zzz),rzuion(ii,zzz),
!     &              qquion(ii,zzz)
!            end do
!         end do

         do 30 ic = 1, ncount
            do 20 zzz = 1,2
! --- zzz = 1: old configuration 
! --- zzz = 2: new configuration

               sumr(zzz) = 0.0d0
               sumi(zzz) = 0.0d0
               imolty = moltion(zzz)
               do ii = 1, nunit(imolty)
!                  if ( lqchg(ntype(imolty,ii)) ) then
                     arg = kx(ic,ibox)*rxuion(ii,zzz) +
     &                    ky(ic,ibox)*ryuion(ii,zzz) +
     &                    kz(ic,ibox)*rzuion(ii,zzz)
                     sumr(zzz) = sumr(zzz) + 
     &                    qquion(ii,zzz)*dcos(arg)
                     sumi(zzz) = sumi(zzz) + 
     &                    qquion(ii,zzz)*dsin(arg)
!                  end if
               end do
 20         continue
            ssumrn(ic,ibox) = ssumr(ic,ibox) - sumr(1)
     &           + sumr(2)
            ssumin(ic,ibox) = ssumi(ic,ibox) - sumi(1)
     &           + sumi(2)
 30      continue
         vrecipnew = 0.0d0
         vrecipold = 0.0d0
         do ic = 1,ncount
            vrecipnew = vrecipnew + (ssumrn(ic,ibox)*
     &           ssumrn(ic,ibox) + ssumin(ic,ibox)*
     &           ssumin(ic,ibox))*prefact(ic,ibox)
            vrecipold = vrecipold + (ssumr(ic,ibox)*
     &           ssumr(ic,ibox) + ssumi(ic,ibox)*
     &           ssumi(ic,ibox))*prefact(ic,ibox)
         end do

         vrecipnew = vrecipnew*qqfact
         vrecipold = vrecipold*qqfact

      elseif (type .eq. 2) then

! *** update the reciprocal space k vectors

         do ic = 1, ncount
            ssumr(ic,ibox) = ssumrn(ic,ibox)
            ssumi(ic,ibox) = ssumin(ic,ibox)
         end do

      elseif (type .eq. 3) then

! *** store the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumro(ic,ibox) = ssumr(ic,ibox)
            ssumio(ic,ibox) = ssumi(ic,ibox)
         end do

      elseif (type .eq. 4) then

! *** restore the reciprocal space k vectors         
         
         do ic = 1, ncount
            ssumr(ic,ibox) = ssumro(ic,ibox)
            ssumi(ic,ibox) = ssumio(ic,ibox)
         end do

      end if

!      write(iou,*) 'in recip:',ssumr(100,ibox),ibox,ssumrn(100,ibox)

      return 
      end

