      subroutine recip_atom(ibox,vrecipnew,vrecipold,type,ii)

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
      integer(KIND=normal_int)::ic,izz,ii,imolty,ibox,ncount,type
      real(KIND=double_precision)::vrecipnew,vrecipold,sumr(2),sumi(2)
     & ,arg

      ncount = numvect(ibox)

      if ( type .eq. 1 ) then
         
! *** recalculate the reciprocal space part for one-particle move, translation,
! *** rotation, swap, flucq, and swatch.
! *** old conformation izz = 1 (which is 0 for swap inserted molecule)
! *** new conformation izz = 2 (which is 0 for swap removed molecule)

!         write(6,*) 'in recip:',moltion(1),moltion(2)
!         do izz = 1,2
!            imolty = moltion(izz)
!            do ii = 1, nunit(imolty)
!               write(6,*) rxuion(ii,izz),ryuion(ii,izz),rzuion(ii,izz),
!     &              qquion(ii,izz)
!            end do
!         end do

         do 30 ic = 1, ncount
            do 20 izz = 1,2
! --- izz = 1: old configuration 
! --- izz = 2: new configuration

               sumr(izz) = 0.0d0
               sumi(izz) = 0.0d0
               imolty = moltion(izz)
                  if ( lqchg(ntype(imolty,ii)) ) then
                     arg = kx(ic,ibox)*rxuion(ii,izz) +
     &                    ky(ic,ibox)*ryuion(ii,izz) +
     &                    kz(ic,ibox)*rzuion(ii,izz)
                     sumr(izz) = sumr(izz) + 
     &                    qquion(ii,izz)*dcos(arg)
                     sumi(izz) = sumi(izz) + 
     &                    qquion(ii,izz)*dsin(arg)
                  end if
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

!      write(6,*) 'in recip:',ssumr(100,ibox),ibox,ssumrn(100,ibox)

      return 
      end

