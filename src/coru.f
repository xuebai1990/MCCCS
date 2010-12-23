      function coru(imolty,jmolty,rho,ibox)

! **********************************
! *** tail-corrections in energy ***
! **********************************

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'poten.inc'
!$$$      include 'conver.inc'
!$$$      include 'system.inc'
!$$$      include 'expsix.inc'
!$$$      include 'merck.inc'
!$$$      include 'nsix.inc'

      real(KIND=double_precision)::coru,rci3,rho,epsilon2,sigma2
      real(KIND=double_precision)::rci1
      integer(KIND=int)::imolty,jmolty,ii,jj, ntii, ntjj, ntij,ibox

      coru = 0.0d0

      do ii = 1, nunit(imolty) 
         ntii = ntype(imolty,ii)

         do jj = 1, nunit(jmolty) 
            ntjj = ntype(jmolty,jj)
            if (lexpsix) then
               ntij = (ntii+ntjj)/2
               coru = coru + rho*consu(ntij)
            elseif (lmmff) then
               ntij = (ntii+ntjj)/2
               coru = coru + rho * epsimmff(ntij) * coru_cons(ntij) *
     &              sigimmff(ntij)**3.0d0*twopi
            elseif (lninesix) then
               ntij = (ntii-1)*nxatom + ntjj
               coru = coru + 8.0d0*onepi*rho*epsnx(ntij)*
     &            rzero(ntij)**3*(rzero(ntij)/rcut(ibox))**3*
     &            ((rzero(ntij)/rcut(ibox))**3/3.0d0 - 1.0d0) 
            elseif (lgenlj) then
               ntij = (ntii-1)*nntype + ntjj
               rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
               rci1 = rci3 **(1.0d0/3.0d0)

               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigma(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)
     &                 *epsilon(jmolty,jj))
               elseif ( lexpand(imolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigi(ntjj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)*epsi(ntjj))
               elseif ( lexpand(jmolty) ) then
                  sigma2 = (sigma(jmolty,jj)+sigi(ntii))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               end if
               coru = coru
     &          + 2.0d0 * onepi * epsilon2 * sigma2 ** (1.50d0) * rho *
     &        (  (( (2.0d0**(4.0d0*n1/n0))/(2.0d0*n1-3.0d0))
     & * rci1 **(2.0d0*n1-3.0d0) ) -
     &       ( (2.0d0**((2.0d0*n1/n0)+1.0d0))/(n1-3.0d0))
     & * rci1 **(n1-3.0d0) )

            else
               ntij = (ntii-1)*nntype + ntjj
               rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigma(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)
     &                 *epsilon(jmolty,jj))
               elseif ( lexpand(imolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigi(ntjj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)*epsi(ntjj))
               elseif ( lexpand(jmolty) ) then
                  sigma2 = (sigma(jmolty,jj)+sigi(ntii))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(jmolty,jj)*epsi(ntii))
               else
                  sigma2 = sig2ij(ntij)
                  epsilon2 = epsij(ntij)
               end if
               coru = coru + 
     &              8.0d0 * onepi * epsilon2 * 
     &              sigma2**(1.5d0) *rho * 
     &              (rci3 * rci3 * rci3 / 9.0d0 - rci3 / 3.0d0)
            end if
         end do
      end do
      return
      end





