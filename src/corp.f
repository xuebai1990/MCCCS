!> \brief Impulsive force and tail corrections to pressure
!>
!> @todo Impulsive force correction only done for LJ 12-6 potential
      function corp(imolty,jmolty,rhosq,ibox)

! ************************************
! *** tail-corrections in pressure ***
! ************************************

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
!$$$      include 'conver.inc'
!$$$      include 'poten.inc'
!$$$      include 'system.inc'
!$$$      include 'expsix.inc'
!$$$      include 'merck.inc'
!$$$      include 'nsix.inc'

      real(KIND=double_precision)::corp,rci3,rhosq, epsilon2, sigma2
      real(KIND=double_precision)::rci1
      integer(KIND=normal_int)::imolty,jmolty,ii,jj, ntii, ntjj, ntij ,ibox
      corp = 0.0d0

      do ii = 1, nunit(imolty) 
         ntii = ntype(imolty,ii)

         do jj = 1, nunit(jmolty) 
            ntjj = ntype(jmolty,jj)

            if (lexpsix .and. ltailc) then
               ntij = (ntii+ntjj)/2
               corp = corp + consp(ntij)
            elseif (lmmff .and. ltailc) then
               ntij = (ntii+ntjj)/2
               corp = corp+((-2.0d0)/3.0d0)*onepi*epsimmff(ntij)*sigimmff(ntij)**3.0d0*corp_cons(ntij)
            elseif (lninesix .and. ltailc) then
               ntij = (ntii-1)*nxatom + ntjj
               corp = corp + 16.0d0 * onepi * epsnx(ntij) * rzero(ntij)**3 * (0.5d0*(rzero(ntij)/rcut(ibox))**6 - (rzero(ntij)/rcut(ibox))**3)
            elseif (lgenlj .and. ltailc) then
               ntij = (ntii-1)*nntype + ntjj
               rci3 = sig2ij(ntij)**(3.0d0/2.0d0) / rcut(ibox)**3
               rci1 = rci3 **(1.0d0/3.0d0)
               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigma(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii) *epsilon(jmolty,jj))
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
                 corp = corp + (2.0d0/3.0d0) * onepi * epsilon2 * sigma2 ** (1.50d0) * n1 * ( (( (2.0d0**((4.0d0*n1/n0)+1.0d0))/(2.0d0*n1-3.0d0)) * rci1 **(2.0d0*n1-3.0d0) ) - ( ((2.0d0**((2.0d0*n1/n0)+1.0d0))/(n1-3.0d0))* rci1 **(n1-3.0d0) )  )
            else if (.not.(lexpsix.or.lmmff.or.lninesix.or.lgenlj)) then
               ntij = (ntii-1)*nntype + ntjj
               
               if ( lexpand(imolty) .and. lexpand(jmolty) ) then
                  sigma2 = (sigma(imolty,ii)+sigma(jmolty,jj))**2/4.0d0
                  epsilon2 = dsqrt(epsilon(imolty,ii)* epsilon(jmolty,jj))
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
               rci3 = (dsqrt(sigma2)/rcut(ibox))**3 
               if (ltailc) then ! both impulsive force and tail corrections 
                  corp = corp+8.0d0*onepi*epsilon2*sigma2**(1.5d0)*(rci3*rci3*rci3*7.0d0/9.0d0-rci3)
               else ! only impulsive force corrections 
                  corp=corp+(8.0d0/3.0d0)*onepi*epsilon2*sigma2**(1.5d0)*(rci3*rci3*rci3 - rci3)
               end if
            end if
         end do
      end do

      corp=corp*rhosq

      return
      end


