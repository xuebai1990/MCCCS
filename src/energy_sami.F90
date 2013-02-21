MODULE energy_sami
  use sim_system
  implicit none
  private
  save
  public::susami,ljsami,exsami,ljmuir,exmuir

! EXTERNALMUIR.INC
  real,public::sigpri,c9ch2,c3ch2,c9ch3,c3ch3,zprmin,v2prmin,v3prmin,beta1,beta2,betac2,betac3
! LJSAMIPARA.INC
  real::alpha1,alpha2
  real::sij(9),eij(9),vsh(9),vsha(9)
  integer::tau1,tau2
! Sami's parameters: to be used with lsami = .true. AND lmuir = .true.
  parameter ( alpha1=21.162d0, alpha2=-21.162d0, beta1=661.6d0, beta2=6616.0d0, tau1=-32, tau2=-16 )
contains
!    *******************************
!    ** Set-Up SAMI's potentials. **
!    *******************************
  subroutine susami

      real::hsig,heps,tsig,teps

      parameter ( hsig=4.220d0, heps=110.68816d0,  tsig=3.527d0, teps=79.982210d0 )

      real::rcsami
      integer::ij

! --------------------------------------------------------------------

      rcsami = 2.5d0*tsig

      sij(1)=hsig
      sij(2)=0.5d0*(hsig+tsig)
      sij(3)=sij(2)
      sij(4)=sij(2)
      sij(5)=tsig
      sij(6)=tsig
      sij(7)=sij(2)
      sij(8)=tsig
      sij(9)=tsig

      eij(1)=heps
      eij(2)=dsqrt(heps*teps)
      eij(3)=eij(2)
      eij(4)=eij(2)
      eij(5)=teps
      eij(6)=teps
      eij(7)=eij(2)
      eij(8)=teps
      eij(9)=teps

      vsh(1)  = eij(1) * ( ( 13.0d0 * (sij(1)/rcsami)**12 ) + (  4.0d0 * (sij(1)/rcsami)**3  ) )
      vsha(1) = eij(1) * ( ( 12.0d0 * sij(1)**12 / rcsami**13 ) + (  3.0d0 * sij(1)**3  / rcsami**4  ) )

      do ij = 2, 9
         vsh(ij)  = 4.0d0 * eij(ij) *  ( ( 13.0d0 * (sij(ij)/rcsami)**12 ) - (  7.0d0 * (sij(ij)/rcsami)**6  ) )
         vsha(ij) = 4.0d0 * eij(ij) * ( ( 12.0d0 * sij(ij)**12 / rcsami**13 ) - (  6.0d0 * sij(ij)**6  / rcsami**7  ) )
      end do

!      do ij = 1,9
!         write(io_output,*) 'ij',ij,'vsh',(vsh(ij)/80.0d0),
!     +                      'vsha',(vsha(ij)/80.0d0),
!     +                      'eij',(eij(ij)/80.0d0)
!      end do

      return

! ----------------------------------------------------------------------------

  end subroutine susami

!    *********************************************************
!    ** calculates SAMI's LJ and 12+3 energy for a bead.    **
!    *********************************************************
  function ljsami( rijsq, ntij )

      real::ljsami, rijsq, rij, sr
      integer::ntij

! --------------------------------------------------------------------

      rij = dsqrt( rijsq )
      sr = sij(ntij) / rij

      if ( ntij .eq. 1 ) then
! *** head-head interaction ( repulsive 12+3 interaction )
         ljsami = ( eij(1) * sr**3 * ( 1.0d0 + sr**9 ) ) - vsh(1) + ( rij * vsha(1) )
      else
! *** head-tail or tail-tail interaction ( LJ 12-6 interaction )
         ljsami = ( 4.0d0 * eij(ntij) * sr**6 * ( sr**6 - 1.0d0 ) ) - vsh(ntij) + ( rij * vsha(ntij) )
      end if

      return
  end function ljsami

!    **********************************************************
!    ** calculates the SAMI's external energy for a bead.    **
!    **********************************************************
  function exsami( z, ntj )

      real::exsami, z
      integer::ntj

! --------------------------------------------------------------------

      if ( ntj .eq. 1 ) then
! --- HEADgroup potential ---
         if ( z .le. alpha2 ) then
            exsami = 0.0d0
         else
            exsami = beta2 / ( 1.0d0 + ( (z/alpha2) - 1.0d0 )**tau2 )
         end if
      else
! --- TAILgroup potential ---
         if ( z .ge. alpha1 ) then
            exsami = 0.0d0
         else
            exsami = beta1 / ( 1.0d0 + ( 1.0d0 - (z/alpha1) )**tau1 )
         end if
      end if

      return

  end function exsami

!    *************************************************
!    ** calculates SAMI's 12+3 energy for headgroup **
!    **            and normal LJ energy for tail    **
!    *************************************************
  function ljmuir ( rijsq, ntij )

      real::ljmuir, rijsq, sr, sr2, sr6, epshead, sighead
      integer::ntij

! --- attention: eps_hh / 4 used, since later multiplied by 4 ---
!      parameter (epshead=27.67204d0,sighead=4.22d0)
      parameter (epshead=27.7204d0,sighead=6.5d0)

! --------------------------------------------------------------------

!       write(io_output,*) 'sig2ij',sig2ij
!       write(io_output,*) 'epsij',epsij

      if ( ntij .eq. 1 ) then
         sr = sighead / dsqrt( rijsq )
!       write(io_output,*) 'sr',sr,'v',4.0d0*epshead*sr**3*(sr**9+1.0d0)
         ljmuir = epshead * sr**3 * ( sr**9 + 1.0d0 )
      else
         sr2 = sig2ij(ntij) / rijsq
         sr6 = sr2 * sr2 * sr2
         ljmuir = epsij(ntij) * sr6 * ( sr6 - 1.0d0)
!         if (ljmuir .gt. 100.0d0)
!     &       write(18,*) sig2ij(ntij),rijsq,'sr',dsqrt(sr2),'v',
!     &            4.0d0 * epsij(ntij) * sr6 * ( sr6 - 1.0d0)
      end if

      return

! ----------------------------------------------------------------------------

  end function ljmuir

!    *********************************************************
!    ** calculates the lmuir external energy for a bead.    **
!    *********************************************************
  function exmuir( z, ntj )

      real::exmuir, z
      integer::ntj

! --------------------------------------------------------------------

      if ( ntj .eq. 1 ) then
! --- HEADgroup potential ---
         if ( z .le. alpha2 ) then
            exmuir = 0.0d0
         else
            exmuir = beta2 / ( 1.0d0 + ( (z/alpha2) - 1.0d0 )**tau2 )
         end if
      else
! --- TAILgroup potential ---
         if ( z .ge. alpha1 ) then
            exmuir = 0.0d0
         else
            if ( z .lt. zprmin ) then
               if ( ntj .eq. 2 ) then
                  exmuir = betac2 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + v2prmin
               else
                  exmuir = betac3 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + v3prmin
               end if
            else
               if ( ntj .eq. 2 ) then
                  exmuir = betac2 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + c9ch2 / z**9 -  c3ch2 / z**3
               else
                  exmuir = betac3 / (1.0d0+(1.0d0-(z/alpha1))**tau1 ) + c9ch3 / z**9 -  c3ch3 / z**3
               end if
            end if
         end if
      end if
      return

  end function exmuir
end MODULE energy_sami
