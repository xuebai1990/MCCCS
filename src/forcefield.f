      subroutine forcefield(rczeo)

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
   
!$$$      include 'zeolite.inc'
!$$$      include 'zeopoten.inc'
!$$$      include 'control.inc'

      real(KIND=double_precision)::epsilon,sigma,rczeo
      integer(KIND=int)::         itype,jtype,idz

!     In atomtype.f atoms O of the zeolite have been given id=1
!     Hopefully in some other place the guestmolecules CH4 have been
!     assigned id=2

!      if (zntype.ne.2) call cleanup('** forcefield: ntype ne 2 not allowed **')

!     Force field data:

!     The interaction parameters are in this case calculated from 
!     atomic values

!     Fill the matrices

! Needed are epsilon and sigma such that in Angstroms and kcal/mol:
! Energy = 4*epsilon*((sigma/r)**12 - (sigma/r)**6)

      read(25,*)  
      read(25,*)
      idz=1
      do itype = 1,zntype
	 read(25,*) epsilon,sigma
         zeps(itype,idz) = epsilon
         zsig2(itype,idz) = sigma*sigma
      end do

! Write epsilons and sigma's:

      write(16,100)
      jtype=idz
      do itype = 1,zntype
         write(16,1001) itype,jtype,zeps(itype,jtype),
     &                  itype,jtype,dsqrt(zsig2(itype,jtype))
      end do
 
! Force field cutoff radius rczeo for interactions between guests
! and zeolite atoms:

      write(16,101) rczeo
      do itype = 1,zntype
        zrc2(itype,idz)=rczeo**2
      end do
!
! === calculate cutt-off of the potential
!
      if (lshift) then
        write(iou,102)
        do itype = 1,zntype
          zencut(itype,idz)=4.*zeps(itype,idz)*
     &         ( (zsig2(itype,idz)/zrc2(itype,idz))**6
     &          -(zsig2(itype,idz)/zrc2(itype,idz))**3)
          write(iou,103) itype,zencut(itype,idz)
        end do
      end if     

  100 format(/,' FORCE FIELD DEFINITION:',/,
     &         ' ------------------------------------------------')
  101 format(/,' Force field cutoffs: rczeo = ',f6.1,/,
     &         '                      rcads = ',f6.1,' Angstrom',/,
     &         ' ------------------------------------------------',//)
 1000 format(/,' Lennard-Jones Parameters (K, Angstrom):')
 1001 format(  ' epsilon(',i2,',',i2,') = ',f6.1,
     &         ' K &  sigma(',i2,',',i2,') = ',f6.1,' Angstrom')
 102  format(' Value at cut-off distance ' )
 103  format('    interaction ',i3, '   : ',f8.2,'[K]')        
      return
      end

