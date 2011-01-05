      subroutine  calcsolpar(pres,Heat_vapor_T,Heat_vapor_LJ,
     &  Heat_vapor_COUL, pdV, CED_T, CED_LJ,CED_COUL,HSP_T,HSP_LJ,
     &  HSP_COUL,ibox,jbox)
! ***********************************************************************************
!  Written by Neeraj Rai. Date 08/22/2006

!  This subroutine calculates heat of vaporization and solubility parameter
!  Separates the Coulombic and LJ part of solubility parameter
! Previos implementation of heat of vaporization was good only for single component
! system. Now we can do any number of components.
!  It is being written for two box gibbs ensemble simulations. If it is more than two boxes it will c return without doing anything.
! or if the box is solid and non rectangular it will return.
! If you have more than two boxes then there are couple changes that need to be done
! in order for this subroutine to work. nprop in blkavg.inc needs to be set properly
! ************************************************************************************
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
!$$$      include 'inputdata.inc'
!$$$      include 'blkavg.inc'
!$$$      include 'coord.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'system.inc'
!$$$      include 'cell.inc'

      integer(KIND=normal_int)::ibox, jbox, ig, il, imolty
      integer(KIND=normal_int),dimension(nbxmax):: temp_nmol,box_volume
      real(KIND=double_precision),dimension(nbxmax)::mol_vol      
      real(KIND=double_precision),dimension(nbxmax)::pres
      real(KIND=double_precision)::enchg1, enchg2,enchg3
      real(KIND=double_precision)::Heat_vapor_T, Heat_vapor_LJ,
     & Heat_vapor_COUL
      real(KIND=double_precision)::CED_T, HSP_T, CED_LJ, HSP_LJ
      real(KIND=double_precision)::CED_COUL,HSP_COUL
      real(KIND=double_precision)::cal2joule, joule2cal
      real(KIND=double_precision)::pdV

      real(KIND=double_precision)::T_Energy_Liq, T_Energy_Gas
      real(KIND=double_precision)::LJ_Energy_Liq, LJ_Energy_Gas
      real(KIND=double_precision)::Coul_Energy_Liq, Coul_Energy_Gas


!     -- This is different than the conventional 4.184
      cal2joule = 4.184d0
      joule2cal = 1.0d0/cal2joule

      temp_nmol(ibox) = 0

      do imolty = 1,nmolty
         temp_nmol(ibox) = temp_nmol(ibox) + ncmt(ibox,imolty)
      end do

      if (lsolid(ibox).and.(.not.lrect(ibox))) then
         box_volume(ibox) = cell_vol(ibox)
      else   
         box_volume(ibox) = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
      end if
!     molar volume in cc/mol    
      if (temp_nmol(ibox).eq.0) then    
!     set the molar volume artibrarily high as it will not affect any since pressure will be zero
         mol_vol(ibox) = 10.0d10    
      else
         mol_vol(ibox)=box_volume(ibox)/temp_nmol(ibox)*0.6022d0 
      end if

      temp_nmol(jbox) = 0

      do imolty = 1,nmolty
         temp_nmol(jbox) = temp_nmol(jbox) + ncmt(jbox,imolty)
      end do

      if (lsolid(jbox).and.(.not.lrect(jbox))) then
         box_volume(jbox) = cell_vol(jbox)
      else
         box_volume(jbox) = boxlx(jbox)*boxly(jbox)*boxlz(jbox)
      end if
!     -- molar volume in cc/mol
      if (temp_nmol(jbox).eq.0) then
         mol_vol(jbox) = 10.0D6
      else
         mol_vol(jbox)=box_volume(jbox)/temp_nmol(jbox)*0.6022d0
      end if

      if ( mol_vol(ibox).gt.mol_vol(jbox) ) then
         ig = ibox
         il = jbox
      else
         ig = jbox
         il = ibox
      end if
!     -- This is total energy
      if (temp_nmol(il).eq.0) then
         T_Energy_Liq = 0.0d0
         LJ_Energy_Liq = 0.0d0
         Coul_energy_Liq = 0.0d0 
      else
         T_Energy_Liq = vbox(il)/temp_nmol(il)
         LJ_Energy_Liq = (vinterb(il)+vintrab(il)+
     &        vtailb(il))/temp_nmol(il)
         Coul_Energy_Liq = velectb(il)/temp_nmol(il)
      end if
      
      if (temp_nmol(ig).eq.0) then
         T_Energy_Gas = 0.0d0
         LJ_Energy_Gas = 0.0d0
         Coul_energy_Gas = 0.0d0
      else
         T_Energy_Gas = vbox(ig)/temp_nmol(ig)
         LJ_Energy_Gas = (vinterb(ig)+vintrab(ig)+
     &        vtailb(ig))/temp_nmol(ig)
         Coul_Energy_Gas = velectb(ig)/temp_nmol(ig)
      end if
      
      
      enchg1 = 0.008314510d0*(T_Energy_Gas -
     &     T_Energy_Liq)

      pdV = pres(ig)*(mol_vol(ig)-mol_vol(il))*1.0d-6
      
      Heat_vapor_T = enchg1 + pres(ig)*
     &     (mol_vol(ig)-mol_vol(il))*1.0d-6
!     write(iou,1505) il,ig,abs(enthchg1)
!     -- This is inter+intra LJ
      enchg2 = 0.008314510d0*(LJ_Energy_Gas-LJ_Energy_Liq)
      Heat_vapor_LJ = enchg2 + pres(ig)*
     &     (mol_vol(ig)-mol_vol(il))*1.0d-6
!     -- This is Coulomb part
      enchg3 = 0.008314510d0*(Coul_Energy_Gas-Coul_Energy_Liq)
      Heat_vapor_COUL = enchg3 + pres(ig)*
     &     (mol_vol(ig)-mol_vol(il))*1.0d-6
!     -- Calculating Hildebrand solubility parameter (total)
      CED_T = (abs(enchg1)*
     &     1000.0d0*joule2cal)/mol_vol(il)
      if (CED_T.lt.0.0d0) then
         HSP_T = 0.0
      else
         HSP_T = sqrt(CED_T)
      end if
!     write(iou,1508) temp,HSP_TOTAL
!     -- Calculating Hildebrand solubility parameter (LJ part)
      CED_LJ = ((abs(enchg2))*1000.0d0*
     &     joule2cal)/mol_vol(il)
      if (CED_LJ.lt.0.0d0) then
         HSP_LJ =0.0d0
      else
         HSP_LJ = sqrt(CED_LJ)
      end if
!     write(iou,1509) HSP_LJ
!     -- Calculating Hildebrand solubility parameter (Coulomb part)
      CED_COUL = (abs(enchg3)*
     &     1000.0d0*joule2cal)/mol_vol(il)
      if (CED_COUL.lt.0.0d0) then
         HSP_COUL=0.0d0
      else
         HSP_COUL = sqrt(CED_COUL)
      end if
      end subroutine calcsolpar
