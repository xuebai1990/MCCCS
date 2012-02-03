      subroutine  calcsolpar(pres,Heat_vapor_T,Heat_vapor_LJ,
     &  Heat_vapor_COUL, pdV, CED_T, CED_LJ,CED_COUL,HSP_T,HSP_LJ,
     &  HSP_COUL,ibox,jbox)
c ***********************************************************************************
c  Written by Neeraj Rai. Date 08/22/2006

c  This subroutine calculates heat of vaporization and solubility parameter
c  Separates the Coulombic and LJ part of solubility parameter
c Previos implementation of heat of vaporization was good only for single component
c system. Now we can do any number of components.
c  It is being written for two box gibbs ensemble simulations. If it is more than two boxes it will c return without doing anything.
c or if the box is solid and non rectangular it will return.
c If you have more than two boxes then there are couple changes that need to be done
c in order for this subroutine to work. nprop in blkavg.inc needs to be set properly
c ************************************************************************************
      implicit none

      include 'control.inc'
      include 'inputdata.inc'
      include 'blkavg.inc'
      include 'coord.inc'
      include 'ensemble.inc'
      include 'system.inc'
      include 'cell.inc'

      integer ibox, jbox, ig, il, imolty
      integer, dimension(nbxmax):: temp_nmol,box_volume
      integer, dimension(nbxmax,ntmax):: acnbox
      double precision, dimension(nbxmax)::mol_vol      
      double precision, dimension(nbxmax)::pres
      double precision enchg1, enchg2,enchg3
      double precision Heat_vapor_T, Heat_vapor_LJ, Heat_vapor_COUL
      double precision CED_T, HSP_T, CED_LJ, HSP_LJ
      double precision CED_COUL,HSP_COUL
      double precision cal2joule, joule2cal
      double precision pdV

      double precision T_Energy_Liq, T_Energy_Gas
      double precision LJ_Energy_Liq, LJ_Energy_Gas
      double precision Coul_Energy_Liq, Coul_Energy_Gas


c     -- This is different than the conventional 4.184
      cal2joule = 4.184d0
      joule2cal = 1.0d0/cal2joule

      temp_nmol(ibox) = 0

      do imolty = 1,nmolty
         temp_nmol(ibox) = temp_nmol(ibox) + ncmt(ibox,imolty)
      enddo

      if (lsolid(ibox).and.(.not.lrect(ibox))) then
         box_volume(ibox) = cell_vol(ibox)
      else   
         box_volume(ibox) = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
      endif
c     molar volume in cc/mol    
      if (temp_nmol(ibox).eq.0) then    
c     set the molar volume artibrarily high as it will not affect any since pressure will be zero
         mol_vol(ibox) = 10.0d10    
      else
         mol_vol(ibox)=box_volume(ibox)/temp_nmol(ibox)*0.6022d0 
      endif

      temp_nmol(jbox) = 0

      do imolty = 1,nmolty
         temp_nmol(jbox) = temp_nmol(jbox) + ncmt(jbox,imolty)
      enddo

      if (lsolid(jbox).and.(.not.lrect(jbox))) then
         box_volume(jbox) = cell_vol(jbox)
      else
         box_volume(jbox) = boxlx(jbox)*boxly(jbox)*boxlz(jbox)
      endif
c     -- molar volume in cc/mol
      if (temp_nmol(jbox).eq.0) then
         mol_vol(jbox) = 10.0D6
      else
         mol_vol(jbox)=box_volume(jbox)/temp_nmol(jbox)*0.6022d0
      endif

      if ( mol_vol(ibox).gt.mol_vol(jbox) ) then
         ig = ibox
         il = jbox
      else
         ig = jbox
         il = ibox
      endif
c     -- This is total energy
      if (temp_nmol(il).eq.0) then
         T_Energy_Liq = 0.0d0
         LJ_Energy_Liq = 0.0d0
         Coul_energy_Liq = 0.0d0 
      else
         T_Energy_Liq = vbox(il)/temp_nmol(il)
         LJ_Energy_Liq = (vinterb(il)+vintrab(il)+
     &        vtailb(il))/temp_nmol(il)
         Coul_Energy_Liq = velectb(il)/temp_nmol(il)
      endif
      
      if (temp_nmol(ig).eq.0) then
         T_Energy_Gas = 0.0d0
         LJ_Energy_Gas = 0.0d0
         Coul_energy_Gas = 0.0d0
      else
         T_Energy_Gas = vbox(ig)/temp_nmol(ig)
         LJ_Energy_Gas = (vinterb(ig)+vintrab(ig)+
     &        vtailb(ig))/temp_nmol(ig)
         Coul_Energy_Gas = velectb(ig)/temp_nmol(ig)
      endif
      
      
      enchg1 = 0.008314510d0*(T_Energy_Gas -
     &     T_Energy_Liq)

      pdV = pres(ig)*(mol_vol(ig)-mol_vol(il))*1.0d-6
      
      Heat_vapor_T = enchg1 + pres(ig)*
     &     (mol_vol(ig)-mol_vol(il))*1.0d-6
c     write(iou,1505) il,ig,abs(enthchg1)
c     -- This is inter+intra LJ
      enchg2 = 0.008314510d0*(LJ_Energy_Gas-LJ_Energy_Liq)
      Heat_vapor_LJ = enchg2 + pres(ig)*
     &     (mol_vol(ig)-mol_vol(il))*1.0d-6
c     -- This is Coulomb part
      enchg3 = 0.008314510d0*(Coul_Energy_Gas-Coul_Energy_Liq)
      Heat_vapor_COUL = enchg3 + pres(ig)*
     &     (mol_vol(ig)-mol_vol(il))*1.0d-6
c     -- Calculating Hildebrand solubility parameter (total)
      CED_T = (abs(enchg1)*
     &     1000.0d0*joule2cal)/mol_vol(il)
      if (CED_T.lt.0.0d0) then
         HSP_T = 0.0
      else
         HSP_T = sqrt(CED_T)
      endif
c     write(iou,1508) temp,HSP_TOTAL
c     -- Calculating Hildebrand solubility parameter (LJ part)
      CED_LJ = ((abs(enchg2))*1000.0d0*
     &     joule2cal)/mol_vol(il)
      if (CED_LJ.lt.0.0d0) then
         HSP_LJ =0.0d0
      else
         HSP_LJ = sqrt(CED_LJ)
      endif
c     write(iou,1509) HSP_LJ
c     -- Calculating Hildebrand solubility parameter (Coulomb part)
      CED_COUL = (abs(enchg3)*
     &     1000.0d0*joule2cal)/mol_vol(il)
      if (CED_COUL.lt.0.0d0) then
         HSP_COUL=0.0d0
      else
         HSP_COUL = sqrt(CED_COUL)
      endif
      end subroutine calcsolpar
