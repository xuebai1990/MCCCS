      function vtorso( thetac, itype )

c vtorso
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
c John Stubbs, and Collin Wick and Ilja Siepmann  
c                     
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to 
c
c Free Software Foundation, Inc. 
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer itype
      double precision vtorso, thetac, theta
      double precision tac2,tac3,tac4,tac5,tac6,tac7

c *** common blocks ***
      include 'conver.inc' 
      include 'contorsion.inc'

c ----------------------------------------------------------------

      if ((itype.ge.200).and.(itype.le.202)) then
        if (thetac.gt.1.0d0) thetac = 1.0d0
         if (thetac.lt.-1.0d0) thetac = -1.0d0
         theta = dacos(thetac)
         vtorso = vtt0(itype)*(1.0d0-dcos(2.0d0*theta))

      elseif ((itype.ge.100).and.(itype.le.140)) then
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)+onepi
         vtorso =  vtt0(itype) + vtt1(itype)*dcos(1.0d0*theta)+
     &             vtt2(itype)*dcos(2.0d0*theta)+
     &             vtt3(itype)*dcos(3.0d0*theta)+
     &             vtt4(itype)*dcos(4.0d0*theta)+
     &             vtt5(itype)*dcos(5.0d0*theta)+
     &             vtt6(itype)*dcos(6.0d0*theta)+
     &             vtt7(itype)*dcos(7.0d0*theta)+
     &             vtt8(itype)*dcos(8.0d0*theta)+
     &             vtt9(itype)*dcos(9.0d0*theta)
 

      elseif ( itype .ge. 1 .and. itype .le. 7) then
c - parameters for linear and branched alkane molecules - ALKANE CURRENTLY USED
c - 5 + 6 parameters for alcohols
c - Jorgensen potential
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)+onepi
c --- remember: 1 + cos( theta+onepi ) = 1 - cos( theta )
         vtorso = vtt0(itype) + vtt1(itype)*(1.0d0-thetac) +
     +            vtt2(itype)*(1.d0-dcos(2.d0*theta)) +
     +            vtt3(itype)*(1.d0+dcos(3.d0*theta))

      elseif ( itype .eq. 8 ) then
c - Cummings torsional potential
c - PERFLUOROCARBON CURRENTLY USED starting 10-1-97
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0

         tac2 = thetac*thetac
         tac3 = tac2*thetac
         tac4 = tac3*thetac
         tac5 = tac4*thetac
         tac6 = tac5*thetac
         tac7 = tac6*thetac
c extra displacement is to get the curve above zero
         vtorso = 595.4d0 + 345.0d0
     &        -(282.7d0)*thetac 
     &        +(1355.2d0)*tac2 
     &        +(6800d0)*tac3
     &        -(7875.3d0)*tac4
     &        -(14168.0d0)*tac5
     &        +(9213.7d0)*tac6
     &        +(4123.7d0)*tac7
c         write(iou,*) 'thetac,vtorso',thetac,vtorsoAK
c - Roethlisberger torsional potential for linear perfluorocarbon
c - PERFLUOROCARBON no longer USED
c         if (thetac.gt.1.d0) thetac=1.0d0
c         if (thetac.lt.-1.d0) thetac=-1.0d0
c         theta=dacos(thetac)
c         vtorso = - 269.5616d0 + 503.6229d0*(1.0d0-thetac)
c     &            + 679.921d0*(1.0d0-(dcos(3.0d0*theta)))
c     &            + 3.0085d0*(1.0d0-thetac)**5
c     &            + 420.5883d0*dexp(-30.0d0*theta**2)

      elseif ( itype .eq. 9 ) then
c  Toxvaerd II Torsion potential parms:  to be used for anisotropic alkanes
c  JCP 94, 5650-54 (1991)
         vtorso = 1037.76d0 + 2426.07d0*thetac + 81.64d0*thetac**2
     &    -3129.46d0*thetac**3 -163.28d0*thetac**4 -252.73d0*thetac**5
c      elseif ( itype .ge. 1 .and. itype .le. 3 ) then
c - parameters for improper torsion in carboxylic headgroup
c - C'-C2-O-O'
c         if (thetac.gt.1.d0) thetac=1.0d0
c         if (thetac.lt.-1.d0) thetac=-1.0d0
c         theta=dacos(thetac)+onepi
c         vtorso = vtt2(itype)*(1.d00-dcos(2.d00*theta))
      elseif ( itype .eq. 10 ) then
c        --- dummy torsion just to set up inclusion table right
         vtorso = 0.0d0

      elseif ( itype .eq. 11 .or. itype .eq. 12) then
c - parameters for trans and conformations of double bonds
c - derived by Marcus Martin using data from pcmodel for butene
c   trans torsion paramter 4-14-99 MGM

         if (thetac .gt. 1.0d0) then
            theta = 0.0d0
         elseif (thetac .lt. -1.0d0) then
            theta = onepi
         else
            theta = dacos(thetac)
         endif

         vtorso = vtt0(itype)*(theta - vtt1(itype) )**2.0d0

      elseif ( itype .eq. 13 ) then
c   acetic acid torsional potential H3C--C--O--H
c   modified from J Phys Chem 94, 1683-1686 1990
c   had to divide the potential by 2 and did a bit of trig.
         vtorso = 630.0d0*(1.0d0 - thetac) 
     &        + 1562.4d0*(1.0d0 - thetac*thetac)

      elseif ( itype .eq. 14 ) then
c   acetic acid torsional potential  O==C--O--H
c   modified from J Phys Chem 94, 1683-1686 1990
c   had to divide the potential by 2 and did a bit of trig.
         vtorso = 630.0d0*(1.0d0 + thetac) 
     &        + 1562.4d0*(1.0d0 - thetac*thetac)

      elseif ( itype .eq. 15 ) then
c - Rice torsional potential for linear perfluorocarbons - OLD-FASHIONED
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)
         vtorso = 1784.2812d0 + 1357.1705d0*thetac - 1444.08d0*thetac**2
     &          + 1176.6605d0*thetac**3 + 2888.1600d0*thetac**4 
     &          - 7166.9209d0*thetac**5
     &          + 1684.7600d0*dexp(-12.7176d0*theta**2)

      elseif ( itype .eq. 16 ) then
c - normal parameters for linear molecules - OLD-FASHIONED ALKANE
c - Ryckaert-Bellemans potential
         vtorso = 1116.0d0 + 1462.0d0*thetac - 1578.0d0*thetac**2
     &    - 368.1d0*thetac**3 + 3156.1d0*thetac**4 - 3788.0d0*thetac**5

      elseif ( itype .eq. 19 ) then
c --- methyl group rotations explicit  H-C-C-H McQuarrie
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)
         vtorso = 716.77d0*(1.0d0-dcos(3.0d0*theta))

      elseif (itype .eq. 20) then
c --   methyl group rotation explicit hydrogen model Scott+Scheraga
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)
         vtorso = 853.93d0*(1.0d0-dcos(3.0d0*theta))
      elseif (itype .eq. 21) then
c --   methyl group rotation explicit hydrogen model Scott+Scheraga
c --   designed to be used with fully flexible (ie divided by 3)
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)
         vtorso = 284.64d0*(1.0d0-dcos(3.0d0*theta))
      elseif ( itype .eq. 22) then
c --   torsional motion about the central C-O in ester
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)+onepi
         vtorso = 1253.07*(1.0d0-thetac) + 
     &        1560.08*(1.0d0-dcos(2.0d0*theta))

      elseif ( itype .eq. 23) then
c - Jorgensen potential for segment containing a (H-)-O-C-(CH3)_3 OPLS
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)+onepi
c --- remember: 1 + cos( theta+onepi ) = 1 - cos( theta )
         vtorso = 163.56d0*(1.d00+dcos(3.d0*theta))
         
      elseif ( itype .eq. 24) then
         vtorso = 0.0d0

      elseif ( itype .eq. 25) then
c *** OPLS C-C-O-C for ether paper JCC 1990
         if (thetac.gt.1.d0) thetac=1.0d0
         if (thetac.lt.-1.d0) thetac=-1.0d0
         theta=dacos(thetac)+onepi
         vtorso = 725.35d0*(1.d0 + dcos(theta)) -
     &        163.75d0*(1.d0 - dcos(2.d0*theta)) +
     &        558.2d0*(1.d0 + dcos(3.d0*theta))

      elseif (itype .eq. 26) then
c --   for sp3 carbon to aromatic bond
         if (thetac.gt.1.0d0) thetac = 1.0d0
         if (thetac.lt.-1.0d0) thetac=-1.0d0
         theta=dacos(thetac)
         vtorso = 1344.0d0*(1.0d0-dcos(2.0d0*theta+dacos(-1.0d0)))
      elseif ( itype .eq. 27 ) then
c *** ethylene glycol O-C-C-O torsional potential from Hayashi et al,
c * J. Chem. Soc. Faraday Trans. 1995, 91(1), 31-39.

         if (thetac .gt. 1.0d0) thetac=1.0d0
         if (thetac .lt. -1.0d0) thetac=-1.0d0

c         theta = dacos(thetac)

         vtorso = -123.4d0 -2341.2d0*thetac +1728.3d0*thetac**2
     &        +10788.3d0*thetac**3 -1155.2d0*thetac**4 
     &        -8896.9d0*thetac**5

      elseif ( itype .eq. 28 ) then
c *** polyethylene O-C-C-O torsional potential from amber, testing
c * JMS 7/21/03

         if (thetac .gt. 1.0d0) thetac=1.0d0
         if (thetac .lt. -1.0d0) thetac=-1.0d0

         theta=dacos(thetac)+onepi

c         vtorso = 251.619d0*(1.d0 + dcos(2.d0*theta)) +
c     &        1006.475d0*(1.d0 + dcos(3.d0*theta))

         vtorso = 503.24d0 -
     &        251.62d0*(1.d0 - dcos(2.d0*theta)) +
     &        1006.47d0*(1.d0 + dcos(3.d0*theta))

      elseif ( itype .eq. 29 ) then
c *** polyethylene O-C-C-O torsional potential from Collin, based on Grant Smith's, testing
c * JMS 11/24/03

         if (thetac .gt. 1.0d0) thetac=1.0d0
         if (thetac .lt. -1.0d0) thetac=-1.0d0

         theta = dacos(thetac)
         vtorso = 0.5d0 * ( 950.0d0 *(1.0d0-dcos(theta)) + 
     &        950.0d0 * (1.0d0 - dcos( 2.0d0 * (theta+0.25d0*twopi))))

      elseif (itype .eq. 30) then
c * formic acid O=C-O-H torsion from llnl 4/6/04 jms 
         if (thetac .gt. 1.0d0) thetac = 1.0d0
         if (thetac .lt. -1.0d0) thetac = -1.0d0

c same convention as topmon
c backwards!         vtorso = 2576.5d0*(1.0d0 - dcos(2.0d0*theta))
         vtorso = 1258.0d0*(1.0d0 + dcos(theta))

      elseif (itype .eq. 31) then
c * formic acid H-C-O-H torsion from llnl 4/6/04 jms 
         if (thetac .gt. 1.0d0) thetac = 1.0d0
         if (thetac .lt. -1.0d0) thetac = -1.0d0

c same convention as topmon
c backwards!         vtorso = 1258.0d0*(1.0d0 + dcos(theta))
         vtorso = 2576.5d0*(1.0d0 - dcos(2.0d0*theta))

cc - added 7/12/06 C-C-N-O torsion for nitro group also #61
         elseif (itype .eq.32) then
            if (thetac .gt. 1.0d0) thetac = 1.0d0
            if (thetac .lt. -1.0d0) thetac = -1.0d0

            theta = dacos(thetac)
            vtorso = 69.2d0 - 41.4d0*dcos(theta)-14.5d0*dcos(2*theta) - 
     &           19.1d0*dcos(3*theta) + 8.03d0*dcos(4*theta) - 
     &           2.91d0*dcos(5*theta) + 0.95d0*dcos(6*theta)

cc - added 1/29/07 for N-C-C-C torsion also #70
            elseif (itype .eq. 33) then
               if (thetac .gt. 1.0d0) thetac = 1.0d0
               if (thetac .lt. -1.0d0) thetac = -1.0d0
               theta = dacos(thetac)
               vtorso = 438.0d0 + 481.0*dcos(theta) + 
     &              150.0d0*dcos(2*theta) - 115*dcos(3*theta) - 
     &              0.57*dcos(4*theta) + 0.8*cos(5*theta) - 
     &              0.01*dcos(6*theta)

cc - added 06/27/07 for acrylates
            elseif (itype .ge. 34 .and. itype.le. 46) then
               if (thetac .gt. 1.0d0) thetac = 1.0d0
               if (thetac .lt. -1.0d0) thetac = -1.0d0

               theta = dacos(thetac) + onepi

               vtorso = vtt0(itype) + vtt1(itype)*dcos(theta) +
     &              vtt2(itype)*dcos(2.0d0*theta) + 
     &              vtt3(itype)*dcos(3.0d0*theta) + 
     &              vtt4(itype)*dcos(4.0d0*theta)


      elseif (itype .ge. 48 .and. itype .le. 50) then
c -- torsion from Neimark DMMP JPCA v108, 1435 (2004)
         if (thetac.gt.1.0d0) thetac = 1.0d0
         if (thetac.lt.-1.0d0) thetac=-1.0d0
c -- pi added to torsion because the topmon code is backwards.  Trans
c -- configuration is defined as 0.	 
	 theta=dacos(thetac)+onepi
         vtorso = vtt0(itype)*(1+cos(theta))
     &        + vtt1(itype)*(1+cos(2.*theta))
     &        + vtt2(itype)*(1+cos(3.*theta))
     &        + vtt3(itype)*(1+cos(4.*theta))
     &        + vtt4(itype)*(1+cos(5.*theta))
     &        + vtt5(itype)*(1+cos(6.*theta))
    
      elseif(((itype.ge.51).and.(itype.le.52)).or.(itype.eq.56)) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*(1.0d0-dcos(theta))
     &             + vtt2(itype)*(1.0d0+dcos(theta*2.0d0))
     &             + vtt3(itype)*(1.0d0-dcos(theta*3.0d0))

      elseif (itype .eq. 53 .or. itype .eq. 55) then
         if (thetac.gt.1.0d0) thetac = 1.0d0
         if (thetac.lt.-1.0d0) thetac=-1.0d0
         theta=dacos(thetac)+onepi
         vtorso = vtt0(itype)
     &        + vtt1(itype)*(cos(theta))
     &        + vtt2(itype)*(cos(2.*theta))
     &        + vtt3(itype)*(cos(3.*theta))
     &        + vtt4(itype)*(cos(4.*theta))
     &        + vtt5(itype)*(cos(5.*theta))
     &        + vtt6(itype)*(cos(6.*theta))
     &        + vtt7(itype)*(cos(7.*theta))
     &        + vtt8(itype)*(cos(8.*theta))
     &        + vtt9(itype)*(cos(9.*theta))

      elseif (itype .eq. 54) then
         if (thetac.gt.1.0d0) thetac = 1.0d0
         if (thetac.lt.-1.0d0) thetac=-1.0d0
         theta=dacos(thetac)+onepi
         vtorso = vtt0(itype)
     &        + vtt1(itype)*(1+cos(theta))
     &        + vtt2(itype)*(1+cos(2.*theta))
     &        + vtt3(itype)*(1+cos(3.*theta))
     &        + vtt4(itype)*(1+cos(4.*theta))
     &        + vtt5(itype)*(1+cos(5.*theta))
     &        + vtt6(itype)*(1+cos(6.*theta))


      elseif(itype.eq.101) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*(1.0d0-dcos(theta))
     &             + vtt2(itype)*(1.0d0+dcos(theta*2.0d0))
     &             + vtt3(itype)*(1.0d0-dcos(theta*3.0d0))
 
      elseif(itype.eq.103) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*dcos(theta)
     &             + vtt2(itype)*dcos(theta*2.0d0)
     &             + vtt3(itype)*dcos(theta*3.0d0)
     &             + vtt4(itype)*dcos(theta*4.0d0)
     &             + vtt5(itype)*dcos(theta*5.0d0)
     &             + vtt6(itype)*dcos(theta*6.0d0)
     &             + vtt7(itype)*dcos(theta*7.0d0)
     &             + vtt8(itype)*dcos(theta*8.0d0)
     &             + vtt9(itype)*dcos(theta*9.0d0)
   


      elseif(itype.eq.144) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*(1.0d0-dcos(theta))
     &             + vtt2(itype)*(1.0d0+dcos(theta*2.0d0))
     &             + vtt3(itype)*(1.0d0-dcos(theta*3.0d0))
     &             + vtt4(itype)*(1.0d0+dcos(theta*4.0d0))              


      elseif(itype.eq.145) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*dcos(theta)
     &             + vtt2(itype)*dcos(theta*2.0d0)
     &             + vtt3(itype)*dcos(theta*3.0d0)
     &             + vtt4(itype)*dcos(theta*4.0d0)
     &             + vtt5(itype)*dcos(theta*5.0d0)
     &             + vtt6(itype)*dcos(theta*6.0d0)
     &             + vtt7(itype)*dcos(theta*7.0d0) 

      elseif(itype.eq.146) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*dcos(theta)
     &             + vtt2(itype)*dcos(theta*2.0d0)
     &             + vtt3(itype)*dcos(theta*3.0d0)
     &             + vtt4(itype)*dcos(theta*4.0d0)
     &             + vtt5(itype)*dcos(theta*5.0d0)
     &             + vtt6(itype)*dcos(theta*6.0d0)


      elseif(itype.eq.60 .or. itype.eq.61) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*dcos(theta)
     &             + vtt2(itype)*dcos(theta*2.0d0)
     &             + vtt3(itype)*dcos(theta*3.0d0)


      elseif((itype.ge.65).and.(itype.le.66)) then
          if (thetac.gt.1.0d0) thetac = 1.0d0
          if (thetac.lt.-1.0d0) thetac=-1.0d0
c -- pi added to torsion because the topmon code is backwards.  Trans
c -- configuration is defined as 0.
          theta=dacos(thetac)+onepi
          vtorso =   vtt0(itype)
     &             + vtt1(itype)*(1.0d0+dcos(theta))
     &             + vtt2(itype)*(1.0d0-dcos(theta*2.0d0))
     &             + vtt3(itype)*(1.0d0+dcos(theta*3.0d0))


      elseif ( itype .ge. 70 .and. itype .le. 80) then

c --- OPLS SEVEN PARAMETER FIT

         if (thetac.gt.1.0d0) thetac = 1.0d0
         if (thetac.lt.-1.0d0) thetac=-1.0d0

         theta = dacos(thetac)

         vtorso = vtt0(itype) + vtt1(itype)*thetac
     &        + vtt2(itype)*dcos(theta*2.0d0)
     &        + vtt3(itype)*dcos(theta*3.0d0)
     &        + vtt4(itype)*dcos(theta*4.0d0)
     &        + vtt5(itype)*dcos(theta*5.0d0)
     &        + vtt6(itype)*dcos(theta*6.0d0)



      else
         write(6,*) 'you picked a non-defined torsional type'
         stop
      endif

      return
      end
