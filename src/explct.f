      subroutine explct(ichain,vmethyl,lcrysl,lswitch)

!     adds H-atoms to a linear carbon chain
!     Potential for methyl-group rotation
!     is: V = 0.5*E0*(1-cos(3*alpha)), where
!     alpha is Ryckaert torsion angle!
!     ------

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
!$$$      include 'connect.inc'

! - arguments
      integer(KIND=int)::ichain,nngrow,negrow,i,iplus,imins,nn,
     &        imolty,iben,iend,ii,jj
      real(KIND=double_precision)::vmethyl,ch,cc,cch,ca,ah,hch2,hk,ck,en0
     &     ,a1,b1,c1,dln1,a2,b2,c2,dln2,a3,b3,c3,x12,y12,z12
     &     ,x32,y32,z32,xa,ya,za
     &     ,r,rx,ry,rz,dr,ven,random,prob,a4,b4,c4,rn
     &     ,hch,coh,ce,ratio
      real(KIND=double_precision)::oa,hoh,hoh2,oh,ok,om,oc,v1,v2,alpha,rr
      logical::lcrysl,londone,lswitch,lalkanol

! - find intramolecular structure 
      imolty = moltyp(ichain)
      nngrow = nugrow(imolty)
      negrow = nunit(imolty) - 3
      vmethyl = 0.0d0
      if ((nunit(imolty)-6) .eq. 3*(nngrow-3) ) then
! *** alkanol cases
         lalkanol = .true.
         iend = 2
      else
         lalkanol = .false.
         iend = 1
      end if


      if ( nunit(imolty) .eq. 3 .or. nunit(imolty) .eq. 4) then

! - Water case
        
         if ( nngrow .eq. 3 ) then

! - TIP-4P geometry with 2H and O as the growing unit
! - M site has been determined by their positions through geometry 
! - constraint

            om = brvib(itvib(imolty,4,1))
            a1 = 0.5d0*(rxu(ichain,2)+rxu(ichain,3)) 
     &           - rxu(ichain,1)
            b1 = 0.5d0*(ryu(ichain,2)+ryu(ichain,3)) 
     &           - ryu(ichain,1)
            c1 = 0.5d0*(rzu(ichain,2)+rzu(ichain,3)) 
     &           - rzu(ichain,1)
            dr = dsqrt(a1*a1+b1*b1+c1*c1)
            rxu(ichain,4) = rxu(ichain,1) + om*a1/dr
            ryu(ichain,4) = ryu(ichain,1) + om*b1/dr
            rzu(ichain,4) = rzu(ichain,1) + om*c1/dr
            return
            
         elseif ( nngrow .eq. 1 ) then

            if ( itvib(imolty,nngrow+1,1) .eq. 
     &           itvib(imolty,nngrow+2,1) ) then
               oh = brvib(itvib(imolty,nngrow+1,1))
               if ( nunit(imolty) .eq. 3 ) then
! - SPC geometry

                  hoh = brben(itben(imolty,nngrow+1,1))
               else
! - TIP-4P geometry
            
                  hoh = brben(itben(imolty,nngrow+1,1))
               end if
               if ( lcrysl ) then
                  hoh2 = hoh / 2.0d0
                  hk = oh*dsin(hoh2)
                  ok = oh*dcos(hoh2)
                  rxu(ichain,2) = rxu(ichain,1)
                  ryu(ichain,2) = ryu(ichain,1) + hk
                  rzu(ichain,2) = rzu(ichain,1) + ok
                  rxu(ichain,3) = rxu(ichain,1)
                  ryu(ichain,3) = ryu(ichain,1) - hk
                  rzu(ichain,3) = rzu(ichain,1) + ok
               else
                  oa = oh*dcos(onepi-hoh)
                  ah = oh*dsin(onepi-hoh)
! --- generate a random vector on a sphere for the first H ---
 111              rx = 2.0d0*random() - 1.0d0
                  ry = 2.0d0*random() - 1.0d0
                  dr = rx*rx + ry*ry
                  if ( dr .gt. 1.0d0 ) goto 111
                  rz = 2.0d0*dsqrt(1-dr)
                  a2 = rx*rz
                  b2 = ry*rz
                  c2 = 1 - 2.0d0*dr

                  rxu(ichain,2) = rxu(ichain,1) - oh*a2
                  ryu(ichain,2) = ryu(ichain,1) - oh*b2
                  rzu(ichain,2) = rzu(ichain,1) - oh*c2
! --- generate another random vector on a sphere for the second H ---
 222              rx = 2.0d0*random() - 1.0d0
                  ry = 2.0d0*random() - 1.0d0
                  dr = rx*rx + ry*ry
                  if ( dr .gt. 1.0d0 ) goto 222
                  rz = 2.0d0*dsqrt(1-dr)
                  rx = rx*rz
                  ry = ry*rz
                  rz = 1 - 2.0d0*dr
                
! --- The two vectors above form a plane identified by n1 ---
                  a1 = b2*rz - c2*ry
                  b1 = -(a2*rz - rx*c2)
                  c1 = a2*ry - rx*b2
                  dln1 = dsqrt(a1*a1 + b1*b1 + c1*c1)
                  a1 = a1/dln1
                  b1 = b1/dln1
                  c1 = c1/dln1

! --- cross product n1 x n2 calc.-> n3 ---
            
                  a3 = b1*c2 - b2*c1
                  b3 = -(a1*c2 - a2*c1)
                  c3 = a1*b2 - a2*b1
                  rxu(ichain,3) = rxu(ichain,1) + a2*oa + ah*a3
                  ryu(ichain,3) = ryu(ichain,1) + b2*oa + ah*b3
                  rzu(ichain,3) = rzu(ichain,1) + c2*oa + ah*c3
               end if
               if ( nunit(imolty) .eq. 4 ) then
                  om = brvib(itvib(imolty,4,1))
                  a1 = 0.5d0*(rxu(ichain,2)+rxu(ichain,3)) 
     &                 - rxu(ichain,1)
                  b1 = 0.5d0*(ryu(ichain,2)+ryu(ichain,3)) 
     &                 - ryu(ichain,1)
                  c1 = 0.5d0*(rzu(ichain,2)+rzu(ichain,3)) 
     &                 - rzu(ichain,1)
                  dr = dsqrt(a1*a1+b1*b1+c1*c1)
                  rxu(ichain,4) = rxu(ichain,1) + om*a1/dr
                  ryu(ichain,4) = ryu(ichain,1) + om*b1/dr
                  rzu(ichain,4) = rzu(ichain,1) + om*c1/dr
               end if
               return
            else
! --- HF model
! --- generate a random vector on a sphere for the H and M sites ---
 333           rx = 2.0d0*random() - 1.0d0
               ry = 2.0d0*random() - 1.0d0
               dr = rx*rx + ry*ry
               if ( dr .gt. 1.0d0 ) goto 333
               rz = 2.0d0*dsqrt(1-dr)
               a2 = rx*rz
               b2 = ry*rz
               c2 = 1 - 2.0d0*dr

               do i = 2,3
                  ch = brvib(itvib(imolty,i,1))
                  rxu(ichain,i) = rxu(ichain,1) + ch*a2
                  ryu(ichain,i) = ryu(ichain,1) + ch*b2
                  rzu(ichain,i) = rzu(ichain,1) + ch*c2
               end do
            end if 
         end if
         
      elseif ( nunit(imolty) .eq. 5 .or. nunit(imolty) .eq. 7
     &        .or. nunit(imolty) .eq. 9 ) then
         

! - Methane case or other rigid molecule with 5 units

         ch = brvib(itvib(imolty,nngrow+1,1))
         hch = brben(itben(imolty,nngrow+1,1))
         if ( nngrow .eq. 3 ) then

! *** for hydro-furan

!            a1 = rxu(ichain,3)-rxu(ichain,1)
!            b1 = ryu(ichain,3)-ryu(ichain,1)
!            c1 = rzu(ichain,3)-rzu(ichain,1)
!            dln1 = dsqrt(a1*a1 + b1*b1 + c1*c1)
!            a1 = a1/dln1
!            b1 = b1/dln1
!            c1 = c1/dln1
!            a2 = 0.5d0*(rxu(ichain,1)+rxu(ichain,3))-
!     &           rxu(ichain,2)
!            b2 = 0.5d0*(ryu(ichain,1)+ryu(ichain,3))-
!     &           ryu(ichain,2)
!            c2 = 0.5d0*(rzu(ichain,1)+rzu(ichain,3))-
!     &           rzu(ichain,2)
!            dln1 = dsqrt(a2*a2 + b2*b2 + c2*c2)
!            a2 = a2/dln1
!            b2 = b2/dln1
!            c2 = c2/dln1
!            rxu(ichain,4) = rxu(ichain,2) + 2.2758059d0*a2 + 
!     &           0.765d0*a1
!            ryu(ichain,4) = ryu(ichain,2) + 2.2758059d0*b2 + 
!     &           0.765d0*b1
!            rzu(ichain,4) = rzu(ichain,2) + 2.2758059d0*c2 + 
!     &           0.765d0*c1
!            rxu(ichain,5) = rxu(ichain,4) - 1.53d0*a1
!            ryu(ichain,5) = ryu(ichain,4) - 1.53d0*b1
!            rzu(ichain,5) = rzu(ichain,4) - 1.53d0*c1            

            a1 = rxu(ichain,2)-rxu(ichain,1)
            b1 = ryu(ichain,2)-ryu(ichain,1)
            c1 = rzu(ichain,2)-rzu(ichain,1)
            
            a2 = rxu(ichain,3)-rxu(ichain,2)
            b2 = ryu(ichain,3)-ryu(ichain,2)
            c2 = rzu(ichain,3)-rzu(ichain,2)

!     -------cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(a1*c2 - a2*c1)
            c3 = a1*b2 - a2*b1
            dln1 = dsqrt(a3*a3 + b3*b3 + c3*c3)
            a2 = a3/dln1
            b2 = b3/dln1
            c2 = c3/dln1

            a1 = 0.5d0*(rxu(ichain,1)+rxu(ichain,3))-
     &           rxu(ichain,2)
            b1 = 0.5d0*(ryu(ichain,1)+ryu(ichain,3))-
     &           ryu(ichain,2)
            c1 = 0.5d0*(rzu(ichain,1)+rzu(ichain,3))-
     &           rzu(ichain,2)
            dln1 = dsqrt(a1*a1 + b1*b1 + c1*c1)
            a1 = a1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
            a3 = b1*c2 - b2*c1
            b3 = -(a1*c2 - a2*c1)
            c3 = a1*b2 - a2*b1
           
            rxu(ichain,4) = rxu(ichain,2) + 2.25d0*a1 + 
     &           0.41d0*a2+0.77d0*a3
            ryu(ichain,4) = ryu(ichain,2) + 2.25d0*b1 + 
     &           0.41d0*b2+0.77d0*b3
            rzu(ichain,4) = rzu(ichain,2) + 2.25d0*c1 + 
     &           0.41d0*c2+0.77d0*c3
            rxu(ichain,5) = rxu(ichain,4) - 1.54d0*a3
            ryu(ichain,5) = ryu(ichain,4) - 1.54d0*b3
            rzu(ichain,5) = rzu(ichain,4) - 1.54d0*c3
         elseif ( nugrow(imolty) .eq. 2) then
            
            ca = ch*dcos(onepi-hch)
            ah = ch*dsin(onepi-hch)

! - only one H participate in the growing
! - this subroutine put on the rest 3 hydrogens
!     -------H-atoms for the first CH3-group---
!     -------first define vector C-H ------
            x12 = rxu(ichain,2)-rxu(ichain,1)
            y12 = ryu(ichain,2)-ryu(ichain,1)
            z12 = rzu(ichain,2)-rzu(ichain,1)
!     ------ generate a random vector ------
 444        rx = 2.0d0*random() - 1.0d0
            ry = 2.0d0*random() - 1.0d0
            dr = rx*rx + ry*ry
            if ( dr .gt. 1.0d0 ) goto 444
            rz = 2.0d0*dsqrt(1-dr)
            rx = rx*rz
            ry = ry*rz
            rz = 1 - 2.0d0*dr

!     ------ the two vectors above form a plane identified by n1 ------
            a1 = y12*rz - z12*ry
            b1 = -(x12*rz - rx*z12)
            c1 = x12*ry - rx*y12
            dln1 = dsqrt(a1**2 + b1**2 + c1**2)
            a1 = a1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
!     --------normalizing C-H vector -> n2 ----
            a2 = x12/ch
            b2 = y12/ch
            c2 = z12/ch
!     -------cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(a1*c2 - a2*c1)
            c3 = a1*b2 - a2*b1
!     ------- point A ------ 
            xa = rxu(ichain,2)+a2*ca
            ya = ryu(ichain,2)+b2*ca
            za = rzu(ichain,2)+c2*ca
!     ------H-atoms for the first methyl group------
            r = 0.0d0
            do i = 1,3
               a4 = ah*a3*dcos(r) + ah*a1*dsin(r)
               b4 = ah*b3*dcos(r) + ah*b1*dsin(r)
               c4 = ah*c3*dcos(r) + ah*c1*dsin(r)
               rxu(ichain,nngrow+i) = xa + a4 
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0d0*onepi/180.0d0
            end do

         else
! *** only carbon has been grown

!     ------ METHANE ------
            if ( lcrysl ) then
               hch2 = brben(itben(imolty,2,1))/2
               hk = ch*dsin(hch2)
               ck = ch*dcos(hch2)
               rxu(ichain,2) = rxu(ichain,1) - ck
               ryu(ichain,2) = ryu(ichain,1)
               rzu(ichain,2) = rzu(ichain,1) + hk               
               rxu(ichain,3) = rxu(ichain,1) - ck
               ryu(ichain,3) = ryu(ichain,1)
               rzu(ichain,3) = rzu(ichain,1) - hk
               hch2 = brben(itben(imolty,4,1))/2
               ch = brvib(itvib(imolty,4,1))
               hk = ch*dsin(hch2)
               ck = ch*dcos(hch2)
               rxu(ichain,4) = rxu(ichain,1) + ck
               ryu(ichain,4) = ryu(ichain,1) + hk
               rzu(ichain,4) = rzu(ichain,1) 
               rxu(ichain,5) = rxu(ichain,1) + ck
               ryu(ichain,5) = ryu(ichain,1) - hk
               rzu(ichain,5) = rzu(ichain,1) 

! - KEEP THE OLD CONFIGURATION in SWITCH MOVE

            elseif ( lswitch ) then
!               if (imolty .eq. 1) then
!                  ratio = brvib(itvib(1,nngrow+1,1))
!     &                 /brvib(itvib(2,nngrow+1,1))
!               else
!                  ratio = brvib(itvib(2,nngrow+1,1))
!     &                 /brvib(itvib(1,nngrow+1,1))
!               end if
               ratio = 1.0d0
               rxu(ichain,2) = rxu(ichain,1) + ratio*(rxu(ichain,2)
     &              -rxu(ichain,1))
               ryu(ichain,2) = ryu(ichain,1) + ratio*(ryu(ichain,2)
     &              -ryu(ichain,1))
               rzu(ichain,2) = rzu(ichain,1) + ratio*(rzu(ichain,2)
     &              -rzu(ichain,1))
               rxu(ichain,3) = rxu(ichain,1) + ratio*(rxu(ichain,3)
     &              -rxu(ichain,1))
               ryu(ichain,3) = ryu(ichain,1) + ratio*(ryu(ichain,3)
     &              -ryu(ichain,1))
               rzu(ichain,3) = rzu(ichain,1) + ratio*(rzu(ichain,3)
     &              -rzu(ichain,1))
               ratio = 0.7d0/0.25d0
               rxu(ichain,4) = rxu(ichain,1) + ratio*(rxu(ichain,4)
     &              -rxu(ichain,1))
               ryu(ichain,4) = ryu(ichain,1) + ratio*(ryu(ichain,4)
     &              -ryu(ichain,1))
               rzu(ichain,4) = rzu(ichain,1) + ratio*(rzu(ichain,4)
     &              -rzu(ichain,1))
               rxu(ichain,5) = rxu(ichain,1) + ratio*(rxu(ichain,5)
     &              -rxu(ichain,1))
               ryu(ichain,5) = ryu(ichain,1) + ratio*(ryu(ichain,5)
     &              -ryu(ichain,1))
               rzu(ichain,5) = rzu(ichain,1) + ratio*(rzu(ichain,5)
     &              -rzu(ichain,1))
            else
               ch = brvib(itvib(imolty,2,1))
               hch2 = brben(itben(imolty,2,1))/2
               ca = ch*dcos(hch2)
               ah = ch*dsin(hch2)         
!     ------ generate a random vector for the first H ------
 555           rx = 2.0d0*random() - 1.0d0
               ry = 2.0d0*random() - 1.0d0
               dr = rx*rx + ry*ry
               if ( dr .gt. 1.0d0 ) goto 555
               rz = 2.0d0*dsqrt(1-dr)
               a2 = rx*rz
               b2 = ry*rz
               c2 = 1 - 2.0d0*dr
!     ------ generate another random vector for the rotation ------
 666           rx = 2.0d0*random() - 1.0d0
               ry = 2.0d0*random() - 1.0d0
               dr = rx*rx + ry*ry
               if ( dr .gt. 1.0d0 ) goto 666
               rz = 2.0d0*dsqrt(1-dr)
               rx = rx*rz
               ry = ry*rz
               rz = 1 - 2.0d0*dr
!     ------ the two vectors above form a plane identified by n1 ------
               a1 = b2*rz - c2*ry
               b1 = -(a2*rz - rx*c2)
               c1 = a2*ry - rx*b2
               dln1 = dsqrt(a1**2 + b1**2 + c1**2)
               a1 = a1/dln1
               b1 = b1/dln1
               c1 = c1/dln1

               rxu(ichain,2) = rxu(ichain,1) + a2*ca + a1*ah 
               ryu(ichain,2) = ryu(ichain,1) + b2*ca + b1*ah
               rzu(ichain,2) = rzu(ichain,1) + c2*ca + c1*ah
               rxu(ichain,3) = rxu(ichain,1) + a2*ca - a1*ah 
               ryu(ichain,3) = ryu(ichain,1) + b2*ca - b1*ah
               rzu(ichain,3) = rzu(ichain,1) + c2*ca - c1*ah

               ch = brvib(itvib(imolty,4,1))
               hch2 = brben(itben(imolty,4,1))/2
               ca = ch*dcos(hch2)
               ah = ch*dsin(hch2)         

            
!     -------cross product n1 x n2 calc.-> n3 ---
               a3 = b1*c2 - b2*c1
               b3 = -(a1*c2 - a2*c1)
               c3 = a1*b2 - a2*b1

               rxu(ichain,4) = rxu(ichain,1) - a2*ca + a3*ah 
               ryu(ichain,4) = ryu(ichain,1) - b2*ca + b3*ah
               rzu(ichain,4) = rzu(ichain,1) - c2*ca + c3*ah
               rxu(ichain,5) = rxu(ichain,1) - a2*ca - a3*ah 
               ryu(ichain,5) = ryu(ichain,1) - b2*ca - b3*ah
               rzu(ichain,5) = rzu(ichain,1) - c2*ca - c3*ah
               
            end if
            if ( nunit(imolty) .eq. 9 ) then
               ce = brvib(itvib(imolty,6,1))
               ratio = ce / ch
               do i = 1,4
                  ii = i + 1
                  jj = i + 5
                  rxu(ichain,jj) = rxu(ichain,1) + ratio*
     &                 (rxu(ichain,ii)-rxu(ichain,1))
                  ryu(ichain,jj) = ryu(ichain,1) + ratio*
     &                 (ryu(ichain,ii)-ryu(ichain,1))
                  rzu(ichain,jj) = rzu(ichain,1) + ratio*
     &                 (rzu(ichain,ii)-rzu(ichain,1))
               end do
            end if
         end if
         if ( nunit(imolty) .eq. 7 ) then
! *** for seven site water model
            rxu(ichain,6) = 
     &           4d0*rxu(ichain,4)-3d0*rxu(ichain,1)
            ryu(ichain,6) = 
     &           4d0*ryu(ichain,4)-3d0*ryu(ichain,1)
            rzu(ichain,6) = 
     &           4d0*rzu(ichain,4)-3d0*rzu(ichain,1)
            rxu(ichain,7) = 
     &           4d0*rxu(ichain,5)-3d0*rxu(ichain,1)
            ryu(ichain,7) = 
     &           4d0*ryu(ichain,5)-3d0*ryu(ichain,1)
            rzu(ichain,7) = 
     &           4d0*rzu(ichain,5)-3d0*rzu(ichain,1)
         end if
      else

! - ALKANES (not for methane)          
! - find intramolecular structure for longer alkanes
! - WARNING: work only for pure hydro- or perfluoro-carbons!!!

         cc = brvib(itvib(imolty,1,1))
         ch = brvib(itvib(imolty,nngrow+1,1))

         cch = brben(itben(imolty,nngrow+1,1))
         ca = ch*dcos(onepi-cch)
         ah = ch*dsin(onepi-cch)
         if ( nngrow .gt. 2) then
            hch = brben(itben(imolty,nngrow+4,1))
            hch2 = hch/2.0d0
            hk = ch*dsin(hch2)
            ck = ch*dcos(hch2)
            en0 = 853.93d0

!     ------ REGULAR N-ALKANE ------   
!     ------ main loop
!     calculates hydrogen positions for methylene groups
!     WARNING only works for linear alkanes
            do i = 2,nngrow-iend
               iplus = i+1
               imins = i-1
               a1 = (ryu(ichain,imins)-ryu(ichain,i))*(rzu(ichain,iplus)
     &              -rzu(ichain,i)) - (rzu(ichain,imins)-rzu(ichain,i))*
     &              (ryu(ichain,iplus)-ryu(ichain,i))
              b1 = -(rxu(ichain,imins)-rxu(ichain,i))*(rzu(ichain,iplus)
     &              -rzu(ichain,i)) + (rzu(ichain,imins)-rzu(ichain,i))*
     &              (rxu(ichain,iplus)-rxu(ichain,i))
               c1 = (rxu(ichain,imins)-rxu(ichain,i))*(ryu(ichain,iplus)
     &              -ryu(ichain,i)) - (ryu(ichain,imins)-ryu(ichain,i))*
     &              (rxu(ichain,iplus)-rxu(ichain,i))
               dln1 = dsqrt(a1**2 + b1**2 + c1**2)
               a1 = a1/dln1
               b1 = b1/dln1
               c1 = c1/dln1
               a2 = (rxu(ichain,iplus)-rxu(ichain,imins))
               b2 = (ryu(ichain,iplus)-ryu(ichain,imins))
               c2 = (rzu(ichain,iplus)-rzu(ichain,imins))
               dln2 = dsqrt(a2**2 + b2**2 + c2**2)
               a2 = a2/dln2
               b2 = b2/dln2
               c2 = c2/dln2
               a3 = b1*c2 - c1*b2
               b3 = -a1*c2 + a2*c1
               c3 = a1*b2 - b1*a2
               londone = .false.
               do nn = nngrow+4+2*(i-2), nngrow+5+2*(i-2)
                  if (.not.londone) then
                     rxu(ichain,nn) = rxu(ichain,i) + a3*ck + a1*hk
                     ryu(ichain,nn) = ryu(ichain,i) + b3*ck + b1*hk
                     rzu(ichain,nn) = rzu(ichain,i) + c3*ck + c1*hk
                     londone = .true.
                  else
                     rxu(ichain,nn) = rxu(ichain,i) + a3*ck - a1*hk
                     ryu(ichain,nn) = ryu(ichain,i) + b3*ck - b1*hk
                     rzu(ichain,nn) = rzu(ichain,i) + c3*ck - c1*hk
                  end if
               end do
            end do

!     -------H-atoms for the first CH3-group---
!     -------define c1c2c3 plane -> n1 ------
            x12 = rxu(ichain,1)-rxu(ichain,2)
            y12 = ryu(ichain,1)-ryu(ichain,2)
            z12 = rzu(ichain,1)-rzu(ichain,2)
            x32 = rxu(ichain,3)-rxu(ichain,2)
            y32 = ryu(ichain,3)-ryu(ichain,2)
            z32 = rzu(ichain,3)-rzu(ichain,2)
            a1 = y12*z32 - z12*y32
            b1 = -(x12*z32 - x32*z12)
            c1 = x12*y32 - x32*y12
            dln1 = dsqrt(a1**2 + b1**2 + c1**2)
            a1 = a1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
!     --------normalizing c2c1 vector -> n2 ----
            a2 = x12/cc
            b2 = y12/cc
            c2 = z12/cc
!     -------cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(a1*c2 - a2*c1)
            c3 = a1*b2 - a2*b1
!     ------- point A ------ 
            xa = rxu(ichain,1)+a2*ca
            ya = ryu(ichain,1)+b2*ca
            za = rzu(ichain,1)+c2*ca
!     ------H-atoms ------
            if ( lcrysl ) then
               r=0.0d0
               ven = 0.0d0
            else
 101           r = 2.0d0*onepi*random()
               ven = en0*(1.0d0-dcos(3.0d0*r))
               prob = dexp(-beta*ven)
               rn = random()
               if (rn.gt.prob) goto 101
            end if
!     rgr = r*180.0d0/onepi
!     write(iou,*) 'final Ryckaert angle: ',rgr
            vmethyl = vmethyl + ven
            do i = 1,3
               a4 = ah*a3*dcos(r+onepi) + ah*a1*dsin(r+onepi)
               b4 = ah*b3*dcos(r+onepi) + ah*b1*dsin(r+onepi)
               c4 = ah*c3*dcos(r+onepi) + ah*c1*dsin(r+onepi)
               
               rxu(ichain,nngrow+i) = xa + a4 
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0d0*onepi/180.0d0
            end do

! --- if it is an alkanol molecule, it did not have the other ending CH3

            if ( lalkanol ) return

!     ------H-atoms for the last CH3------
!     -------define c1c2c3 plane -> n1 ------
            x12 = rxu(ichain,nngrow)-rxu(ichain,nngrow-1)
            y12 = ryu(ichain,nngrow)-ryu(ichain,nngrow-1)
            z12 = rzu(ichain,nngrow)-rzu(ichain,nngrow-1)
            x32 = rxu(ichain,nngrow-2)-rxu(ichain,nngrow-1)
            y32 = ryu(ichain,nngrow-2)-ryu(ichain,nngrow-1)
            z32 = rzu(ichain,nngrow-2)-rzu(ichain,nngrow-1)
            a1 = y12*z32 - z12*y32
            b1 = -(x12*z32 - x32*z12)
            c1 = x12*y32 - x32*y12
            dln1 = dsqrt(a1**2 + b1**2 + c1**2)
            a1 = a1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
!     --------normalizing c2c1 vector -> n2 ----
            a2 = x12/cc
            b2 = y12/cc
            c2 = z12/cc
!     -------cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(a1*c2 - a2*c1)
            c3 = a1*b2 - a2*b1
!     ------- point A ------ 
            xa = rxu(ichain,nngrow)+a2*ca
            ya = ryu(ichain,nngrow)+b2*ca
            za = rzu(ichain,nngrow)+c2*ca
!     ---- H-atoms -------
            if (lcrysl) then
               r = 0.0d0
               ven = 0.0d0
            else
 200           r = 2.0d0*onepi*random()
               ven = en0*(1.0d0-dcos(3.0d0*r))
               prob = dexp(-beta*ven)
               rn = random()
               if (rn.gt.prob) goto 200
            end if
            vmethyl = vmethyl + ven
            do i = 1,3
               a4 = ah*a3*dcos(r+onepi) + ah*a1*dsin(r+onepi)
               b4 = ah*b3*dcos(r+onepi) + ah*b1*dsin(r+onepi)
               c4 = ah*c3*dcos(r+onepi) + ah*c1*dsin(r+onepi)
               rxu(ichain,negrow+i) = xa + a4 
               ryu(ichain,negrow+i) = ya + b4
               rzu(ichain,negrow+i) = za + c4
               r = r + 120.0d0*onepi/180.0d0
            end do
            
         else if (nngrow .eq. 2 .and. nunit(imolty) .eq. 8) then
!     ------ ETHANE  -----
!     ------ define the new type en0 ------
!     use en0 for torsion type 19
            en0 = 716.77d0 
!     ------ for ethane case ------
!     -------H-atoms for the first CH3-group---
!     -------first define vector C2C1 ------
            x12 = rxu(ichain,1)-rxu(ichain,2)
            y12 = ryu(ichain,1)-ryu(ichain,2)
            z12 = rzu(ichain,1)-rzu(ichain,2)
!     ------ generate a random vector ------
 777        rx = 2.0d0*random() - 1.0d0
            ry = 2.0d0*random() - 1.0d0
            dr = rx*rx + ry*ry
            if ( dr .gt. 1.0d0 ) goto 777
            rz = 2.0d0*dsqrt(1-dr)
            rx = rx*rz
            ry = ry*rz
            rz = 1 - 2.0d0*dr

!     ------ the two vectors above form a plane identified by n1 ------
            a1 = y12*rz - z12*ry
            b1 = -(x12*rz - rx*z12)
            c1 = x12*ry - rx*y12
            dln1 = dsqrt(a1**2 + b1**2 + c1**2)
            a1 = a1/dln1
            b1 = b1/dln1
            c1 = c1/dln1
!     --------normalizing c2c1 vector -> n2 ----
            a2 = x12/cc
            b2 = y12/cc
            c2 = z12/cc
!     -------cross product n1 x n2 calc.-> n3 ---
            a3 = b1*c2 - b2*c1
            b3 = -(a1*c2 - a2*c1)
            c3 = a1*b2 - a2*b1
!     ------- point A ------ 
            xa = rxu(ichain,1)+a2*ca
            ya = ryu(ichain,1)+b2*ca
            za = rzu(ichain,1)+c2*ca
!     ------H-atoms for the first methyl group------
            r = 0.0d0
            do i = 1,3
               a4 = ah*a3*dcos(r) + ah*a1*dsin(r)
               b4 = ah*b3*dcos(r) + ah*b1*dsin(r)
               c4 = ah*c3*dcos(r) + ah*c1*dsin(r)
               rxu(ichain,nngrow+i) = xa + a4 
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0d0*onepi/180.0d0
            end do
!     ------ H-atoms for the second methyl group ------
!     ------ point A' is opposite to point A ------ 
            xa = rxu(ichain,2)-a2*ca
            ya = ryu(ichain,2)-b2*ca
            za = rzu(ichain,2)-c2*ca
            if ( lcrysl ) then
               r=0.0d0
               ven = 0.0d0
            else
 100           r = 2.0d0*onepi*random()
               ven = en0*(1.0d0-dcos(3.0d0*r))
               prob = dexp(-beta*ven)
               rn = random()
               if (rn.gt.prob) goto 100
            end if
!     rgr = r*180.0d0/onepi
!      write(iou,*) 'final Ryckaert angle: ',rgr
            vmethyl = vmethyl + ven

            do i = 4,6
               a4 = ah*a3*dcos(r+onepi) + ah*a1*dsin(r+onepi)
               b4 = ah*b3*dcos(r+onepi) + ah*b1*dsin(r+onepi)
               c4 = ah*c3*dcos(r+onepi) + ah*c1*dsin(r+onepi)        
               rxu(ichain,nngrow+i) = xa + a4 
               ryu(ichain,nngrow+i) = ya + b4
               rzu(ichain,nngrow+i) = za + c4
               r = r + 120.0d0*onepi/180.0d0
            end do
     
         end if
      end if
      return
      end








