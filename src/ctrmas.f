       subroutine ctrmas ( lall, ibox, j, mtype)

! ctrmas
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
! John Stubbs, and Collin Wick and Ilja Siepmann  
!                     
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to 
!
! Free Software Foundation, Inc. 
! 59 Temple Place - Suite 330
! Boston, MA  02111-1307, USA.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     ************************************************************
!     ** finds the center of mass of a chain and returns it to  **
!     ** the periodic box if it has left the box.  lall controls**
!     ** whether this is done for just one chain or for the     **
!     ** entire box. ibox is the box of the particle and j is   **
!     ** the particle number.                                   **
!     ************************************************************

      implicit none

      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'
      include 'cell.inc'

      integer::ibox, i,ii,imolty,iunit,j,stt,edd,mtype,iwarn
     &        ,inboxx,inboxy,inboxz,iadjust,itype
      logical::lall,ldx,ldy,ldz,lintbx
      real(8)::bx,by,bz,boxlen,dx,dy,dz,nxcm,nycm,nzcm
      real(8)::dmaxsq,rxuij,ryuij,rzuij,rijsq

      real(8)::sx,sy,sz,nxcm2,nycm2,nzcm2

      bx = boxlx(ibox)
      by = boxly(ibox)
      bz = boxlz(ibox)

      if (lall) then
         stt = 1
         edd = nchain
      else
         stt = j
         edd = j
      end if

      if ((mtype .eq. 1) .or.(mtype .eq. 2).or.(mtype .eq. 5)) then 

         if (mtype .eq. 1 .and. 
     &        (lsolid(ibox) .and. .not. lrect(ibox))) then
            iwarn = 2
         elseif(mtype.eq.2) then
! kea 6/3/09 --- necessary for non-COM rotations
            iwarn = 2
         else
            iwarn = 1
         end if
      elseif ( mtype .eq. 6 ) then
         iwarn = 0
      else
         iwarn = 2
      end if

      do i = stt, edd
         imolty = moltyp(i)
         iunit = nunit(imolty)
         iadjust = 1            

! ----- Check if the chain i is in the correct box
         if (nboxi(i) .eq. ibox) then
! ----- Determine new center of mass for chain i 
            lintbx = .false.
 25         nxcm = 0.0d0
            nycm = 0.0d0
            nzcm = 0.0d0
            do ii = 1, iunit
               itype = ntype(imolty,ii)
               nxcm = nxcm + rxu(i,ii) * mass(itype)
               nycm = nycm + ryu(i,ii) * mass(itype)
               nzcm = nzcm + rzu(i,ii) * mass(itype)
            end do
            nxcm = nxcm / masst(imolty)
            nycm = nycm / masst(imolty)
            nzcm = nzcm / masst(imolty)

            ldx = .false.
            ldy = .false.
            ldz = .false.
            dx = 0.0d0
            dy = 0.0d0
            dz = 0.0d0

            if (lsolid(ibox) .and. .not. lrect(ibox)) then

               sx = nxcm*hmati(ibox,1)+nycm*hmati(ibox,4)
     &              +nzcm*hmati(ibox,7)
               sy = nxcm*hmati(ibox,2)+nycm*hmati(ibox,5)
     &              +nzcm*hmati(ibox,8)
               sz = nxcm*hmati(ibox,3)+nycm*hmati(ibox,6)
     &              +nzcm*hmati(ibox,9)

               if ( sx .lt. -1.0d-10 ) then
                  sx = sx + 1.0d0
                  ldx = .true.
               elseif( sx .gt. 1d0 ) then
                  sx = sx - 1.0d0
                  ldx = .true.
               end if
               if ( sy .lt. -1.0d-10 ) then
                  sy = sy + 1.0d0
                  ldy = .true.
               elseif( sy .gt. 1d0 ) then
                  sy = sy - 1.0d0
                  ldy = .true.
               end if
               if ( sz .lt. -1.0d-10 ) then
                  sz = sz + 1.0d0
                  ldz = .true.
               elseif( sz .gt. 1d0 ) then
                  sz = sz - 1.0d0
                  ldz = .true.
               end if

               sxcm(i) = sx
               sycm(i) = sy
               szcm(i) = sz
               if ( ldx .or. ldy .or. ldz ) then
                  if ( mtype .eq. 5 ) then
                     write(iou,*) 'sx, sy, sz:',sx,sy,sz
                  end if
                  nxcm2 = sx*hmat(ibox,1)+sy*hmat(ibox,4)
     &                 +sz*hmat(ibox,7)
                  nycm2 = sx*hmat(ibox,2)+sy*hmat(ibox,5)
     &                 +sz*hmat(ibox,8)
                  nzcm2 = sx*hmat(ibox,3)+sy*hmat(ibox,6)
     &                 +sz*hmat(ibox,9)
                  dx = nxcm2-nxcm
                  dy = nycm2-nycm
                  dz = nzcm2-nzcm
                  do ii = 1, iunit
                     rxu(i,ii) = rxu(i,ii) + dx
                     ryu(i,ii) = ryu(i,ii) + dy
                     rzu(i,ii) = rzu(i,ii) + dz
                  end do
                  nxcm = nxcm2
                  nycm = nycm2
                  nzcm = nzcm2
               end if               

            else

               if( lintbx ) then
                  if ( nxcm .gt. bx ) then
                     inboxx = idint(nxcm/bx)
                     ldx = .true.
                  elseif ( nxcm .lt. -1.0d-10 ) then
                     inboxx = idint(nxcm/bx) - 1
                     ldx = .true.
                  end if
                  if ( nycm .gt. by ) then
                     inboxy = idint(nycm/by)
                     ldy = .true.
                  elseif( nycm .lt. -1.0d-10 )then
                     inboxy = idint(nycm/by) - 1
                     ldy = .true.
                  end if
                  if (lpbcz) then
                     if ( nzcm .gt. bz ) then
                        inboxz = idint(nzcm/bz)
                        ldz = .true.
                     elseif( nzcm .lt. -1.0d-10 ) then
                        inboxz = idint(nzcm/bz) - 1
                        ldz = .true.
                     end if
                  end if
                  if ( ldx ) then
                     dx = -(dble(inboxx)*bx)
                  end if
                  if ( ldy ) then
                     dy = -(dble(inboxy)*by)
                  end if
                  if ( ldz ) then
                     dz = -(dble(inboxz)*bz)
                  end if
               else
                  if (nxcm .lt. -1.0d-10) then
                     dx = bx
                     ldx = .true.
                  elseif (nxcm .gt. bx) then
                     dx = -bx
                     ldx = .true.
                  end if
                  
                  if (nycm .lt. -1.0d-10) then
                     dy = by
                     ldy = .true.
                  elseif (nycm .gt. by) then
                     dy = -by
                     ldy = .true.
                  end if
                  
                  if(lpbcz) then
                     if (nzcm .lt. -1.0d-10) then
                        dz = bz
                        ldz = .true.
                     elseif (nzcm .gt. bz) then
                        dz = -bz
                        ldz = .true.
                     end if
                  end if
               end if
               
               if (ldx) then
                  do ii = 1, iunit
                     rxu(i,ii) = rxu(i,ii) + dx
                  end do
               end if
               
               if (ldy) then
                  do ii = 1, iunit
                     ryu(i,ii) = ryu(i,ii) + dy
                  end do
               end if
               
               if (ldz) then
                  do ii = 1, iunit
                     rzu(i,ii) = rzu(i,ii) + dz
                  end do
               end if
               
            end if

            if( ldx .or. ldy .or. ldz ) then
               if ( (iadjust .ge. iwarn) ) then
                  if (mtype .eq. 1) write(iou,*) 'translational move'
                  if (mtype .eq. 2) write(iou,*) 'rotational move'
                  if (mtype .eq. 3) write(iou,*) 'swap move'
                  if (mtype .eq. 4) write(iou,*) 'switch move'
                  if (mtype .eq. 5) write(iou,*) 'volume move'
                  if (mtype .eq. 6) write(iou,*) 'readdat move'
                  if (mtype .eq. 7) write(iou,*) 'config move'
                  if (mtype .eq. 8) write(iou,*) 'swatch move'
                  if (mtype .eq. 9) write(iou,*) 'energy call'
                  write(iou,*) 'ibox,i,iunit,boxlen',ibox,i,iunit,bx,
     &                 by,bz
                  lintbx = .true.
                  write(iou,*) 'nxcm,nycm,nzcm',nxcm,nycm,nzcm
                  write(iou,*) 'dx,dy,dz',dx,dy,dz
                  if (iwarn .ne. 0) call cleanup('')
               end if
               iadjust = iadjust + 1
               goto 25
            end if

!     assign the new center of mass
            xcm(i) = nxcm
            ycm(i) = nycm
            zcm(i) = nzcm
            if ( (lcutcm .or. ldual) .and. iwarn .ne. 1 ) then
               dmaxsq = 0.0d0 
               do ii=1,iunit
                  rxuij = rxu(i,ii)-nxcm
                  ryuij = ryu(i,ii)-nycm
                  rzuij = rzu(i,ii)-nzcm
                  
!     --- minimum image the ctrmas pair separations ---
                  if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )
                  
                  rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij
                  
                  if ( rijsq .gt. dmaxsq ) dmaxsq = rijsq
               end do                  
               rcmu(i) = dsqrt(dmaxsq)+ 1.0d-10
!               write(iou,*) 'rcmu(i)',rcmu(i)
            end if
                  
         else
            if (.not. lall) write(iou,*)'prob with box in ctrmas'
         end if
      end do

      return
      end


