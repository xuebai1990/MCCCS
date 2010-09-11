      subroutine mimage ( rxuij,ryuij,rzuij,ibox)

      implicit none
      include 'peboco.inc'
      include 'control.inc'
      include 'system.inc'
      include 'cell.inc'

      integer::ibox
      real(8)::rxuij,ryuij,rzuij,hsx,hsy,hsz,sx,sy,sz

! ----------------------------------------------------------------

      if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
! *** non-hexagonal box
         hsx = 0.5d0*hmat(ibox,1)
         hsy = 0.5d0*hmat(ibox,5)
         hsz = 0.5d0*hmat(ibox,9)
         sx = rxuij*hmati(ibox,1)+ryuij*hmati(ibox,4)
     &        +rzuij*hmati(ibox,7)
         sy = rxuij*hmati(ibox,2)+ryuij*hmati(ibox,5)
     &        +rzuij*hmati(ibox,8)
         sz = rxuij*hmati(ibox,3)+ryuij*hmati(ibox,6)
     &        +rzuij*hmati(ibox,9)



!         if ( sx .gt. 0.5d0 ) then
!            sx = sx-1d0
!         elseif ( sx .lt. -0.5d0 ) then
!            sx = sx+1d0
!         end if
!         if ( sy .gt. 0.5d0 ) then
!            sy = sy-1d0
!         elseif ( sy .lt. -0.5d0 ) then
!            sy = sy+1d0
!         end if
!         if ( sz .gt. 0.5d0 ) then
!            sz = sz-1d0
!         elseif ( sz .lt. -0.5d0 ) then
!            sz = sz+1d0
!         end if
!         sx = sx-nint(sx)
!         sy = sy-nint(sy)
!         sz = sz-nint(sz)
!         rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
!         ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)+sz*hmat(ibox,8)
!         rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)+sz*hmat(ibox,9)



!         print*,rxuij,ryuij,rzuij

         if ( rzuij .gt. hsz ) then
            rzuij=rzuij-hmat(ibox,9)
            sz=sz-1d0
            if ( rzuij .gt. hsz ) then
               rzuij=rzuij-hmat(ibox,9)
               sz=sz-1d0
            end if
         elseif ( rzuij .lt. -hsz ) then
            rzuij=rzuij+hmat(ibox,9)
            sz=sz+1d0
            if ( rzuij .lt. -hsz ) then
              rzuij=rzuij+hmat(ibox,9)
               sz=sz+1d0
            end if
         end if

         ryuij=sy*hmat(ibox,5)+sz*hmat(ibox,8)
         if ( ryuij .gt. hsy ) then
            ryuij=ryuij-hmat(ibox,5)
            sy=sy-1d0
            if ( ryuij .gt. hsy ) then
               ryuij=ryuij-hmat(ibox,5)
               sy=sy-1d0
            end if
         elseif ( ryuij .lt. -hsy ) then
            ryuij=ryuij+hmat(ibox,5)
            sy=sy+1d0
            if ( ryuij .lt. -hsy ) then
               ryuij=ryuij+hmat(ibox,5)
              sy=sy+1d0
            end if
        end if

         rxuij=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
         if ( rxuij .gt. hsx ) then
            rxuij=rxuij-hmat(ibox,1)
            if ( rxuij .gt. hsx ) then
               rxuij=rxuij-hmat(ibox,1)
            end if
         elseif ( rxuij .lt. -hsx ) then
            rxuij=rxuij+hmat(ibox,1)
            if ( rxuij .lt. -hsx ) then
               rxuij=rxuij+hmat(ibox,1)
            end if
         end if


!         print*, sqrt(rxuij**2+ryuij**2+rzuij**2)
!         print*

      else

         if ( lpbcx ) then
            if ( lfold ) then
               if ( rxuij .gt. hbx ) then
                  rxuij=rxuij-bx
               else
                  if (rxuij.lt.-hbx) rxuij=rxuij+bx
               end if
            else
!            rxuij = rxuij - bx*anint(rxuij*bxi)
               rxuij = rxuij - bx*dint(rxuij*bxi+dsign(0.5d0,rxuij))
            end if
         end if

         if ( lpbcy ) then
            if ( lfold ) then
               if ( ryuij .gt. hby ) then
                  ryuij=ryuij-by
               else
                  if (ryuij.lt.-hby) ryuij=ryuij+by
               end if
            else
!           ryuij  = ryuij - by*anint(ryuij*byi)
               ryuij = ryuij - by*dint(ryuij*byi+dsign(0.5d0,ryuij))
            end if
         end if

         if ( lpbcz ) then
            if ( lfold ) then
               if (rzuij.gt.hbz) then
                  rzuij=rzuij-bz
               else
                  if (rzuij.lt.-hbz) rzuij=rzuij+bz
               end if
            else
!            rzuij  = rzuij - bz*anint(rzuij*bzi)
               rzuij = rzuij - bz*dint(rzuij*bzi+dsign(0.5d0,rzuij))
            end if
         end if

! --- JLR 11-17-09
! --- for orginal RPLC setup, mimage must be applied twice or energy errors can result 
! --- this may be removed in future if larger system size is used!!!
         if (ltwice(ibox)) then !everthing is folded again
            if ( rxuij .gt. hbx ) then
               rxuij=rxuij-bx
            else
               if (rxuij.lt.-hbx) rxuij=rxuij+bx
            end if
            if ( ryuij .gt. hby ) then
               ryuij=ryuij-by
            else
               if (ryuij.lt.-hby) ryuij=ryuij+by
            end if
         end if
! --- END JLR 11-17-09

      end if


      return
      end

