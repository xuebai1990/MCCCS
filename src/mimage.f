      subroutine mimage ( rxuij,ryuij,rzuij,ibox)

      implicit none
      include 'peboco.inc'
      include 'control.inc'
      include 'system.inc'
      include 'cell.inc'

      integer ibox
      double precision rxuij,ryuij,rzuij,hsx,hsy,hsz,sx,sy,sz

c ----------------------------------------------------------------

      if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
c *** non-hexagonal box
         hsx = 0.5d0*hmat(ibox,1)
         hsy = 0.5d0*hmat(ibox,5)
         hsz = 0.5d0*hmat(ibox,9)
         sx = rxuij*hmati(ibox,1)+ryuij*hmati(ibox,4)
     &        +rzuij*hmati(ibox,7)
         sy = rxuij*hmati(ibox,2)+ryuij*hmati(ibox,5)
     &        +rzuij*hmati(ibox,8)
         sz = rxuij*hmati(ibox,3)+ryuij*hmati(ibox,6)
     &        +rzuij*hmati(ibox,9)



c         if ( sx .gt. 0.5d0 ) then
c            sx = sx-1d0
c         elseif ( sx .lt. -0.5d0 ) then
c            sx = sx+1d0
c         endif
c         if ( sy .gt. 0.5d0 ) then
c            sy = sy-1d0
c         elseif ( sy .lt. -0.5d0 ) then
c            sy = sy+1d0
c         endif
c         if ( sz .gt. 0.5d0 ) then
c            sz = sz-1d0
c         elseif ( sz .lt. -0.5d0 ) then
c            sz = sz+1d0
c         endif
c         sx = sx-nint(sx)
c         sy = sy-nint(sy)
c         sz = sz-nint(sz)
c         rxuij = sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
c         ryuij = sx*hmat(ibox,2)+sy*hmat(ibox,5)+sz*hmat(ibox,8)
c         rzuij = sx*hmat(ibox,3)+sy*hmat(ibox,6)+sz*hmat(ibox,9)



c         print*,rxuij,ryuij,rzuij

         if ( rzuij .gt. hsz ) then
            rzuij=rzuij-hmat(ibox,9)
            sz=sz-1d0
            if ( rzuij .gt. hsz ) then
               rzuij=rzuij-hmat(ibox,9)
               sz=sz-1d0
            endif
         elseif ( rzuij .lt. -hsz ) then
            rzuij=rzuij+hmat(ibox,9)
            sz=sz+1d0
            if ( rzuij .lt. -hsz ) then
              rzuij=rzuij+hmat(ibox,9)
               sz=sz+1d0
            endif
         endif

         ryuij=sy*hmat(ibox,5)+sz*hmat(ibox,8)
         if ( ryuij .gt. hsy ) then
            ryuij=ryuij-hmat(ibox,5)
            sy=sy-1d0
            if ( ryuij .gt. hsy ) then
               ryuij=ryuij-hmat(ibox,5)
               sy=sy-1d0
            endif
         elseif ( ryuij .lt. -hsy ) then
            ryuij=ryuij+hmat(ibox,5)
            sy=sy+1d0
            if ( ryuij .lt. -hsy ) then
               ryuij=ryuij+hmat(ibox,5)
              sy=sy+1d0
            endif
        endif

         rxuij=sx*hmat(ibox,1)+sy*hmat(ibox,4)+sz*hmat(ibox,7)
         if ( rxuij .gt. hsx ) then
            rxuij=rxuij-hmat(ibox,1)
            if ( rxuij .gt. hsx ) then
               rxuij=rxuij-hmat(ibox,1)
            endif
         elseif ( rxuij .lt. -hsx ) then
            rxuij=rxuij+hmat(ibox,1)
            if ( rxuij .lt. -hsx ) then
               rxuij=rxuij+hmat(ibox,1)
            endif
         endif


c         print*, sqrt(rxuij**2+ryuij**2+rzuij**2)
c         print*

      else

         if ( lpbcx ) then
            if ( lfold ) then
               if ( rxuij .gt. hbx ) then
                  rxuij=rxuij-bx
               else
                  if (rxuij.lt.-hbx) rxuij=rxuij+bx
               endif
            else
c            rxuij = rxuij - bx*anint(rxuij*bxi)
               rxuij = rxuij - bx*dint(rxuij*bxi+dsign(0.5d0,rxuij))
            endif
         endif

         if ( lpbcy ) then
            if ( lfold ) then
               if ( ryuij .gt. hby ) then
                  ryuij=ryuij-by
               else
                  if (ryuij.lt.-hby) ryuij=ryuij+by
               endif
            else
c           ryuij  = ryuij - by*anint(ryuij*byi)
               ryuij = ryuij - by*dint(ryuij*byi+dsign(0.5d0,ryuij))
            endif
         endif

         if ( lpbcz ) then
            if ( lfold ) then
               if (rzuij.gt.hbz) then
                  rzuij=rzuij-bz
               else
                  if (rzuij.lt.-hbz) rzuij=rzuij+bz
               endif
            else
c            rzuij  = rzuij - bz*anint(rzuij*bzi)
               rzuij = rzuij - bz*dint(rzuij*bzi+dsign(0.5d0,rzuij))
            endif
         endif

c --- JLR 11-17-09
c --- for orginal RPLC setup, mimage must be applied twice or energy errors can result 
c --- this may be removed in future if larger system size is used!!!
         if (ltwice(ibox)) then !everthing is folded again
            if ( rxuij .gt. hbx ) then
               rxuij=rxuij-bx
            else
               if (rxuij.lt.-hbx) rxuij=rxuij+bx
            endif
            if ( ryuij .gt. hby ) then
               ryuij=ryuij-by
            else
               if (ryuij.lt.-hby) ryuij=ryuij+by
            endif
         endif
c --- END JLR 11-17-09

      endif


      return
      end

