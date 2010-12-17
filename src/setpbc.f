      subroutine setpbc (ibox)

      implicit none
      include 'control.inc'
      include 'system.inc'
      include 'peboco.inc'

      integer::ibox

! ----------------------------------------------------------------

      if ( lpbcx ) then
         bx = boxlx(ibox)
         if ( lfold ) then
            hbx = 0.5d0 * bx
         else
            bxi = 1.0d0 / bx
         end if
      end if

      if ( lpbcy ) then
         by = boxly(ibox)
         if ( lfold ) then
            hby = 0.5d0 * by
         else
            byi = 1.0d0 / by
         end if
      end if

      if ( lpbcz ) then
         bz = boxlz(ibox)
         if ( lfold ) then
            hbz = 0.5d0 * bz
         else
            bzi = 1.0d0 / bz
         end if
      end if

      return
      end


