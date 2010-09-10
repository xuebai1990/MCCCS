      subroutine deriv(ibox)
c
c calculate integrand in maginns interphase switch
c
      implicit none
      include 'control.inc'
      include 'system.inc'
      include 'coord.inc'
      include 'ipswpar.inc'
      include 'cell.inc'

      integer::ibox,i,j,k

      real(8)::vol,hmats(3,3),hmatsi(3,3)

      if (lstagea) then
         dvdl = -(1.0d0-etais)*vipsw
      elseif (lstageb) then
         if (.not.lsolid(ibox)) then
            vol = boxlx(ibox)*boxly(ibox)*boxlz(ibox)
            dvdl = -3.0d0*vol/boxlx(ibox)*(lenc-lena)*(etais*pipsw
     &              +lambdais*pwellipsw)+vwellipsw
         elseif (lsolid(ibox).and.(.not.lrect(ibox))) then
            vol = cell_vol(ibox) 
            hmats(1,1) = hmat(ibox,1)
            hmats(1,2) = hmat(ibox,4)
            hmats(1,3) = hmat(ibox,7)
            hmats(2,1) = hmat(ibox,2)
            hmats(2,2) = hmat(ibox,5)
            hmats(2,3) = hmat(ibox,8)
            hmats(3,1) = hmat(ibox,3)
            hmats(3,2) = hmat(ibox,6)
            hmats(3,3) = hmat(ibox,9)
            hmatsi(1,1) = hmati(ibox,1)
            hmatsi(1,2) = hmati(ibox,4)
            hmatsi(1,3) = hmati(ibox,7)
            hmatsi(2,1) = hmati(ibox,2)
            hmatsi(2,2) = hmati(ibox,5)
            hmatsi(2,3) = hmati(ibox,8)
            hmatsi(3,1) = hmati(ibox,3)
            hmatsi(3,2) = hmati(ibox,6)
            hmatsi(3,3) = hmati(ibox,9)
            dvdl = 0.0d0
            do i = 1, 3
               do j = 1, 3
                  do k = 1, 3
                     dvdl = dvdl - hmatsi(j,k)*dhmat(i,j)*
     &                    (etais*pips(i,k)+lambdais*pwellips(i,k))
c	write(6,*) i,j,pips(i,j),pwellips(i,j)
                  enddo
c	write(6,*) i,j,pips(i,j)
               enddo
            enddo
            dvdl = vol*dvdl+vwellipsw
         endif
      else
         dvdl = (1.0d0-etais)*vipsw-vwellipsw
      endif

c	write(6,*) 'deriv', dvdl,vwellipsw,pipsw,pwellipsw,vol

      return
      end
