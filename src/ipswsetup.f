      subroutine ipswsetup
c
c setup maginns interphase switch
c
      implicit none
      include 'control.inc'
      include 'coord.inc'
      include 'ipswpar.inc'
      include 'inputdata.inc'
      include 'system.inc'
      include 'cell.inc'

      integer i,j,k,tnw,ibox

      double precision lx,hmata(9),hmatc(9),hm(9)

      logical lhm,llwell

      ibox = 1
      if (lmipsw.and.(nbox.gt.1)) stop 'ipsw only for 1 box'
      if (lmipsw.and.lnpt.and.pmvol.gt.1.0d-7) stop 'ipsw only for NVT'
      read(35,*)
      read(35,*) (lwell(i),i=1,nmolty)
      read(35,*)
      tnw = 0
      llwell = .false.
      do i = 1, nmolty
         if (lwell(i)) then
            llwell = .true.
            nwell(i) = temtyp(i)
            tnw = tnw+nwell(i)*nunit(i)
            do j = 1, nunit(i)
               read(35,*) (awell(j,k,i), k = 1,nunit(i))
            enddo
         endif
      enddo 
      if (nw.lt.tnw) then
         write(iou,*) 'increase nw in ipswpar to ', tnw
         stop
      endif
      read(35,*)
      read(35,*) bwell
      read(35,*)
      read(35,*) lstagea, lstageb, lstagec
      if (lmipsw) then
         if ((.not.lstagea).and.(.not.lstageb).and.(.not.lstagec))
     &      stop 'one stage must be true'
      endif
      if (llwell.and.lstagea) stop 'ipsw well NOT for stage a'
      if (lstagea) then
         if (lstageb.or.lstagec) stop 'only one lstage must be true'
      endif
      if (lstageb) then
         if (lstagea.or.lstagec) stop 'only one lstage must be true'
      endif
      if (lstagec) then
         if (lstagea.or.lstageb) stop 'only one lstage must be true'
      endif
      read(35,*)
      read(35,*) etais, lambdais
      if ((lambdais.le.-1.0d-6).or.(lambdais.ge.1.000001))
     &   stop 'lambdais must be between 0 and 1'
      if ((etais.le.-1.0d-6).or.(etais.ge.1.000001))
     &   stop 'etais must be between 0 and 1'
      read(35,*)
      if (lsolid(ibox).and.(.not.lrect(ibox))) then
         read(35,*) hmata(1),hmata(2),hmata(3)
         read(35,*) hmata(4),hmata(5),hmata(6)
         read(35,*) hmata(7),hmata(8),hmata(9)
         read(35,*) hmatc(1),hmatc(2),hmatc(3)
         read(35,*) hmatc(4),hmatc(5),hmatc(6)
         read(35,*) hmatc(7),hmatc(8),hmatc(9)
      elseif (.not.lsolid(ibox).and.(.not.lrect(ibox))) then
         read(35,*) lena, lenc
      endif
      if (lmipsw.and.lstageb) then
         if (lsolid(ibox).and.(.not.lrect(ibox))) then
            lhm = .false.
            do i = 1, 9
               hm(i) = (1.0d0-lambdais)*hmata(i)+lambdais*hmatc(i)
               if (dabs(hmat(ibox,i)-hm(i)).gt.1.0d-6) lhm = .true.
            enddo
            if (lhm) then
               write(iou,*) 'enter correct hmat'
               write(iou,*) hm(1),hm(2),hm(3)
               write(iou,*) hm(4),hm(5),hm(6)
               write(iou,*) hm(7),hm(8),hm(9)
               stop
            endif
         elseif (.not.lsolid(ibox)) then
            lx = (1.0d0-lambdais)*lena+lambdais*lenc
            if (dabs(boxlx(1)-lx).gt.1.0d-6) then
               write(iou,*) 'input correct boxl', lx
               stop
            endif
         endif
      endif
      if (lmipsw.and.lstagea) then
         if (lsolid(ibox).and.(.not.lrect(ibox))) then
            lhm = .false.
            do i = 1, 9
               if (dabs(hmat(ibox,i)-hmata(i)).gt.1.0d-6) lhm = .true.
            enddo
            if (lhm) then
               write(iou,*) 'enter correct hmat'
               write(iou,*) hmata(1),hmata(2),hmata(3)
               write(iou,*) hmata(4),hmata(5),hmata(6)
               write(iou,*) hmata(7),hmata(8),hmata(9)
               stop
            endif
c	write(iou,*) boxlx(1),lena
         elseif (.not.lsolid(ibox)) then
            if (dabs(boxlx(1)-lena).gt.1.0d-6) then
               write(iou,*) 'input correct boxl', lena
               stop
            endif
         endif
      endif
      if (lmipsw.and.lstagec) then
         if (lsolid(ibox).and.(.not.lrect(ibox))) then
            lhm = .false.
            do i = 1, 9
               if (dabs(hmat(ibox,i)-hmatc(i)).gt.1.0d-6) lhm = .true.
            enddo
            if (lhm) then
               write(iou,*) 'enter correct hmat'
               write(iou,*) hmatc(1),hmatc(2),hmatc(3)
               write(iou,*) hmatc(4),hmatc(5),hmatc(6)
               write(iou,*) hmatc(7),hmatc(8),hmatc(9)
               stop
            endif
         elseif (.not.lsolid(ibox)) then
            if (dabs(boxlx(1)-lenc).gt.1.0d-6) then
               write(iou,*) 'input correct boxl', lenc
               stop
            endif
         endif
      endif
      dhmat(1,1) = hmatc(1)-hmata(1)
      dhmat(1,2) = hmatc(4)-hmata(4)
      dhmat(1,3) = hmatc(7)-hmata(7)
      dhmat(2,1) = hmatc(2)-hmata(2)
      dhmat(2,2) = hmatc(5)-hmata(5)
      dhmat(2,3) = hmatc(8)-hmata(8)
      dhmat(3,1) = hmatc(3)-hmata(3)
      dhmat(3,2) = hmatc(6)-hmata(6)
      dhmat(3,3) = hmatc(9)-hmata(9)
      read(35,*)
      read(35,*) iratipsw
      if (lmipsw.and.lstageb) then
         if (mod(iratipsw,iratp).ne.0) stop 'iratipsw must be integer
     &       multiple of iratp if lstageb is true'
      endif
      do i = 1, nmolty
         if (lwell(i)) then
            do j = 1, nwell(i)*nunit(i)
               read(35,*) sxwell(j,i),sywell(j,i),szwell(j,i)
            enddo
         endif
      enddo
      if (lmipsw) then
         write(iou,*) '*****************************************'
         write(iou,*) 'Some parameters used for interphase switch'
         write(iou,*) 'moltyp, lwell, and nwell'
         do i = 1, nmolty
            if (lwell(i)) then
               write(iou,*) i, lwell(i), nwell(i)
            else
               write(iou,*) i, lwell(i), 'not defined'
            endif 
         enddo
         write(iou,*) 'lstagea,lstageb,lstagec', lstagea,lstageb,lstagec
         write(iou,*) 'etais, lambdais', etais,lambdais
         write(iou,*) '*****************************************'
      else
         lstagea = .false.
         lstageb = .false.
         lstagec = .false.
         do i = 1, nmolty
            lwell(i) = .false.
         enddo
         lambdais = 0.0d0
      endif
      acdvdl = 0.0d0
      acipsw = 0.0d0

      if (lsolid(ibox).and.(.not.lrect(ibox))) then
         do i = 1, nmolty
            if (lwell(i)) then
               do j = 1, nwell(i)*nunit(i)
                  rxwell(j,i) = sxwell(j,i)*hmat(ibox,1)+
     &               sywell(j,i)*hmat(ibox,4)+szwell(j,i)*hmat(ibox,7)
                  rywell(j,i) = sxwell(j,i)*hmat(ibox,2)+
     &               sywell(j,i)*hmat(ibox,5)+szwell(j,i)*hmat(ibox,8)
                  rzwell(j,i) = sxwell(j,i)*hmat(ibox,3)+
     &               sywell(j,i)*hmat(ibox,6)+szwell(j,i)*hmat(ibox,9)
               enddo
            endif
         enddo
      else
         do i = 1, nmolty
            if (lwell(i)) then
               do j = 1, nwell(i)*nunit(i)
                  rxwell(j,i) = sxwell(j,i)*boxlx(ibox)
                  rywell(j,i) = sywell(j,i)*boxly(ibox)
                  rzwell(j,i) = szwell(j,i)*boxlz(ibox)
               enddo
            endif
         enddo
      endif

      return
      end
