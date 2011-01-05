      subroutine config

!    *******************************************************************
!    ** performs a lengthconserving configurational bias move         **
!    ** for linear, branched, anisotropic, and explicit atom          **
!    ** molecules                                                     **
!    ** rewritten from old config and branch subroutines by           **
!    ** M.G. Martin 9-19-97                                           **
!    ** number of trial attempts starting at unit inb is stored in    **
!    **    bncb ( inb ).                                              **
!    ** number of successful generations of trial configuration is in **
!    **    bscb ( 1,inb ).                                            **
!    ** number of accepted trial configurations is in                 **
!    **    bscb ( 2,inb ).                                            **
!    *******************************************************************
 
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
!$$$      include 'coord2.inc'
!$$$      include 'system.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'rosen.inc' 
!$$$      include 'inputdata.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'poten.inc'
!$$$      include 'neigh.inc'
!$$$      include 'ipswpar.inc'
!$$$      include 'eepar.inc'
      
      logical::lterm, ovrlap, ltors, lneighij,lfixnow

      integer(KIND=normal_int)::i,j,k,iii,ibox,iunit,igrow,icbu,islen
     & ,imolty,iutry
      integer(KIND=normal_int)::istt,iett,nchp1,ic,total,bin
     & ,count,findex,iw
      integer(KIND=normal_int)::ddum,ip

      dimension ddum(27)

      real(KIND=double_precision)::v,vintra,vinter,vext,velect,vewald
     & ,vtorold,vtornew,delen,deleo,vdum,wplace,wrig
      real(KIND=double_precision)::dchain,random,rchain,wnlog,wolog
     & ,wdlog,wratio
      real(KIND=double_precision)::vrecipn,vrecipo,cwtorfo,cwtorfn,x,y,z

! ------------------------------------------------------------------

!      write(iou,*) 'start CONFIG'
! ***    select a chain at random ***
      rchain  = random()
      do icbu = 1,nmolty
         if ( rchain .lt. pmcbmt(icbu) ) then
            imolty = icbu
            exit
         end if
      end do

      if ((lexpee).and.(imolty.ge.nmolty1))
     &   imolty = ee_moltyp(mstate)
                                                                                
      if (temtyp(imolty).eq.0) return

!     *** determine whether to use fecbmc or not ***
      if (random().lt.pmfix(imolty)) then
         lfixnow = .true.
      else
         lfixnow = .false.
      end if

      if (lgrand) then
! ---    select a chain in box 1
!         write(iou,*) 'counters not implemented properly for grand'
         if (ncmt(1,imolty).eq.0) return
         i = idint( dble(ncmt(1,imolty))*random() ) + 1
         i = parbox(i,1,imolty)
         if ( moltyp(i) .ne. imolty ) write(iou,*) 'screwup config'
         ibox=1
      else 
         dchain = dble(temtyp(imolty))
         i = int( dchain*random() + 1 )
         i = parall(imolty,i)
         ibox = nboxi(i)
         if ( moltyp(i) .ne. imolty ) call cleanup('screwup config')
      end if

! *** store number of units in iunit and # to be grown in igrow ***
      iunit = nunit(imolty)
      igrow = nugrow(imolty)
 
! *** store position of trial chain in r x/y/z cbu ***
      do icbu = 1, igrow
         rxnew(icbu) = rxu(i,icbu)
         rynew(icbu) = ryu(i,icbu)
         rznew(icbu) = rzu(i,icbu)
      end do

      if (lfixnow) then
         call safeschedule(igrow,imolty,islen,iutry,findex,1)
      else
         call schedule(igrow,imolty,islen,iutry,0,1)
      end if

!     --- determine how many beads are being regrown
      total = 0
      do icbu = 1,igrow
         if ( .not. lexshed(icbu) ) total = total + 1
      end do

      if (lfixnow) then
         fbncb(imolty,findex-1) = fbncb(imolty,findex-1) + 1.0d0
      else
         bncb(imolty,total) = bncb(imolty,total) + 1.0d0
      end if
 
!      if ( lelect(imolty) ) then
!        ---  Call qqcheck to setup the group based qq cutoff
!         call qqcheck(i,ibox,rxnew(1),rynew(1),rznew(1))
!      end if

!     --- grow new chain conformation
      call rosenbluth ( .true., lterm,i,i,imolty,islen,ibox,igrow
     &     ,vdum,lfixnow,cwtorfn,1 )
 
! --- termination of cbmc attempt due to walk termination ---
      if ( lterm ) then 
!        write(iou,*) 'termination of growth',i
        return
      end if

      if (llrig) then
         call rigfix(.true.,i,ibox,imolty,lterm,wrig)
         if ( lterm ) return
         weight = weight * wrig
      end if

      if (llplace(imolty).and.lfixnow) then
         call place(.true.,lterm,i,imolty,ibox,islen,wplace)
         if ( lterm ) return
         weight = weight * wplace 
      end if


!     --- grow old chain conformation
      call rosenbluth ( .false.,lterm,i,i,imolty,islen,ibox,igrow
     &     ,vdum,lfixnow,cwtorfo,1)

!     --- termination of old walk due to problems generating orientations
      if ( lterm ) then
         write(iou,*) 'CONFIG:old growth rejected in box',ibox
     &    ,' for moltyp',imolty
         return
      end if

      if (llrig) then
         call rigfix(.false.,i,ibox,imolty,lterm,wrig)
         if ( lterm ) then
            write(iou,*) 'CONFIG: old rigid fix rejected'
            return
         end if
         weiold = weiold * wrig
      end if

      if (llplace(imolty).and.lfixnow) then
         call place(.false.,lterm,i,imolty,ibox,islen,wplace)
         
         if ( lterm ) then
            write(iou,*) 'CONFIG: old hydrogen placement rejected'
            return
         end if
         weiold = weiold * wplace 
      end if

! -----------------------------------------------------------------------------
!     Begin DC-CBMC, Explicit Atom and Ewald-sum Corrections 

      if ( ldual .or. lewald .or. iunit .ne. igrow 
     &     .or. ((.not. lchgall) .and. lelect(imolty)) ) then
!     --- Put on hydrogens for explicit AA model for calculation of COM
!     --- and assign all of the grown new and old beads to rxuion
!     --- with old = 1, new = 2
         do j=1,igrow
            rxuion(j,1)=rxu(i,j)
            ryuion(j,1)=ryu(i,j)
            rzuion(j,1)=rzu(i,j)
            qquion(j,1)=qqu(i,j)
         end do	 
         do j = 1,igrow
            rxuion(j,2) = rxnew(j)
            ryuion(j,2) = rynew(j)
            rzuion(j,2) = rznew(j)
            qquion(j,2) = qquion(j,1)
         end do
         nchp1=nchain+1
         nboxi(nchp1) = ibox
         moltyp(nchp1) = imolty
         moltion(1) = imolty
         moltion(2) = imolty

         if ( igrow .ne. iunit ) then
!           -- iii = 1 old conformation
            do j = igrow+1, iunit
               rxuion(j,1) = rxu(i,j)
               ryuion(j,1) = ryu(i,j)
               rzuion(j,1) = rzu(i,j)
               qquion(j,1) = qqu(i,j)
            end do
!           -- iii = 2 new conformation
            do j=1, igrow
               rxu(nchp1,j) = rxnew(j)
               ryu(nchp1,j) = rynew(j)
               rzu(nchp1,j) = rznew(j)
            end do
            call explct(nchp1,vtornew,.false.,.false.)
            do j=igrow+1, iunit
               rxuion(j,2) = rxu(nchp1,j)
               ryuion(j,2) = ryu(nchp1,j)
               rzuion(j,2) = rzu(nchp1,j)
               qquion(j,2) = qquion(j,1)
            end do
         end if
      end if

      if (ldual .or. ((.not. lchgall) .and. lelect(imolty))
     &     .or. (lchgall .and. lewald .and. (.not. ldual))) then
         istt = 1
         iett = igrow

!        -- check new before old
         do iii = 2,1,-1
!          calculate the Full rcut Lennard-Jones energy for the grown beads
!          iii = 1 old conformation
!          iii = 2 new conformation
            
            call energy (i,imolty,v,vintra,vinter,vext,velect
     &           ,vewald,iii,ibox,istt,iett,.true.,ovrlap,
     &           .false.,vdum,.false.,.false.)

            if (ovrlap .and. (iii .eq. 1)) then
!            if (ovrlap) then
               write(iou,*) 'disaster: overlap in old conf config',i
               call cleanup('')
            end if

            if (iii .eq. 2) then
               delen = ( vnewinter + vnewext + vnewelect +
     &              vnewewald + vnewintra) 
               if (lstagea) then
                  delen = (1.0d0-(1.0d0-etais)*lambdais)*delen
               elseif (lstageb) then
                  delen = etais*delen
               elseif (lstagec) then
                  delen = (etais+(1.0d0-etais)*lambdais)*delen
               end if
               delen = v - delen
! --- JLR 11-19-09 Commenting this out, it makes no sense and gives me an energy error!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM   this may not exactly right, but for
!!!!! numerical reasons I think we need it 
!               IF(delen*beta .LT. -2.3d0*softcut) THEN
!                  delen=-2.3d0*softcut/beta
!               ELSEIF(delen*beta .GT. 2.3d0*softcut) THEN
!                  delen=2.3d0*softcut
!               end if
! --- END JLR 11-19-09
               weight    = weight*dexp(-(beta*delen))
               vnewt     = vnewt + delen
               vnewinter = vinter
               vnewext   = vext
               vnewelect = velect
               vnewintra = vintra
               vnewewald = vewald
            else
               deleo = ( voldinter + voldext + voldelect +
     &              voldewald + voldintra) 
               if (lstagea) then
                  deleo = (1.0d0-(1.0d0-etais)*lambdais)*deleo
               elseif (lstageb) then
                  deleo = etais*deleo
               elseif (lstagec) then
                  deleo = (etais+(1.0d0-etais)*lambdais)*deleo
               end if
               deleo = v - deleo
! --- JLR 11-19-09 Commenting this out, it gives me energy error!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM   this may not exactly right, but for
!!!!! numerical reasons I think we need it 
!               IF(deleo*beta .LT. -2.3d0*softcut) THEN
!                  deleo=-2.3d0*softcut/beta
!               ELSEIF(deleo*beta .GT. 2.3d0*softcut) THEN
!                  deleo=2.3d0*softcut
!               end if
! --- END JLR 11-19-09
               weiold    = weiold*dexp(-(beta*deleo))
               voldt     = voldt + deleo
               voldinter = vinter
               voldext   = vext
               voldelect = velect
               voldintra = vintra
               voldewald = vewald
            end if
         end do
         
      end if

      if ( iunit .ne. igrow ) then
         istt = igrow+1
         iett = iunit

!        -- check new before old
         do iii = 2,1,-1
!           calculate the true Lennard-Jones energy for the hydrogens
!           iii=1 new conformation
!           iii=2 old conformation
!           hydrogens were placed and rxuion was assigned above

            if (iii .eq. 1) then
               ltors = .true.
            else
               ltors = .false.
            end if

! Calculate the energy of the non-backbone beads 
            call energy (i,imolty,v,vintra,vinter,vext,velect
     &           ,vewald,iii,ibox,istt,iett, .true.,ovrlap
     &           ,ltors,vtorold,.true.,.false.)

            if (iii .eq. 2) then
               if (ovrlap) return
               delen = v + vtornew
               if ( delen*beta .gt. (2.3d0*softcut) ) then
!                  write(iou,*) '##softcut in config caught explicit atoms'
                  return
               end if
               weight = weight*dexp(-(beta*delen))
               vnewt  = vnewt + delen
               vnewintra = vnewintra + vintra 
               vnewinter = vnewinter + vinter 
               vnewext   = vnewext + vext 
               vnewtg = vnewtg + vtornew
               vnewelect = vnewelect + velect
               vnewewald = vnewewald + vewald
            else
               if (ovrlap) then
                  write(iou,*) 'ovrlap problem in old confomation',
     &                 ' - CONFIG'
                  return
               end if
               deleo = v + vtorold
               weiold = weiold*dexp(-(beta*deleo))
               if ( weiold .lt. softlog ) then
                  write(iou,*) '##old weight for explicit too low'
               end if
               voldt     = voldt + deleo
               voldintra = voldintra + vintra
               voldinter = voldinter + vinter
               voldext   = voldext + vext
               voldtg    = voldtg + vtorold
               voldelect = voldelect + velect
               voldewald = voldewald + vewald
            end if
         end do
      end if

      if ( lewald .and. lelect(imolty) ) then
!        --- reciprocal space sum ---
!        --- rxuion: 1= old configuration; 2= new configuration
         call recip(ibox,vrecipn,vrecipo,1)
         delen = vrecipn
         deleo = vrecipo
         vnewelect = vnewelect + vrecipn
         voldelect = voldelect + vrecipo
         vipswn = vipswn + vrecipn
         vipswo = vipswo + vrecipo
         if (lstagea) then
            vrecipn = (1.0d0-(1.0d0-etais)*lambdais)*vrecipn
            vrecipo = (1.0d0-(1.0d0-etais)*lambdais)*vrecipo
         elseif (lstageb) then
            vrecipn = etais*vrecipn
            vrecipo = etais*vrecipo
         elseif (lstagec) then
            vrecipn = (etais+(1.0d0-etais)*lambdais)*vrecipn
            vrecipo = (etais+(1.0d0-etais)*lambdais)*vrecipo
         end if
         weight = weight * dexp(-(beta*vrecipn))
         weiold = weiold * dexp(-(beta*vrecipo))
         vnewt = vnewt + vrecipn
         voldt = voldt + vrecipo
      end if

!     End of DC-CBMC, Explicit Atom and Ewald-sum Corrections

! *** check for acceptance of trial configuration ***
      wnlog = dlog10 ( weight )
      wolog = dlog10 ( weiold )
!      write(iou,*) 'weight:',weight
!      write(iou,*) 'weiold:',weiold
      wdlog = wnlog - wolog
      if ( wdlog .lt. -softcut ) then
!         write(iou,*) 'cbmc softcut',i
         return
      end if
 
      if (lfixnow) then
         wratio = weight * cwtorfo / ( weiold * cwtorfn)
         fbscb(imolty,1,findex-1) = 
     &        fbscb(imolty,1,findex-1) + 1.0d0
      else
         wratio = weight / weiold
         bscb(imolty,1,total) = bscb(imolty,1,total) + 1.0d0
      end if

      if ( random() .le. wratio ) then
!         write(iou,*) 'CONFIG accepted',i,ibox
!        --- we can now accept !!!!! ***
         if (lfixnow) then
            fbscb(imolty,2,findex-1) = fbscb(imolty,2,findex-1) 
     &           + 1.0d0
         else
            bscb(imolty,2,total) = bscb(imolty,2,total) + 1.0d0
         end if


         vbox(ibox)    = vbox(ibox)    + ( vnewt - voldt )
         vinterb(ibox) = vinterb(ibox) + (vnewinter - voldinter)
         vintrab(ibox) = vintrab(ibox) + (vnewintra- voldintra)
         vvibb(ibox)   =  vvibb(ibox)  + (vnewbvib- voldbvib)
         vtgb(ibox)    = vtgb(ibox)    + (vnewtg- voldtg)
         vextb(ibox)   = vextb(ibox)   + (vnewext - voldext)
         vbendb(ibox)  = vbendb(ibox)  + (vnewbb - voldbb)
         velectb(ibox) = velectb(ibox) + (vnewelect - voldelect)
     &        + (vnewewald - voldewald)
         vipswb(ibox) = vipswb(ibox) + (vipswn-vipswo)
         vwellipswb(ibox) = vwellipswb(ibox) + (vwellipswn-vwellipswo)
         vipsw = vipswb(ibox)
         vwellipsw = vwellipswb(ibox)
         do ic = 1, igrow
            rxu(i,ic) = rxnew(ic)
            ryu(i,ic) = rynew(ic)
            rzu(i,ic) = rznew(ic)
         end do
         do ic = igrow+1, iunit
            rxu(i,ic)  = rxuion(ic,2)
            ryu(i,ic)  = ryuion(ic,2)
            rzu(i,ic)  = rzuion(ic,2)
         end do

         if ( lewald .and. lelect(imolty) ) then
! *** update reciprocal-space sum
            call recip(ibox,vdum,vdum,2)
         end if

         if (ldielect) then
            call dipole(ibox,1)
         end if  

! ***    update center of mass
         call ctrmas(.false.,ibox,i,7)
! *** update linkcell, if applicable
         if ( licell .and. (ibox.eq.boxlink)) then
            call linkcell(2,i,vdum,vdum,vdum,ddum)
         end if

! ***    update the neighbour map ***
         if ( lneigh ) call updnn( i )

         if ( lneighbor ) then
            do ic = 1, neigh_cnt(i)
               j = neighbor(ic,i)
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. i ) then
                     neighbor(ip,j)=neighbor(neigh_cnt(j),j)
                     neigh_cnt(j) = neigh_cnt(j)-1
                     exit
                  end if
               end do
            end do
            neigh_cnt(i) = neigh_icnt
            do ic = 1,neigh_icnt
               j = neighi(ic)
               neighbor(ic,i)=j
               lneighij = .false.
               do ip = 1,neigh_cnt(j)
                  if ( neighbor(ip,j) .eq. i ) then
                     lneighij = .true.
                  end if
               end do
               if ( .not. lneighij ) then
                  neigh_cnt(j) = neigh_cnt(j)+1
                  neighbor(neigh_cnt(j),j) = i
               end if
            end do
         end if

       end if

       if (lpresim.or.lfixnow) then
!     --- record bond distances for presimulation and reweighting
          counthist = counthist + 1
          do iw = 1, islen
             do count = 1, grownum(iw)
                k = growlist(iw,count)
                do j = 1, nunit(imolty)
                   if (k.eq.j) cycle
                   x = rxu(i,j) - rxu(i,k)
                   y = ryu(i,j) - ryu(i,k)
                   z = rzu(i,j) - rzu(i,k)
                   
                   bin = anint(10.0d0*dsqrt(x**2+y**2+z**2))
                   
                   if (bin.gt.maxbin) cycle
                   
                   hist(j,k,bin) = hist(j,k,bin) + 1
                   hist(k,j,bin) = hist(k,j,bin) + 1
                end do
             end do
          end do
       end if
         
! -----------------------------------------------------------------
!       write(iou,*) 'end CONFIG'
       return
       end

