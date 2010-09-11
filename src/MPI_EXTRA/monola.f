      subroutine monola  

! monola
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

! -----------------------------------------------------------------
! subroutine monola
! - reads the control-data from unit 4
! - reads a starting configuration from unit 7
! - calculates interaction table
! - starts and controls the simulation
! -----------------------------------------------------------------

      implicit none

      include 'control.inc'
      include 'coord.inc'      
      include 'poten.inc'
      include 'system.inc'
      include 'ensemble.inc'
      include 'cbmc.inc'
      include 'conver.inc'
      include 'inpar.inc' 
      include 'external.inc'
      include 'zeolite.inc'
      include 'inputdata.inc'
      include 'blkavg.inc'
      include 'bnbsma.inc'
      include 'swtcmove.inc'
      include 'ewaldsum.inc'
      include 'neigh.inc'
      include 'clusterbias.inc'
      include 'cell.inc'
      include 'ipswpar.inc'
      include 'eepar.inc'
!   KEA -- spline torsion inc file
      include 'torsion.inc'
      include 'garofalini.inc'
!   RP added for MPI
      include 'mpif.h'
      include 'mpi.inc'
! -----------------------

! - variables added for GCMC histogram reweighting
      integer::fmax,idum,nummol
      parameter (fmax=1e6)
      character*20 file_flt,file_hist,file_ndis(ntmax)
      character*20 file_config,file_cell
      character*20 fname2,fname3,fname4,ftemp
      character*50 fileout
      integer::fname,ntii,findpos
      integer::imax,itmax
      integer::n,nconfig,nentry,nminp(ntmax),nmaxp(ntmax)
      integer::ncmt_list(fmax,ntmax),ndist(0:nmax,ntmax)
      real(8)::eng_list(fmax)
      real(8)::vhist

      real(8)::Temp_Energy, Temp_Mol_Vol

      integer::point_of_start, point_to_end

      integer::im,mnbox,i,j,inb,nblock,ibox,jbox,nend,nnn,ii,itemp
     &  ,itype,itype2,intg,imolty,ilunit,nbl,itel,ig,il,ucheck
     &  ,jbox_max,k,histtot,Temp_nmol
      integer::nvirial,zz,steps,igrow,ddum,total
      real(8)::starvir,stepvir,starviro
      real(8)::acv,acvsq,aflv,acpres,acnp,acmove,acsurf
     &  ,acboxl,acboxa,asetel,acdens,acnbox,dsq,v,vinter,vtail,vend
     &  ,vintra,vvib,vbend,vtg,vext,vstart,press1,press2,dsq1
     &  ,velect,vflucq,boxlen,acnbox2,pres(nbxmax),surf,acvolume
      real(8)::rm,random,temvol,setx,sety,setz,setel
     &  ,pscb1,pscb2,ratvol,avv,temacd,temspd,dblock,dbl1
     &  ,sterr,stdev,errme,qelect
      real(8)::bsswap,bnswap,ostwald,stdost,dummy,debroglie
     &  ,bnswap_in,bnswap_out, acvkjmol(nener,nbxmax)

      double precision, dimension(nprop1,nbxmax,nbxmax)::acsolpar
 
      double precision, dimension(nbxmax)::acEnthalpy,acEnthalpy1

      real(8):::: enthalpy,enthalpy2,sigma2Hsimulation
      real(8):::: inst_enth, inst_enth2, tmp,sigma2H,Cp

      real(8):::: ennergy,ennergy2,sigma2Esimulation
      real(8):::: inst_energy, inst_energy2, sigma2E,Cv

      
  
      real(8)::binvir,binvir2,inside,bvirial
      real(8)::enchg1,enthchg1,srand
      real(8)::enchg2,enthchg2
      real(8)::enchg3,enthchg3
      real(8)::cal2joule, joule2cal
      real(8)::HSP_T, HSP_LJ, HSP_COUL 
      real(8)::CED_T, CED_LJ, CED_COUL
      real(8)::Heat_vapor_T,Heat_vapor_LJ,Heat_vapor_COUL   
      real(8)::Heat_vapor_EXT, Ext_Energy_Liq, Ext_Energy_Gas
      real(8)::DeltaU_Ext, pdV 

      double precision, dimension(nprop1,nbxmax,nbxmax)::stdev1,
     &                     sterr1,errme1       
 
      dimension binvir(maxvir,maxntemp),binvir2(maxvir,maxntemp)
      dimension vstart(nbxmax),vend(nbxmax),avv(nener,nbxmax),
     &   acv(nener,nbxmax),acvsq(nener,nbxmax),aflv(nbxmax)
      dimension acboxl(nbxmax,3),acboxa(nbxmax,3),acpres(nbxmax)
      dimension acsurf(nbxmax),acvolume(nbxmax)
      dimension acnbox(nbxmax,ntmax),acnbox2(nbxmax,ntmax,20)
      dimension bsswap(ntmax,npabmax,nbxmax*2),
     &     bnswap(ntmax,npabmax,nbxmax*2),bnswap_in(ntmax,2),
     &     bnswap_out(ntmax,2)
      real(8)::molfra,molfrac,gconst,vdum
      dimension mnbox(nbxmax,ntmax),asetel(nbxmax,ntmax),
     &     acdens(nbxmax,ntmax),molfra(nbxmax,ntmax)
      real(8)::molvol(nbxmax),speden(nbxmax)
      real(8)::flucmom(nbxmax),flucmom2(nbxmax),flucev(nbxmax)
     &     ,flucv(nbxmax),dielect,acvol(nbxmax),acvolsq(nbxmax)
      dimension qelect(nntype)
! --- dimension statements for block averages ---
      character::*15 vname(nener)
      dimension dsq(nprop,nbxmax), stdev(nprop,nbxmax),
     &          dsq1(nprop1,nbxmax,nbxmax),  
     &          sterr(nprop,nbxmax),errme(nprop,nbxmax)
      dimension ucheck(ntmax),ddum(27)

      logical::ovrlap,lratio,lratv,lprint,lmv,lrsave,lblock,lucall
     &     ,lvirial2,ltfix,lratfix,ltsolute,lsolute,lpr

      dimension lratfix(ntmax),lsolute(ntmax)
      character::*25 enth
      character::*25 enth1

      integer::bin,cnt_wf1(0:6,0:6,4),cnt_wf2(0:6,0:6,4),
     &     cnt_wra1(1000,4),cnt_wra2(1000,4)
      real(8)::binstep,profile(1000)

! KEA
      integer::ttor

! -------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
      call init_vars

      cal2joule = 4.184d0
      joule2cal = 1.0d0/cal2joule 

      do bin = 1,1000
         profile(bin) = 0.0d0
      end do
      do i = 0,6
         do j = 0,6
            do nnn = 1,4
               cnt_wf1(i,j,nnn) = 0
               cnt_wf2(i,j,nnn) = 0
            end do
         end do
      end do
      do i = 1,1000
         do j = 1,4
            cnt_wra1(i,j) = 0
            cnt_wra2(i,j) = 0
         end do
      end do
      binstep = 0.05d0


      lvirial2 = .false.

      liswatch = .false.

! RP added for MPI       
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
! --------------------------------------

! KM for MPI
! only one processor at a time reads and writes data from files
      do i=1,numprocs
         if (myid.eq.i-1) then
            call readdat(lucall,ucheck,nvirial,starvir
     &           ,stepvir,qelect)
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         close(4)
      end do

! KM for MPI
! program will hang if call cleanup('') called from readdat
! set ldie in readdat and have all processors call cleanup('') here
      if (ldie) call cleanup('')


! KM for MPI
! check that if lneigh or lgaro numprocs .eq. 1
      if (lneigh.and.numprocs.ne.1) then
         write(iou,*) 'Cannot run on more than 1 processor with 
     &        neighbor list!!'
         call cleanup('')
      end if
      if (lgaro.and.numprocs.ne.1) then
         write(iou,*) 'Cannot run on more than 1 processor with 
     &        lgaro = .true.!!'
         call cleanup('')
      end if

! kea don't call cleanup('') for lgaro
!      if (lchgall .and. (.not. lewald).and.(.not.lgaro)) then
!         write(iou,*) 'lchgall is true and lewald is false.',
!     &        ' Not checked for accuracy!'
!         call cleanup('')
!      end if

      vname(1)  = ' Total energy'
      vname(2)  = ' Inter LJ    '
      vname(3)  = ' Bond bending'
      vname(4)  = ' Torsion     '
      vname(5)  = ' Intra LJ    '
      vname(6)  = ' External pot'
      vname(7)  = ' Stretch     '
      vname(8)  = ' Coulomb     '
      vname(9)  = ' Tail  LJ    '
      vname(10) = ' Fluc Q      '
      enth      = ' Enthalpy Inst. Press '
      enth1     = ' Enthalpy Ext.  Press '

      fname = run_num
!      write(6,*) 'fname ', fname
! --  SETTING UP ARRAYS  FOR ANALYSYS PURPOSE

!--- JLR 11-11-09
!--- do not call analysis if you set ianalyze to be greater than number of cycles 
! KM 01/10 remove analysis
!       if (ianalyze.le.nstep) then
!          call analysis(0)
!       end if
!--- END JLR 11-11-09

! --  set up initial linkcell
      if (licell) then
         call linkcell(1,0,vdum,vdum,vdum,ddum)
      end if

! --- set up thermodynamic integration stuff
      if (lmipsw) call ipswsetup
 
      if (.not.lmipsw) then
          lstagea = .false.
          lstageb = .false.
          lstagec = .false.
      end if 
      
!      write(6,*) 'lexpee ', lexpee

! --- set up expanded ensemble stuff
      if (lexpee) call eesetup(qelect)

      if (lexpee.and.lmipsw) call cleanup('not for BOTH lexpee AND lmipsw')
      
! - use internal read/write to get integer::number in character::format
      write(ftemp,*) fname
      read(ftemp,*) fname2

! KM for MPI
! only processor 0 opens files     
      if (myid.eq.0) then 
         file_cell =
     &        'cell_param'//fname2(1:len_trim(fname2))//suffix//'.dat'
         open(unit=13,file=file_cell, status='unknown')
         close(unit=13)
      end if
      
! - setup files for histogram reweighting
! KM fom MPI
! will need to check this file I/O if want to run grand canonical in parallel
      if(lgrand) then
         if (myid.eq.0) then
            file_flt =
     &           'nfl'//fname2(1:len_trim(fname2))//suffix//'.dat'
            file_hist =
     &           'his'//fname2(1:len_trim(fname2))//suffix//'.dat'

            do i=1,nmolty
               write(ftemp,*) i
               read(ftemp,*) fname3
               file_ndis(i) = 'n'//fname3(1:len_trim(fname3))
     &              //'dis'//fname2(1:len_trim(fname2)) 
     &              //suffix//'.dat'
            end do

            open(unit=50, file = file_flt, status='unknown')  
            close(unit=50)
         
            open(unit=51,file=file_hist,status='unknown')
            write(51,'(f8.4,2x,i5,2x,g15.5,3f12.3)') 
     &           temp, nmolty, (temp*log(B(i)),i=1,nmolty), boxlx(1),
     &           boxly(1), boxlz(1)         
            close(unit=51)
         end if

! --- extra zero accumulator for grand-canonical ensemble

         nentry = 0
         nconfig = 0
         do itmax = 1,ntmax
            nminp(itmax) = 1e6
            nmaxp(itmax) = -1e6
            do imax=0,nmax
               ndist(imax,itmax) = 0
            end do
         end do

      end if

! *** zero accumulators ***
      do i = 1,11
        do ibox = 1,nbox-1
           do jbox = ibox+1, nbox
                acsolpar(i,ibox,jbox)=0.0d0
           end do
        end do
      end do  

      do i=1,nbox 
         do j=1,nener
            acv(j,i) = 0.0d0
            acvsq(j,i) = 0.0d0
            acvkjmol(j,i) = 0.0d0
	 end do
         aflv(i) = 0.0d0
         do j = 1, nmolty
            acchem(i,j) = 0.0d0
            bnchem(i,j) = 0.0d0

            solcount(i,j) = 0
            avsolinter(i,j) = 0.0d0
            avsolintra(i,j) = 0.0d0
            avsolbend(i,j) = 0.0d0
            avsoltor(i,j) = 0.0d0
            avsolelc(i,j) = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM     someone reversed these indices
!!!!!!!!!!!!!!!!!!!!!!!!!!  acntrax(ntmax,nbxmax) in blkavg.inc
            acntrax(j,i) = 0.d0
            acntray(j,i) = 0.d0
            acntraz(j,i) = 0.d0
            acnrotx(j,i) = 0.d0
            acnroty(j,i) = 0.d0
            acnrotz(j,i) = 0.d0
            acstrax(j,i) = 0.d0
            acstray(j,i) = 0.d0
            acstraz(j,i) = 0.d0
            acsrotx(j,i) = 0.d0
            acsroty(j,i) = 0.d0
            acsrotz(j,i) = 0.d0
            
         end do
         acpres(i) = 0.0d0
         acsurf(i) = 0.0d0
         acEnthalpy(i) = 0.0d0
         acEnthalpy1(i) = 0.0d0
      end do

      inst_enth = 0.0d0
      inst_enth2 = 0.0d0
      inst_energy = 0.0d0
      inst_energy2 = 0.0d0     

      acnp = 0.0d0
      do ibox = 1, nbox

         acsvol(ibox) = 0.d0
         acnvol(ibox) = 0.d0
         acvol(ibox) = 0.0d0
         acvolsq(ibox) = 0.0d0

         acvolume(ibox) = 0.0d0

         if (lsolid(ibox) .and. .not. lrect(ibox)) then
            do j = 1,9
               acshmat(ibox,j) = 0.0d0
               acnhmat(ibox,j) = 0.0d0
            end do
         end if
         
      end do

      acmove = 0.0d0

! *** 2nd viral coefficient
      if (lvirial) then
         do i=1,maxvir
            do j = 1, ntemp
               binvir(i,j) = 0.0d0
               binvir2(i,j) = 0.0d0
            end do
         end do
!         if ( lvirial2 ) then
!            call virial2(binvir,binvir2,nvirial,starvir,stepvir)
!            goto 2000
!         end if
      end if
! *** permanent accumulators for box properties ***  
      do ibox = 1, nbox
         do i = 1,3
            acboxl(ibox,i) = 0.0d0
            acboxa(ibox,i) = 0.0d0
         end do
      end do
      
      do i = 1, nmolty
         do ibox = 1, nbox
            asetel(ibox,i) = 0.0d0
            mnbox(ibox,i)  = 0
            acdens(ibox,i) = 0.0d0
            molfra(ibox,i) = 0.0d0
            acnbox(ibox,i) = 0.0d0
         end do
         if ( lexpand(i) ) then
            do itype = 1, numcoeff(i)
               do ibox = 1,2
                  acnbox2(1,i,itype) = 0.0d0
                  acnbox2(2,i,itype) = 0.0d0
               end do
            end do
         end if
      end do


! *** temporary accumulators for max. displacement updates ***
      do im=1,nbox
         do imolty = 1,nmolty
            bstrax(imolty,im) = 0.0d0
            bstray(imolty,im) = 0.0d0
            bstraz(imolty,im) = 0.0d0
            bsrotx(imolty,im) = 0.0d0
            bsroty(imolty,im) = 0.0d0
            bsrotz(imolty,im) = 0.0d0
            bsexpc(imolty,im) = 0.0d0
            bntrax(imolty,im) = 0.0d0
            bntray(imolty,im) = 0.0d0
            bntraz(imolty,im) = 0.0d0
            bnrotx(imolty,im) = 0.0d0
            bnroty(imolty,im) = 0.0d0
            bnrotz(imolty,im) = 0.0d0
            bnexpc(imolty,im) = 0.0d0
         end do
      end do

! *** For the atom displacements
      Abstrax = 0.0d0
      Abstray = 0.0d0
      Abstraz = 0.0d0
      Abntrax = 0.0d0
      Abntray = 0.0d0
      Abntraz = 0.0d0

      do ibox = 1, nbox
         bsvol(ibox) = 0.0d0
      end do

      do ibox = 1, nbox
         bnvol(ibox) = 0.0d0 
      end do

      do ibox = 1,nbox
         if (lsolid(ibox) .and. .not. lrect(ibox)) then
            do j = 1,9
               bshmat(ibox,j) = 0.0d0
               bnhmat(ibox,j) = 0.0d0
            end do
         end if
      end do

! *** temporary accumulators for conf.bias performance ***
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
      DO imolty=1,ntmax
         DO j=1,npabmax
            DO ibox=1,nbxmax*2
               bnswap(imolty,j,ibox)=0.0d0
               bsswap(imolty,j,ibox)=0.0d0
            end do
         end do
      end do
      do i = 1, nmolty
         do j = 1,nswapb(i)
            do ibox = 1, 2*nbox
               bsswap(i,j,ibox) = 0.0d0
               bnswap(i,j,ibox) = 0.0d0
            end do
         end do
         bnswap_in(i,1) = 0.0d0
         bnswap_out(i,1) = 0.0d0
         bnswap_in(i,2) = 0.0d0
         bnswap_out(i,2) = 0.0d0
         bsregr(i,1) = 0.0d0
         bsregr(i,2) = 0.0d0
         bnregr(i) = 0.0d0 
         do inb = 1, numax
            bncb ( i,inb ) = 0.0d0
            bscb ( i,1,inb ) = 0.0d0
            bscb ( i,2,inb ) = 0.0d0
         end do
      end do
!     --- accumulators for swatch performance
      do i = 1, nswaty
         do j = 1,nswtcb(i)
            bnswat(i,j) = 0.0d0
            bsswat(i,j) = 0.0d0
            bnswat_empty(i,j) = 0.0d0
         end do
      end do
!     --- accumulators for fluctuating charge performance
      do i = 1, nmolty
         do j = 1,nbox
            bnflcq(i,j) = 0.0d0
            bnflcq2(i,j) = 0.0d0
            bsflcq(i,j) = 0.0d0
            bsflcq2(i,j) = 0.0d0
         end do
      end do
! --- accumulators for block averages ---
      do i = 1, nprop
         do j = 1, nbox
            naccu(i,j) = 0.0d0
            accum(i,j) = 0.0d0
            nccold(i,j) = 0.0d0
            bccold(i,j) = 0.0d0
            dsq(i,j) = 0.0d0
         end do
      end do
      do i = 1, nprop1
         do ibox = 1,nbox-1
            do jbox = ibox+1,nbox       
               naccu1(i,ibox,jbox) = 0.0d0
               accum1(i,ibox,jbox) = 0.0d0
               nccold1(i,ibox,jbox) = 0.0d0
               bccold1(i,ibox,jbox) = 0.0d0
               dsq1(i,ibox,jbox) = 0.0d0
            end do
         end do 
      end do
      
      nblock = 0
      
! -----------------------------------------------------------------

      if (lneighbor) then
         do i = 1,maxneigh
            do ii = 1,nmax
               neighbor(i,ii) = 0
            end do
         end do
      end if
      
! *** calculate initial energy and check for overlaps ***
      do ibox=1,nbox
         call sumup( ovrlap, v, vinter,vtail, vintra,vvib,
     &                  vbend,vtg,vext,velect,vflucq, ibox, .false.)

         vbox(ibox) = v
         vinterb(ibox)  = vinter
         vtailb(ibox)   = vtail
         vintrab(ibox)  = vintra
         vvibb(ibox)    = vvib  
         vbendb(ibox)   = vbend
         vtgb(ibox)     = vtg  
         vextb(ibox)    = vext 
         velectb(ibox)  = velect
         vflucqb(ibox)  = vflucq
         vipswb(ibox) = vipsw
         vwellipswb(ibox) = vwellipsw
!     kea
         v3garob(ibox) = v3garo
         
         if( ovrlap ) then
            write(iou,*) ' overlap in initial configuration '
            call cleanup('')
         end if
         vstart(ibox) = vbox(ibox)
         if (myid.eq.0) then
            write(iou,*)
            write(iou,*) 'box  ',ibox,' initial v   = ', vbox(ibox)
         end if
         if ( lneigh ) then
! ***        call for initial set-up of the near-neighbour bitmap ***
            call setnn (ibox)
         end if
! *** calculate initial pressure ***
         call pressure ( press1, surf, ibox )
         if (myid.eq.0) then
            write(iou,74) ibox, surf
            write(iou,64) ibox, press1
         end if
      end do

      if (myid.eq.0) then
         write(iou,*)
         write(iou,*) '+++++ start of markov chain +++++'
         write(iou,*)
         write(iou,*) 
     &      'Cycle   Total   Energy    Boxlength   Pressure  Molecules'
!     set up info at beginning of fort.12 for analysis
         write(12,*) nstep,nmolty,(masst(i),i=1,nmolty)
      end if
! *******************************************************************
! ** loops over all cycles and all molecules                       **
! *******************************************************************
 
      nend = nnstep + nstep
      if (lneighbor) then
         write(21,*) 'ii:',ii,(neigh_cnt(i),i=1,nchain)
      end if
        
      do 100 nnn = 1, nstep
         do 99 ii = 1, nchain 

            tmcc = nnstep + nnn - 1
            
!            write(iou,*) 'nstep',(nnn-1)*nchain+ii
! ***       select a move-type at random ***
            rm = random()

! ###       special ensemble dependent moves ###
            if  (rm .le. pmvol) then
!           ---  volume move ---
               if ( lnpt ) then
                  call prvolume  
               else    
                  call volume
               end if

            elseif (rm .le. pmswat) then
!           --- CBMC switch move ---
               call swatch
            elseif ( rm .le. pmswap ) then
!           --- swap move for linear and branched molecules ---
               call swap(bsswap,bnswap,bnswap_in,bnswap_out,
     &              cnt_wf1,cnt_wf2,cnt_wra1,cnt_wra2,qelect)
            elseif ( rm .le. pmcb ) then
!           --- configurational bias move ---
               call config
            elseif ( rm .le. pmflcq ) then
!           --- displacement of fluctuating charges ---
               call flucq(2,0)

            elseif (rm .le. pmexpc ) then
!           --- expanded-ensemble move ---
               call expand

            elseif (rm .le. pmexpc1 ) then
!           --- new expanded-ensemble move ---
!               call expand
               if (random().le.eeratio) then
                  call ee_index_swap
               else
                  call eemove
               end if
            elseif ( rm .le. pm_atom_tra) then
               rm = 3.0d0 * random()
                 if ( rm .le. 1.0d0 ) then
                     call Atom_traxyz (.true.,.false.,.false.)
                 elseif ( rm .le. 2.0d0 ) then
                     call Atom_traxyz (.false.,.true.,.false.)
                 else
                     call Atom_traxyz (.false.,.false.,.true.)
                 end if
            elseif ( rm .le. pmtra ) then
!           --- translational move in x,y, or z direction ---
                 rm = 3.0d0 * random()
                 if ( rm .le. 1.0d0 ) then
                    call traxyz (.true.,.false.,.false.)	
                 elseif ( rm .le. 2.0d0 ) then
                     call traxyz (.false.,.true.,.false.)
                 else
                     call traxyz (.false.,.false.,.true.)
                 end if
            else
!           --- rotation around x,y, or z axis move --
               rm = 3.0d0 * random()
               if ( rm .le. 1.0d0 ) then
                  call rotxyz(.true.,.false.,.false.)
               elseif ( rm .le. 2.0d0 ) then
                  call rotxyz(.false.,.true.,.false.)
               else
                  call rotxyz(.false.,.false.,.true.)
               end if
            end if

! RP added for MPI
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                   
            acmove = acmove + 1.0d0
! ***       accumulate probability of being in an expanded ensemble
! ***       state

            if (lexpee) then
                ee_prob(mstate) = ee_prob(mstate)+1
            end if

! ***       calculate instantaneous values ***
! ***       accumulate averages ***
            
	    do ibox=1,nbox
               do itype = 1, nmolty
                  acnbox(ibox,itype) = acnbox(ibox,itype) + 
     &                 dble(ncmt(ibox,itype))
                  if ( lexpand(itype) ) then
                     do itype2 = 1, numcoeff(itype)
                        acnbox2(ibox,itype,itype2) = 
     &                       acnbox2(ibox,itype,itype2) +
     &                       dble(ncmt2(ibox,itype,itype2))
!                        write(iou,*) '1:',acnbox2(ibox,itype,itype2)
                     end do
!                     write(iou,*) '2:', acnbox(ibox,itype)
                  end if
               end do

               acv(1,ibox)    = acv(1,ibox)   + vbox(ibox)
               acvsq(1,ibox)  = acvsq(1,ibox) + vbox(ibox)**2
               acv(2,ibox)    = acv(2,ibox)   + vinterb(ibox)
               acvsq(2,ibox)  = acvsq(2,ibox) + vinterb(ibox)**2
               acv(3,ibox)    = acv(3,ibox)   + vbendb(ibox)
               acvsq(3,ibox)  = acvsq(3,ibox) + vbendb(ibox)**2
               acv(4,ibox)    = acv(4,ibox)   + vtgb(ibox)
               acvsq(4,ibox)  = acvsq(4,ibox) + vtgb(ibox)**2
               acv(5,ibox)    = acv(5,ibox)   + vintrab(ibox)
               acvsq(5,ibox)  = acvsq(5,ibox) + vintrab(ibox)**2
               acv(6,ibox)    = acv(6,ibox)   + vextb(ibox)
               acvsq(6,ibox)  = acvsq(6,ibox) + vextb(ibox)**2
               acv(7,ibox)    = acv(7,ibox)   + vvibb(ibox)
               acvsq(7,ibox)  = acvsq(7,ibox) + vvibb(ibox)**2
               acv(8,ibox)    = acv(8,ibox)   + velectb(ibox)
               acvsq(8,ibox)  = acvsq(8,ibox) + velectb(ibox)**2
               acv(9,ibox)    = acv(9,ibox)   + vtailb(ibox)
               acvsq(9,ibox)  = acvsq(9,ibox) + vtailb(ibox)**2
               acv(10,ibox)    = acv(10,ibox)   + vflucqb(ibox)
               acvsq(10,ibox)  = acvsq(10,ibox) + vflucqb(ibox)**2
! KEA added 17 for v3garo
               acv(17,ibox)    = acv(17,ibox) + v3garob(ibox)
               acvsq(17,ibox)  = acvsq(17,ibox) + v3garob(ibox)**2

! leftover from Bin, not currently used
                if ( ldielect ) then
                  acv(11,ibox) = acv(11,ibox)+dipolex(ibox)
                  acvsq(11,ibox) = acvsq(11,ibox)+dipolex(ibox)**2
                  acv(12,ibox) = acv(12,ibox)+dipoley(ibox)
                  acvsq(12,ibox) = acvsq(12,ibox)+dipoley(ibox)**2
                  acv(13,ibox) = acv(13,ibox)+dipolez(ibox)
                  acvsq(13,ibox) = acvsq(13,ibox)+dipolez(ibox)**2
                  acvsq(14,ibox) = acvsq(11,ibox) + acvsq(12,ibox)
     &                 + acvsq(13,ibox)
                  acv(15,ibox) = acv(15,ibox)+dsqrt(dipolex(ibox)*
     &                 dipolex(ibox)+dipoley(ibox)*dipoley(ibox)+
     &                 dipolez(ibox)*dipolez(ibox))
               end if

               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  temvol = cell_vol(ibox) 
               else
                  if ( lpbcz ) then
                     temvol = boxlx(ibox)*boxly(ibox)*boxlz(ibox) 
                  else
                     temvol = boxlx(ibox)*boxly(ibox)
                  end if
               end if

! KMB/KEA Energy in kJ/mol
               Temp_nmol = 0
               do itype=1,nmolty
                  Temp_nmol =   Temp_nmol + ncmt(ibox,itype)
               end do
               Temp_Energy  = (vbox(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(1,ibox) = acvkjmol(1,ibox) + Temp_Energy
               Temp_Energy  = (vinterb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(2,ibox) = acvkjmol(2,ibox) + Temp_Energy
               Temp_Energy  = (vbendb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(3,ibox) = acvkjmol(3,ibox) + Temp_Energy
               Temp_Energy  = (vtgb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(4,ibox) = acvkjmol(4,ibox) + Temp_Energy
               Temp_Energy  = (vintrab(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(5,ibox) = acvkjmol(5,ibox) + Temp_Energy
               Temp_Energy  = (vextb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(6,ibox) = acvkjmol(6,ibox) + Temp_Energy
               Temp_Energy  = (vvibb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(7,ibox) = acvkjmol(7,ibox) + Temp_Energy
               Temp_Energy  = (velectb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(8,ibox) = acvkjmol(8,ibox) + Temp_Energy
               Temp_Energy  = (vtailb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(9,ibox) = acvkjmol(9,ibox) + Temp_Energy
               Temp_Energy  = (vflucqb(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(10,ibox) = acvkjmol(10,ibox) + Temp_Energy
               Temp_Energy  = (v3garob(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acvkjmol(17,ibox) = acvkjmol(17,ibox) + Temp_Energy

               if ( lnpt ) then
!                  acv(16,ibox) = acv(16,ibox) + vbox(ibox)*boxlx(ibox)
!     &                 *boxly(ibox)*boxlz(ibox)
                  acv(16,ibox) = acv(16,ibox) + vbox(ibox)*temvol
               end if


               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  boxlx(ibox) = cell_length(ibox,1) 
                  boxly(ibox) = cell_length(ibox,2)
                  boxlz(ibox) = cell_length(ibox,3)

                  acboxa(ibox,1) = acboxa(ibox,1) + cell_ang(ibox,1)
                  acboxa(ibox,2) = acboxa(ibox,2) + cell_ang(ibox,2)
                  acboxa(ibox,3) = acboxa(ibox,3) + cell_ang(ibox,3)

               end if

               acboxl(ibox,1) = acboxl(ibox,1) + boxlx(ibox)
               acboxl(ibox,2) = acboxl(ibox,2) + boxly(ibox)
               acboxl(ibox,3) = acboxl(ibox,3) + boxlz(ibox)

               acvol(ibox) = acvol(ibox) + temvol
               acvolsq(ibox) = acvolsq(ibox) + temvol*temvol
               acvolume(ibox) = acvolume(ibox) + temvol

               do itype = 1, nmolty
                  acdens(ibox,itype) = acdens(ibox,itype) +
     &                 ncmt(ibox,itype) / temvol
                  if ( nchbox(ibox) .gt. 0 ) then
                     molfra(ibox,itype) = molfra(ibox,itype) +
     &                    dble(ncmt(ibox,itype)) / dble(nchbox(ibox))
                  end if
               end do
	    end do


            if (lstop) then
	       if ( int(acmove) .ge. nstep ) goto 101
	    end if

! *************************************************************
! ** ends loop over chains                                   **
! *************************************************************

        If  ( lnpt.and..not.lgibbs ) then
            ibox = 1
            tmp= vbox(ibox) + express(ibox) * ( boxlx(ibox)*boxly(ibox)
     &        * boxlz(ibox))
!         write(49,*) nnn, tmp
         inst_enth=inst_enth+ tmp
         inst_enth2=inst_enth2+(tmp*tmp)
         end if

         If  (.not. lnpt .and..not.lgibbs ) then
            ibox = 1
            tmp= vbox(ibox) 
!         write(49,*) nnn, tmp
         inst_energy=inst_energy + tmp
         inst_energy2=inst_energy2+ (tmp*tmp)
         end if


! - collect histogram data (added 8/30/99)
            if(lgrand) then
               ibox = 1
               vhist = vinterb(ibox) +
     &              velectb(ibox) + vflucqb(ibox)	      
               nconfig = nconfig + 1
! KM for MPI
! check this if want to run grand canonical
               if (mod(nconfig,ninstf).eq.0.and.myid.eq.0) then 

	         open(unit=50, file = file_flt, status='old', 
     &					position='append')
     
!	         write(50, '(i10,5x,i7,5x,g15.6)') nconfig, 
!     +			  (ncmt(ibox,i), i=1,nmolty),vhist

! Formatting removed until I can figure out how to get <n>i7
! to work
	         write(50, * ) nconfig, 
     &			  (ncmt(ibox,i), i=1,nmolty),vhist
	         close(unit=50)
	      end if
	   	    
	      if(mod(nconfig,ninsth).eq.0.and.nconfig.gt.nequil) then

	      do imolty = 1, nmolty
	         nminp(imolty) = min(nminp(imolty),ncmt(ibox,imolty))
		 nmaxp(imolty) = max(nmaxp(imolty),ncmt(ibox,imolty))
	      end do
                 nentry = nentry + 1

	         do imolty=1,nmolty
                    ncmt_list(nentry,imolty) = ncmt(ibox,imolty)
		    ndist(ncmt(ibox,imolty),imolty) = 
     &	   	    ndist(ncmt(ibox,imolty),imolty) + 1
		 end do

		 eng_list(nentry) = vhist
			 			
	      end if

	      if (mod(nconfig,ndumph).eq.0.and.myid.eq.0) then
	         open(unit=51,file = file_hist,status='old',
     &	    	      position='append') 
                 do i=1,nentry 
		    write(51, * ) 
     &		      (ncmt_list(i,imolty), imolty=1,nmolty), 
     &					eng_list(i)     
		 end do
		 close(unit=51)
		 nentry = 0

		 do imolty=1,nmolty
		   open(unit=52, file=file_ndis(imolty), status='unknown')
		      do n=nminp(imolty),nmaxp(imolty)
		         write(52,*) n,ndist(n,imolty)
		      end do
		      close(unit=52)
		 end do
	      end if
	 end if

!9       continue

! *** perform periodic operations  ***
 
         if ( lgibbs .or. lnpt .and. (.not. lvirial) ) then
            do ibox = 1,nbox
               if ( lpbcz ) then
                  if (lsolid(ibox) .and. .not. lrect(ibox).and.
     &                 myid.eq.0) then
                     write(12,'(7e13.5,15i5)') hmat(ibox,1)
     &                    ,hmat(ibox,4),hmat(ibox,5)
     &                    ,hmat(ibox,7),hmat(ibox,8)
     &                    ,hmat(ibox,9),vbox(ibox),
     &                    (ncmt(ibox,itype),itype=1,nmolty)
                     
                     open(unit=13,file = file_cell,status='old',
     &                    position='append')
                     write(13,'(i8,6f12.4)') nnn+nnstep,
     &                    cell_length(ibox,1)/Num_cell_a,
     &                    cell_length(ibox,2)/Num_cell_b,
     &                    cell_length(ibox,3)/Num_cell_c,
     &                    cell_ang(ibox,1)*180.0d0/onepi,
     &                    cell_ang(ibox,2)*180.0d0/onepi,
     &                    cell_ang(ibox,3)*180.0d0/onepi
                     close(unit=13)

!                     write(13,'(i8,3f12.4)') nnn,cell_ang(ibox,1)
!     +                        , cell_ang(ibox,2),cell_ang(ibox,3)
                  else
!                  do ibox = 1, nbox
                     if (myid.eq.0) then
                        write(12,'(4e13.5,15i5)')boxlx(ibox),boxly(ibox)
     &                       ,boxlz(ibox),vbox(ibox),
     &                       (ncmt(ibox,itype),itype=1,nmolty)
                     end if
!                  end do
                  end if
               else
!               do ibox = 1, nbox
                  if (myid.eq.0) then
                     write(12,'(2e12.5,15i4)') boxlx(ibox)*boxly(ibox)
     &                    ,vbox(ibox),(ncmt(ibox,itype),itype=1,nmolty)
                  end if
!               end do
               end if
            end do
         end if

         if (lucall) then
            write(iou,*) 'not recently checked for accuracy'
            call cleanup('')
!            do j = 1,nmolty
!               if ( ucheck(j) .gt. 0 ) then
!                  call chempt(bsswap,j,ucheck(j),qelect)
!               end if
!            end do
         end if

         if ( mod(nnn,iratp) .eq. 0 ) then
! *** calculate pressure ***
            acnp = acnp + 1.0d0
            do ibox = 1, nbox
               call pressure ( press1, surf, ibox )
!              write(iou,*) 'control pressure', press1
               pres(ibox) = press1
               acpres(ibox) = acpres(ibox) + press1
               acsurf(ibox) = acsurf(ibox) + surf
            end do

! Enthalpy calculation

            do ibox = 1,nbox  
               Temp_nmol = 0
               do itype=1,nmolty
                  Temp_nmol =   Temp_nmol + ncmt(ibox,itype)
               end do
               Temp_Mol_Vol = temvol/Temp_nmol*0.6022d-06 ! m3/mol
               Temp_Energy  = (vbox(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acEnthalpy(ibox) = acEnthalpy(ibox) + Temp_Energy +
     &                       pres(ibox)*Temp_Mol_Vol     !kJ/mol
               acEnthalpy1(ibox) = acEnthalpy1(ibox) + Temp_Energy +
     &                       (express(ibox)/7.2429d-5)*Temp_Mol_Vol     !kJ/mol   
            end do

! --- cannot calculate a heat of vaporization for only one box,
! --- and some compilers choke because Heat_vapor_T will not be
! --- defined if nbox == 1
            if(lgibbs) then
               do ibox = 1,nbox-1
                  do jbox = ibox+1,nbox
!                     WRITE(6,*) 'ieouwfe ',ibox,jbox
                     call calcsolpar(pres,Heat_vapor_T,Heat_vapor_LJ,
     &                    Heat_vapor_COUL,pdV, CED_T,CED_LJ,CED_COUL,
     &                    HSP_T,HSP_LJ, HSP_COUL,ibox,jbox)

! --- Heat of vaporization            
                     acsolpar(1,ibox,jbox)=
     &                    acsolpar(1,ibox,jbox)+Heat_vapor_T
                     acsolpar(2,ibox,jbox)=
     &                    acsolpar(2,ibox,jbox)+Heat_vapor_LJ
                     acsolpar(3,ibox,jbox)=
     &                    acsolpar(3,ibox,jbox)+Heat_vapor_COUL
                     acsolpar(4,ibox,jbox)=
     &                    acsolpar(4,ibox,jbox)+CED_T
                     acsolpar(5,ibox,jbox)=
     &                    acsolpar(5,ibox,jbox)+CED_LJ
                     acsolpar(6,ibox,jbox)=
     &                    acsolpar(6,ibox,jbox)+CED_COUL
                     acsolpar(7,ibox,jbox)=
     &                    acsolpar(7,ibox,jbox)+HSP_T
                     acsolpar(8,ibox,jbox)=
     &                    acsolpar(8,ibox,jbox)+HSP_LJ
                     acsolpar(9,ibox,jbox)=
     &                    acsolpar(9,ibox,jbox)+HSP_COUL 
!                     acsolpar(10,ibox,jbox) = 
!     &                    acsolpar(10,ibox,jbox)+DeltaU_Ext
                     acsolpar(11,ibox,jbox) =
     &                    acsolpar(11,ibox,jbox)+pdV
                  end do  
               end do
            end if
         end if

! calculate the integrand of thermosynamic integration
         if (lmipsw.and.(mod(nnn,iratipsw).eq.0)) then
            acipsw = acipsw+1.0d0
            call deriv(1)
            acdvdl = acdvdl+dvdl
         end if

 64      format(' pressure check:   box',i2,' =',f14.2)
 74      format(' surf. tension :   box',i2,' =',f14.5)


! *** Add a call for subroutine to compute the 

         lratio = .false.
         lratv = .false.
         lprint = .false.
         lmv = .false.
         lrsave = .false.
         lblock = .false.

         if ( mod(nnn,iratio) .eq. 0 ) then
            lratio = .true.
         end if

         if ( lgibbs .or. lnpt ) then
            if ( mod(nnn,iratv) .eq. 0 .and. pmvol .gt. 0.0d0 ) then
               lratv = .true.
            end if
         end if

         if ( mod(nnn,iprint) .eq. 0 ) then
            lprint = .true.
         end if
 
         if ( mod(nnn,imv) .eq. 0 ) then
            if ( lvirial ) then
               call virial(binvir,binvir2,nvirial,starvir,stepvir)
            else
               lmv = .true.
            end if
         end if
         if ( mod(nnn,irsave) .eq. 0 ) then
            lrsave = .true.
         end if
 
!--- JLR 11-11-09
!--- do not call analysis if ianalyze is greater than number of cycles
! KM 01/10 remove analysis
!	 if(mod(nnn,ianalyze).eq.0) then
!	    call analysis(1)
!         end if  
!--- END JLR 11-11-09

         do intg = 1, nchain
            ibox = nboxi(intg)
            imolty = moltyp(intg)
! *** accumulate m-n-box and m-s-e-t-e-l ***
!     only count the main chain - not the hydrogens
            ilunit = nugrow(imolty)
            setx = rxu(intg,1) - rxu(intg,ilunit)
            sety = ryu(intg,1) - ryu(intg,ilunit)
            setz = rzu(intg,1) - rzu(intg,ilunit)
            setel = setx*setx + sety*sety + setz*setz
!            if ( imolty .eq. 2 ) then
!               write(??,*) imolty,setel
!            end if
            mnbox( ibox, imolty ) = mnbox( ibox, imolty ) + 1
            asetel( ibox, imolty ) = asetel( ibox, imolty ) + setel
         end do

         if ( mod(nnn,iblock) .eq. 0 ) then
            lblock = .true.
         end if
         
         ltsolute = .false.
         ltfix    = .false.
         
         do imolty = 1, nmolty
            if (pmfix(imolty).gt.0.0001d0.and
     &           .mod(nnn,iupdatefix).eq.0) then
               ltfix = .true.
               lratfix(imolty) = .true.
            else
               lratfix(imolty) = .false.
            end if

            if (mod(nnn,isolute(imolty)).eq.0) then
               ltsolute = .true.
               lsolute(imolty) = .true.
            else
               lsolute(imolty) = .false.
            end if
         end do

         if (lratio .or. lratv .or. lprint .or. lmv .or. lrsave
     &        .or. lblock .or. ltfix .or. ltsolute) then
            call monper(acv,acpres,acsurf,acvolume,molfra,mnbox,asetel
     &       ,acdens,acmove,acnp,pres,nbox,nnn,nblock,lratio,lratv
     &       ,lprint,lmv,lrsave,lblock,lratfix,lsolute,acsolpar,
     &        acEnthalpy,acEnthalpy1)
         end if
!     not currently used
         if (ldielect.and.(mod(nnn,idiele).eq.0).and.myid.eq.0) then
            dielect = acvsq(14,ibox)/acmove

! *If you really want this quantity comment should be taken out**

!            write(14,*) nnn,6.9994685465110493d5*dielect*beta/
!     &           (boxlx(ibox)*boxly(ibox)*boxlz(ibox))

            dielect = acvsq(14,ibox)/acmove  
     &           -(acv(11,ibox)/acmove)**2
     &           -(acv(12,ibox)/acmove)**2 - (acv(13,ibox)/acmove)**2
! ** use fort.27 to calculate dielectric constant
!            write(15,*) nnn,6.9994685465110493d5*dielect*beta/
!     &           (boxlx(ibox)*boxly(ibox)*boxlz(ibox))
!            write(16,*) nnn,acv(11,ibox)/acmove, acv(12,ibox)/acmove,
!     &                acv(13,ibox)/acmove
         end if
         if ( lnpt .and. nmolty .eq. 1 ) then
! *** output the fluctuation information
            if (mod(nnn,idiele) .eq. 0.and.myid.eq.0) then
               write(14,*) nnn,acvol(ibox)/acmove
               write(15,*) nnn,acvolsq(ibox)/acmove
               write(16,*) nnn,acv(1,ibox)/acmove
               write(17,*) nnn,acvsq(1,ibox)/acmove
               write(18,*) nnn,acvol(ibox)*acv(1,ibox)/(acmove*acmove)
               write(19,*) nnn,acv(16,ibox)/acmove - acvol(ibox)
     &              *acv(1,ibox)/(acmove*acmove)
            end if
         end if
         
! - set idiele = 1 to print every cycle
         if ( ldielect .and. mod(nnn,idiele).eq. 0.and.myid.eq.0) then
            do ibox = 1,nbox
               write(27,*) dipolex(ibox),dipoley(ibox),dipolez(ibox)
            end do
         end if

!         if ( mod(nnn,idiele) .eq. 0 ) write(25,*) nnn+nnstep,vbox(1)

!         ibox = 1
!         imolty = 1
!         do i = 1,nchain
!            if ( nboxi(i) .eq. ibox ) then
!               if ( moltyp(i) .eq. imolty ) then
!                  bin = dint(zcm(i)/binstep) + 1
!                  temvol = boxlx(ibox)*boxly(ibox)*binstep
!                  profile(bin) = profile(bin)+1.0d0/temvol
!               end if
!            end if
!         end do

!--Residual Heat capacity ---
        if(mod(nnn,iheatcapacity) .eq. 0) then

          if(lnpt.and..not.lgibbs) then
             enthalpy= inst_enth/(dble(nchain*nnn))
             enthalpy2= inst_enth2/(dble(nchain*nnn))
             sigma2Hsimulation=(enthalpy2)-(enthalpy*enthalpy)
             sigma2H=sigma2Hsimulation*(6.022d23)*((1.38066d-23)**2) /
     &            (dble(nchain)) !(J2/mol)
             Cp=sigma2H/((1.38066d-23)*(temp**2))
             if (myid.eq.0) then
                write(56,'(I12,F18.6,F18.2,F18.6)')nnn,Cp,enthalpy2,
     &               enthalpy
             end if
       
          elseif( .not.lnpt .and..not.lgibbs) then
             ennergy= inst_energy/(dble(nchain*nnn))
             ennergy2= inst_energy2/(dble(nchain*nnn))
             sigma2Esimulation=(ennergy2)-(ennergy*ennergy)
             sigma2E=sigma2Esimulation*(6.022d23)*((1.38066d-23)**2) /
     &            (dble(nchain)) !(J2/mol)
             Cv=sigma2E/((1.38066d-23)*(temp**2))
             if (myid.eq.0) then
                write(55,'(I12,F18.6,F18.2,F18.6)')nnn,Cv,ennergy2,
     &               ennergy
             end if
          end if
       end if



!00   continue

      lpr = .false.
      do i = 1,ntmax
         if (lbias(i)) then
            lpr = .true.
         end if
      end do

      if (lpr.and.myid.eq.0) then
         do nnn = 1,4
            write(31,*)
            write(31,*) 'nnn:',nnn
            do i = 0,6
               do j = 0,6
                  if (cnt_wf1(i,j,nnn) .gt. 0 ) then
                     write(31,*) i,j,cnt_wf1(i,j,nnn),cnt_wf2(i,j,nnn)
                  end if
               end do
            end do
            write(32,*) 
            write(32,*) 'nnn:', nnn
            write(33,*)
            write(33,*) 'nnn:', nnn
            
            do i = 1,1000
               if ( cnt_wra1(i,nnn) .gt. 0 ) 
     &              write(32,*) (dble(i)-0.5d0)*
     &              0.1d0-95.0d0,cnt_wra1(i,nnn)
               if ( cnt_wra2(i,nnn) .gt. 0 ) 
     &              write(33,*) (dble(i)-0.5d0)*
     &              0.1d0-95.0d0,cnt_wra2(i,nnn)
            end do
         end do
      end if

!      do bin = 1,1000
!         write(26,*) binstep*(dble(bin)-0.5d0),profile(bin)/nstep
!      end do

! *******************************************************************
! ** ends the loop over cycles                                     **
! *******************************************************************
  101 continue 

      if (lneighbor) then
         write(21,*) 'ii:',ii,(neigh_cnt(i),i=1,nchain)
      end if
      if (myid.eq.0) then
         write(iou,*)
         write(iou,*) '+++++ end of markov chain +++++'
 
! *** write some information about translations and rotations
         write(iou,*)
         write(iou,*) '### Translations ###'
         write(iou,*)
         do ibox = 1,nbox
            do i=1,nmolty
               write(iou,*) 'molecule typ =',i,' in box',ibox
               acntrax(i,ibox) = acntrax(i,ibox) + bntrax(i,ibox)
               acstrax(i,ibox) = acstrax(i,ibox) + bstrax(i,ibox)
               if ( acntrax(i,ibox) .ne. 0.0d0 ) then
                  ratvol = acstrax(i,ibox) / acntrax(i,ibox)
               else
                  ratvol = 0.0d0
               end if
               write(iou,71) acntrax(i,ibox),ratvol,rmtrax(i,ibox)
               
               acntray(i,ibox) = acntray(i,ibox) + bntray(i,ibox)
               acstray(i,ibox) = acstray(i,ibox) + bstray(i,ibox)
               if ( acntray(i,ibox) .ne. 0.0d0 ) then
                  ratvol = acstray(i,ibox) / acntray(i,ibox)
               else
                  ratvol = 0.0d0
               end if
               write(iou,72) acntray(i,ibox),ratvol,rmtray(i,ibox)
               
               acntraz(i,ibox) = acntraz(i,ibox) + bntraz(i,ibox)
               acstraz(i,ibox) = acstraz(i,ibox) + bstraz(i,ibox)
               if ( acntraz(i,ibox) .ne. 0.0d0 ) then
                  ratvol = acstraz(i,ibox) / acntraz(i,ibox)
               else
                  ratvol = 0.0d0
               end if
               write(iou,73) acntraz(i,ibox),ratvol,rmtraz(i,ibox)
               write(iou,*)
               
            end do
         end do
         
         write(iou,*) '### Rotations ###'
         write(iou,*)
         do ibox = 1,nbox
            do i=1,nmolty
               write(iou,*) 'molecule typ =',i,' in box',ibox
               acnrotx(i,ibox) = acnrotx(i,ibox) + bnrotx(i,ibox)
               acsrotx(i,ibox) = acsrotx(i,ibox) + bsrotx(i,ibox)
               if ( acnrotx(i,ibox) .ne. 0.0d0 ) then
                  ratvol = acsrotx(i,ibox) / acnrotx(i,ibox)
               else
                  ratvol = 0.0d0
               end if
               write(iou,71) acnrotx(i,ibox),ratvol,rmrotx(i,ibox)
               
               acnroty(i,ibox) = acnroty(i,ibox) + bnroty(i,ibox)
               acsroty(i,ibox) = acsroty(i,ibox) + bsroty(i,ibox)
               if ( acnroty(i,ibox) .ne. 0.0d0 ) then
                  ratvol = acsroty(i,ibox) / acnroty(i,ibox)
               else
                  ratvol = 0.0d0
               end if
               write(iou,72) acnroty(i,ibox),ratvol,rmroty(i,ibox)
               
               acnrotz(i,ibox) = acnrotz(i,ibox) + bnrotz(i,ibox)
               acsrotz(i,ibox) = acsrotz(i,ibox) + bsrotz(i,ibox)
               if ( acnrotz(i,ibox) .ne. 0.0d0 ) then
                  ratvol = acsrotz(i,ibox) / acnrotz(i,ibox)
               else
                  ratvol = 0.0d0
               end if
               write(iou,73) acnrotz(i,ibox),ratvol,rmrotz(i,ibox)
               write(iou,*)
            end do
         end do
 71      format(' x-dir: attempts =',F10.1,'   ratio =',f6.3,
     &        '   max.displ. =',e10.4)
 72      format(' y-dir: attempts =',F10.1,'   ratio =',f6.3,
     &      '   max.displ. =',e10.4)
 73      format(' z-dir: attempts =',F10.1,'   ratio =',f6.3,
     &        '   max.displ. =',e10.4)
! *** write some information about config performance ***
         if ( pmcb .gt. 0.0d0 ) then
            write(iou,*)
            write(iou,*) '### Configurational-bias ###'
            write(iou,*)
            do i = 1, nmolty
               write(iou,*) 'molecule typ =',i
              write(iou,*) '    length  attempts  succ.growth  accepted'
     &              ,'   %su.gr.    %accep.'
              do inb = 1, nunit(i)
                 if ( bncb(i,inb) .gt. 0.0d0 ) then
                    pscb1 = bscb(i,1,inb) * 100.0d0 / bncb(i,inb)
                    pscb2 = bscb(i,2,inb) * 100.0d0 / bncb(i,inb)
                    write(iou,'(i9,3f10.1,2f10.2)') inb, bncb(i,inb),
     &                   bscb(i,1,inb), bscb(i,2,inb), pscb1, pscb2
                 end if
              end do
              if (pmfix(i).gt.0.0d0) then
                 write(iou,*) ' SAFE-CBMC move '
                 write(iou,*) '    length  attempts  succ.growth  ',
     &                'accepted   %su.gr.    %accep.'               
                 do inb = 1, nunit(i)
                    if (fbncb(i,inb) .gt. 0.0d0 ) then
                       pscb1 = fbscb(i,1,inb) * 100.0d0 
     &                      / fbncb(i,inb)
                       pscb2 = fbscb(i,2,inb) * 100.0d0 
     &                      / fbncb(i,inb)
                       write(iou,'(i9,3f10.1,2f10.2)') inb,fbncb(i,inb),
     &                      fbscb(i,1,inb), fbscb(i,2,inb)
     &                      , pscb1, pscb2
                    end if
                 end do
              end if
           end do
           write(iou,*)
        end if
! *** write some information about volume performance ***
        if ( lgibbs .or. lnpt) then
           write(iou,*)
           write(iou,*) '### Volume change       ###'
           do ibox = 1,nbox
              if (lsolid(ibox) .and. .not. lrect(ibox)) then
                 do j = 1,9
                    acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
                    acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
                    if ( acshmat(ibox,j) .gt. 0.5d0) then
                       write(iou,70) acnhmat(ibox,j),
     &                    acshmat(ibox,j)/acnhmat(ibox,j),rmhmat(ibox,j)
                    else
                       write(iou,70)acnhmat(ibox,j),0.0d0,rmhmat(ibox,j)
                    end if
                 end do
              else
                 acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
                 acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
                 if ( acnvol(ibox) .ne. 0.0d0 ) then
                    ratvol = acsvol(ibox) / acnvol(ibox)
                 else
                    ratvol = 0.0d0
                 end if
                 write(iou,61) acnvol(ibox),ratvol,rmvol(ibox)
              end if
           end do
        end if
        write(iou,*)  
        write(iou,*) '### Molecule swap       ###'
        write(iou,*)
        do i = 1, nmolty
           write(iou,*) 'molecule typ =',i
           do j=1,nswapb(i)
              if ( box1(i,j) .eq. box2(i,j) ) then
                 jbox_max = 1
              else
                 jbox_max = 2
              end if
              do jbox = 1,jbox_max
                 if ( jbox .eq. 1 ) ibox = box1(i,j)
                 if ( jbox .eq. 2 ) ibox = box2(i,j)
                 write(iou,66) box1(i,j),box2(i,j),ibox,
     &                bsswap(i,j,ibox),bnswap(i,j,ibox),
     &                bnswap(i,j,ibox+nbox)
                 if (bnswap(i,j,ibox) .gt. 0.5d0) then
                    bsswap(i,j,ibox) = bsswap(i,j,ibox+nbox)*
     &                   100.0d0/bsswap(i,j,ibox)
                    bnswap(i,j,ibox) = bnswap(i,j,ibox+nbox)*
     &                   100.0d0/bnswap(i,j,ibox)
                    write(iou,63) bsswap(i,j,ibox),bnswap(i,j,ibox)
                 end if
              end do
           end do
           write(iou,68) bnswap_in(i,1), bnswap_in(i,2)
           write(iou,69) bnswap_out(i,1), bnswap_out(i,2)
        end do

        write(iou,*)
        write(iou,*) '### Molecule swatch     ###'
        write(iou,*)
        do i = 1, nswaty
           write(iou,*) 'pair typ =',i
           write(iou,*) 'moltyps = ',nswatb(i,1),' and',nswatb(i,2)
           do j = 1, nswtcb(i)
! --- JLR 12-1-09 changing to exclude empty box attempts from swatch rate 
              write(iou,62) box3(i,j),box4(i,j),
     &             bnswat(i,j),bnswat(i,j)-bnswat_empty(i,j),bsswat(i,j)
              if (bnswat(i,j) .gt. 0.5d0 ) then
                 write(iou,65) 100.0d0 * bsswat(i,j)/
     &                (bnswat(i,j)-bnswat_empty(i,j))
              end if
! --- EN JLR 12-1-09
           end do
        end do

        write(iou,*)
        write(iou,*)    '### Charge Fluctuation  ###'
        write(iou,*)
      
        do i = 1, nmolty
           do j = 1,nbox
              bnflcq2(i,j) = bnflcq2(i,j) + bnflcq(i,j) 
              bsflcq2(i,j) = bsflcq2(i,j) + bsflcq(i,j) 
              if (bnflcq2(i,j) .gt. 0.5d0) then
                 write(iou,*) 'molecule typ =',i,'  box =',j
                 bsflcq2(i,j) = bsflcq2(i,j)/bnflcq2(i,j)
                 write(iou,61) bnflcq2(i,j),bsflcq2(i,j),rmflcq(i,j)
              end if
           end do
        end do

        write(iou,*) 
        write(iou,*)    '### Expanded Ensemble Move  ###'
        write(iou,*) 
        do i = 1, nmolty
           do j = 1,nbox
              if (lexpand(i) .and. bnexpc(i,j) .gt. 0.5) then
                 write(iou,*) 'molecule typ =',i,'  box =',j
                 write(iou,67) bnexpc(i,j),bsexpc(i,j),
     &                bsexpc(i,j)/bnexpc(i,j)
              end if
           end do
        end do
        
 61     format(' attempts =',f8.1,'   ratio =',f6.3,
     &      '   max.displ. =',e10.4)
! --- JLR 12-1-09
!62   format('between box ',i2,'and ',i2,
!     +     ' attempts =',f9.1,'   accepted =',f8.1) 
 62     format('between box ',i2,' and ',i2,
     &       '   uattempts =',f12.1,   '  attempts =',f9.1,
     &       '  accepted =',f8.1)
! --- END JLR 12-1-09
 63     format(' suc.growth % =',f7.3,'   accepted % =',f7.3)
 65     format(' accepted % =',f7.3)
 66     format('between box ',i2,' and ',i2,' into box',i2,
     &       '   uattempts =',f12.1,' attempts =',f9.1
     &       ,'   accepted =',f8.1)
 67     format(' attempts =',f8.1,'   accepted =',f8.1,
     &       ' accepted % =',f7.3)
        
 68     format('number of times move in: ', f12.1, 
     &       '  accepted=',f8.1)
 69     format('number of times move out: ', f12.1, 
     &       '  accepted=',f8.1)
 70     format(' h-matrix attempts =',f8.1,'   ratio =',f6.3,
     &       '   max.displ. =',e10.4)
        if (lexzeo) then
!        --- write end-conf to picture file
           call writepdb(ncmt(1,1),nunit(1),1)
        end if
! *** checks final value of the potential energy is consistent ***

      end if  ! end if (myid.eq.0)
      
      do ibox=1,nbox
         if ( ldielect ) then
! *** store old dipole moment
            call dipole(ibox,2)
         end if

         call sumup( ovrlap, v, vinter,vtail,vintra,vvib,
     &        vbend,vtg,vext,velect,vflucq,ibox, .false.)
         vend(ibox) = v


!---need to check
         if (myid.eq.0) then
            if ( abs(v - vbox(ibox)) .gt. 0.0001) then
               write(iou,*) '### problem with energy ###  box ',ibox
               write(iou,*) ' Total energy: ',v,vbox(ibox),v-vbox(ibox)
            end if
            if ( abs(vinter - vinterb(ibox)) .gt. 0.000001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Inter mol.en.: ',vinter,vinterb(ibox)
               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                 write(iou,*)'You might check the cutoff wrt box widths'
                 write(iou,*) 'Normal PBC might be failing'
               end if
            end if
            if ( abs(vtail - vtailb(ibox)) .gt. 0.000001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Tail corr.en.: ',vtail,vtailb(ibox)
            end if
            if ( abs(vintra - vintrab(ibox)) .gt. 0.000001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Intra mol.en.: ',vintra,vintrab(ibox)
            end if
            if ( abs(vvib - vvibb(ibox)) .gt. 0.001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' bond vib. en.: ',vvib,vvibb(ibox)
            end if
            if ( abs(vbend - vbendb(ibox)) .gt. 0.001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Bond ben.en.: ',vbend,vbendb(ibox)
            end if
            if ( abs(vtg - vtgb(ibox)) .gt. 0.001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Torsion.en.: ',vtg,vtgb(ibox)
            end if
            if ( abs(vext - vextb(ibox)) .gt. 0.0001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Externa.en.: ',vext,vextb(ibox)
            end if
            if ( abs(velect - velectb(ibox)) .gt. 0.000001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Coulomb.en.: ',velect,velectb(ibox)
            end if
            if ( abs(vflucq - vflucqb(ibox)) .gt. 0.0001) then
               write(iou,*) '### problem  ###'
               write(iou,*) ' Fluc Q en.: ',vflucq,vflucqb(ibox)
            end if
            if ( abs(v3garo - v3garob(ibox) ) .gt.0.001) then
               write(iou,*) '### problem ###'
               write(iou,*) ' 3-body en.: ',v3garo,v3garob(ibox)
            end if
            if ( ldielect ) then
               if ( abs(dipolexo - dipolex(ibox)) .gt. 0.0001) then
                  write(iou,*) '### problem  ###'
                  write(iou,*) ' Dipole X: ',dipolexo,dipolex(ibox)
               end if
            end if
            if (lmipsw) then
               if (abs(vwellipsw-vwellipswb(ibox)).gt.0.001) then
                  write(iou,*) '### problem  ###'
                  write(iou,*) ' well en.: ',vwellipsw,vwellipswb(ibox)
               end if
            end if
         end if  ! end if myid.eq.0
      end do

! KM for MPI
! only processor 0 needs to calculate and write out averages, final config, etc

      if (myid.eq.0) then
         write(iou,*)
         write(iou,1501) (vstart(i) ,i=1,nbox)
         write(iou,1502) (vend(i)   ,i=1,nbox)
         write(iou,1504) (vbox(i)   ,i=1,nbox)
         write(iou,*)
            
!     ** normalize and write out presim results in fort.22 **
          
         if (lpresim) then
            if (counttot.eq.0) then
               write(21,*) counthist
            else
               write(21,*) counttot
            end if
            
            
            do j = 1, iring(1) 
               do k = 1, iring(1)
                  if (j.eq.k) goto 150
                  histtot = 0
                  do bin = 1, maxbin
                     hist(j,k,bin) = hist(j,k,bin) + 1.0d0
                     histtot = histtot + hist(j,k,bin)
                  end do
               
                  do bin = 1, maxbin
                     hist(j,k,bin) = hist(j,k,bin) / histtot
                     write(21,*) bin,hist(j,k,bin)
                  end do
 150              continue
               end do
            end do
         end if
         
!     *** put new distribution back into a file
         do imolty = 1, nmolty
            if (pmfix(imolty).gt.0) then
               if (counttot.eq.0) then
                  write(21,*) counthist
               else
                  write(21,*) counttot
               end if
               do j = 1, iring(1) 
                  do k = 1, iring(1)
                     if (j.eq.k) goto 160
                     do bin = 1, maxbin
                        write(21,*) bin,probf(j,k,bin)
                     end do
 160                 continue
                  end do
               end do
            end if
         end do
         

! ** write out the final configuration for each box, Added by Neeraj 06/26/2006 3M ***
         do ibox = 1,nbox
            write(ftemp,*) ibox
            read(ftemp,*) fname4
            fileout = 'box'//fname4(1:len_trim(fname4))//'config'//
     &           fname2(1:len_trim(fname2))
     &           //suffix//'.xyz'
            open (unit=200+ibox,FILE=fileout,status="unknown")
            
            nummol = 0
            do i = 1,nchain
               if (nboxi(i).eq.ibox) then
                  nummol = nummol + nunit(moltyp(i))
               end if
            end do 
            write(200+ibox,*) nummol 
            write(200+ibox,*) 
            do i = 1,nchain
               if(nboxi(i).eq.ibox) then 
                  imolty = moltyp(i) 
                  do ii = 1,nunit(imolty)
                     ntii = ntype(imolty,ii)               
                     write(200+ibox,'(a4,5x,3f15.4)') chemid(ntii),
     &                    rxu(i,ii), ryu(i,ii), rzu(i,ii)
                  end do
               end if
            end do   
            close(200+ibox)
         end do  
         
         
! ** write out the final configuration from the run ***
         
         file_config = 'config'//fname2(1:len_trim(fname2))
     &        //suffix//'.dat'
         open(8, file=file_config,status='unknown')
         
         write(8,*) nend
         if ( nend .gt. 0 ) then
            write(8,*) Armtrax, Armtray, Armtraz 
            do im=1,nbox
               do imolty=1,nmolty
                  write(8,*) rmtrax(imolty,im), rmtray(imolty,im)
     &                 , rmtraz(imolty,im)
                  write(8,*) rmrotx(imolty,im), rmroty(imolty,im)
     &                 , rmrotz(imolty,im)
               end do
            end do
            do im=1, nbox
               write (8,*) (rmflcq(i,im),i=1,nmolty)
            end do
! -- changed formatting so fort.77 same for all ensembles
! -- 06/08/09 KM
            write(8,*) (rmvol(ibox),ibox=1,nbox)
            do ibox = 1,nbox
               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  write(8,*) (rmhmat(ibox,i),i=1,9)
                  write(8,*) (hmat(ibox,i),i=1,9)
               else
                  write(8,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
               end if
            end do
         end if
         
         if (L_add) then
            do i = nchain+1, nchain+N_add    
               moltyp(i)=N_moltyp2add
               nboxi(i) = N_box2add
               rxu(i,1) = random()*boxlx(N_box2add)
               ryu(i,1) = random()*boxly(N_box2add)
               rzu(i,1) = random()*boxlz(N_box2add)
               qqu(i,1) = 0.0
            end do
            nchain = nchain+N_add
            write(8,*) nchain
            write(8,*) nmolty
            write(8,*) (nunit(i),i=1,nmolty)
            write(8,*) (moltyp(i),i=1,nchain)
            write(8,*) (nboxi(i),i=1,nchain)
            do i = 1, nmolty
               if ( lexpand(i) ) write(8,*) eetype(i)
            end do
            do i = 1, nmolty
               if ( lexpand(i) ) write(8,*) rmexpc(i)
            end do
            do  i = 1, nchain
               imolty = moltyp(i)
               do  j = 1, nunit(imolty)
                  write(8,'(4f15.6)') rxu(i,j), ryu(i,j), rzu(i,j), 
     &                 qqu(i,j)
               end do
            end do 
            
         elseif(L_sub) then
            point_of_start = 0
            do i =1,N_moltyp2sub
               point_of_start=point_of_start+temtyp(i)
            end do
            point_of_start = point_of_start-N_sub+1
            
            point_to_end = nchain-N_sub
            
            do i = point_of_start,point_to_end
               nboxi(i) = nboxi(i+N_sub)
               moltyp(i) = moltyp(i+N_sub) 
            end do
            
            write(8,*) nchain-N_sub
            write(8,*) nmolty
            write(8,*) (nunit(i),i=1,nmolty)
            write(8,*) (moltyp(i),i=1,(nchain-N_sub))
            write(8,*) (nboxi(i),i=1,(nchain-N_sub))
            do i = 1, nmolty
               if ( lexpand(i) ) write(8,*) eetype(i)
            end do
            do i = 1, nmolty
               if ( lexpand(i) ) write(8,*) rmexpc(i)
            end do
            do  i = 1, nchain
               if(i.lt.(point_of_start).or.i.gt.
     &              (point_of_start+N_sub-1)) then
                  imolty = moltyp(i)
                  do  j = 1, nunit(imolty)
                     write(8,'(4f15.6)') rxu(i,j), ryu(i,j), rzu(i,j), 
     &                    qqu(i,j)
                  end do
               end if
            end do
         else      
            write(8,*) nchain
            write(8,*) nmolty
            write(8,*) (nunit(i),i=1,nmolty)
            write(8,*) (moltyp(i),i=1,nchain)
            write(8,*) (nboxi(i),i=1,nchain)
            do i = 1, nmolty
               if ( lexpand(i) ) write(8,*) eetype(i)
            end do
            do i = 1, nmolty
               if ( lexpand(i) ) write(8,*) rmexpc(i)
            end do
            do i = 1, nchain
               imolty = moltyp(i)
               do j = 1, nunit(imolty)
                  write(8,'(4f15.6)') rxu(i,j), ryu(i,j), rzu(i,j), 
     &                 qqu(i,j)
               end do 
            end do
         end if
         
         close(8)
         
! *** calculate and write out running averages ***
         do ibox=1,nbox
! - energies
            do j=1,nener
               avv(j,ibox)   = acv(j,ibox) / acmove
               acvsq(j,ibox) = (acvsq(j,ibox)/acmove) - avv(j,ibox) ** 2
               acvkjmol(j,ibox) = acvkjmol(j,ibox)/acmove
            end do
            if ( ldielect ) then
               flucmom(ibox) = acvsq(14,ibox)-avv(15,ibox)*avv(15,ibox)
!            momconst = 6.9994685465110493d5
               flucmom(ibox) = 6.9994685465110493d5*flucmom(ibox)*beta/
     &              (boxlx(ibox)*boxly(ibox)*boxlz(ibox))
               flucmom2(ibox) = acvsq(14,ibox)-avv(11,ibox)*avv(11,ibox)
     &              -avv(12,ibox)*avv(12,ibox) 
     &              - avv(13,ibox)*avv(13,ibox)
!            momconst = 6.9994685465110493d5
               flucmom2(ibox) = 6.9994685465110493d5*flucmom2(ibox)*beta
     &              /(boxlx(ibox)*boxly(ibox)*boxlz(ibox))            
            end if
! - boxlength
            acboxl(ibox,1) = acboxl(ibox,1) / acmove
            acboxl(ibox,2) = acboxl(ibox,2) / acmove
            acboxl(ibox,3) = acboxl(ibox,3) / acmove
            
            if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
               acboxa(ibox,1) = acboxa(ibox,1) / acmove
               acboxa(ibox,2) = acboxa(ibox,2) / acmove
               acboxa(ibox,3) = acboxa(ibox,3) / acmove
            end if

            acvol(ibox) = acvol(ibox) / acmove
            acvolsq(ibox) = acvolsq(ibox) / acmove
            
            do itype = 1, nmolty
! - number of molecules
               acnbox(ibox,itype) = acnbox(ibox,itype) / acmove
! - molfraction
               molfra(ibox,itype) = molfra(ibox,itype) / acmove
! - square end-to-end length
               if ( mnbox(ibox,itype) .gt. 0 ) then 
                  asetel(ibox,itype) = 
     &                 asetel(ibox,itype) / dble(mnbox(ibox,itype))
               end if
            end do
            
            if ( lpbcz ) then
               do itype = 1, nmolty
! - number density
                  acdens(ibox,itype)=1000.0d0*acdens(ibox,itype)/acmove
               end do
! - sum over all types of molecules
               temacd = 0.0d0
               do itype = 1, nmolty
                  temacd = temacd + acdens(ibox,itype)
               end do

! - molar volume
               molvol(ibox) = 602.2045d0 / temacd
               temspd = 0.0d0
               do itype = 1, nmolty
                  temspd = temspd + 
     &                 ( acdens(ibox,itype) * masst(itype)/602.2045d0)
               end do
! - specific density
               speden(ibox) = temspd
            else
               do itype = 1, nmolty
! - number density
                  acdens(ibox,itype)=100.0d0*acdens(ibox,itype)/acmove
               end do
               temacd = 0.0d0
               do itype = 1, nmolty
                  temacd = temacd + acdens(ibox,itype)
               end do
! - molar volume
               molvol(ibox) = 100.0d0 / temacd
            end if

! - system volume- convert to average
            acvolume(ibox) = acvolume(ibox) / acmove

! - pressure and surface tension
            if ( acnp .gt. 0.5d0 ) then
               acpres(ibox) = acpres(ibox) / acnp
               acsurf(ibox) = acsurf(ibox) / acnp
            end if
            
! thermodynamic integration stuff
            if (acipsw.gt.0.5d0) acdvdl = acdvdl/acipsw

! - chemical potential
            do itype = 1, nmolty
               if (.not. lrigid(itype)) then
                  if( bnchem(ibox,itype) .gt. 0.5d0 ) then
!              --- determine how many steps it takes to grow molecule
!              --- not counting the first inserted bead
                     igrow = nugrow(itype)
                     debroglie = 17.458d0/( dsqrt(masst(itype)/beta ))
                     if (lrigid(itype)) then 
                        call schedule(igrow,itype,steps,1,0,4)
                     else
                        call schedule(igrow,itype,steps,1,0,2)
                     end if
                     acchem(ibox,itype) = ((-1.0d0)/beta) * 
     &                    dlog(acchem(ibox,itype) /
     &                    ( dble( nchoi1(itype) ) 
     &                    * dble( nchoi(itype))**steps
     &                    * dble( nchoih(itype) ) 
     &                    * bnchem(ibox,itype) 
     &                    * debroglie*debroglie*debroglie ) )
                  end if
               else
                  if( bnchem(ibox,itype) .gt. 0.5d0 ) then
!              --- determine how many steps it takes to grow molecule
!              --- not counting the first inserted bead
                     debroglie = 17.458d0/( dsqrt(masst(itype)/beta ))
                     acchem(ibox,itype) = ((-1.0d0)/beta) *
     &                    dlog(acchem(ibox,itype) /
     &                    ( dble( nchoi1(itype) )
     &                    * dble( nchoir(itype))
     &                    * dble( nchoih(itype) )
     &                    * bnchem(ibox,itype)
     &                    * debroglie*debroglie*debroglie ) )
                  end if
               end if
            end do
            if (acvsq(1,ibox).gt.0.0d0) aflv(ibox)=dsqrt(acvsq(1,ibox))
         end do

         write(iou,1215) ('       Box ',i,i=1,nbox) 
         write(iou,*)
         write(iou,1209) (acpres(i) ,i=1,nbox)
         write(iou,1212) ((acpres(i)*7.2429d-5),i=1,nbox)
         write(iou,1216) (acsurf(i) ,i=1,nbox)
         do itype = 1, nmolty
            write(iou,1210) itype, (acchem(i,itype) ,i=1,nbox)
         end do
         write(iou,*)
         
         do i = 1,3
            write(iou,1202) (acboxl(ibox,i) ,ibox=1,nbox)
         end do

         do ibox = 1, nbox
            if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
               do i = 1,3
                  write(iou,1200) acboxa(ibox,i)*180.0d0/onepi
               end do
            end if
         end do
      
         do itype = 1, nmolty
            write(iou,1201) itype, (acnbox(i,itype) ,i=1,nbox)
         end do
         if ( lpbcz ) then
            write(iou,1204) (molvol(i) ,i=1,nbox)
            write(iou,1205) (speden(i) ,i=1,nbox)
            do itype = 1, nmolty
               write(iou,1203) itype, (acdens(i,itype) ,i=1,nbox)
               if ( lexpand(itype) ) then
                  do itype2 = 1, numcoeff(itype) 
                     write(iou,1503) itype, itype2,acdens(itype,itype)*
     &                    acnbox2(itype,itype,itype2)/
     &                    (acnbox(itype,itype)*acmove),
     &                    acnbox2(itype,itype,itype2)/
     &                    (acnbox(itype,itype)*acmove)
                  end do
               end if
            end do
         else
            write(iou,1214) (molvol(i), i=1,nbox)
            do itype = 1, nmolty
               write(iou,1213) itype, (acdens(i,itype), i=1,nbox)
            end do
         end if
         do itype = 1, nmolty
            write(iou,1211) itype, (molfra(i,itype), i=1,nbox)
         end do
         do itype = 1, nmolty
            write(iou,1208) itype, (asetel(i,itype) ,i=1,nbox)
         end do
         write(iou,*)
         do j=1,10
! *** only 1 to 10 is the energy information         
            write(iou,1206) vname(j),avv(j,1:nbox),acvkjmol(j,1:nbox) 
         end do
         
         write(iou,*)
         write(iou,1207) (aflv(i) ,i=1,nbox)
         write(iou,*)
      
! ---   Output 2nd virial coefficient data
 2000    if (lvirial) then
            starviro = starvir
            dummy = dble(nstep/imv)
            do itemp = 1,ntemp
               starvir = starviro
               binvir(1,itemp) = binvir(1,itemp)/dummy
!            write(45,*) starvir,binvir(1,itemp)
               inside = starvir*starvir*binvir(1,itemp)
!            write(46,*) starvir,inside
               bvirial = 0.5d0*inside
!            write(47,*) starvir,bvirial
               starvir = starvir + stepvir
            
               do i = 2,nvirial-1
                  binvir(i,itemp) = binvir(i,itemp)/dummy
!               write(45,*) starvir,binvir(i,itemp)
                  inside = starvir*starvir*binvir(i,itemp)
!               write(46,*) starvir,inside
                  bvirial = bvirial + inside
!               write(47,*) starvir,bvirial
                  starvir = starvir + stepvir
               end do
            
               binvir(nvirial,itemp) = binvir(nvirial,itemp)/dummy
!            write(45,*) starvir,binvir(nvirial,itemp)
               inside = starvir*starvir*binvir(nvirial,itemp)
!            write(46,*) starvir,inside
!            write(47,*) starvir,bvirial
               starvir = starvir + stepvir
               bvirial = bvirial + 0.5d0*inside
          
               write(iou,*) 'At temperature of',virtemp(itemp)
               write(iou,*) 'bvirial ',
     &              -(twopi*stepvir*bvirial),' [A^3 / molecule]'
               write(iou,*) 'bvirial ',-0.602d0*twopi*
     &              stepvir*bvirial,' [cm^3 / mole]'

!            if ( lvirial2 ) then
               starvir = starviro + 0.5d0*stepvir
               do i = 2, nvirial
                  binvir2(i,itemp) = 
     &                 binvir2(i,itemp)/dummy
                  inside = starvir*starvir*binvir2(i,itemp)
                  bvirial = bvirial + inside
                  starvir = starvir + stepvir
               end do
               bvirial = -(twopi*stepvir*bvirial)
               write(iou,*) 'With quantum correction:'
               write(iou,*) 'bvirial ',bvirial,' [A^3 / molecule]'
               write(iou,*) 'bvirial ',0.602d0*bvirial,' [cm^3 / mole]'
            end do
         end if

! - solute values
         write(iou,*) 'type  box     vinter      vintra      vtor',
     &        '        vbend       vtail'
            
         do itype = 1, nmolty
            do ibox = 1, nbox
               if (solcount(ibox,itype).gt.0) then
                  write(iou,1372) itype,ibox,avsolinter(ibox,itype)
     &                 /solcount(ibox,itype),avsolintra(ibox,itype)
     &                 /solcount(ibox,itype),avsoltor(ibox,itype)
     &                 /solcount(ibox,itype),avsolbend(ibox,itype)
     &                 /solcount(ibox,itype),avsolelc(ibox,itype)
     &                 /solcount(ibox,itype)
               else
                  write(iou,1372) itype,ibox,0.0,0.0,0.0,0.0,0.0
               end if              
            end do
         end do

! --- calculate statistical errors ---
         if ( nblock .ge. 2 ) then
            dblock = dble(nblock)
            dbl1 = dblock - 1.0d0
! -      global averages -
            do i = 1,nprop
               do j = 1,nbox
                  if ( naccu(i,j) .lt. 0.5d-5 ) then
                     aver(i,j) = 0.0d0
                  else
                     aver(i,j) = accum(i,j) / naccu(i,j)
                 
                  end if
               end do
            end do
            do i = 1,nprop
               do j = 1,nbox
                  do nbl = 1, nblock
                     dsq(i,j) = dsq(i,j) + 
     &                    ( baver(i,j,nbl) - aver(i,j) )**2
                  end do
                  stdev(i,j) = dsqrt( dsq(i,j) / dblock )
                  sterr(i,j) = dsqrt( dsq(i,j) / dbl1 )
                  errme(i,j) = sterr(i,j) / dsqrt(dblock)
               end do
            end do


            do i = 1,nprop1
               do ibox = 1,nbox-1
                  do jbox = ibox+1,nbox    
                     if ( naccu1(i,ibox,jbox) .lt. 0.5d-5 ) then
                        aver1(i,ibox,jbox) = 0.0d0
                     else
                        aver1(i,ibox,jbox) = accum1(i,ibox,jbox) / 
     &                       naccu1(i,ibox,jbox)
                     end if
                  end do
               end do 
            end do
            
            do i = 1,nprop1
               do ibox = 1,nbox-1 
                  do jbox = ibox+1,nbox
                     do nbl = 1, nblock
                        dsq1(i,ibox,jbox) = dsq1(i,ibox,jbox) +
     &                       ( baver1(i,ibox,jbox,nbl)
     &                       - aver1(i,ibox,jbox) )**2
                     end do
                     stdev1(i,ibox,jbox) = dsqrt( dsq1(i,ibox,jbox) /
     &                    dblock )
                     sterr1(i,ibox,jbox) = dsqrt( dsq1(i,ibox,jbox) / 
     &                    dbl1 )
                     errme1(i,ibox,jbox) = sterr1(i,ibox,jbox) / 
     &                    dsqrt(dblock)
                  end do
               end do
            end do
 
! - write out the heat of vaporization and solubility parameters
            do ibox = 1,nbox-1
               do jbox = ibox+1,nbox 
                  write(iou,1508) ibox,jbox,aver1(1,ibox,jbox),
     &                 stdev1(1,ibox,jbox),errme1(1,ibox,jbox)
                  write(iou,1509) ibox,jbox,aver1(2,ibox,jbox),
     &                 stdev1(2,ibox,jbox),errme1(2,ibox,jbox)
                  write(iou,1510) ibox,jbox,aver1(3,ibox,jbox),
     &                 stdev1(3,ibox,jbox),errme1(3,ibox,jbox)
!               write(iou,1518) ibox,jbox,aver1(10,ibox,jbox),
!     &                 stdev1(10,ibox,jbox), errme1(10,ibox,jbox)
                  write(iou,1519) ibox,jbox,aver1(11,ibox,jbox),
     &                 stdev1(11,ibox,jbox), errme1(11,ibox,jbox) 

                  write(iou,1511) ibox,jbox,aver1(4,ibox,jbox),
     &                 stdev1(4,ibox,jbox),errme1(4,ibox,jbox)
                  write(iou,1512) ibox,jbox,aver1(5,ibox,jbox),
     &                 stdev1(5,ibox,jbox),errme1(5,ibox,jbox)
                  write(iou,1513) ibox,jbox,aver1(6,ibox,jbox),
     &                 stdev1(6,ibox,jbox),errme1(6,ibox,jbox)
                  write(iou,1514) ibox,jbox,aver1(7,ibox,jbox),
     &                 stdev1(7,ibox,jbox),errme1(7,ibox,jbox)
                  write(iou,1515) ibox,jbox,aver1(8,ibox,jbox),
     &                 stdev1(8,ibox,jbox),errme1(8,ibox,jbox)
                  write(iou,1516) ibox,jbox,aver1(9,ibox,jbox),
     &                 stdev1(9,ibox,jbox),errme1(9,ibox,jbox)
               end do
            end do
 

! - specific density
            do ibox = 1, nbox
               write(iou,1331) ibox,aver(1,ibox),stdev(1,ibox),
     &              errme(1,ibox)
            end do
            
! * system volume
            itel = 4 + nener + 4*nmolty
            do ibox = 1, nbox
               write(iou,1343) ibox,
     &              aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
            end do
            
! - pressure
            do ibox = 1, nbox
               write(iou,1341) ibox,aver(2,ibox),stdev(2,ibox),
     &              errme(2,ibox)
            end do

! - surface tension
            itel = 2+nener+ 4*nmolty+1
            do ibox = 1, nbox
               write(iou,1342) ibox,
     &              aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
            end do

            write(iou,*) 
! - energies
!         write(iou,*) 'average value', 'STD', 'SEM'
            do ibox = 1, nbox
               do j=3,2+10
! *** only 1 to 10 is the energy information
                  write(iou,1311) vname(j-2),ibox,aver(j,ibox),
     &                 stdev(j,ibox),errme(j,ibox)
               end do
            end do
      
            write(iou,*)

!-- Enthalpy
            do ibox = 1,nbox
               j = 4+nener +4*nmolty + 1
               write(iou, 1517) enth, ibox, aver(j,ibox),
     &              stdev(j,ibox),errme(j,ibox)
               j = 4+nener + 4*nmolty + 2
               write(iou,1517) enth1,ibox,aver(j,ibox),
     &              stdev(j,ibox),errme(j,ibox)
            end do 
            write(iou,*)

!--Residual Heat capacity --- 

            if(lnpt.and..not.lgibbs) then
               inst_enth= inst_enth/(dble(nchain*nstep))
               inst_enth2= inst_enth2/(dble(nchain*nstep))
               sigma2Hsimulation=((inst_enth2)- (inst_enth*inst_enth))
               sigma2H=sigma2Hsimulation*(6.022d23)*((1.38066d-23)**2) /
     &              (dble(nchain)) !(J2/mol) 
               Cp=sigma2H/((1.38066d-23)*(temp**2))
               write(iou,*) 'Cp residual(J/Kmol) =', Cp 
               write(iou,*) ' H2=', inst_enth2
               write(iou,*)  ' H=', inst_enth
            end if

            if( .not. lnpt .and..not.lgibbs) then
               inst_energy= inst_energy/(dble(nchain*nstep))
               inst_energy2= inst_energy2/(dble(nchain*nstep))
               sigma2Esimulation=((inst_energy2)- 
     &              (inst_energy*inst_energy))
               sigma2E=sigma2Esimulation*(6.022d23)*((1.38066d-23)**2) /
     &              (dble(nchain)) !(J2/mol) 
               Cv=sigma2E/((1.38066d-23)*(temp**2))
               write(iou,*) 'Cv residual(J/Kmol) =', Cv
               write(iou,*) ' E2=', inst_energy2
               write(iou,*) ' E=', inst_energy
            end if


! - chemical potential
            do itype = 1, nmolty
               itel = (2+nener) + itype
               do ibox = 1, nbox
                  if ( aver(itel,ibox) .ne. 0.0d0 ) then
                     write(iou,1351) itype,ibox,
     &                    ((-1.0d0)/beta)*dlog(aver(itel,ibox)),
     &                    (1.0d0/beta)*stdev(itel,ibox)/aver(itel,ibox),
     &                    (1.0d0/beta)*errme(itel,ibox)/aver(itel,ibox)
                  else
                     write(iou,1351) itype,ibox,0.0d0,0.0d0,0.0d0
                  end if
               end do
            end do
            
! - square end-to-end length
            do itype = 1, nmolty
               itel = (2+nener) + nmolty + itype
               do ibox = 1, nbox
                  write(iou,1321) itype,ibox,aver(itel,ibox)
     &                 ,stdev(itel,ibox),errme(itel,ibox)
               end do
            end do
            
! - number density
            do itype = 1, nmolty
               itel = (2+nener) + 2 * nmolty + itype
               do ibox = 1, nbox
                  if ( lpbcz ) then
                     write(iou,1361) itype,ibox,1.0d3*aver(itel,ibox)
     &                    ,1.0d3*stdev(itel,ibox),1.0d3*errme(itel,ibox)
                  else
                     write(iou,1361) itype,ibox,1.0d2*aver(itel,ibox)
     &                    ,1.0d2*stdev(itel,ibox),1.0d2*errme(itel,ibox)
                  end if
                  if ( lexpand(itype) .and. 
     &                 acnbox(ibox,itype) .gt. 0.5) then
                     do itype2 = 1, numcoeff(itype) 
                        molfrac = acnbox2(ibox,itype,itype2)
     &                       /(acmove*acnbox(ibox,itype))
                        write(iou,1366) itype,itype2,1.0d3*
     &                       aver(itel,ibox)*molfrac,molfrac
                     end do
                  end if
               end do
            end do
            
! - molfraction
            do itype = 1, nmolty
               itel = (2+nener) + 3 * nmolty + itype
               do ibox = 1, nbox
                  write(iou,1371) itype,ibox,aver(itel,ibox),
     &                 stdev(itel,ibox),errme(itel,ibox)
               end do
            end do
            
            if (lgibbs) then  
! --- write density results in fitting format ---
               do ibox = 1, nbox-1
                  do jbox = ibox+1,nbox
                     if (speden(ibox).lt.speden(jbox)) then
                        ig = ibox
                        il = jbox
                     else
                        ig = jbox
                        il = ibox
                     end if
!                 write(41,1401) temp,aver(1,ig),stdev(1,ig)
!     &                ,aver(1,il),stdev(1,il)
! --- write ostwald values for each moltyp
                     gconst = 8.314/(1000*beta)
                     do itype = 1,nmolty
                        itel = (2+nener) + 2 * nmolty + itype
                        ostwald = aver(itel,il)/aver(itel,ig)
                        stdost  = ostwald
     &                       * dsqrt( (stdev(itel,il)/aver(itel,il))**2
     &                       + (stdev(itel,ig)/aver(itel,ig))**2 )
!                    write(42,*) nunit(itype),ostwald,stdost
!                    write(43,*) nunit(itype),
!     &                   -(gconst*log(ostwald)) + (eta2(ig,itype) 
!     &                   - eta2(il,itype)) / 120.27167
!     &                   ,gconst*stdost/ostwald
                        write(iou,1506) itype,ig,il,ostwald,stdost
                        write(iou,1507) itype,ig,il,
     &                       -(gconst*log(ostwald)) + (eta2(ig,itype) 
     &                       - eta2(il,itype)) / 120.27167
     &                       ,gconst*stdost/ostwald
                     end do
                  end do
               end do
            end if

            write(iou,*)

! ---    write block averages  ---
            write(iou,*)
            write(iou,*) '-----block averages ------'
            do ibox=1,nbox
               write(iou,1403) ibox
               do nbl = 1, nblock
! -- changed so output the same for all ensembles
! -- 06/08/09 KM
                  write(iou,1402) nbl,baver(3,ibox,nbl),
     &                 baver(1,ibox,nbl),baver(2,ibox,nbl),
     &                 baver(3+nener+4*nmolty,ibox,nbl),
     &                 (baver(2+nener+3*nmolty+zz,ibox,nbl),zz=1,nmolty)
               end do
               if (lmipsw) then
                  write(iou,*) 'lambdais', lambdais
                  write(iou,*) 'maginn interphase switch integrand'
                  do nbl = 1, nblock
                     write(iou,*) nbl,baver(nprop,ibox,nbl)
                  end do
               end if
            end do
        
         end if

! KM 01/10 remove analysis
!      if (ianalyze.le.nstep) then
!         call analysis(2) 	
!      end if

! --- ee prob
         IF(lexpee) then                                
            write(iou,*)
            write(iou,*) 'probability of each mstate in ee'
            do nnn = 1, fmstate
               write(iou,1601) nnn,ee_prob(nnn)
            end do
         end if
         if (L_movie_xyz) then
            do ibox=1,nbox
               close (unit=210+ibox)
            end do
         end if

      end if  ! end if myid. eq. 0

! RP added for MPI
      call MPI_FINALIZE(fierr)

 1012 format(3(1x,f10.6),2i5)

 1101 format(' max trans. displacement:        ',3f10.6)
 1102 format(' max rot. displacement:          ',3f10.6)
 1103 format(' max volume displacement:        ',e12.6)
 1104 format(' dimension box 1:                ',3f10.6)
 1105 format(' dimension box 2:                ',3f10.6)

 1200 format(' box angle                                deg =',3f12.3)
 1201 format(' no. of chains of type      ',i4,'              =',3f12.3)
 1202 format(' boxlength                                [A] =',3f12.3)
 1203 format(' number density of type     ',i4,' [chain/nm^3] =',3f12.5)
 1503 format(' number density of type     ',i4,' eetype ',i4,'  =',
     &     2f12.5)
 1204 format(' molar volume                      [cm^3/mol] =',3f12.3)
 1213 format(' number density of type     ',i4,' [chain/nm^2] =',3f12.6)
 1214 format(' area per chain                     [A^2/chain] =',3f12.4)
 1205 format(' specific density                    [g/cm^3] =',3f12.6)
 1206 format(a15,'[K per system and kJ/mol per chain] =',3(f14.2,f12.2))
 1207 format(' fluctuation in <vtot>                          =',3f12.2)
 1208 format(' mean sete length of type   ',i4,'        [A^2] =',3f12.3)
 1209 format(' pressure                               [kPa] =',3f12.2)
 1210 format(' chem. potential of type    ',i4,'          [K] =',3f12.3)
 1211 format(' molfraction of type        ',i4,'              =',3f12.7)
 1212 format(' pressure                  [simulation units] =',3f12.6)
 1216 format(' surface tension                       [mN/m] =',3f12.4)
 1215 format(' Averages and fluctuations',21x,4(a11,i1))
! 1215 format(' Averages and fluctuations',21x,
!     & 3('       Box ',i1))

 1311 format(a15 ,' box ',i3, ' = ',3e14.5)
 1321 format(' mean sete length    itype ',i3,' box ',i3, ' = ',3f12.3)
 1331 format(' specific density    box ',i3, ' = ',3e12.5)
 1341 format(' pressure            box ',i3, ' = ',3f12.2)
 1342 format(' surface tension     box ',i3, ' = ',3f12.5)
 1343 format(' system volume       box ',i3, ' = ',3e12.5)
 1351 format(' chemical potential  itype ',i3,' box ',i3, ' = ',3f12.3)
 1361 format(' number density      itype ',i3,' box ',i3, ' = ',3e12.5)
 1366 format(' number density      itype ',i3,' typ ',i3, ' = ',2e12.5)
 1371 format(' mole fraction       itype ',i3,' box ',i3, ' = ',3f12.7)
 1372 format(i5,i5,3f12.5,3f12.5,3f12.5,3f12.5,3f12.5)


 1401 format(2x,f6.1,4(1x,e10.5))
 1400 format('  ------------ box: ' ,i4,/,
     & ' block   mu    msetl  density  press  '
     &       ,' dens(av/dif)  energies ... ') 
 1402 format(2x,i2,15(2x,e9.3))
 1403 format('  ------------ box: ' ,i4,/,
     & ' block   energy    density   pressure   surf ten  mol fracs') 
 1501 format(' vstart       =',3f24.10)
 1502 format(' vend         =',3f24.10)
 1504 format(' vbox         =',3f24.10)
 1505 format('Heat of Vaporization [kJ/mol] between box ',i2,' and ',i2,
     & '      ',f24.10)
 1506 format('Ostwald Coefficient  itype ',i3,' between box ',i2,
     &     ' and ',i2,f18.6,f18.6)
 1507 format('Free Enrgy of Transf itype ',i3,' between box ',i2,
     &     ' and ',i2,f18.6,f18.6,' kJ/mol')

 1508 format(' H_vap      [kJ/mol] btwn box   ',i4,' and',i4, ' =',
     & 3f15.4)
 1509 format(' H_vap LJ  [kJ/mol] btwn box   ',i4,' and',i4, ' =',
     & 3f15.4)
 1510 format(' H_vap Coul [kJ/mol] btwn box  ',i4,' and',i4, ' =',
     & 3f15.4)
 1511 format(' CED [cal/cc]   btwn box        ',i4,' and',i4, ' =',
     & 3f15.4)
 1512 format(' CED_LJ[cal/cc] btwn box        ',i4,' and',i4, ' =',
     & 3f15.4)
 1513 format(' CED_Coul[cal/cc] btwn box      ',i4,' and',i4, ' =',
     & 3f15.4)
 1514 format(' HSP [(cal/cc)^1/2]  btwn box   ',i4,' and',i4, ' =',
     & 3f15.4)
 1515 format(' HSP_LJ[(cal/cc)^1/2] btwn box  ',i4,' and',i4, ' =',
     & 3f15.4)
 1516 format(' HSP_Cou[(cal/cc)^1/2] btwn box ',i4,' and',i4, ' =',
     & 3f15.4)

! 1518 format(' DeltaU Ext [kJ/mol] btwn box   ',i4,' and',i4, ' =',
!     & 3f15.4)
 1519 format(' pdV        [kJ/mol] btwn box   ',i4,' and',i4, ' =',
     & 3f15.4)

 1517 format(a15,'[kJ/mol] for box',i3,' =',3(f12.4))
 1601 format(i5,1x,i10)


      call cleanup('')
      end



!=====================================================================72
      subroutine stopwatch_start(ctimer)
      include "mpif.h"
      character*20 ctimer
      parameter(NMAX=1000)
      real*8       starttime(NMAX), sumtime(NMAX), t0
      integer::     ntimes(NMAX,2)
      character*20 ctimers(NMAX)
      common /timelap/ starttime,sumtime,t0,ntimers,ntimes,ctimers
      data ntimers/0/

      if(ntimers .eq. 0) t0 = MPI_WTIME()

      ! Is this a pre existing timer?
      isum = -1
      do i=1,ntimers
        if(ctimer .eq. ctimers(i)) isum = i
      end do

      if(isum .gt. 0) then
        starttime(isum) = MPI_WTIME()-t0
        ntimes(isum,1)  = ntimes(isum,1) + 1
      else
        if(ntimers .ge. NMAX) return
        ntimers = ntimers + 1
        ctimers(ntimers)   = ctimer
        starttime(ntimers) = MPI_WTIME()-t0
        sumtime(ntimers)   = 0.0d0
        ntimes(ntimers,1)  = 1
        ntimes(ntimers,2)  = 0
      end if

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc MPI timing subroutines from David Porter at MSI (porter@msi.umn.edu)
!cc call start with a character::variable, then call call cleanup('') with the
!cc same character::variable
!cc at the end of monola call stopwatch_write to get the timing info
!cc this has not been personally tested
!cc KM 02/08/10
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!=====================================================================72
      subroutine stopwatch_stop(ctimer)
      include "mpif.h"
      character*20 ctimer
      parameter(NMAX=1000)
      real*8       starttime(NMAX), sumtime(NMAX), t0
      integer::     ntimes(NMAX,2)
      character*20 ctimers(NMAX)
      common /timelap/ starttime,sumtime,t0,ntimers,ntimes,ctimers
      real*8 time

      ! Is this a pre existing timer?
      isum = -1
      do i=1,ntimers
        if(ctimer .eq. ctimers(i)) isum = i
      end do
      if(isum .lt. 0) return

      time = MPI_WTIME()-t0
      sumtime(isum) = sumtime(isum) + (time-starttime(isum))
      starttime(isum) = time
        ntimes(isum,2)  = ntimes(isum,2) + 1

      return
      end

!=====================================================================72
      subroutine stopwatch_write(cfile)
      include "mpif.h"
      character*20 cfile
      parameter(NMAX=1000)
      real*8       starttime(NMAX), sumtime(NMAX), t0
      integer::     ntimes(NMAX,2)
      character*20 ctimers(NMAX)
      common /timelap/ starttime,sumtime,t0,ntimers,ntimes,ctimers

      if(ntimers .le. 0) return

      if(cfile(1:6) .eq. "stdout") then
        write (6,*) "Timer name         Total time [sec]",
     &              "    # starts     # stops"
        do i=1,ntimers
          write (6,999) ctimers(i), sumtime(i), (ntimes(i,j),j=1,2)
!99       format(a20, f16.6, 2i12)
          sumtime(i) = 0.0d0
        end do
      else
        open(unit=11, file=cfile, form="formatted")
        write (11,*) "Timer name         Total time [sec]",
     &               "    # starts     # stops"
        do i=1,ntimers
          write (11,999) ctimers(i), sumtime(i), (ntimes(i,j),j=1,2)
          sumtime(i) = 0.0d0
        end do
        close(11)
      end if

      return
      end

