      subroutine monola  

c monola
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

C -----------------------------------------------------------------
C subroutine monola
C - reads the control-data from unit 4
C - reads a starting configuration from unit 7
C - calculates interaction table
C - starts and controls the simulation
C -----------------------------------------------------------------

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
c   KEA -- spline torsion inc file
      include 'torsion.inc'
      include 'garofalini.inc'


c - variables added for GCMC histogram reweighting
      integer fmax,idum,nummol
      parameter (fmax=1e6)
      character*20 file_flt,file_hist,file_ndis(ntmax)
      character*20 file_config,file_cell
      character*20 fname2,fname3,fname4,ftemp
      character*50 fileout
      integer fname,ntii,findpos
      integer imax,itmax
      integer n,nconfig,nentry,nminp(ntmax),nmaxp(ntmax)
      integer ncmt_list(fmax,ntmax),ndist(0:nmax,ntmax)
      double precision eng_list(fmax)
      double precision vhist

      double precision Temp_Energy, Temp_Mol_Vol

      integer point_of_start, point_to_end

      integer im,mnbox,i,j,inb,nblock,ibox,jbox,nend,nnn,ii,itemp
     &  ,itype,itype2,intg,imolty,ilunit,nbl,itel,ig,il,ucheck
     &  ,jbox_max,k,histtot,Temp_nmol
      integer nvirial,zz,steps,igrow,ddum,total
      double precision starvir,stepvir,starviro
      double precision acv,acvsq,aflv,acpres,acnp,acmove,acsurf
     &  ,acboxl,acboxa,asetel,acdens,acnbox,dsq,v,vinter,vtail,vend
     &  ,vintra,vvib,vbend,vtg,vext,vstart,press1,press2,dsq1
     &  ,velect,vflucq,boxlen,acnbox2,pres(nbxmax),surf,acvolume
      double precision rm,random,temvol,setx,sety,setz,setel
     &  ,pscb1,pscb2,ratvol,avv,temacd,temspd,dblock,dbl1
     &  ,sterr,stdev,errme,qelect
      double precision bsswap,bnswap,ostwald,stdost,dummy,debroglie
     &  ,bnswap_in,bnswap_out

      double precision, dimension(nprop1,nbxmax,nbxmax)::acsolpar
 
      double precision, dimension(nbxmax)::acEnthalpy,acEnthalpy1  

      double precision binvir,binvir2,inside,bvirial
      double precision enchg1,enthchg1,srand
      double precision enchg2,enthchg2
      double precision enchg3,enthchg3
      double precision cal2joule, joule2cal
      double precision HSP_T, HSP_LJ, HSP_COUL 
      double precision CED_T, CED_LJ, CED_COUL
      double precision Heat_vapor_T,Heat_vapor_LJ,Heat_vapor_COUL    

      double precision, dimension(nprop1,nbxmax,nbxmax)::stdev1,
     &                     sterr1,errme1       
 
      dimension binvir(maxvir,maxntemp),binvir2(maxvir,maxntemp)
      dimension vstart(nbxmax),vend(nbxmax),avv(nener,nbxmax),
     +   acv(nener,nbxmax),acvsq(nener,nbxmax),aflv(nbxmax)
      dimension acboxl(nbxmax,3),acboxa(nbxmax,3),acpres(nbxmax)
      dimension acsurf(nbxmax),acvolume(nbxmax)
      dimension acnbox(nbxmax,ntmax),acnbox2(nbxmax,ntmax,20)
      dimension bsswap(ntmax,npabmax,nbxmax*2),
     &     bnswap(ntmax,npabmax,nbxmax*2),bnswap_in(ntmax,2),
     &     bnswap_out(ntmax,2)
      double precision molfra,molfrac,gconst,vdum
      dimension mnbox(nbxmax,ntmax),asetel(nbxmax,ntmax),
     &     acdens(nbxmax,ntmax),molfra(nbxmax,ntmax)
      double precision molvol(nbxmax),speden(nbxmax)
      double precision flucmom(nbxmax),flucmom2(nbxmax),flucev(nbxmax)
     &     ,flucv(nbxmax),dielect,acvol(nbxmax),acvolsq(nbxmax)
      dimension qelect(nntype)
c --- dimension statements for block averages ---
      character *15 vname(nener)
      dimension dsq(nprop,nbxmax), stdev(nprop,nbxmax),
     +          dsq1(nprop1,nbxmax,nbxmax),  
     +          sterr(nprop,nbxmax),errme(nprop,nbxmax)
      dimension ucheck(ntmax),ddum(27)

      logical ovrlap,lratio,lratv,lprint,lmv,lrsave,lblock,lucall
     &     ,lvirial2,ltfix,lratfix,ltsolute,lsolute,lpr

      dimension lratfix(ntmax),lsolute(ntmax)
      character *25 enth
      character *25 enth1

      integer bin,cnt_wf1(0:6,0:6,4),cnt_wf2(0:6,0:6,4),
     &     cnt_wra1(1000,4),cnt_wra2(1000,4)
      double precision binstep,profile(1000)

c KEA
      integer ttor
        

C -------------------------------------------------------------------

      cal2joule = 4.184d0
      joule2cal = 1.0d0/cal2joule 

      do bin = 1,1000
         profile(bin) = 0.0d0
      enddo
      do i = 0,6
         do j = 0,6
            do nnn = 1,4
               cnt_wf1(i,j,nnn) = 0
               cnt_wf2(i,j,nnn) = 0
            enddo
         enddo
      enddo
      do i = 1,1000
         do j = 1,4
            cnt_wra1(i,j) = 0
            cnt_wra2(i,j) = 0
         enddo
      enddo
      binstep = 0.05d0


      lvirial2 = .false.

      liswatch = .false.

c kea don't stop for lgaro
      if (lchgall .and. (.not. lewald).and.(.not.lgaro)) then
         write(2,*) 'lchgall is true and lewald is false.',
     &        ' Not checked for accuracy!'
         stop
      endif

      call readdat(lucall,ucheck,nvirial,starvir
     &     ,stepvir,qelect)

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
c --  SETTING UP ARRAYS  FOR ANALYSYS PURPOSE

       call analysis(0)


c --  set up initial linkcell
      if (licell) then
         call linkcell(1,0,vdum,vdum,vdum,ddum)
      endif

c --- set up thermodynamic integration stuff
      if (lmipsw) call ipswsetup
 
      if (.not.lmipsw) then
          lstagea = .false.
          lstageb = .false.
          lstagec = .false.
      endif 
      
c --- set up expanded ensemble stuff
      if (lexpee) call eesetup(qelect)

      if (lexpee.and.lmipsw) stop 'not for BOTH lexpee AND lmipsw'
      
c - use internal read/write to get integer number in character format
        write(ftemp,*) fname
        read(ftemp,*) fname2


        file_cell =
     &         'cell_param'//fname2(1:len_trim(fname2))//suffix//'.dat'


        open(unit=13,file=file_cell, status='unknown')
        close(unit=13)


	
c - setup files for histogram reweighting
      if(lgrand) then
         file_flt =
     +        'nfl'//fname2(1:len_trim(fname2))//suffix//'.dat'
         file_hist =
     +        'his'//fname2(1:len_trim(fname2))//suffix//'.dat'

         do i=1,nmolty
            write(ftemp,*) i
            read(ftemp,*) fname3
            file_ndis(i) = 'n'//fname3(1:len_trim(fname3))
     +           //'dis'//fname2(1:len_trim(fname2)) 
     +           //suffix//'.dat'
         enddo

         open(unit=50, file = file_flt, status='unknown')  
         close(unit=50)
         
         open(unit=51,file=file_hist,status='unknown')
         write(51,'(f8.4,2x,i5,2x,g15.5,3f12.3)') 
     +        temp, nmolty, (temp*log(B(i)),i=1,nmolty), boxlx(1),
     +        boxly(1), boxlz(1)         
         close(unit=51)

c --- extra zero accumulator for grand-canonical ensemble

         nentry = 0
         nconfig = 0
         do itmax = 1,ntmax
            nminp(itmax) = 1e6
            nmaxp(itmax) = -1e6
            do imax=0,nmax
               ndist(imax,itmax) = 0
            enddo
         enddo

      endif

c *** zero accumulators ***

      do i = 1,9
        do ibox = 1,nbox-1
           do jbox = ibox+1, nbox
                acsolpar(i,ibox,jbox)=0.0d0
           enddo
        enddo
      enddo  

      do i=1,nbox 
         do j=1,nener
            acv(j,i) = 0.0d0
            acvsq(j,i) = 0.0d0
	 enddo
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
            
         enddo
         acpres(i) = 0.0d0
         acsurf(i) = 0.0d0
         acEnthalpy(i) = 0.0d0
         acEnthalpy1(i) = 0.0d0
      enddo

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
            enddo
         endif
         
      enddo

      acmove = 0.0d0

c *** 2nd viral coefficient
      if (lvirial) then
         do i=1,maxvir
            do j = 1, ntemp
               binvir(i,j) = 0.0d0
               binvir2(i,j) = 0.0d0
            enddo
         enddo
c         if ( lvirial2 ) then
c            call virial2(binvir,binvir2,nvirial,starvir,stepvir)
c            goto 2000
c         endif
      endif
c *** permanent accumulators for box properties ***      
      do ibox = 1, nbox
         do i = 1,3
            acboxl(ibox,i) = 0.0d0
            acboxa(ibox,i) = 0.0d0
         enddo
      enddo
      
      do i = 1, nmolty
         do ibox = 1, nbox
            asetel(ibox,i) = 0.0d0
            mnbox(ibox,i)  = 0
            acdens(ibox,i) = 0.0d0
            molfra(ibox,i) = 0.0d0
            acnbox(ibox,i) = 0.0d0
         enddo
         if ( lexpand(i) ) then
            do itype = 1, numcoeff(i)
               do ibox = 1,2
                  acnbox2(1,i,itype) = 0.0d0
                  acnbox2(2,i,itype) = 0.0d0
               enddo
            enddo
         endif
      enddo


c *** temporary accumulators for max. displacement updates ***
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
         enddo
      enddo

c *** For the atom displacements
      Abstrax = 0.0d0
      Abstray = 0.0d0
      Abstraz = 0.0d0
      Abntrax = 0.0d0
      Abntray = 0.0d0
      Abntraz = 0.0d0


       
      do ibox = 1, nbox
         bsvol(ibox) = 0.0d0
      enddo

      do ibox = 1, nbox
         bnvol(ibox) = 0.0d0 
      enddo

      do ibox = 1,nbox
         if (lsolid(ibox) .and. .not. lrect(ibox)) then
            do j = 1,9
               bshmat(ibox,j) = 0.0d0
               bnhmat(ibox,j) = 0.0d0
            enddo
         endif
      enddo

c *** temporary accumulators for conf.bias performance ***
      do i = 1, nmolty
         do j = 1,nswapb(i)
            do ibox = 1, 2*nbox
               bsswap(i,j,ibox) = 0.0d0
               bnswap(i,j,ibox) = 0.0d0
            enddo
         enddo
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
         enddo
      enddo
c     --- accumulators for swatch performance
      do i = 1, nswaty
         do j = 1,nswtcb(i)
            bnswat(i,j) = 0.0d0
            bsswat(i,j) = 0.0d0
         enddo
      enddo
c     --- accumulators for fluctuating charge performance
      do i = 1, nmolty
         do j = 1,nbox
            bnflcq(i,j) = 0.0d0
            bnflcq2(i,j) = 0.0d0
            bsflcq(i,j) = 0.0d0
            bsflcq2(i,j) = 0.0d0
         enddo
      enddo
c --- accumulators for block averages ---
      do i = 1, nprop
         do j = 1, nbox
            naccu(i,j) = 0.0d0
            accum(i,j) = 0.0d0
            nccold(i,j) = 0.0d0
            bccold(i,j) = 0.0d0
            dsq(i,j) = 0.0d0
         enddo
      enddo
      do i = 1, nprop1
         do ibox = 1,nbox-1
            do jbox = ibox+1,nbox       
               naccu1(i,ibox,jbox) = 0.0d0
               accum1(i,ibox,jbox) = 0.0d0
               nccold1(i,ibox,jbox) = 0.0d0
               bccold1(i,ibox,jbox) = 0.0d0
               dsq1(i,ibox,jbox) = 0.0d0
            enddo
         enddo 
      enddo

      nblock = 0
 
C -----------------------------------------------------------------

      if (lneighbor) then
         do i = 1,100
            do ii = 1,nmax
               neighbor(i,ii) = 0
            enddo
         enddo
      endif
 
c *** calculate initial energy and check for overlaps ***
      do ibox=1,nbox
          call sumup( ovrlap, v, vinter,vtail, vintra,vvib,
     +                  vbend,vtg,vext,velect,vflucq, ibox, .false.)
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
ckea
          v3garob(ibox) = v3garo

          if( ovrlap ) then
             write(2,*) ' overlap in initial configuration '
             stop
          endif
	  vstart(ibox) = vbox(ibox)
          write(2,*)
          write(2,*) 'box  ',ibox,' initial v   = ', vbox(ibox)
	  if ( lneigh ) then
c ***        call for initial set-up of the near-neighbour bitmap ***
             call setnn (ibox)
	   endif
c *** calculate initial pressure ***
         call pressure ( press1, surf, ibox )
         write(2,74) ibox, surf
         write(2,64) ibox, press1
      enddo

      write(2,*)
      write(2,*) '+++++ start of markov chain +++++'
      write(2,*)
      write(2,*) 
     & 'Cycle   Total   Energy    Boxlength   Pressure  Molecules'
c     set up info at beginning of fort.12 for analysis
      write(12,*) nstep,nmolty,(masst(i),i=1,nmolty)

c *******************************************************************
c ** loops over all cycles and all molecules                       **
c *******************************************************************
 
      nend = nnstep + nstep
      if (lneighbor) then
         write(21,*) 'ii:',ii,(neigh_cnt(i),i=1,nchain)
      endif
      do 100 nnn = 1, nstep
         do 99 ii = 1, nchain 

            tmcc = nnstep + nnn - 1

c           write(2,*) 'nstep',(nnn-1)*nchain+ii
c ***       select a move-type at random ***
            rm = random()

c ###       special ensemble dependent moves ###
            if  (rm .le. pmvol) then
c           ---  volume move ---
               if ( lnpt ) then
                  call prvolume  
               else    
                  call volume
               endif

            elseif (rm .le. pmswat) then
c           --- CBMC switch move ---
               call swatch

            elseif ( rm .le. pmswap ) then
c           --- swap move for linear and branched molecules ---
               call swap(bsswap,bnswap,bnswap_in,bnswap_out,
     &              cnt_wf1,cnt_wf2,cnt_wra1,cnt_wra2,qelect)
               
            elseif ( rm .le. pmcb ) then
c           --- configurational bias move ---
               call config

            elseif ( rm .le. pmflcq ) then
c           --- displacement of fluctuating charges ---
               call flucq(2,0)

            elseif (rm .le. pmexpc ) then
c           --- expanded-ensemble move ---
               call expand

            elseif (rm .le. pmexpc1 ) then
c           --- new expanded-ensemble move ---
c               call expand
               if (random().le.eeratio) then
                  call ee_index_swap
               else
                  call eemove
               endif
            elseif ( rm .le. pm_atom_tra) then
               rm = 3.0d0 * random()
                 if ( rm .le. 1.0d0 ) then
                     call Atom_traxyz (.true.,.false.,.false.)
                 elseif ( rm .le. 2.0d0 ) then
                     call Atom_traxyz (.false.,.true.,.false.)
                 else
                     call Atom_traxyz (.false.,.false.,.true.)
                 endif
            elseif ( rm .le. pmtra ) then
c           --- translational move in x,y, or z direction ---
                 rm = 3.0d0 * random()
                 if ( rm .le. 1.0d0 ) then
                    call traxyz (.true.,.false.,.false.)	
                 elseif ( rm .le. 2.0d0 ) then
                     call traxyz (.false.,.true.,.false.)
                 else
                     call traxyz (.false.,.false.,.true.)
                 endif
            else
c           --- rotation around x,y, or z axis move --
               rm = 3.0d0 * random()
               if ( rm .le. 1.0d0 ) then
                  call rotxyz(.true.,.false.,.false.)
               elseif ( rm .le. 2.0d0 ) then
                  call rotxyz(.false.,.true.,.false.)
               else
                  call rotxyz(.false.,.false.,.true.)
               endif
            endif
            
            acmove = acmove + 1.0d0
c ***       accumulate probability of being in an expanded ensemble
c ***       state

            if (lexpee) then
                ee_prob(mstate) = ee_prob(mstate)+1
            endif

c ***       calculate instantaneous values ***
c ***       accumulate averages ***
            
	    do ibox=1,nbox
               do itype = 1, nmolty
                  acnbox(ibox,itype) = acnbox(ibox,itype) + 
     +                 dble(ncmt(ibox,itype))
                  if ( lexpand(itype) ) then
                     do itype2 = 1, numcoeff(itype)
                        acnbox2(ibox,itype,itype2) = 
     &                       acnbox2(ibox,itype,itype2) +
     &                       dble(ncmt2(ibox,itype,itype2))
c                        write(2,*) '1:',acnbox2(ibox,itype,itype2)
                     enddo
c                     write(2,*) '2:', acnbox(ibox,itype)
                  endif
               enddo

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
c KEA added 17 for v3garo
               acv(17,ibox)    = acv(17,ibox) + v3garob(ibox)
               acvsq(17,ibox)  = acvsq(17,ibox) + v3garob(ibox)**2

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
               endif

               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  temvol = cell_vol(ibox) 
               else
                  if ( lpbcz ) then
                     temvol = boxlx(ibox)*boxly(ibox)*boxlz(ibox) 
                  else
                     temvol = boxlx(ibox)*boxly(ibox)
                  endif
               endif

               if ( lnpt ) then
c                  acv(16,ibox) = acv(16,ibox) + vbox(ibox)*boxlx(ibox)
c     &                 *boxly(ibox)*boxlz(ibox)
                  acv(16,ibox) = acv(16,ibox) + vbox(ibox)*temvol
               endif


               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  boxlx(ibox) = cell_length(ibox,1) 
                  boxly(ibox) = cell_length(ibox,2)
                  boxlz(ibox) = cell_length(ibox,3)

                  acboxa(ibox,1) = acboxa(ibox,1) + cell_ang(ibox,1)
                  acboxa(ibox,2) = acboxa(ibox,2) + cell_ang(ibox,2)
                  acboxa(ibox,3) = acboxa(ibox,3) + cell_ang(ibox,3)

               endif

               acboxl(ibox,1) = acboxl(ibox,1) + boxlx(ibox)
               acboxl(ibox,2) = acboxl(ibox,2) + boxly(ibox)
               acboxl(ibox,3) = acboxl(ibox,3) + boxlz(ibox)

               acvol(ibox) = acvol(ibox) + temvol
               acvolsq(ibox) = acvolsq(ibox) + temvol*temvol
               acvolume(ibox) = acvolume(ibox) + temvol

               do itype = 1, nmolty
                  acdens(ibox,itype) = acdens(ibox,itype) +
     +                 ncmt(ibox,itype) / temvol
                  if ( nchbox(ibox) .gt. 0 ) then
                     molfra(ibox,itype) = molfra(ibox,itype) +
     +                    dble(ncmt(ibox,itype)) / dble(nchbox(ibox))
                  endif
               enddo
	    enddo


            if (lstop) then
	       if ( int(acmove) .ge. nstep ) goto 101
	    endif

c *************************************************************
c ** ends loop over chains                                   **
c *************************************************************

c - collect histogram data (added 8/30/99)
            if(lgrand) then
               ibox = 1
               vhist = vinterb(ibox) +
     +              velectb(ibox) + vflucqb(ibox)	      
               nconfig = nconfig + 1
               if (mod(nconfig,ninstf).eq.0) then 

	         open(unit=50, file = file_flt, status='old', 
     +					position='append')
     
c	         write(50, '(i10,5x,i7,5x,g15.6)') nconfig, 
c     +			  (ncmt(ibox,i), i=1,nmolty),vhist

c Formatting removed until I can figure out how to get <n>i7
c to work
	         write(50, * ) nconfig, 
     +			  (ncmt(ibox,i), i=1,nmolty),vhist
	         close(unit=50)
	      endif
	   	    
	      if(mod(nconfig,ninsth).eq.0.and.nconfig.gt.nequil) then

	      do imolty = 1, nmolty
	         nminp(imolty) = min(nminp(imolty),ncmt(ibox,imolty))
		 nmaxp(imolty) = max(nmaxp(imolty),ncmt(ibox,imolty))
	      enddo
                 nentry = nentry + 1

	         do imolty=1,nmolty
                    ncmt_list(nentry,imolty) = ncmt(ibox,imolty)
		    ndist(ncmt(ibox,imolty),imolty) = 
     +	   	    ndist(ncmt(ibox,imolty),imolty) + 1
		 enddo

		 eng_list(nentry) = vhist
			 			
	      endif

	      if (mod(nconfig,ndumph).eq.0) then
	         open(unit=51,file = file_hist,status='old',
     +	    	      position='append') 
                 do i=1,nentry 
		    write(51, * ) 
     +		      (ncmt_list(i,imolty), imolty=1,nmolty), 
     +					eng_list(i)     
		 enddo
		 close(unit=51)
		 nentry = 0

		 do imolty=1,nmolty
		   open(unit=52, file=file_ndis(imolty), status='unknown')
		      do n=nminp(imolty),nmaxp(imolty)
		         write(52,*) n,ndist(n,imolty)
		      enddo
		      close(unit=52)
		 enddo
	      endif
	 endif

99       continue

c *** perform periodic operations  ***
 
         if ( lgibbs .or. lnpt .and. (.not. lvirial) ) then
            do ibox = 1,nbox
               if ( lpbcz ) then
                  if (lsolid(ibox) .and. .not. lrect(ibox)) then
                     write(12,'(7e13.5,15i5)') hmat(ibox,1)
     &                    ,hmat(ibox,4),hmat(ibox,5)
     &                    ,hmat(ibox,7),hmat(ibox,8)
     &                    ,hmat(ibox,9),vbox(ibox),
     &                    (ncmt(ibox,itype),itype=1,nmolty)

                 open(unit=13,file = file_cell,status='old',
     +                position='append')
                     write(13,'(i8,6f12.4)') nnn+nnstep,
     +               cell_length(ibox,1)/Num_cell_a,
     +               cell_length(ibox,2)/Num_cell_b,
     +               cell_length(ibox,3)/Num_cell_c,
     +               cell_ang(ibox,1)*180.0d0/onepi,
     +               cell_ang(ibox,2)*180.0d0/onepi,
     +               cell_ang(ibox,3)*180.0d0/onepi
                 close(unit=13)

c                     write(13,'(i8,3f12.4)') nnn,cell_ang(ibox,1)
c     +                        , cell_ang(ibox,2),cell_ang(ibox,3)
                  else
c                  do ibox = 1, nbox
                     write(12,'(4e13.5,15i5)') boxlx(ibox),boxly(ibox)
     +                    ,boxlz(ibox),vbox(ibox),
     +                    (ncmt(ibox,itype),itype=1,nmolty)
c                  enddo
                  endif
               else
c               do ibox = 1, nbox
                  write(12,'(2e12.5,15i4)') boxlx(ibox)*boxly(ibox)
     +                 ,vbox(ibox),(ncmt(ibox,itype),itype=1,nmolty)
c               enddo
               endif
            enddo
         endif

         if (lucall) then
            stop 'not recently checked for accuracy'
c            do j = 1,nmolty
c               if ( ucheck(j) .gt. 0 ) then
c                  call chempt(bsswap,j,ucheck(j),qelect)
c               endif
c            enddo
         endif

         if ( mod(nnn,iratp) .eq. 0 ) then
c *** calculate pressure ***
            acnp = acnp + 1.0d0
            do ibox = 1, nbox
               call pressure ( press1, surf, ibox )
c              write(2,*) 'control pressure', press1
               pres(ibox) = press1
               acpres(ibox) = acpres(ibox) + press1
               acsurf(ibox) = acsurf(ibox) + surf
            enddo

c Enthalpy calculation

            do ibox = 1,nbox  
               Temp_nmol = 0
               do itype=1,nmolty
                  Temp_nmol =   Temp_nmol + ncmt(ibox,itype)
               enddo
               Temp_Mol_Vol = temvol/Temp_nmol*0.6022d-06 ! m3/mol
               Temp_Energy  = (vbox(ibox)/Temp_nmol)*0.00831451d0 ! kJ/mol
               acEnthalpy(ibox) = acEnthalpy(ibox) + Temp_Energy +
     +                       pres(ibox)*Temp_Mol_Vol     !kJ/mol
               acEnthalpy1(ibox) = acEnthalpy1(ibox) + Temp_Energy +
     +                       (express(ibox)/7.2429d-5)*Temp_Mol_Vol     !kJ/mol   
            enddo

c --- cannot calculate a heat of vaporization for only one box,
c --- and some compilers choke because Heat_vapor_T will not be
c --- defined if nbox == 1
            if(lgibbs) then
               do ibox = 1,nbox-1
                  do jbox = ibox+1,nbox
                       call calcsolpar(pres,Heat_vapor_T,Heat_vapor_LJ,
     &                      Heat_vapor_COUL,CED_T,CED_LJ,CED_COUL,
     &                      HSP_T,HSP_LJ,
     &                      HSP_COUL,ibox,jbox)
c --- Heat of vaporization            
                            acsolpar(1,ibox,jbox)=
     &                            acsolpar(1,ibox,jbox)+Heat_vapor_T
                            acsolpar(2,ibox,jbox)=
     &                            acsolpar(2,ibox,jbox)+Heat_vapor_LJ
                            acsolpar(3,ibox,jbox)=
     &                            acsolpar(3,ibox,jbox)+Heat_vapor_COUL
                            acsolpar(4,ibox,jbox)=
     &                            acsolpar(4,ibox,jbox)+CED_T
                            acsolpar(5,ibox,jbox)=
     &                            acsolpar(5,ibox,jbox)+CED_LJ
                            acsolpar(6,ibox,jbox)=
     &                            acsolpar(6,ibox,jbox)+CED_COUL
                            acsolpar(7,ibox,jbox)=
     &                            acsolpar(7,ibox,jbox)+HSP_T
                            acsolpar(8,ibox,jbox)=
     &                            acsolpar(8,ibox,jbox)+HSP_LJ
                            acsolpar(9,ibox,jbox)=
     &                            acsolpar(9,ibox,jbox)+HSP_COUL 
                  enddo  
               enddo
            endif
         endif

c calculate the integrand of thermosynamic integration
         if (lmipsw.and.(mod(nnn,iratipsw).eq.0)) then
            acipsw = acipsw+1.0d0
            call deriv(1)
            acdvdl = acdvdl+dvdl
         endif

 64      format(' pressure check:   box',i2,' =',f14.2)
 74      format(' surf. tension :   box',i2,' =',f14.5)


c *** Add a call for subroutine to compute the 

         lratio = .false.
         lratv = .false.
         lprint = .false.
         lmv = .false.
         lrsave = .false.
         lblock = .false.

         if ( mod(nnn,iratio) .eq. 0 ) then
            lratio = .true.
         endif

         if ( lgibbs .or. lnpt ) then
            if ( mod(nnn,iratv) .eq. 0 .and. pmvol .gt. 0.0d0 ) then
               lratv = .true.
            endif
         endif

         if ( mod(nnn,iprint) .eq. 0 ) then
            lprint = .true.
         endif
 
         if ( mod(nnn,imv) .eq. 0 ) then
            if ( lvirial ) then
               call virial(binvir,binvir2,nvirial,starvir,stepvir)
            else
               lmv = .true.
            endif
         endif
         if ( mod(nnn,irsave) .eq. 0 ) then
            lrsave = .true.
         endif
 
	 if(mod(nnn,ianalyze).eq.0) then
	    call analysis(1)
         endif  

         do intg = 1, nchain
            ibox = nboxi(intg)
            imolty = moltyp(intg)
c *** accumulate m-n-box and m-s-e-t-e-l ***
c     only count the main chain - not the hydrogens
            ilunit = nugrow(imolty)
            setx = rxu(intg,1) - rxu(intg,ilunit)
            sety = ryu(intg,1) - ryu(intg,ilunit)
            setz = rzu(intg,1) - rzu(intg,ilunit)
            setel = setx*setx + sety*sety + setz*setz
c            if ( imolty .eq. 2 ) then
c               write(??,*) imolty,setel
c            endif
            mnbox( ibox, imolty ) = mnbox( ibox, imolty ) + 1
            asetel( ibox, imolty ) = asetel( ibox, imolty ) + setel
         enddo

         if ( mod(nnn,iblock) .eq. 0 ) then
            lblock = .true.
         endif
         
         ltsolute = .false.
         ltfix    = .false.
         
         do imolty = 1, nmolty
            if (pmfix(imolty).gt.0.0001d0.and
     &           .mod(nnn,iupdatefix).eq.0) then
               ltfix = .true.
               lratfix(imolty) = .true.
            else
               lratfix(imolty) = .false.
            endif

            if (mod(nnn,isolute(imolty)).eq.0) then
               ltsolute = .true.
               lsolute(imolty) = .true.
            else
               lsolute(imolty) = .false.
            endif
         enddo

         if (lratio .or. lratv .or. lprint .or. lmv .or. lrsave
     &        .or. lblock .or. ltfix .or. ltsolute) then
            call monper(acv,acpres,acsurf,acvolume,molfra,mnbox,asetel
     &       ,acdens,acmove,acnp,pres,nbox,nnn,nblock,lratio,lratv
     &       ,lprint,lmv,lrsave,lblock,lratfix,lsolute,acsolpar,
     &        acEnthalpy,acEnthalpy1)
         endif
         if ( ldielect .and. (mod(nnn,idiele) .eq. 0)) then
            dielect = acvsq(14,ibox)/acmove

! *If you really want this quantity comment should be taken out**

!            write(14,*) nnn,6.9994685465110493d5*dielect*beta/
!     &           (boxlx(ibox)*boxly(ibox)*boxlz(ibox))

            dielect = acvsq(14,ibox)/acmove  
     &           -(acv(11,ibox)/acmove)**2
     &           -(acv(12,ibox)/acmove)**2 - (acv(13,ibox)/acmove)**2
            write(15,*) nnn,6.9994685465110493d5*dielect*beta/
     &           (boxlx(ibox)*boxly(ibox)*boxlz(ibox))
            write(16,*) nnn,acv(11,ibox)/acmove, acv(12,ibox)/acmove,
     &                acv(13,ibox)/acmove
         endif
         if ( lnpt .and. nmolty .eq. 1 ) then
c *** output the fluctuation information
            if (mod(nnn,idiele) .eq. 0) then
               write(14,*) nnn,acvol(ibox)/acmove
               write(15,*) nnn,acvolsq(ibox)/acmove
               write(16,*) nnn,acv(1,ibox)/acmove
               write(17,*) nnn,acvsq(1,ibox)/acmove
               write(18,*) nnn,acvol(ibox)*acv(1,ibox)/(acmove*acmove)
               write(19,*) nnn,acv(16,ibox)/acmove - acvol(ibox)
     &              *acv(1,ibox)/(acmove*acmove)
            endif
         endif
         
         if ( ldielect .and. mod(nnn,idiele) .eq. 0 ) then
            do ibox = 1,nbox
                write(17,'(2I8,3f16.8)') nnn,ibox,dipolex(ibox),
     &        dipoley(ibox),dipolez(ibox)
            enddo
         endif

c         if ( mod(nnn,idiele) .eq. 0 ) write(25,*) nnn+nnstep,vbox(1)

c         ibox = 1
c         imolty = 1
c         do i = 1,nchain
c            if ( nboxi(i) .eq. ibox ) then
c               if ( moltyp(i) .eq. imolty ) then
c                  bin = dint(zcm(i)/binstep) + 1
c                  temvol = boxlx(ibox)*boxly(ibox)*binstep
c                  profile(bin) = profile(bin)+1.0d0/temvol
c               endif
c            endif
c         enddo

100   continue

      lpr = .false.
      do i = 1,ntmax
         if (lbias(i)) then
            lpr = .true.
         endif
      enddo

      if (lpr) then
         do nnn = 1,4
            write(31,*)
            write(31,*) 'nnn:',nnn
            do i = 0,6
               do j = 0,6
                  if (cnt_wf1(i,j,nnn) .gt. 0 ) then
                     write(31,*) i,j,cnt_wf1(i,j,nnn),cnt_wf2(i,j,nnn)
                  endif
               enddo
            enddo
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
            enddo
         enddo
      endif

c      do bin = 1,1000
c         write(26,*) binstep*(dble(bin)-0.5d0),profile(bin)/nstep
c      enddo

c *******************************************************************
c ** ends the loop over cycles                                     **
c *******************************************************************
  101 continue 

 
      if (lneighbor) then
         write(21,*) 'ii:',ii,(neigh_cnt(i),i=1,nchain)
      endif
      write(2,*)
      write(2,*) '+++++ end of markov chain +++++'
 
c *** write some information about config performance ***
      if ( pmcb .gt. 0.0d0 ) then
         write(2,*)
         write(2,*) '### Configurational-bias ###'
         write(2,*)
         do i = 1, nmolty
            write(2,*) 'molecule typ =',i
            write(2,*) '    length  attempts  succ.growth  accepted   ',
     +           '%su.gr.    %accep.'
            do inb = 1, nunit(i)
               if ( bncb(i,inb) .gt. 0.0d0 ) then
                  pscb1 = bscb(i,1,inb) * 100.0d0 / bncb(i,inb)
                  pscb2 = bscb(i,2,inb) * 100.0d0 / bncb(i,inb)
                  write(2,'(i9,3f10.1,2f10.2)') inb, bncb(i,inb),
     +                 bscb(i,1,inb), bscb(i,2,inb), pscb1, pscb2
               endif
            enddo
            if (pmfix(i).gt.0.0d0) then
               write(2,*) ' SAFE-CBMC move '
               write(2,*) '    length  attempts  succ.growth  ',
     +              'accepted   %su.gr.    %accep.'               
               do inb = 1, nunit(i)
                  if (fbncb(i,inb) .gt. 0.0d0 ) then
                     pscb1 = fbscb(i,1,inb) * 100.0d0 
     +                    / fbncb(i,inb)
                     pscb2 = fbscb(i,2,inb) * 100.0d0 
     +                    / fbncb(i,inb)
                     write(2,'(i9,3f10.1,2f10.2)') inb, fbncb(i,inb),
     +                    fbscb(i,1,inb), fbscb(i,2,inb)
     +                    , pscb1, pscb2
                  endif
               enddo
            endif
         enddo
         write(2,*)
      endif
c *** write some information about volume performance ***
      if ( lgibbs .or. lnpt) then
         write(2,*)
         write(2,*) '### Volume change       ###'
         do ibox = 1,nbox
            if (lsolid(ibox) .and. .not. lrect(ibox)) then
               do j = 1,9
                  acnhmat(ibox,j) = acnhmat(ibox,j) + bnhmat(ibox,j)
                  acshmat(ibox,j) = acshmat(ibox,j) + bshmat(ibox,j)
                  if ( acshmat(ibox,j) .gt. 0.5d0) then
                     write(2,70) acnhmat(ibox,j),
     &                    acshmat(ibox,j)/acnhmat(ibox,j),rmhmat(ibox,j)
                  else
                     write(2,70) acnhmat(ibox,j),0.0d0,rmhmat(ibox,j)
                  endif
               enddo
            else
               acnvol(ibox) = acnvol(ibox) + bnvol(ibox)
               acsvol(ibox) = acsvol(ibox) + bsvol(ibox)
               if ( acnvol(ibox) .ne. 0.0d0 ) then
                  ratvol = acsvol(ibox) / acnvol(ibox)
               else
                  ratvol = 0.0d0
               endif
               write(2,61) acnvol(ibox),ratvol,rmvol(ibox)
            endif
         enddo
      endif
      write(2,*)  
      write(2,*) '### Molecule swap       ###'
      write(2,*)
      do i = 1, nmolty
         write(2,*) 'molecule typ =',i
         do j=1,nswapb(i)
            if ( box1(i,j) .eq. box2(i,j) ) then
               jbox_max = 1
            else
               jbox_max = 2
            endif
            do jbox = 1,jbox_max
               if ( jbox .eq. 1 ) ibox = box1(i,j)
               if ( jbox .eq. 2 ) ibox = box2(i,j)
               write(2,66) box1(i,j),box2(i,j),ibox,
     &              bsswap(i,j,ibox),bnswap(i,j,ibox),
     &              bnswap(i,j,ibox+nbox)
               if (bnswap(i,j,ibox) .gt. 0.5d0) then
                  bsswap(i,j,ibox) = bsswap(i,j,ibox+nbox)*
     &                 100.0d0/bsswap(i,j,ibox)
                  bnswap(i,j,ibox) = bnswap(i,j,ibox+nbox)*
     &                 100.0d0/bnswap(i,j,ibox)
                  write(2,63) bsswap(i,j,ibox),bnswap(i,j,ibox)
               endif
            enddo
         enddo
         write(2,68) bnswap_in(i,1), bnswap_in(i,2)
         write(2,69) bnswap_out(i,1), bnswap_out(i,2)
      enddo

      write(2,*)
      write(2,*) '### Molecule swatch     ###'
      write(2,*)
      do i = 1, nswaty
         write(2,*) 'pair typ =',i
         write(2,*) 'moltyps = ',nswatb(i,1),' and',nswatb(i,2)
         do j = 1, nswtcb(i)
            write(2,62) box3(i,j),box4(i,j),
     &           bnswat(i,j),bsswat(i,j)
            if (bnswat(i,j) .gt. 0.5d0 ) then
               write(2,65) 100.0d0 * bsswat(i,j)/bnswat(i,j)
            endif
         enddo
      enddo

      write(2,*)
      write(2,*)    '### Charge Fluctuation  ###'
      write(2,*)
      
      do i = 1, nmolty
         do j = 1,nbox
            bnflcq2(i,j) = bnflcq2(i,j) + bnflcq(i,j) 
            bsflcq2(i,j) = bsflcq2(i,j) + bsflcq(i,j) 
            if (bnflcq2(i,j) .gt. 0.5d0) then
               write(2,*) 'molecule typ =',i,'  box =',j
               bsflcq2(i,j) = bsflcq2(i,j)/bnflcq2(i,j)
               write(2,61) bnflcq2(i,j),bsflcq2(i,j),rmflcq(i,j)
            endif
         enddo
      enddo

      write(2,*) 
      write(2,*)    '### Expanded Ensemble Move  ###'
      write(2,*) 
      do i = 1, nmolty
         do j = 1,nbox
            if (lexpand(i) .and. bnexpc(i,j) .gt. 0.5) then
               write(2,*) 'molecule typ =',i,'  box =',j
               write(2,67) bnexpc(i,j),bsexpc(i,j),
     &              bsexpc(i,j)/bnexpc(i,j)
            endif
         enddo
      enddo

 61   format(' attempts =',f8.1,'   ratio =',f6.3,
     +      '   max.displ. =',e10.4)
 62   format('between box ',i2,'and ',i2,
     +     ' attempts =',f9.1,'   accepted =',f8.1) 
 63   format(' suc.growth % =',f7.3,'   accepted % =',f7.3)
 65   format(' accepted % =',f7.3)
 66   format('between box ',i2,' and ',i2,' into box',i2,
     +     '   uattempts =',f12.1,' attempts =',f9.1
     &     ,'   accepted =',f8.1)
 67   format(' attempts =',f8.1,'   accepted =',f8.1,
     &     ' accepted % =',f7.3)

 68   format('number of times move in: ', f12.1, 
     &     '  accepted=',f8.1)
 69   format('number of times move out: ', f12.1, 
     &     '  accepted=',f8.1)
 70   format(' h-matrix attempts =',f8.1,'   ratio =',f6.3,
     +      '   max.displ. =',e10.4)
      if (lexzeo) then
c        --- write end-conf to picture file
         call writepdb(ncmt(1,1),nunit(1),1)
      endif
c *** checks final value of the potential energy is consistent ***
      
      do ibox=1,nbox
         
         if ( ldielect ) then
c *** store old dipole moment
            call dipole(ibox,2)
         endif
         call sumup( ovrlap, v, vinter,vtail,vintra,vvib,
     +                  vbend,vtg,vext,velect,vflucq,ibox, .false.)
         vend(ibox) = v

c---need to check

         if ( abs(v - vbox(ibox)) .gt. 0.0001) then
            write(2,*) '### problem with energy ###  box ',ibox
	    write(2,*) ' Total energy: ',v,vbox(ibox),v-vbox(ibox)
	 endif
         if ( abs(vinter - vinterb(ibox)) .gt. 0.000001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Inter mol.en.: ',vinter,vinterb(ibox)
             if (lsolid(ibox) .and. .not. lrect(ibox)) then
                write(2,*) 'You might check the cutoff wrt box widths'
                write(2,*) 'Normal PBC might be failing'
             endif
	 endif
         if ( abs(vtail - vtailb(ibox)) .gt. 0.000001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Tail corr.en.: ',vtail,vtailb(ibox)
	 endif
         if ( abs(vintra - vintrab(ibox)) .gt. 0.000001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Intra mol.en.: ',vintra,vintrab(ibox)
	 endif
         if ( abs(vvib - vvibb(ibox)) .gt. 0.001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' bond vib. en.: ',vvib,vvibb(ibox)
	 endif
         if ( abs(vbend - vbendb(ibox)) .gt. 0.001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Bond ben.en.: ',vbend,vbendb(ibox)
	 endif
         if ( abs(vtg - vtgb(ibox)) .gt. 0.001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Torsion.en.: ',vtg,vtgb(ibox)
	 endif
         if ( abs(vext - vextb(ibox)) .gt. 0.0001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Externa.en.: ',vext,vextb(ibox)
	 endif
         if ( abs(velect - velectb(ibox)) .gt. 0.000001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Coulomb.en.: ',velect,velectb(ibox)
	 endif
         if ( abs(vflucq - vflucqb(ibox)) .gt. 0.0001) then
             write(2,*) '### problem  ###'
	     write(2,*) ' Fluc Q en.: ',vflucq,vflucqb(ibox)
	 endif
         if ( abs(v3garo - v3garob(ibox) ) .gt.0.001) then
            write(6,*) '### problem ###'
            write(6,*) ' 3-body en.: ',v3garo,v3garob(ibox)
         endif
         if ( ldielect ) then
            if ( abs(dipolexo - dipolex(ibox)) .gt. 0.0001) then
               write(2,*) '### problem  ###'
               write(2,*) ' Dipole X: ',dipolexo,dipolex(ibox)
            endif
         endif
         if (lmipsw) then
            if (abs(vwellipsw-vwellipswb(ibox)).gt.0.001) then
                write(2,*) '### problem  ###'
                write(2,*) ' well en.: ',vwellipsw,vwellipswb(ibox)
            endif
         endif
      enddo
      write(2,*)
      write(2,1501) (vstart(i) ,i=1,nbox)
      write(2,1502) (vend(i)   ,i=1,nbox)
      write(2,1504) (vbox(i)   ,i=1,nbox)
      write(2,*)

c     ** normalize and write out presim results in fort.22 **

      if (lpresim) then
         if (counttot.eq.0) then
            write(21,*) counthist
         else
            write(21,*) counttot
         endif

         do j = 1, iring(1) 
            do k = 1, iring(1)
               if (j.eq.k) goto 150
               histtot = 0
               do bin = 1, maxbin
                  hist(j,k,bin) = hist(j,k,bin) + 1.0d0
                  histtot = histtot + hist(j,k,bin)
               enddo
               
               do bin = 1, maxbin
                  hist(j,k,bin) = hist(j,k,bin) / histtot
                  write(21,*) bin,hist(j,k,bin)
               enddo
 150           continue
            enddo
         enddo
      endif

c     *** put new distribution back into a file
      do imolty = 1, nmolty
         if (pmfix(imolty).gt.0) then
            if (counttot.eq.0) then
               write(21,*) counthist
            else
               write(21,*) counttot
            endif
            do j = 1, iring(1) 
               do k = 1, iring(1)
                  if (j.eq.k) goto 160
                  do bin = 1, maxbin
                     write(21,*) bin,probf(j,k,bin)
                  enddo
 160              continue
               enddo
            enddo
         endif
      enddo
 

c ** write out the final configuration for each box, Added by Neeraj 06/26/2006 3M ***

      do ibox = 1,nbox
         write(ftemp,*) ibox
         read(ftemp,*) fname4
         fileout = 'box'//fname4(1:len_trim(fname4))//'config'//
     &               fname2(1:len_trim(fname2))
     &                         //suffix//'.xyz'
         open (unit=200+ibox,FILE=fileout,status="unknown")
         
         nummol = 0
         do i = 1,nchain
            if (nboxi(i).eq.ibox) then
               nummol = nummol + nunit(moltyp(i))
            endif
         enddo 
         write(200+ibox,*) nummol 
         write(200+ibox,*) 
         do i = 1,nchain
          if(nboxi(i).eq.ibox) then 
            imolty = moltyp(i) 
            do ii = 1,nunit(imolty)
               ntii = ntype(imolty,ii)               
               write(200+ibox,'(a4,5x,3f15.4)') chemid(ntii),rxu(i,ii),
     &                  ryu(i,ii), rzu(i,ii)
            enddo
          endif
         enddo   
         close(200+ibox)
      enddo  

      
c ** write out the final configuration from the run ***



        file_config = 'config'//fname2(1:len_trim(fname2))
     +                        //suffix//'.dat'
        open(8, file=file_config,status='unknown')
      
      write(8,*) nend
      if ( nend .gt. 0 ) then
         write(8,*) Armtrax, Armtray, Armtraz 
         do im=1,nbox
            do imolty=1,nmolty
               write(8,*) rmtrax(imolty,im), rmtray(imolty,im)
     &              , rmtraz(imolty,im)
               write(8,*) rmrotx(imolty,im), rmroty(imolty,im)
     &              , rmrotz(imolty,im)
            enddo
         enddo
         do im=1, nbox
            write (8,*) (rmflcq(i,im),i=1,nmolty)
         enddo
         if ( lgibbs .or. lgrand .or. lnpt ) then
            write(8,*) (rmvol(ibox),ibox=1,nbox)
            do ibox = 1,nbox
               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  write(8,*) (rmhmat(ibox,i),i=1,9)
                  write(8,*) (hmat(ibox,i),i=1,9)
               else
                  write(8,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
               endif
            enddo
         endif
      endif
  
      if (L_add) then
         do i = nchain+1, nchain+N_add    
            moltyp(i)=N_moltyp2add
            nboxi(i) = N_box2add
            rxu(i,1) = random()*boxlx(N_box2add)
            ryu(i,1) = random()*boxly(N_box2add)
            rzu(i,1) = random()*boxlz(N_box2add)
            qqu(i,1) = 0.0
         enddo
         nchain = nchain+N_add
         write(8,*) nchain
         write(8,*) nmolty
         write(8,*) (nunit(i),i=1,nmolty)
         write(8,*) (moltyp(i),i=1,nchain)
         write(8,*) (nboxi(i),i=1,nchain)
         do i = 1, nmolty
           if ( lexpand(i) ) write(8,*) eetype(i)
         enddo
         do i = 1, nmolty
           if ( lexpand(i) ) write(8,*) rmexpc(i)
         enddo
         do  i = 1, nchain
           imolty = moltyp(i)
           do  j = 1, nunit(imolty)
             write(8,'(4f15.6)') rxu(i,j), ryu(i,j), rzu(i,j), qqu(i,j)
           enddo
         enddo 

      elseif(L_sub) then
         point_of_start = 0
         do i =1,N_moltyp2sub
            point_of_start=point_of_start+temtyp(i)
         enddo
         point_of_start = point_of_start-N_sub+1

         point_to_end = nchain-N_sub
            
         do i = point_of_start,point_to_end
            nboxi(i) = nboxi(i+N_sub)
            moltyp(i) = moltyp(i+N_sub) 
         enddo

         write(8,*) nchain-N_sub
         write(8,*) nmolty
         write(8,*) (nunit(i),i=1,nmolty)
         write(8,*) (moltyp(i),i=1,(nchain-N_sub))
         write(8,*) (nboxi(i),i=1,(nchain-N_sub))
         do i = 1, nmolty
           if ( lexpand(i) ) write(8,*) eetype(i)
         enddo
         do i = 1, nmolty
           if ( lexpand(i) ) write(8,*) rmexpc(i)
         enddo
         do  i = 1, nchain
          if(i.lt.(point_of_start).or.i.gt.
     &              (point_of_start+N_sub-1)) then
            imolty = moltyp(i)
            do  j = 1, nunit(imolty)
              write(8,'(4f15.6)') rxu(i,j), ryu(i,j), rzu(i,j), qqu(i,j)
            enddo
          endif
         enddo
      else      
        write(8,*) nchain
        write(8,*) nmolty
        write(8,*) (nunit(i),i=1,nmolty)
        write(8,*) (moltyp(i),i=1,nchain)
        write(8,*) (nboxi(i),i=1,nchain)
        do i = 1, nmolty
           if ( lexpand(i) ) write(8,*) eetype(i)
        enddo
        do i = 1, nmolty
           if ( lexpand(i) ) write(8,*) rmexpc(i)
        enddo
        do i = 1, nchain
           imolty = moltyp(i)
           do j = 1, nunit(imolty)
             write(8,'(4f15.6)') rxu(i,j), ryu(i,j), rzu(i,j), qqu(i,j)
           enddo 
         enddo
      endif

      close(8)



c *** calculate and write out running averages ***
      do ibox=1,nbox
c - energies
         do j=1,nener
            avv(j,ibox)   = acv(j,ibox) / acmove
            acvsq(j,ibox) = (acvsq(j,ibox)/acmove) - avv(j,ibox) ** 2
	 enddo
         if ( ldielect ) then
            flucmom(ibox) = acvsq(14,ibox)-avv(15,ibox)*avv(15,ibox)
c            momconst = 6.9994685465110493d5
            flucmom(ibox) = 6.9994685465110493d5*flucmom(ibox)*beta/
     &           (boxlx(ibox)*boxly(ibox)*boxlz(ibox))
            flucmom2(ibox) = acvsq(14,ibox)-avv(11,ibox)*avv(11,ibox)
     &           -avv(12,ibox)*avv(12,ibox) 
     &           - avv(13,ibox)*avv(13,ibox)
c            momconst = 6.9994685465110493d5
            flucmom2(ibox) = 6.9994685465110493d5*flucmom2(ibox)*beta/
     &           (boxlx(ibox)*boxly(ibox)*boxlz(ibox))            
         endif
c - boxlength
         acboxl(ibox,1) = acboxl(ibox,1) / acmove
         acboxl(ibox,2) = acboxl(ibox,2) / acmove
         acboxl(ibox,3) = acboxl(ibox,3) / acmove

         if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
            acboxa(ibox,1) = acboxa(ibox,1) / acmove
            acboxa(ibox,2) = acboxa(ibox,2) / acmove
            acboxa(ibox,3) = acboxa(ibox,3) / acmove
         endif

         acvol(ibox) = acvol(ibox) / acmove
         acvolsq(ibox) = acvolsq(ibox) / acmove

         do itype = 1, nmolty
c - number of molecules
            acnbox(ibox,itype) = acnbox(ibox,itype) / acmove
c - molfraction
            molfra(ibox,itype) = molfra(ibox,itype) / acmove
c - square end-to-end length
            if ( mnbox(ibox,itype) .gt. 0 ) then 
               asetel(ibox,itype) = 
     &              asetel(ibox,itype) / dble(mnbox(ibox,itype))
            endif
         enddo

         if ( lpbcz ) then
            do itype = 1, nmolty
c - number density
               acdens(ibox,itype) = 1000.0d0*acdens(ibox,itype) / acmove
            enddo
c - sum over all types of molecules
            temacd = 0.0d0
            do itype = 1, nmolty
               temacd = temacd + acdens(ibox,itype)
            enddo

c - molar volume
            molvol(ibox) = 602.2045d0 / temacd
            temspd = 0.0d0
            do itype = 1, nmolty
               temspd = temspd + 
     &              ( acdens(ibox,itype) * masst(itype) / 602.2045d0 )
            enddo
c - specific density
            speden(ibox) = temspd
         else
            do itype = 1, nmolty
c - number density
               acdens(ibox,itype) = 100.0d0*acdens(ibox,itype) / acmove
            enddo
            temacd = 0.0d0
            do itype = 1, nmolty
               temacd = temacd + acdens(ibox,itype)
            enddo
c - molar volume
            molvol(ibox) = 100.0d0 / temacd
         endif

c - system volume- convert to average
         acvolume(ibox) = acvolume(ibox) / acmove

c - pressure and surface tension
         if ( acnp .gt. 0.5d0 ) then
            acpres(ibox) = acpres(ibox) / acnp
            acsurf(ibox) = acsurf(ibox) / acnp
         endif

c thermodynamic integration stuff
          if (acipsw.gt.0.5d0) acdvdl = acdvdl/acipsw

c - chemical potential
         do itype = 1, nmolty
            if (.not. lrigid(itype)) then
               if( bnchem(ibox,itype) .gt. 0.5d0 ) then
c              --- determine how many steps it takes to grow molecule
c              --- not counting the first inserted bead
                  igrow = nugrow(itype)
                  debroglie = 17.458d0/( dsqrt(masst(itype)/beta ))
                  if (lrigid(itype)) then 
 	             call schedule(igrow,itype,steps,1,0,4)
                  else
                     call schedule(igrow,itype,steps,1,0,2)
                  endif
                  acchem(ibox,itype) = ((-1.0d0)/beta) * 
     &                 dlog(acchem(ibox,itype) /
     &                 ( dble( nchoi1(itype) ) 
     &                 * dble( nchoi(itype))**steps
     &                 * dble( nchoih(itype) ) 
     &                 * bnchem(ibox,itype) 
     &                 * debroglie*debroglie*debroglie ) )
               endif
            else
               if( bnchem(ibox,itype) .gt. 0.5d0 ) then
c              --- determine how many steps it takes to grow molecule
c              --- not counting the first inserted bead
                  debroglie = 17.458d0/( dsqrt(masst(itype)/beta ))
                  acchem(ibox,itype) = ((-1.0d0)/beta) *
     &                 dlog(acchem(ibox,itype) /
     &                 ( dble( nchoi1(itype) )
     &                 * dble( nchoir(itype))
     &                 * dble( nchoih(itype) )
     &                 * bnchem(ibox,itype)
     &                 * debroglie*debroglie*debroglie ) )
               endif
            endif
         enddo
         if (acvsq(1,ibox) .gt. 0.0d0) aflv(ibox) =dsqrt(acvsq(1,ibox))
      enddo

      write(2,1215) ('       Box ',i,i=1,nbox) 
      write(2,*)
      write(2,1209) (acpres(i) ,i=1,nbox)
      write(2,1212) ((acpres(i)*7.2429d-5),i=1,nbox)
      write(2,1216) (acsurf(i) ,i=1,nbox)
      do itype = 1, nmolty
         write(2,1210) itype, (acchem(i,itype) ,i=1,nbox)
      enddo
      write(2,*)

      do i = 1,3
         write(2,1202) (acboxl(ibox,i) ,ibox=1,nbox)
      enddo

      do ibox = 1, nbox
         if ( lsolid(ibox) .and. .not. lrect(ibox) ) then
            do i = 1,3
               write(2,1200) acboxa(ibox,i)*180.0d0/onepi
            enddo
         endif
      enddo

      do itype = 1, nmolty
         write(2,1201) itype, (acnbox(i,itype) ,i=1,nbox)
      enddo
      if ( lpbcz ) then
         write(2,1204) (molvol(i) ,i=1,nbox)
         write(2,1205) (speden(i) ,i=1,nbox)
         do itype = 1, nmolty
            write(2,1203) itype, (acdens(i,itype) ,i=1,nbox)
            if ( lexpand(itype) ) then
               do itype2 = 1, numcoeff(itype) 
                  write(2,1503) itype, itype2,acdens(itype,itype)*
     &                 acnbox2(itype,itype,itype2)/(acnbox(itype,itype)
     &                 *acmove),acnbox2(itype,itype,itype2)/
     &                 (acnbox(itype,itype)*acmove)
               enddo
            endif
         enddo
      else
         write(2,1214) (molvol(i), i=1,nbox)
         do itype = 1, nmolty
            write(2,1213) itype, (acdens(i,itype), i=1,nbox)
         enddo
      endif
      do itype = 1, nmolty
         write(2,1211) itype, (molfra(i,itype), i=1,nbox)
      enddo
      do itype = 1, nmolty
         write(2,1208) itype, (asetel(i,itype) ,i=1,nbox)
      enddo
      write(2,*)
      do j=1,10
c *** only 1 to 10 is the energy information         
         if ( lgibbs .and. .not. lnpt ) then
            write(2,1206) vname(j),(avv(j,i), 
c     +           0.00831451d0*avv(j,i)/acnbox(i,1),i=1,nbox)
     +           0.00831451d0*avv(j,i),i=1,nbox)
         else
            write(2,1206) vname(j),(avv(j,i), 
     +           0.008314510*avv(j,i),i=1,nbox)
         endif
      enddo
      
      write(2,*)
      write(2,1207) (aflv(i) ,i=1,nbox)
      write(2,*)

c ---   Output 2nd virial coefficient data
 2000 if (lvirial) then
         starviro = starvir
         dummy = dble(nstep/imv)
         do itemp = 1,ntemp
            starvir = starviro
            binvir(1,itemp) = binvir(1,itemp)/dummy
c            write(45,*) starvir,binvir(1,itemp)
            inside = starvir*starvir*binvir(1,itemp)
c            write(46,*) starvir,inside
            bvirial = 0.5d0*inside
c            write(47,*) starvir,bvirial
            starvir = starvir + stepvir
            
            do i = 2,nvirial-1
               binvir(i,itemp) = binvir(i,itemp)/dummy
c               write(45,*) starvir,binvir(i,itemp)
               inside = starvir*starvir*binvir(i,itemp)
c               write(46,*) starvir,inside
               bvirial = bvirial + inside
c               write(47,*) starvir,bvirial
               starvir = starvir + stepvir
            enddo
            
            binvir(nvirial,itemp) = binvir(nvirial,itemp)/dummy
c            write(45,*) starvir,binvir(nvirial,itemp)
            inside = starvir*starvir*binvir(nvirial,itemp)
c            write(46,*) starvir,inside
c            write(47,*) starvir,bvirial
            starvir = starvir + stepvir
            bvirial = bvirial + 0.5d0*inside
            
            write(2,*) 'At temperature of',virtemp(itemp)
            write(2,*) 'bvirial ',
     &           -(twopi*stepvir*bvirial),' [A^3 / molecule]'
            write(2,*) 'bvirial ',-0.602d0*twopi*
     &           stepvir*bvirial,' [cm^3 / mole]'

c            if ( lvirial2 ) then
               starvir = starviro + 0.5d0*stepvir
               do i = 2, nvirial
                  binvir2(i,itemp) = 
     &                 binvir2(i,itemp)/dummy
                  inside = starvir*starvir*binvir2(i,itemp)
                  bvirial = bvirial + inside
                  starvir = starvir + stepvir
               enddo
               bvirial = -(twopi*stepvir*bvirial)
               write(2,*) 'With quantum correction:'
               write(2,*) 'bvirial ',bvirial,' [A^3 / molecule]'
               write(2,*) 'bvirial ',0.602d0*bvirial,' [cm^3 / mole]' 
c            endif
         enddo
         
      endif


c - solute values
      write(2,*) 'type  box     vinter      vintra      vtor',
     &     '        vbend       vtail'
            
      do itype = 1, nmolty
         do ibox = 1, nbox
            if (solcount(ibox,itype).gt.0) then
               write(2,1372) itype,ibox,avsolinter(ibox,itype)
     &              /solcount(ibox,itype),avsolintra(ibox,itype)
     &              /solcount(ibox,itype),avsoltor(ibox,itype)
     &              /solcount(ibox,itype),avsolbend(ibox,itype)
     &              /solcount(ibox,itype),avsolelc(ibox,itype)
     &              /solcount(ibox,itype)
            else
               write(2,1372) itype,ibox,0.0,0.0,0.0,0.0,0.0
            endif              
         enddo
      enddo

c --- calculate statistical errors ---
      if ( nblock .ge. 2 ) then
         dblock = dble(nblock)
         dbl1 = dblock - 1.0d0
c -      global averages -
         do i = 1,nprop
            do j = 1,nbox
               if ( naccu(i,j) .lt. 0.5d-5 ) then
                  aver(i,j) = 0.0d0
               else
                  aver(i,j) = accum(i,j) / naccu(i,j)
                 
                  write(2,*) aver(i,j), accum(i,j) ,naccu(i,j) 

               endif
            enddo
         enddo
         do i = 1,nprop
            do j = 1,nbox
               do nbl = 1, nblock
                  dsq(i,j) = dsq(i,j) + 
     &                   ( baver(i,j,nbl) - aver(i,j) )**2
               enddo
               stdev(i,j) = dsqrt( dsq(i,j) / dblock )
               sterr(i,j) = dsqrt( dsq(i,j) / dbl1 )
               errme(i,j) = sterr(i,j) / dsqrt(dblock)
             enddo
         enddo


         do i = 1,nprop1
            do ibox = 1,nbox-1
               do jbox = ibox+1,nbox    
                  if ( naccu1(i,ibox,jbox) .lt. 0.5d-5 ) then
                     aver1(i,ibox,jbox) = 0.0d0
                  else
                     aver1(i,ibox,jbox) = accum1(i,ibox,jbox) / 
     &                                    naccu1(i,ibox,jbox)
                  endif
               enddo
            enddo 
         enddo

         do i = 1,nprop1
            do ibox = 1,nbox-1 
               do jbox = ibox+1,nbox
                  do nbl = 1, nblock
                      dsq1(i,ibox,jbox) = dsq1(i,ibox,jbox) +
     &                                  ( baver1(i,ibox,jbox,nbl)
     &                                   - aver1(i,ibox,jbox) )**2
                  enddo
                  stdev1(i,ibox,jbox) = dsqrt( dsq1(i,ibox,jbox) /
     &                                                 dblock )
                  sterr1(i,ibox,jbox) = dsqrt( dsq1(i,ibox,jbox) / 
     &                                                 dbl1 )
                  errme1(i,ibox,jbox) = sterr1(i,ibox,jbox) / 
     &                                              dsqrt(dblock)
               enddo
            enddo
         enddo
 
c - write out the heat of vaporization and solubility parameters
         do ibox = 1,nbox-1
            do jbox = ibox+1,nbox 
               write(2,1508) ibox,jbox,aver1(1,ibox,jbox),
     &                 stdev1(1,ibox,jbox),errme1(1,ibox,jbox)
               write(2,1509) ibox,jbox,aver1(2,ibox,jbox),
     &                 stdev1(2,ibox,jbox),errme1(2,ibox,jbox)
               write(2,1510) ibox,jbox,aver1(3,ibox,jbox),
     &                 stdev1(3,ibox,jbox),errme1(3,ibox,jbox)
               write(2,1511) ibox,jbox,aver1(4,ibox,jbox),
     &                 stdev1(4,ibox,jbox),errme1(4,ibox,jbox)
               write(2,1512) ibox,jbox,aver1(5,ibox,jbox),
     &                 stdev1(5,ibox,jbox),errme1(5,ibox,jbox)
               write(2,1513) ibox,jbox,aver1(6,ibox,jbox),
     &                 stdev1(6,ibox,jbox),errme1(6,ibox,jbox)
               write(2,1514) ibox,jbox,aver1(7,ibox,jbox),
     &                 stdev1(7,ibox,jbox),errme1(7,ibox,jbox)
               write(2,1515) ibox,jbox,aver1(8,ibox,jbox),
     &                 stdev1(8,ibox,jbox),errme1(8,ibox,jbox)
               write(2,1516) ibox,jbox,aver1(9,ibox,jbox),
     &                 stdev1(9,ibox,jbox),errme1(9,ibox,jbox)
            enddo
         enddo
 

c - specific density
         do ibox = 1, nbox
            write(2,1331) ibox,aver(1,ibox),stdev(1,ibox),errme(1,ibox)
         enddo

c * system volume
         itel = 4 + nener + 4*nmolty
         do ibox = 1, nbox
            write(2,1343) ibox,
     &           aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
         enddo

c - pressure
         do ibox = 1, nbox
            write(2,1341) ibox,aver(2,ibox),stdev(2,ibox),errme(2,ibox)
         enddo

c - surface tension
         itel = 2+nener+ 4*nmolty+1
         do ibox = 1, nbox
            write(2,1342) ibox,
     &           aver(itel,ibox),stdev(itel,ibox),errme(itel,ibox)
         enddo

         write(2,*)
c - energies
         do ibox = 1, nbox
           do j=3,2+10
c *** only 1 to 10 is the energy information
	      write(2,1311) vname(j-2),ibox,aver(j,ibox),
     +              stdev(j,ibox),errme(j,ibox)
           enddo
        enddo
      
        write(2,*)

c-- Enthalpy
        do ibox = 1,nbox
          j = 4+nener +4*nmolty + 1
          write(2, 1517) enth, ibox, aver(j,ibox),
     +          stdev(j,ibox),errme(j,ibox)
          j = 4+nener + 4*nmolty + 2
          write(2,1517) enth1,ibox,aver(j,ibox),
     +          stdev(j,ibox),errme(j,ibox)
        enddo 
        write(2,*)

c - chemical potential
        do itype = 1, nmolty
           itel = (2+nener) + itype
           do ibox = 1, nbox
              if ( aver(itel,ibox) .ne. 0.0d0 ) then
                 write(2,1351) itype,ibox,
     +                ((-1.0d0)/beta)*dlog(aver(itel,ibox)),
     +                (1.0d0/beta)*stdev(itel,ibox)/aver(itel,ibox),
     +                (1.0d0/beta)*errme(itel,ibox)/aver(itel,ibox)
              else
                 write(2,1351) itype,ibox,0.0d0,0.0d0,0.0d0
              endif
           enddo
        enddo

c - square end-to-end length
        do itype = 1, nmolty
           itel = (2+nener) + nmolty + itype
           do ibox = 1, nbox
              write(2,1321) itype,ibox,aver(itel,ibox)
     &         ,stdev(itel,ibox),errme(itel,ibox)
           enddo
        enddo

c - number density
        do itype = 1, nmolty
           itel = (2+nener) + 2 * nmolty + itype
           do ibox = 1, nbox
              if ( lpbcz ) then
                 write(2,1361) itype,ibox,1.0d3*aver(itel,ibox)
     +                ,1.0d3*stdev(itel,ibox),1.0d3*errme(itel,ibox)
              else
                 write(2,1361) itype,ibox,1.0d2*aver(itel,ibox)
     +                ,1.0d2*stdev(itel,ibox),1.0d2*errme(itel,ibox)
              endif
              if ( lexpand(itype) .and. 
     &             acnbox(ibox,itype) .gt. 0.5) then
                 do itype2 = 1, numcoeff(itype) 
                    molfrac = acnbox2(ibox,itype,itype2)
     &                   /(acmove*acnbox(ibox,itype))
                    write(2,1366) itype,itype2,1.0d3*aver(itel,ibox)*
     &                   molfrac,molfrac
                 enddo
              endif
           enddo
        enddo

c - molfraction
        do itype = 1, nmolty
           itel = (2+nener) + 3 * nmolty + itype
           do ibox = 1, nbox
              write(2,1371) itype,ibox,aver(itel,ibox),stdev(itel,ibox),
     +             errme(itel,ibox)
           enddo
        enddo

        if (lgibbs) then  
c --- write density results in fitting format ---
           do ibox = 1, nbox-1
              do jbox = ibox+1,nbox
                 if (speden(ibox).lt.speden(jbox)) then
                    ig = ibox
                    il = jbox
                 else
                    ig = jbox
                    il = ibox
                 endif
c                 write(41,1401) temp,aver(1,ig),stdev(1,ig)
c     &                ,aver(1,il),stdev(1,il)
c --- write ostwald values for each moltyp
                 gconst = 8.314/(1000*beta)
                 do itype = 1,nmolty
                    itel = (2+nener) + 2 * nmolty + itype
                    ostwald = aver(itel,il)/aver(itel,ig)
                    stdost  = ostwald
     &                   * dsqrt( (stdev(itel,il)/aver(itel,il))**2
     &                   + (stdev(itel,ig)/aver(itel,ig))**2 )
c                    write(42,*) nunit(itype),ostwald,stdost
c                    write(43,*) nunit(itype),
c     &                   -(gconst*log(ostwald)) + (eta2(ig,itype) 
c     &                   - eta2(il,itype)) / 120.27167
c     &                   ,gconst*stdost/ostwald
                    write(2,1506) itype,ig,il,ostwald,stdost
                    write(2,1507) itype,ig,il,
     &                   -(gconst*log(ostwald)) + (eta2(ig,itype) 
     &                   - eta2(il,itype)) / 120.27167
     &                   ,gconst*stdost/ostwald
                 enddo
              enddo
           enddo
        endif

        write(2,*)

c ---    write block averages  ---
        write(2,*)
        write(2,*) '-----block averages ------'
        do ibox=1,nbox
           write(2,1403) ibox
           do nbl = 1, nblock
              if (lgibbs .or. lgrand .or. lnpt ) then
                 write(2,1402) nbl,baver(3,ibox,nbl),baver(1,ibox,nbl),
     &                baver(2,ibox,nbl),baver(3+nener+4*nmolty,ibox,nbl)
     &                ,(baver(2+nener+3*nmolty+zz,ibox,nbl),zz=1,nmolty)
              else
                 write(2,1402) nbl,baver(3,ibox,nbl)
              endif
           enddo
           if (lmipsw) then
              write(2,*) 'lambdais', lambdais
              write(2,*) 'maginn interphase switch integrand'
              do nbl = 1, nblock
                 write(2,*) nbl,baver(nprop,ibox,nbl)
              enddo
           endif
        enddo
        
      endif

      call analysis(2) 	

c --- ee prob
      IF(lexpee) then                                
         write(2,*)
         write(2,*) 'probability of each mstate in ee'
         do nnn = 1, fmstate
            write(2,1601) nnn,ee_prob(nnn)
         enddo
      endif
      if (L_movie_xyz) then
        do ibox=1,nbox
           close (unit=210+ibox)
        enddo
      endif

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
 1206 format(a15,'[K and (kJ/mol)] =',3(f14.2,f12.2))
 1207 format(' fluctuation in <vtot>                          =',3f12.2)
 1208 format(' mean sete length of type   ',i4,'        [A^2] =',3f12.3)
 1209 format(' pressure                               [kPa] =',3f12.2)
 1210 format(' chem. potential of type    ',i4,'          [K] =',3f12.3)
 1211 format(' molfraction of type        ',i4,'              =',3f12.7)
 1212 format(' pressure                  [simulation units] =',3f12.6)
 1216 format(' surface tension                       [mN/m] =',3f12.4)
 1215 format(' Averages and fluctuations',21x,4(a11,i1))
c 1215 format(' Averages and fluctuations',21x,
c     & 3('       Box ',i1))

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
     + ' block   mu    msetl  density  press  '
     +       ,' dens(av/dif)  energies ... ') 
 1402 format(2x,i2,15(2x,e9.3))
 1403 format('  ------------ box: ' ,i4,/,
     + ' block   energy    density   pressure   surf ten  mol fracs') 
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
 1509 format(' H_vap LJ   [kJ/mol] btwn box   ',i4,' and',i4, ' =',
     & 3f15.4)
 1510 format(' H_vap Coul [kJ/mol] btwn box   ',i4,' and',i4, ' =',
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


 1517 format(a15,'[kJ/mol] for box',i3,' =',3(f12.4))
 1601 format(i5,1x,i10)


      stop
      end



