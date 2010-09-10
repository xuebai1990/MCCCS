      subroutine readdat(lucall,ucheck,nvirial,starvir
     &     ,stepvir,qelect)

c readdat
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

      implicit none
 
      include 'control.inc'
      include 'coord.inc'
      include 'cbmc.inc'
      include 'conver.inc'
      include 'system.inc'
      include 'poten.inc'
      include 'inpar.inc'
      include 'external.inc'
      include 'externalmuir.inc'
      include 'zeolite.inc'
      include 'nrtab.inc'
      include 'connect.inc'
      include 'inputdata.inc'
      include 'ewaldsum.inc'
      include 'swtcmove.inc'
      include 'fepsi.inc'
      include 'expand.inc'
      include 'qqlist.inc'
      include 'clusterbias.inc'
      include 'neigh.inc'
      include 'cell.inc'
      include 'nsix.inc'
c KM 01/10 remove analysis
c      include 'gor.inc'
      include 'torsion.inc'
      include 'tabulated.inc'
      include 'mpi.inc'
 
c -- variables for histograms	
      integer::fname
      character*20 file_movie,file_run,ftemp,fname2,file_dipole 
      character*20 fname4
      character*50 fileout
 
      integer::temnc, imol, iutemp, imolty, itype,ipair,bdum,bin,histtot
      integer::idummy(ntmax), atemp 

      integer::i,j,k,ncres, nmtres, iensem, inpbc, nmcount
      integer::im,nures, ibox,  ij, tcount,ucheck,nnframe
      integer::nijspecial,ispecial,jspecial,ji,ii,jj,nexclu,ndum,ntii

      integer::zz,temphe,z,itemp,zzz
      integer::nvirial,k_max_l,k_max_m,k_max_n
      integer::inclnum,inclmol,inclbead,inclsign,ncarbon
      dimension inclmol(ntmax*numax*numax),inclsign(ntmax*numax*numax)
      dimension inclbead(ntmax*numax*numax,2)

      integer::ainclnum,ainclmol,ainclbead,a15t
      dimension ainclmol(ntmax*numax*numax)
      dimension ainclbead(ntmax*numax*numax,2)
      dimension a15t(ntmax*numax*numax)

      real(8)::starvir,stepvir,fqtemp,qelect,qbox,vol,v(3),w(3)
      real(8)::debroglie, qtot,min_boxl

      real(8)::pie2,rcnnsq,umatch,aspecd,bspecd,dum,pm,pcumu
      logical::newtab,lnrtab,lucall,lpolar,lqqelect,lee,lratfix,lreadq
      logical:: linit, lecho, lmixlb, lmixjo, lhere,lsetup,lsolute
      logical::lprint,lverbose,lxyz,lfound

      dimension lratfix(ntmax)
      dimension qbox(nbxmax)
      
c -- variables added (3/24/05) for scaling of 1-4 interactions      
      real(8)::ofscale,ofscale2
      dimension ofscale(ntmax*numax*numax),ofscale2(ntmax*numax*numax)
      
c -- Variables added (6/30/2006) for fort.4 consistency check
  
      integer::numvib,numbend,numtor,vib1,bend2,bend3,tor2,tor3,tor4
      integer::vibtype,bendtype,tortype

c      real(8)::temx,temy,temz
      
      dimension nures(ntmax)
      dimension ncarbon(ntmax)
      dimension lhere(nntype)
      dimension lsolute(ntmax)
      dimension ucheck(ntmax)
      dimension qelect(nntype)
      dimension temphe(nntype)
c      dimension temx(nmax,numax),temy(nmax,numax),temz(nmax,numax)
      dimension k_max_l(nbxmax),k_max_m(nbxmax),k_max_n(nbxmax)

c Conversion factor for Mpa to simulation unit
      real(8)::MPa2SimUnits

      character(100) line

c KEA torsion variables
      logical::Lttor,lspline,linter
      integer::mmm,ttor
c KM tabulated potential variables
      integer::tvib, tbend, iivdW,jjvdW, iielect, jjelect
c KM variable added when analysis removed
      integer::nhere

c -- reads input data and initializes the positions
c
C --------------------------------------------------------------------
 
      MPa2SimUnits = 7.2429d-02

c *** set input arrays to zero ***
      do j=1, ntmax
         lratfix(j) = .false.
         nugrow(j) = 0
         nunit(j) = 0
         do i=1, numax
            ntype(j,i) = 0
         enddo
      enddo
      do i = 1,nntype**2
         lspecial(i) = .false.
      enddo
      do i = 1, nntype
         lhere(i) = .false.
      enddo

      lee = .false.
      qtot = 0.0d0
      ldie = .false.
 
c *** Output unit (if 2, write to runXX.dat file; if 6, write to stdout/screen; outherwise, 
c *** user designate a file to which to write) KEA 6/3/09 (defined in control.inc)

c KM for MPI
c must start at the top of the file, first two lines read in topmon.f for processor 0
      if (myid.ne.0) then
         open(4)
         read(4,*)
         read(4,*)
      endif

      read(4,*)
      read(4,*) iou

c *** To add or remove helium atoms
      read(4,*)
      read(4,*) L_add,N_add,N_box2add,N_moltyp2add
      read(4,*) 
      read(4,*) L_sub,N_sub,N_box2sub,N_moltyp2sub


c *** read echoing and long output flags
      read(4,*)
      read(4,*) lecho,lverbose

c - read whether to compute electrostatic interaction or not during CBMC/SWAP
      read(4,*)
      read(4,*) L_Coul_CBMC 
c-- read the number of unitcell replicated in each directions (a, b, and c)
      read(4,*)
      read(4,*) Num_cell_a,Num_cell_b,Num_cell_c 
c - read run information
      read(4,*) 
      read(4,*) run_num, suffix	
      read(4,*)
      read(4,*) nstep, lstop, lpresim, iupdatefix
c - read torsion, decide whether to use torsion in function form or Table
      read(4,*)
      read(4,*) L_tor_table,L_spline,L_linear
      read(4,*)
      read(4,*) L_vib_table, L_bend_table, L_vdW_table, 
     &     L_elect_table


c -- generate output file name
c - create output file name
      fname = run_num

c - use internal read/write to get integer::number in character format
      write(ftemp,*) fname
      read(ftemp,*) fname2
      file_run = 'run'//fname2(1:len_trim(fname2))//suffix//'.dat' 
      file_movie = 'movie'//fname2(1:len_trim(fname2))//suffix//'.dat' 
      file_dipole='dipole'//fname2(1:len_trim(fname2))//suffix//'.dat'

c KM ldielect writes to fort.27      
c      if (ldielect) then
c        open(unit=17,file=file_dipole,status='unknown')  
!        write(17,*) '# step  ibox   dipole_x   dipole_y   dipole_z'
c      endif

! kea 6/3/09 -- only open runXX.dat if io=2
      if(iou.eq.2.and.myid.eq.0) then
         open(unit=2,file=file_run,status='unknown')  
      endif  
c KM for MPI
c only processor 0 writes data  
      if ( lecho.and.myid.eq.0) then
         if (lverbose) then
            write(iou,*) 'Number of processors: ', numprocs
            write(iou,*) 'L_Coul_CBMC:',L_Coul_CBMC 
            write(iou,*) 'Number of unit cells in a dir = ', Num_cell_a
            write(iou,*) 'Number of unit cells in b dir = ', Num_cell_b
            write(iou,*) 'Number of unit cells in c dir = ', Num_cell_c
   
            if (lstop) then
               write(iou,*) 'number of steps:',nstep
            else
               write(iou,*) 'number of cycles:',nstep
            endif
            write(iou,*) 'lstep:',lstop
            write(iou,*) 'lpresim:',lpresim
            write(iou,*) 'iupdatefix:',iupdatefix
            write(iou,*) 'L_tor_table:',L_tor_table
            write(iou,*) 'L_spline:',L_spline
            write(iou,*) 'L_linear:',L_linear 
            write(iou,*) 'L_vib_table:', L_vib_table
            write(iou,*) 'L_bend_table:', L_bend_table
            write(iou,*) 'L_vdW_table:', L_vdW_table
            write(iou,*) 'L_elect_table:', L_elect_table
         else
            write(iou,*) L_Coul_CBMC,nstep, lstop, lpresim, iupdatefix
            write(iou,*) L_tor_table, L_spline,  L_linear
            write(iou,*) L_vib_table, L_bend_table, L_vdW_table, 
     &           L_elect_table
         endif
      endif

      
      
      read(4,*)
      read(4,*) nbox
      if ( lecho.and.myid.eq.0 ) write(iou,*) 'number of boxes 
     &     in the system:',nbox

      read(4,*)
      read(4,*) (express(ibox),ibox=1,nbox)
 
      read(4,*) 
      read(4,*) (ghost_particles(ibox),ibox=1,nbox) 

      read(4,*)
      read(4,*) L_Ewald_Auto

      read(4,*)
      read(4,*) temp, fqtemp,(Elect_field(ibox),ibox=1,nbox)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'temperature:',temp,' K'
            write(iou,*) 'external pressure:',express(1:nbox),' MPa'
            write(iou,*) 'Ghost particles: ',ghost_particles(1:nbox)
            write(iou,*) 'fluctuating charge temperature:',fqtemp,' K'
            write(iou,*) 'Electric field in z direction:',
     &           Elect_field(1:nbox),' V/A'
         else 
            write(iou,*) temp, (express(ibox),ibox=1,nbox), fqtemp,
     &           Elect_field
         endif
      endif
 
      do ibox = 1,nbox
         express(ibox)  = express(ibox)*MPa2SimUnits
      enddo

C$$   read the analysis information
      read(4,*)
      read(4,*)ianalyze,nbin,lrdf,lintra,lstretch,lgvst,lbend,lete,
     &           lrhoz,bin_width
      if (lecho.and.myid.eq.0) then
         if(lverbose) then
            write(iou,*) 'ianalyze:', ianalyze
	    write(iou,*) 'nbin', nbin
            write(iou,*) 'lrdf:',lrdf
            write(iou,*) 'lintra:',lintra
            write(iou,*) 'lstretch:', lstretch
            write(iou,*) 'lgvst:' ,lgvst
            write(iou,*) 'lbend:' ,lbend
            write(iou,*) 'lete:' ,lete
            write(iou,*) 'lrhoz:', lrhoz 
            write(iou,*) 'bin_width:' ,bin_width
         else
            write(iou,*) ianalyze,nbin,lrdf,lintra,lstretch,lgvst,lbend
     &    ,lete, lrhoz, bin_width  
         endif
      endif

C -------------------------------------------------------------------

c *** set up constants and conversion factors ***     
      pie2 = 8.0d0 * datan(1.0d0)
      raddeg = 360.0d0 / pie2
      degrad = pie2 / 360.0d0
      twopi=pie2
      onepi=pie2/2.d00
c  e times V to Kelvin
      eXV_to_K = 11600.0d0  

      
      beta = 1.0d0 / temp
      fqbeta = 1.0d0 / fqtemp
 
C -------------------------------------------------------------------

      read(4,*)
      read(4,*) iprint, imv, iratio, iblock, idiele, L_movie_xyz,
     & iheatcapacity

      if (L_movie_xyz) then
        do ibox = 1,nbox
          write(ftemp,*) ibox
          read(ftemp,*) fname4
          fileout = 'box'//fname4(1:len_trim(fname4))//'movie'//
     &               fname2(1:len_trim(fname2))
     &                         //suffix//'.xyz'
	if (myid.eq.0) then
          open (unit=210+ibox,FILE=fileout,status='unknown')
        endif
        enddo
      endif

      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'iprint:',iprint
            write(iou,*) 'imv:',imv
            write(iou,*) 'iratio:',iratio
            write(iou,*) 'iblock:',iblock
            write(iou,*) 'idiele:',idiele
            write(iou,*) 'L_movie_xyz',L_movie_xyz
         else
            write(iou,*) iprint, imv, iratio, iblock, idiele,L_movie_xyz
         endif
      endif
c - read information for histogram output (added 8/30/99)
      if (lgrand) then
         read(4,*)
         read(4,*) nequil,ninstf, ninsth, ndumph 
         if ( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               write(iou,*) 'nequil:',nequil
               write(iou,*) 'ninstf:',ninstf
               write(iou,*) 'ninsth:',ninsth
               write(iou,*) 'ndumph:',ndumph
               write(iou,*) 'run_num:',run_num
               write(iou,*) 'suffix:',suffix
            else
               write(iou,*) nequil,ninstf,ninsth,ndumph,run_num,suffix
            endif
         endif
      end if
      
c - create output file name
      fname = run_num

c - use internal read/write to get integer::number in character format
c      write(ftemp,*) fname
c      read(ftemp,*) fname2
c      file_run = 'run'//fname2(1:len_trim(fname2))//suffix//'.dat' 
c      file_movie = 'movie'//fname2(1:len_trim(fname2))//suffix//'.dat' 

c KM for MPI
c jobs stop in monola so that all processors die
      if (dint(dble(nstep)/dble(iblock)) .gt. 100) then
         write(iou,*) 'too many blocks'
         ldie = .true.
         return
      endif

c - read system information

      do i = 1,nbox
         read(4,*)
         read(4,*) boxlx(i),boxly(i),boxlz(i),lsolid(i),lrect(i),
     +        kalp(i),rcut(i),rcutnn(i)
         if ( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               write(iou,*) 'box:',i
               write(iou,*) '   boxlx:',boxlx(i),' A'
               write(iou,*) '   boxly:',boxly(i),' A'
               write(iou,*) '   boxlz:',boxlz(i),' A'
               write(iou,*) '   lsolid:',lsolid(i)
               write(iou,*) '   lrect:',lrect(i)
               write(iou,*) 'neighbor list cutoff (rcutnn):',rcutnn(i),
     &              'A'
               write(iou,*) '   rcut:',rcut(i),'A'
               if(.not.L_Ewald_Auto) then
                   write(iou,*) '   kalp:',kalp(i)
               endif 
            else
               write(iou,*) boxlx(i),boxly(i),boxlz(i),
     +        lsolid(i),lrect(i),rcut(i),rcutnn(i)
              if (.not.L_Ewald_Auto) then  
                 write(iou,*) kalp(i)
              endif
            endif
         endif
      enddo

      read(4,*) 
      read(4,*) nchain, nmolty
      if ( nmolty .gt. ntmax ) then
         write(iou,*) 'nmolty gt ntmax'
         ldie = .true.
         return
      endif
      if ( nchain .gt. nmax-2 ) then
         write(iou,*) 'nchain gt nmax-2'
         ldie = .true.
         return
      endif
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'number of chains:',nchain
            write(iou,*) 'number of molecule types:',nmolty
         else 
            write(iou,*) nchain, nmolty
         endif
      endif
      read(4,*)
      read(4,*) (temtyp(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            do i = 1,nmolty
               write(iou,*) 'number of chains of molecule type',i,':',
     &              temtyp(i)
            enddo
         else 
            write(iou,*) (temtyp(i),i=1,nmolty)
         endif
      endif
      temnc = 0
      do i = 1, nmolty
         do j = 1, temtyp(i)
            temnc = temnc + 1
            moltyp(temnc) = i
         enddo
      enddo

      if (lgrand) then
c - read chemical potentials (added 8/30/99 by jpotoff)
         read(4,*) 
         read(4,*) (B(i),i=1,nmolty)
         if (lecho.and.myid.eq.0) then
            if (lverbose) then
               do i = 1,nmolty
                  write(iou,*) 'chemical potential for molecule type',
     &                 i,':',B(i)
               enddo
            else 
               write(iou,*) "B ", (B(i),i=1,nmolty)
            endif
         endif
C --- removing from here. It will be calculated once we have molecular mass 
C --- to calculate debroglie wavelength.

c - convert chemical potentials to activities 
c         do i=1,nmolty
c            B(i) = exp(B(i)/temp) 
c         enddo
         
      end if

c      if ( lecho.and.myid.eq.0 ) write(iou,*) 'moltyp',(moltyp(i),i=1,nchain)

      if (lgrand) then
         nchain=nmax
         write(iou,*)'in GCMC total number of chains set by NMAX!'
      endif

      read(4,*)
      read(4,*) lmixlb, lmixjo
      if (lmixlb .and. lmixjo) then
         write(iou,*) 'cant use both combining rules!'
         ldie = .true.
         return
      endif
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            if (lmixlb) then
               write(iou,*) 'Lorentz-Berthelot combining rules apply'
            else
               write(iou,*) 'Jorgensen combining rules apply'
            endif
            write(iou,*) '   lmixlb:',lmixlb,' lmixjo:',lmixjo
         else 
            write(iou,*) lmixlb,lmixjo
         endif
      endif
c --- read special combining rule information
      read(4,*)
      read(4,*) nijspecial
      if (lecho.and.myid.eq.0) then
         if (lverbose) then
            write(iou,*) 'number of special combining parameters:',
     &           nijspecial
         else
            write(iou,*) nijspecial
         endif
      endif
      read(4,*)
      if ( nijspecial .eq. 0 ) then
         read(4,*)
      else
         do i = 1,nijspecial
            read(4,*) ispecial,jspecial,aspecd,bspecd
            ij = (ispecial-1)*nntype + jspecial
            ji = (jspecial-1)*nntype + ispecial
            aspecial(ij) = aspecd
            bspecial(ij) = bspecd
            aspecial(ji) = aspecd
            bspecial(ji) = bspecd
            lspecial(ij) = .true. 
            lspecial(ji) = .true.
            if (lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(iou,*) 'special parameter number',i
                  write(iou,*) '   ispecial:',ispecial
                  write(iou,*) '   jspecial:',jspecial
                  write(iou,*) '   aspecd:',aspecd
                  write(iou,*) '   bspecd:',bspecd
               else 
                  write(iou,*) ispecial,jspecial,aspecd,bspecd
               endif
            endif
         enddo
      endif

      read(4,*) 
      read(4,*) rmin, softcut,rcutin,
     &     rbsmax,rbsmin
      if ( lecho .and.myid.eq.0) then
         if (lverbose) then
            write(iou,*) 'minimum cutoff (rmin):',rmin,' A'
            write(iou,*) 'softcut:',softcut
            write(iou,*) 'CBMC inner cutoff (rcutin):',rcutin,' A'
            write(iou,*) 'AVBMC outer cutoff (rbsmax):',rbsmax,' A'
            write(iou,*) 'AVBMC inner cutoff (rbsmin):',rbsmin,' A'
         else 
            write(iou,*) rmin, softcut
     &           ,rcutin,rbsmax,rbsmin 
         endif
      endif
      do i = 1, nbox
         if( rcut(i)/boxlx(i) .gt. 0.5d0) then
            write(iou,*) 'rcut > 0.5*boxlx'
            ldie = .true.
            return
         endif
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
         call setpbc(i)
      enddo

      softlog = 10.0d0**(-softcut)
      vol_eff = (4.0d0/3.0d0)*onepi*
     &     (rbsmax*rbsmax*rbsmax-rbsmin*rbsmin*rbsmin)

c - set up the strectching and bending constants
      call suvibe
c - set up the forcefield and the masses
      call suijtab( lmixlb,lmixjo,qelect )      
     
c - read bead potential information
      do imol = 1, nmolty
         read(4,*) 
         read(4,*) nunit(imol),nugrow(imol),ncarbon(imol),nmaxcbmc(imol)
     &        , iurot(imol),lelect(imol),lflucq(imol),lqtrans(imol)
     &        ,lexpand(imol),lavbmc1(imol),lavbmc2(imol),lavbmc3(imol)
     &        ,fqegp(imol)
         read(4,*)
         read(4,*) maxgrow(imol),lring(imol),lrigid(imol)
     &        ,lrig(imol),lsetup,isolute(imol),(eta2(i,imol), i=1,nbox)

         read(4,*) 
         read(4,*) lq14scale(imol),qscale(imol)

         if (isolute(imol).lt.nstep) then
            lsolute(imol) = .true.
         else
            lsolute(imol) = .false.
         endif

         if (lring(imol)) then
            read(4,*)
            read(4,*) iring(imol)
         else
            iring(imol) = nunit(imol)
         endif

         do i = 1, nunit(imol)
            lrigi(imol,i) = .false.
         enddo

c     *** irig is the site rigid sites will be grown from
c     *** and frig will be the previous site (not kept rigid)

         if (lrig(imol)) then
            read(4,*) 
            read(4,*) nrig(imol)

            if (nrig(imol).gt.0) then
c     --- read in specific points to keep rigid in growth
               read(4,*) 
               do i = 1, nrig(imol)
                  read(4,*) irig(imol,i),frig(imol,i)
                  lrigi(imol,irig(imol,i)) = .true.
               enddo
            else
c     --- we will pick irig at random in each case if nrig = 0
               read(4,*) 
               read(4,*) nrigmin(imol),nrigmax(imol)
               
c     --- nrigmin is the minimum amount of the chain to keep rigid
c     --- nrigmax is the maximum

            endif
         endif

         if (lrigid(imol)) then
            read(4,*) 
c     - number of flexible parts
            read(4,*) rindex(imol)
            if ( rindex(imol).gt.0) then
               do i = 1, rindex(imol)
                  read(4,*) riutry(imol,i)
               enddo
            else
               riutry(imol,1) = 1
            endif
         endif   
         

         if ( nunit(imol) .gt. numax ) then
            write(iou,*) 'nunit gt numax'
            ldie = .true.
            return
         endif
         if ( lflucq(imol)
     &        .and. (.not. lelect(imol) ) ) then
            write(iou,*) 'lelect must be true if flucq is true' 
            ldie = .true.
            return
         endif
         if ( lqtrans(imol) ) then
            if (.not. lflucq(imol) ) then
               write(iou,*) 'lflucq must be true if interm. 
     &              CT is allowed'
               ldie = .true.
               return
            endif
            write(iou,*) 'Intermolecular Charge Transfer is allowed'
         endif

         lbias(imol) = .false.
         if (lavbmc1(imol) .or. lavbmc2(imol) .or. lavbmc3(imol)) then
            lbias(imol) = .true.
         endif

c *** choose only one from the three AVBMC algorithms

         if ( lavbmc1(imol) ) then
            lavbmc2(imol) = .false.
            lavbmc3(imol) = .false.
         elseif ( lavbmc2(imol) ) then
            lavbmc3(imol) = .false.
         endif

c         lneighbor = .true.
         lneighbor = .false.
         if ( (lavbmc2(imol) .or. lavbmc3(imol)).and.(.not.lgaro) ) 
     &                    lneighbor = .true.

         if ( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               write(iou,*) 'molecule type:',imol
               write(iou,*) '   number of units:',nunit(imol)
               write(iou,*) '   number of units for CBMC growth:',
     &              nugrow(imol)
               write(iou,*) '   number of carbons for EH alkane:',
     &              ncarbon(imol)
               write(iou,*) '   maximum number of units for CBMC:',
     &              nmaxcbmc(imol)
               write(iou,*) '   iurot:',iurot(imol)
               write(iou,*) '   lelect:',lelect(imol)
               write(iou,*) '   lflucq:',lflucq(imol)
               write(iou,*) '   lqtrans:',lqtrans(imol)
               write(iou,*) '   lexpand:',lexpand(imol)
               write(iou,*) '   lavbmc1:',lavbmc1(imol)
               write(iou,*) '   lavbmc2:',lavbmc2(imol)
               write(iou,*) '   lavbmc3:',lavbmc3(imol)
               write(iou,*) '   fqegp:',fqegp(imol)
               write(iou,*) '   lsetup:',lsetup
               write(iou,*) '   lq14scale:',lq14scale(imol)
               write(iou,*) '   qscale:',qscale(imol)
               do i = 1,nbox
                  write(iou,*) '   energy offset for box',
     &                 i,':',eta2(i,imol),' K'
               enddo
            else 
               write(iou,*) nunit(imol),nugrow(imol),ncarbon(imol)
     &        ,nmaxcbmc(imol),iurot(imol) ,lelect(imol),lflucq(imol) 
     &        ,lqtrans(imol),lexpand(imol),lavbmc1(imol),lavbmc2(imol)
     &        ,lavbmc3(imol),fqegp(imol)
     &        ,lsetup,(eta2(i,imol), i=1,nbox)
            endif
         endif
         masst(imol) = 0.0d0

         if (lsetup) then
            call molsetup(imol)

            do i = 1,nunit(imol)
               lhere(ntype(imol,i)) = .true.
            enddo

            goto 112
         endif

         do i = 1, nunit(imol)
c - linear/branched chain with connectivity table -
            read(4,*)
            if ( lelect(imol) .and. .not. lchgall ) then
               read(4,*) j, ntype(imol,i), leaderq(imol,i)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
c for other potentials, sigi isn't defined, so I added this if/else statement
c maybe should have a switch for every potential flag as a check?
               IF(llj)THEN
               if(sigi(ntype(imol,i)).lt.1d-06.and.epsi(ntype(imol,i))
     &          .lt.1d-06.and.abs(qelect(ntype(imol,i))).lt.1d-06) then
                  write(iou,*)
                  write(iou,*) '****PROBLEM IN SUIJTAB****'
                  write(iou,*) 'check if the beadtyp',ntype(imol,i),
     &           'is defined'
                 ldie = .true.
                 return
               endif
c KM this isn't necessary
c               ELSE
c                  if (myid.eq.0)  write(iou,*) 
c     &                 'Confirm that your parameters are defined.'
               ENDIF

               if ( lecho.and.myid.eq.0 ) 
     &              write(iou,*) '   bead ',j,' beadtype ',ntype(imol,i)
     &              ,' charge leader ',leaderq(imol,i)
               if ( leaderq(imol,i) .gt. j .and. .not. lchgall) then
                  write(iou,*) 'group-based cut-off screwed for qq'
                  ldie = .true.
                  return
               endif
            else
               read(4,*) j, ntype(imol,i)
               if ( lecho.and.myid.eq.0 ) then
                  if (lverbose) then
                     write(iou,*) '   bead ',j,' beadtype ',
     &                    ntype(imol,i),chname(ntype(imol,i))
                  else 
                     write(iou,*) '   bead ',j,' beadtype ',
     &                    ntype(imol,i)
                  endif
               endif
            endif
            iutemp = ntype(imol,i)

            if (lpl(iutemp)) then
               lplace(imol,i) = .true.
               llplace(imol) = .true.
            else
               lplace(imol,i) = .false.
            endif

            masst(imol)=masst(imol)+mass(iutemp)
            lhere(iutemp) = .true.
            
c - bond vibration -
            read(4,*)
            read(4,*) invib(imol,i)
            if ( invib(imol,i) .gt. 6 ) then
               write(iou,*) 'imol',imol,'   i',i,'  invib',invib(imol,i)
               write(iou,*) 'too many vibrations'
               ldie = .true.
               return
            endif
            do j = 1, invib(imol,i)
               read(4,*) ijvib(imol,i,j),itvib(imol,i,j)
                if(brvib(itvib(imol,i,j)).lt.1d-06) then
                  write(iou,*)
                  write(iou,*) '****PROBLEM IN SUVIBE****'
                  write(iou,*) 'check if the vib type ',itvib(imol,i,j),
     &           'is defined'
                 ldie = .true.
                 return
               endif            
               
               if((ijvib(imol,i,j).eq.i).or.(ijvib(imol,i,j).gt.
     &   nunit(imol))) then
                 write(iou,*) 'check vibrations for mol type',imol,
     &                        'and bead',i
                 ldie = .true.
                 return
               endif              
               if (lverbose.and.myid.eq.0) then
c                  write(iou,*) '      bead',i,' bonded to bead',
c     &                 ijvib(imol,i,j),' with bond type:',
c     &                 itvib(imol,i,j)
                  write(iou,*) '      bead',i,' bonded to bead',
     &                 ijvib(imol,i,j)
                  write(iou,1013) '          bond type:',
     &                 itvib(imol,i,j),' bond length:',
     &                 brvib(itvib(imol,i,j)),' k/2:',
     &                 brvibk(itvib(imol,i,j))
               endif
            enddo
c - bond bending -
            read(4,*)
            read(4,*) inben(imol,i)
c            write(iou,*) inben(imol,i)
            if ( inben(imol,i) .gt. 12 ) then
               write(iou,*) 'too many bends'
               ldie = .true.
               return
            endif
            do j = 1, inben(imol,i)
               read(4,*) ijben2(imol,i,j),ijben3(imol,i,j)
     &              ,itben(imol,i,j)
               if(brben(itben(imol,i,j)).lt.1d-06) then
                 write(iou,*)
                 write(iou,*) '****PROBLEM IN SUVIBE****'
                 write(iou,*) 'check if thE bend type',itben(imol,i,j),
     &                 'is defined'
                 ldie = .true.
                 return
               endif
                         
               
               if ((ijben2(imol,i,j).gt.nunit(imol)).or.(
     &                ijben3(imol,i,j).gt.nunit(imol))) then
                   write(iou,*) 'check bending for the mol type',imol,
     &                      'bead',i
                   ldie = .true.
                   return
               endif
               if ((ijben2(imol,i,j).eq.i).or.(
     &                ijben3(imol,i,j).eq.i).or.(ijben2(imol,i,j)
     &                  .eq.ijben3(imol,i,j))) then
                   write(iou,*) 'check bending for the mol type',imol,
     &                      'bead',i
                   ldie = .true.
                   return
               endif
 
   
               if (lverbose.and.myid.eq.0) then
c                  write(iou,*) '      bead',i, ' bending interaction',
c     &                 ' through',ijben2(imol,i,j),' with bead',
c     &                 ijben3(imol,i,j),' of bend type:',itben(imol,i,j)
                  write(iou,1017) '      bead',i, ' bending interaction'
     &                 ,' through',ijben2(imol,i,j),' with bead',
     &                 ijben3(imol,i,j)
                  write(iou,1013) '          bend type:',itben(imol,i,j)
     &                 ,' bend angle :',
     &                 brben(itben(imol,i,j))*180.0d0/onepi,
     &                 ' k/2:',brbenk(itben(imol,i,j))
               endif
            enddo
c - bond torsion -
            read(4,*)
            read(4,*) intor(imol,i)
            if ( intor(imol,i) .gt. 12 ) then
               write(iou,*) 'too many torsions'
               ldie = .true.
               return
            endif
            do j = 1, intor(imol,i)
               read(4,*) ijtor2(imol,i,j),ijtor3(imol,i,j),
     &              ijtor4(imol,i,j),ittor(imol,i,j)

               if(ijtor2(imol,i,j).gt.nunit(imol).or.ijtor3(imol,i,j)
     &         .gt.nunit(imol).or.ijtor4(imol,i,j).gt.nunit(imol)) then 
                   write(iou,*) 'check torsion for the mol type',imol,
     &                      'bead',i
                   ldie = .true.
                   return
               endif
     
               if((ijtor2(imol,i,j).eq.i.or.ijtor3(imol,i,j).eq.i
     &          .or.ijtor4(imol,i,j).eq.i).or.(ijtor2(imol,i,j).eq.
     &          ijtor3(imol,i,j).or.ijtor2(imol,i,j).eq.
     &           (ijtor4(imol,i,j)).or.(ijtor3(imol,i,j).eq.
     &              ijtor4(imol,i,j)))) then
                   write(iou,*) 'check torsion for the mol type',imol,
     &                      'bead',i
                   ldie = .true.
                   return
               endif

               if (lverbose.and.myid.eq.0) then
                  write(iou,1018) '      bead',i, ' torsional '
     &                 ,'interaction through',ijtor2(imol,i,j),' and',
     &                 ijtor3(imol,i,j),' with bead',ijtor4(imol,i,j),
     &                 ' of torsional type:',ittor(imol,i,j)
               endif
            enddo
         enddo

 112     continue


c -- starting the self consistency check for the bond vibrations, bending, and torsions
c -- this would help in catching errors in fort.4 connectivity. Starting after continue
c -- so that if we use molsetup subroutine it will provide extra checking. (Neeraj)
      if (.not.lrigid(imol)) then 
      do i = 1,nunit(imol)
         numvib=invib(imol,i)
         numbend=inben(imol,i)
         numtor=intor(imol,i)      
         lfound = .false.       
         do j = 1,numvib
            vib1 = ijvib(imol,i,j)
            vibtype  = itvib(imol,i,j)
            if(invib(imol,vib1).eq.0) then
               write(iou,*) 'ERROR IN FORT.4 VIBRATIONS'
               write(iou,*) 'Check vibration for mol. type:',imol,
     &                   'bead',vib1,'with',i
               ldie = .true.
               return
            endif              
            do k =1,invib(imol,vib1)
               if(ijvib(imol,vib1,k).eq.i) then
                  lfound = .true. 
                  if(vibtype.ne.itvib(imol,vib1,k)) then
                     write(iou,*) 'Error in fort.4 vibration',
     &                  ' specifications'
                     write(iou,*) 'check vibration type of bead',i,
     &                  'with',vib1,'molecule type',imol,'vice versa'
                     ldie = .true.
                     return
                   endif
               endif
            enddo
            if(.not.lfound) then
              write(iou,*) 'Error in fort.4 vibration iformation'
              write(iou,*) 'Check vibration for mol. type:',imol,
     &                   'bead ',vib1,'with ',i
              ldie = .true.
              return
            endif
         enddo 
         lfound= .false.
         do j = 1,numbend
            bend2 = ijben2(imol,i,j)
            bend3 = ijben3(imol,i,j)
c            if ( i .eq. 17) then
c              write(iou,*) bend2, bend3, numbend
c            endif
            bendtype = itben(imol,i,j)
            if(inben(imol,bend3).eq.0) then
               write(iou,*) 'ERROR IN FORT.4 BENDING'
               write(iou,*) 'Check bending for mol. type:',imol,
     &                   'bead ',bend3,'with ',i
               ldie = .true.
               return
            endif
            do k = 1,inben(imol,bend3)
               if((ijben2(imol,bend3,k).eq.bend2).and.
     &                 (ijben3(imol,bend3,k).eq.i)) then
                 lfound = .true.
                 if(itben(imol,bend3,k).ne.bendtype) then
                   write(iou,*) 'Error in fort.4 bending',
     &                  ' specifications'
                   write(iou,*) 'check bending type of bead',i,
     &                  'with',bend3,'mol. typ.',imol,'and vice versa'
                   ldie = .true.
                   return
                 endif
               endif
            enddo   
            if(.not.lfound) then
              write(iou,*) 'Error in fort.4 bending iformation'
              write(iou,*) 'Check bending for mol. type:',imol,
     &                 'bead ',bend3,'with ',i
              ldie = .true.
              return
            endif
         enddo
         lfound = .false. 
         do j = 1,numtor
            tor2 = ijtor2(imol,i,j)
            tor3 = ijtor3(imol,i,j)
            tor4 = ijtor4(imol,i,j)
            tortype = ittor(imol,i,j) 
            if(intor(imol,tor4).eq.0) then
               write(iou,*) 'ERROR IN FORT.4 TORSION'
               write(iou,*) 'Check torsion for mol. type:',imol,
     &                   'bead ',tor4,'with ',i,'and vice versa'
               ldie = .true.
               return
            endif
            do k = 1,intor(imol,tor4)
               if((ijtor2(imol,tor4,k).eq.tor3).and.(ijtor3(imol,tor4
     &            ,k).eq.tor2).and.(ijtor4(imol,tor4,k).eq.i)) then
                 lfound=.true. 
                 if(ittor(imol,tor4,k).ne.tortype) then
                   write(iou,*) 'Error in fort.4 torsion',
     &                  ' specifications'
                   write(iou,*) 'check torsion type of bead',i,
     &                  'with',tor4,'mol. typ.',imol,'and vice versa'
                   ldie = .true.
                   return
                 endif 
               endif
            enddo
            if(.not.lfound) then
               write(iou,*) 'Error in fort.4 torsion iformation'
               write(iou,*) 'Check torsion for mol. type:',imol,
     &                'bead ',tor4,'with ',i
               ldie = .true.
               return
            endif
         enddo
      enddo    
      endif   

c -- Neeraj Adding molecule neutrality check
ckea skip if lgaro
      if(.not.(lgaro .or.lionic)) then
         do i=1,nmolty
            qtot =0.0d0
            do j = 1,nunit(i)
               qtot = qtot+qelect(ntype(i,j))
            enddo
            if(dabs(qtot).gt.1d-7) then
               write(iou,*)'molecule type',i,'not neutral check charges'
               ldie = .true.
               return
            endif
         enddo
      endif
        

         if ( lexpand(imol) ) then
            if ( temtyp(imol) .gt. 1 ) then
               write(iou,*) 'Only one molecule of this type is allowed!'
               ldie = .true.
               return
            endif
            lee = .true.
            read(7,*) 
            read(7,*) numcoeff(imol)
            do j = 1,numcoeff(imol)
               read(7,*)
               read(7,*) (epsil(imol,ii,j),ii=1,nunit(imol))
               write(iou,*) 'itype:',j
               write(iou,*) (epsil(imol,ii,j),ii=1,nunit(imol))
               read(7,*) (sigm(imol,ii,j),ii=1,nunit(imol))
               write(iou,*) (sigm(imol,ii,j),ii=1,nunit(imol))
               read(7,*) (qcharge(imol,ii,j),ii=1,nunit(imol))
               write(iou,*) (qcharge(imol,ii,j),ii=1,nunit(imol))
               read(7,*)
               read(7,*) (eta(ii,imol,j),ii=1,2)
               write(iou,*) 'eta:',(eta(ii,imol,j),ii=1,2)
            enddo
         endif
         if ( lbias(imol) ) then
            read(4,*)
            read(4,*) pmbias(imol),(pmbsmt(ii),ii=1,nmolty)
     &           ,pmbias2(imol)
            if ( lecho .and.myid.eq.0) then
               if (lverbose) then
                  write(iou,*) '   AVBMC pmbias',pmbias(imol)
                  do ii = 1,nmolty
                     write(iou,*) '   AVBMC2 and 3 probability for',
     &                    ' molecule',' type',ii,':',pmbsmt(ii)
                  enddo
                  write(iou,*) '   AVBMC3 pmbias2:',pmbias2(imol)
               else
                  write(iou,*) '   AVBMC bias for cluster formation and'
     &                 ,' destruction'
                  write(iou,*) pmbias(imol),(pmbsmt(ii),ii=1,nmolty)
     &                 ,pmbias2(imol)
               endif
            endif
            if (rbsmax .lt. rbsmin) then
               write(iou,*)'rbsmax should be greater than rbsmin'
               ldie = .true.
               return
            endif
         endif
ckea 6/4/09 -- added for multiple rotation centers
c -- To assign multiple rotation centers, set iurot(imol) < 0
c -- Add line after molecule specification, avbmc parameters
c -- First, number of rotation centers
c -- Second, identity of centers (0=COM,integer::> 0 = bead number)
c -- Third, give probability to rotate around different centers
         if(iurot(imol).lt.0) then
            read(4,*)
            read(4,*) nrotbd(imol),irotbd(1:nrotbd(imol),imol),
     &           pmrotbd(1:nrotbd(imol),imol)
            if( lecho.and.myid.eq.0 ) then
               if ( lverbose ) then
                  write(iou,*) ' Multiple rotation centers',nrotbd(imol)
                  do ii=1,nrotbd(imol)
                     write(iou,*) 'Rotation center',ii,':',
     &                    irotbd(ii,imol),'    probability to select:',
     &                    pmrotbd(ii,imol)
                  enddo
               endif
            endif
         endif
      enddo

c   KEA - read in tabulated torsion potentials and set up derivatives
c         for spline interpolation if L_tor_table 
       if (L_tor_table) then
          read(40,*) nttor
          write(iou,*)
c          write(6,*) 'in spline read loop, nttor',nttor
          do mmm = 1,nttor
             read(40,*) ttor
             i=1
 10          read(40,*,end=11) deg(i,ttor),tabtorso(i,ttor)
             if(deg(i,ttor).eq.1000) goto 11
             if (myid.eq.0) then
                write(56,*) i,deg(i,ttor),tabtorso(i,ttor)
             endif
             i=i+1
             goto 10
 11          splpnts(ttor)=i-1
             if (L_spline) then
                if (myid.eq.0) write(iou,*) 'using spline interpolation'
                call spline(1.0d31,1.0d32,ttor)
             elseif (L_linear) then
                if (myid.eq.0) write(iou,*) 'using linear interpolation'
                do i=1,splpnts(ttor)-1
                   tordif(i,ttor) = tabtorso(i+1,ttor)-tabtorso(i,ttor)
                enddo
             endif
          enddo
          close(40)
       endif

c   KM 12/02/08 - read in tabulated vibrational potential
c   and set up linear interpolation

       if (L_vib_table) then
          read(41,*) ntabvib
          write(iou,*)
          mmm=1
          do mmm=1,ntabvib
             read(41,*) tvib
             read(41,*) num_int_vib(tvib)
             i=1
 12          read(41,*,end=13) vib(i,tvib), tabvib(i,tvib)
             if (vib(i,tvib).eq.1000) goto 13
c             write(57,*) i, vib(i,tvib),tabvib(i,tvib)
             i=i+1
             goto 12
 13          vibsplits(tvib)=i-1
             if (myid.eq.0)  write(iou,*) 'using linear interpolation 
     &            for vibrations'
             do i=1,vibsplits(tvib)-1
                vibdiff(i,tvib)=tabvib(i+1,tvib)-tabvib(i,tvib)
             enddo
          enddo
          close(41)
       endif

c   KM 12/02/08 - read in tabulated 1-3 nonbonded 'bending' potential
c   and set up linear interpolation

       if (L_bend_table) then
          read(42,*) ntabbend
          write(iou,*)
          mmm=1
          do mmm=1,ntabbend
             read(42,*) tbend
             read(42,*) num_int_bend(tbend)
             i=1
 14          read(42,*,end=15) bend(i,tbend), tabbend(i,tbend)
             if (bend(i,tbend).eq.1000) goto 15
c             write(58,*) i, bend(i,tbend),tabbend(i,tbend)
             i=i+1
             goto 14
 15          bendsplits(tbend)=i-1
             if (myid.eq.0) write(iou,*) 'using linear interpolation for 
     &            1-3 nonbonded bending'
             do i=1,bendsplits(tbend)-1
                benddiff(i,tbend)=tabbend(i+1,tbend)-tabbend(i,tbend)
             enddo
          enddo
          close(42)
       endif

c   KM 12/02/08 - read in tabulated nonbonded potential
c   and set up linear interpolation

       if (L_vdW_table) then
          read(43,*) ntabvdW
          write(iou,*)
          mmm=1
          do mmm=1,ntabvdW
c            iinvdW and jjvdW are bead types
             read(43,*) iivdW, jjvdW
             read(43,*) num_int_vdW(iivdW,jjvdW)
             i=1
 16          read(43,*,end=17) rvdW(i,iivdW,jjvdW), 
     &            tabvdW(i,iivdW,jjvdW)
             if (rvdW(i,iivdW,jjvdW).eq.1000) goto 17
c             write(59,*) i, rvdW(i,iivdW,jjvdW),
c     &            tabvdW(i,iivdW,jjvdW)
             i=i+1
             goto 16
 17          vdWsplits(iivdW,jjvdW)=i-1
             if (myid.eq.0) write(iou,*) 'using linear interpolation for nonbonded 
     &            van der Waals interactions'
             do i=1,vdWsplits(iivdW,jjvdW)-1
                vdWdiff(i,iivdW,jjvdW)=
     &               tabvdW(i+1,iivdW,jjvdW)-
     &               tabvdW(i,iivdW,jjvdW)
             enddo
          enddo
          close(43)
       endif

c   KM 04/23/09 - read in tabulated electrostatic potential
c   and set up linear interpolation

       if (L_elect_table) then
          read(44,*) ntabelect
          write(iou,*)
          mmm=1
          do mmm=1,ntabelect
c            iielect and jjelect are bead types
             read(44,*) iielect, jjelect
c             write(59,*) iielect,jjelect
             read(44,*) num_int_elect(iielect,jjelect)
c             write(59,*) num_int_elect(iielect,jjelect)
             i=1
 18          read(44,*,end=19) relect(i,iielect,jjelect), 
     &            tabelect(i,iielect,jjelect)
c             write(59,*) relect(i,iielect,jjelect), 
c     &            tabelect(i,iielect,jjelect)
             if (relect(i,iielect,jjelect).eq.1000) goto 19
             i=i+1
             goto 18
 19          electsplits(iielect,jjelect)=i-1
             if (myid.eq.0) write(iou,*) 'using linear interpolation for
     &            electrostatic interactions'
             do i=1,electsplits(iielect,jjelect)-1
                electdiff(i,iielect,jjelect)=
     &               tabelect(i+1,iielect,jjelect)-
     &               tabelect(i,iielect,jjelect)
             enddo
          enddo
          close(44)
       endif

c -- check whether there is a polarizable molecule

      lpolar = .false.
      lqqelect = .false.

      do imol = 1, nmolty
         if (lflucq(imol)) lpolar = .true.
         if (lelect(imol)) lqqelect = .true.
      enddo

      if ( .not. lqqelect ) then
         if ( lewald .or. lchgall ) then
            write(iou,*) 'no charges in the system and turn off lewald'
            ldie = .true.
            return
         endif
      endif

      if ( .not. lpolar ) then
         if ( lanes  ) then
            write(iou,*) 'lanes should be false for nonpolarizable 
     &           systems!'
            ldie = .true.
            return
         endif
         if ( lfepsi ) then 
            write(iou,*) 'lfepsi should be false for nonpolarizable 
     &           systems!'
            ldie = .true.
            return
         endif
      endif

C This B(i) goes in the acceptance rules

      if (lgrand) then
        do i=1,nmolty
            debroglie = 17.458d0/( dsqrt(masst(i)/beta ))
            B(i) = exp(B(i)/temp)/(debroglie*debroglie*debroglie)
         enddo
      endif
      
c - read linkcell information
      read(4,*)
      read(4,*) licell,rintramax,boxlink

      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'licell:',licell
            write(iou,*) 'rintramax:',rintramax,' A'
            write(iou,*) 'boxlink:',boxlink
         else 
            write(iou,*) licell,rintramax,boxlink
         endif
      endif

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
      IF(boxlink .LE. nbxmax)THEN
         if (lsolid(boxlink).and.(.not.lrect(boxlink))) then

             write(iou,*) 
     &        'Linkcell not implemented for nonrectangular boxes'  
             ldie = .true.
             return

         endif
      ENDIF


c -- read the atomic displacements

      read(4,*)
      read(4,*) Armtrax, Armtray, Armtraz

      if(lecho.and.myid.eq.0) then
         if(lverbose) then
           write(iou,1019) 'initial maximum displacements for atoms:',
     &          Armtrax, Armtray, Armtraz
         else
           write(iou,*) Armtrax, Armtray, Armtraz
         endif
       endif


c - read displacement information
      read(4,*) 
      read(4,*) rmtrax(1,1),rmtray(1,1),rmtraz(1,1)
      do im = 1,nbox
         do imol = 1,nmolty
            rmtrax(imol,im) = rmtrax(1,1)
            rmtray(imol,im) = rmtray(1,1)
            rmtraz(imol,im) = rmtraz(1,1)
         enddo
      enddo
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,1019) ' initial maximum x, y and z displacement:',
     &           rmtrax(1,1),rmtray(1,1),rmtraz(1,1)
         else 
            write(iou,*) rmtrax(1,1), rmtray(1,1), rmtraz(1,1)
         endif
      endif


      read(4,*) 
      read(4,*) rmrotx(1,1),rmroty(1,1),rmrotz(1,1)
      do im = 1,nbox
         do imol = 1,nmolty
            rmrotx(imol,im) = rmrotx(1,1)
            rmroty(imol,im) = rmroty(1,1)
            rmrotz(imol,im) = rmrotz(1,1)
         enddo
      enddo
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,1019) ' initial maximum x, y and z rotation:    ',
     &           rmrotx(1,1),rmroty(1,1),rmrotz(1,1)
         else 
            write(iou,*) rmrotx(1,1), rmroty(1,1), rmrotz(1,1)
         endif
      endif
      read(4,*) 
      read(4,*) tatra,tarot
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'target translational acceptance ratio:',tatra
            write(iou,*) 'target rotational acceptance ratio:',tarot
         else 
            write(iou,*) tatra, tarot
         endif
      endif

c - read initial setup information
      read(4,*)
      read(4,*) linit,newtab,lreadq
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'linit:',linit
            write(iou,*) 'lnewtab:',newtab
            write(iou,*) 'lreadq:',lreadq
         else 
            write(iou,*) linit, newtab, lreadq
         endif
      endif
      read(4,*)
      read(4,*) (lbranch(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            do i = 1,nmolty
               write(iou,*) 'lbranch for molecule type',i,':',lbranch(i)
            enddo
         else 
            write(iou,*) (lbranch(i),i=1,nmolty)
         endif
      endif
      do i = 1, nbox
         read(4,*)
         read(4,*) (ininch(j,i),j=1,nmolty)
         if ( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               write(iou,*) 'box:',i
               do j = 1,nmolty
                  write(iou,*) '   initial number of chains of type',
     &                 j,':',ininch(j,i)
               enddo
            else 
               write(iou,*) 'box:',i,(ininch(j,i),j=1,nmolty)
            endif
         endif
         read(4,*)
         read(4,*) inix(i),iniy(i),iniz(i),inirot(i),inimix(i),
     &        zshift(i),dshift(i),nchoiq(i)
         if ( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               write(iou,1020) '    initial number of chains in x, y',
     &              ' and z directions:',inix(i),iniy(i),iniz(i)
               write(iou,*) '   initial rotational displacement:',
     &              inirot(i)
               write(iou,*) '   inimix:',inimix(i)
               write(iou,*) '   zshift:',zshift(i)
               write(iou,*) '   dshift:',dshift(i)
               write(iou,*) '   nchoiq:',nchoiq(i)
            else 
               write(iou,*) inix(i),iniy(i),iniz(i),inirot(i),
     &        inimix(i),zshift(i),dshift(i),nchoiq(i)
            endif
         endif
      enddo

c - read ensemble specific information
      read(4,*)
      read(4,*) rmvol(1), tavol, iratv, iratp, rmflcq(1,1), taflcq
      do zz = 1,nmolty
         do zzz = 1,nbox
            rmflcq(zz,zzz) = rmflcq(1,1)
         enddo
      enddo
      if ( lgibbs ) then
         do zzz = 2, nbox
            rmvol(zzz) = rmvol(1)
         enddo
      endif
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            do zz = 1,nbox
               write(iou,*) 'initial maximum volume displacement '
     &              ,'(rmvol) in box',zz,':',rmvol(zz)
            enddo
            write(iou,*) 'target volume acceptance ratio (tavol):',tavol
            write(iou,*) 'iratv:',iratv
            write(iou,*) 'iratp:',iratp
            do zz = 1,nmolty
               do zzz = 1,nbox
                  write(iou,*) 'initial maximum fluct. charge',
     &                 ' displ. for chain type',zz,' box',
     &                 zzz,':',rmflcq(zz,zzz)
               enddo
            enddo
            write(iou,*) 'target fluctuating charge acceptance ratio',
     &           ' (taflcq):',taflcq
         else
            write(iou,*) rmvol, tavol, iratv, iratp
            write(iou,*) 'rmflcq',
     &           ((rmflcq(zz,zzz),zzz=1,nbox),zz=1,nmolty),taflcq
         endif
      endif

      read(4,*)
      read(4,*) pmvol,(pmvlmt(j),j=1,nbox)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmvol:',pmvol
            do j = 1,nbox
               write(iou,*) '   pmvlmt for box',j,':',pmvlmt(j)
            enddo
         else 
            write(iou,*) 'pmvol',pmvol,(pmvlmt(j),j=1,nbox)
         endif
      endif
      read(4,*)
      read(4,*) nvolb,(pmvolb(j),j=1,nvolb)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'nvolb:',nvolb
            do j = 1,nvolb
               write(iou,*) '   pmvolb:',pmvolb(j)
            enddo
         else 
            write(iou,*) '   nvolb',nvolb,(pmvolb(j),j=1,nvolb)
         endif
      endif
      read(4,*)
      do j = 1,nvolb
         read(4,*) box5(j),box6(j)
         if ( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               write(iou,*) '   box pair for volume move number',j,':',
     &              box5(j),box6(j)
            else 
               write(iou,*) box5(j),box6(j)
            endif
         endif
      enddo
      
      lxyz = .false.
      do j = 1,nbox
         if (lsolid(j) .and. .not. lxyz) then
            lxyz = .true.
            read(4,*)
            read(4,*) pmvolx,pmvoly
            if (lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(iou,*) 'pmvolx:',pmvolx
                  write(iou,*) 'pmvoly:',pmvoly
               else
                  write(iou,*) pmvolx,pmvoly
               endif
            endif
         endif
      enddo
    
c     --- read swatch information
      read(4,*)
      read(4,*) pmswat,nswaty
      if ( nswaty .gt. npamax ) then
         write(iou,*) 'nswaty gt npamax'
         ldie = .true.
         return
      endif

      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmswat:',pmswat
            write(iou,*) '   number of swatch pairs (nswaty):',nswaty
         else 
            write(iou,*) 'pmswat, nswaty',pmswat,nswaty
         endif
      endif

      if (nswaty .gt. 0) then
         read(4,*)
         read(4,*) ((nswatb(i,j),j=1,2),i=1,nswaty)
         read(4,*)
         read(4,*) (pmsatc(i),i=1,nswaty)

c        --- safety checks on swatch
         do i = 1,nswaty
            if ( nswatb(i,1) .eq. nswatb(i,2) ) then
               write(iou,*) 'nswaty ',i,' has identical moltyp'
               write(iou,*) 'cannot swatch identical moltyp'
               ldie = .true.
               return
            endif
         enddo

         if( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               do i = 1,nswaty
                  write(iou,*) '   swatch molecule type pairs:',
     &                 (nswatb(i,j),j=1,2)
               enddo
               do i = 1,nswaty
                  write(iou,*) '   probability of each swatch pair:',
     &                 pmsatc(i)
               enddo
            else
               do i=1,nswaty
                  write(iou,*) 'nswatb', (nswatb(i,j),j=1,2)
                  write(iou,*) 'pmsatc',pmsatc(i)
               enddo
            endif
         endif
         
         do i=1,nswaty
c           --- number of beads that remain in the same position
            read(4,*)
            read(4,*) nsampos(i),(ncut(i,j),j=1,2)
            if (lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(iou,*) '   nsampos:',nsampos(i)
                  write(iou,*) '   ncut:',(ncut(i,j),j=1,2)
               else 
                  write(iou,*) 'nsampos',nsampos(i),' ncut',
     &                 (ncut(i,j),j=1,2)
               endif
            endif
c           --- bead number
            read(4,*)
            do j = 1,nsampos(i)
               read(4,*) (splist(i,j,k),k=1,2)
               if (lecho.and.myid.eq.0) then
                  if (lverbose) then
                     write(iou,*) '   splist:',
     &                    (splist(i,j,k),k=1,2)
                  else 
                     write(iou,*) 'splist',(splist(i,j,k),k=1,2)
                  endif
               endif
            enddo
            
            read(4,*)
            read(4,*) (( gswatc(i,j,k), k=1,2*ncut(i,j) ), j=1,2 )
!!            if (lecho.and.myid.eq.0) then
!!               if (lverbose) then
!!                  do zz = 1,ncut(i,j)
!!                     write(iou,*) '   grow from and prev for ncut',zz,':',
!!     &                    (gswatc(i,j,zz),j=1,2)
!!                  enddo
!!               else 
!!                  write(iou,*) 'gswatc',(( gswatc(i,j,k), 
!!     &                 k=1,2*ncut(i,j) ), j=1,2 )
!!               endif
!!            endif

            read(4,*) 
            read(4,*) nswtcb(i), (pmswtcb(i,ipair), ipair=1,nswtcb(i))
            read(4,*)
            do ipair = 1,nswtcb(i)
               read(4,*) box3(i,ipair),box4(i,ipair)
               if (lecho.and.myid.eq.0) then
                  if (lverbose) then
                     write(iou,*) '   box pair:',
     &                    box3(i,ipair),box4(i,ipair)
                  else
                     write(iou,*) box3(i,ipair),box4(i,ipair)
                  endif
               endif
            enddo
         enddo

      else
c        --- skip past all of the swatch info
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
         read(4,*)
      endif

c --- read swap info
      read(4,*)
      read(4,*) pmswap, (pmswmt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmswap:',pmswap
            do i = 1,nmolty
               write(iou,1021) '   swap probability for molecule type '
     &              ,'        ',i,' (pmswmt):',pmswmt(i)
            enddo
         else 
            write(iou,*) 'pmswap',pmswap,(pmswmt(i),i=1,nmolty)
         endif
      endif
      do i = 1, nmolty
         read(4,*)
         read(4,*) nswapb(i), (pmswapb(i,ipair),ipair=1,nswapb(i))
         if ( lecho.and.myid.eq.0 ) then
            if (lverbose) then
               write(iou,*) '   number of swap moves for molecule type',
     &              i,':',nswapb(i)
               do ipair = 1,nswapb(i)
                  write(iou,*) '      pmswapb:',pmswapb(i,ipair)
               enddo
            else 
               write(iou,*) nswapb(i), 
     &              (pmswapb(i,ipair),ipair=1,nswapb(i))
            endif
         endif
         read(4,*)
         do ipair = 1, nswapb(i)
            read(4,*) box1(i,ipair), box2(i,ipair)
            if ( lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(iou,*) '      box pair:',
     &                 box1(i,ipair), box2(i,ipair)
               else 
                  write(iou,*) 'box pair:',
     &                 box1(i,ipair), box2(i,ipair)
               endif
            endif
         enddo
      enddo

c --- read cbmc info
      read(4,*)
      read(4,*) pmcb, (pmcbmt(i),i=1,nmolty)
      read(4,*)
      read(4,*) (pmall(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmcb:',pmcb
            do i = 1,nmolty
               write(iou,1021) '   CBMC probability for molecule type '
     &              ,'        ',i,' (pmcbmt):',pmcbmt(i)
            enddo
            do i = 1,nmolty
               write(iou,*) '   pmall for molecule type',i,':',pmall(i)
            enddo
         else 
            write(iou,*) 'pmcb',pmcb,(pmcbmt(i),i=1,nmolty),
     &           'pmall',(pmall(i),i=1,nmolty)
         endif
      endif
      read(4,*)
      read(4,*) (pmfix(i),i=1,nmolty)
      if (lecho.and.myid.eq.0) then
         if (lverbose) then
            do i = 1,nmolty
               write(iou,1021) '   probability for SAFE-CBMC for '
     &              ,'molecule type',i,'  (pmfix):',pmfix(i)
            enddo
         else
            write(iou,*) (pmfix(i),i=1,nmolty)
         endif
      endif
      do i = 1, nmolty
         if (lring(i).and.pmfix(i).lt.0.999.and.
     &        .not.lrig(i)) then
            write(iou,*) 'a ring can only be used with safe-cbmc'
            ldie = .true.
            return
         endif
      enddo
c --- read fluctuating charge info
      read(4,*)
      read(4,*) pmflcq, (pmfqmt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmflcq:',pmflcq
            do i = 1,nmolty
               write(iou,1021) '   flcq probability for molecule type '
     &              ,'        ',i,' (pmfqmt):',pmfqmt(i)
            enddo
         else 
            write(iou,*) 'pmflcq',pmflcq,(pmfqmt(i),i=1,nmolty)
         endif
      endif
c --- read expanded-coefficient move info
      read(4,*)
      read(4,*) pmexpc, (pmeemt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmexpc:',pmexpc
            do i = 1,nmolty
               write(iou,1021) '   expanded ens. prob. for molecule '
     &              ,'type      ',i,' (pmeemt):',pmeemt(i)
            enddo
         else 
            write(iou,*) 'pmexpc',pmexpc,(pmeemt(i),i=1,nmolty)
         endif
      endif
c read new expanded ensemble info
      read(4,*)
      read(4,*) pmexpc1
c -- read atom translation probability
      read(4,*)
      read(4,*) pm_atom_tra        
c --- read translation info
      read(4,*)
      read(4,*) pmtra,(pmtrmt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmtra:',pmtra
            do i = 1,nmolty
               write(iou,1021) '   translation probability for molecule'
     &              ,'   type',i,' (pmtrmt):',pmtrmt(i)
            enddo
         else 
            write(iou,*) 'pmtra',pmtra,(pmtrmt(i),i=1,nmolty)
         endif
      endif

c --- read rotation info
      read(4,*)
      read(4,*) (pmromt(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'pmrot:',1.0d0
            do i = 1,nmolty
               write(iou,1021) '   rotational probability for molecule'
     &              ,'    type',i,' (pmromt):',pmromt(i)
            enddo
         else 
            write(iou,*) (pmromt(i),i=1,nmolty)
         endif
      endif

c *** writeout probabilities, in percent
      if (lverbose.and.myid.eq.0) then
         write(iou,*)
         write(iou,*) 'percentage move probabilities:'
         write(iou,'(1x,a19,f8.2,a2)') 'volume move       :',
     &      100.0d0*pmvol,' %'
         pcumu = pmvol
         if (pmswat .gt. pmvol) then
            pm = pmswat - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(iou,'(1x,a19,f8.2,a2)') 'swatch move       :',
     &      100.0d0*pm,' %'
         if (pmswap .gt. pmswat) then
            pm = pmswap - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(iou,'(1x,a19,f8.2,a2)') 'swap move         :',
     &      100.0d0*pm,' %'
         if (pmcb .gt. pmswap) then
            pm = pmcb - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(iou,'(1x,a19,f8.2,a2)') 'CBMC move         :',
     &      100.0d0*pm,' %'
         if (pmflcq .gt. pmcb) then
            pm = pmflcq - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(iou,'(1x,a19,f8.2,a2)') 'fluct charge move :',
     &      100.0d0*pm,' %'
         if (pmexpc .gt. pmflcq) then
            pm = pmexpc - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(iou,'(1x,a19,f8.2,a2)') 'expanded ens move :',
     &      100.0d0*pm,' %'
         if (pmtra .gt. pmexpc) then
            pm = pmtra - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(iou,'(1x,a19,f8.2,a2)') 'translation move  :',
     &      100.0d0*pm,' %'
         pm = 1.0d0 - pmtra
         write(iou,'(1x,a19,f8.2,a2)') 'rotation move     :',
     &      100.0d0*pm,' %'

         write(iou,*)
         write(iou,*) 'Fraction of atom translations move', pm_atom_tra            

      endif
      
 

c --- read growth details
      read(4,*)
      read(4,*) (nchoi1(i),i=1,nmolty),(nchoi(i),i=1,nmolty)
     &     ,(nchoir(i),i=1,nmolty)
     &     ,(nchoih(i),i=1,nmolty),(nchtor(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*)
            write(iou,*) 'molecule type :',(i,'  ',i=1,nmolty)
            write(iou,*) '     nchoi1   :',(nchoi1(i),' ',i=1,nmolty)
            write(iou,*) '     nchoi    :',(nchoi(i),' ',i=1,nmolty)
            write(iou,*) '     nchoir   :',(nchoir(i),' ',i=1,nmolty)
            write(iou,*) '     nchoih   :',(nchoih(i),' ',i=1,nmolty)
            write(iou,*) '     nchtor   :',(nchtor(i),i=1,nmolty)
         else
            write(iou,*) 'nchoi1',(nchoi1(i),i=1,nmolty)
            write(iou,*) 'nchoi',(nchoi(i),i=1,nmolty)
            write(iou,*) 'nchoir',(nchoir(i),i=1,nmolty)
            write(iou,*) 'nchoih',(nchoih(i),i=1,nmolty)
            write(iou,*) 'nchtor',(nchtor(i),i=1,nmolty)
         endif
      endif

      do imol = 1, nmolty
         if (pmfix(imol).gt.0) then
            read(23,*) counttot
c     --- read in from fort.23 the bead-bead distribution
            do i = 1, iring(imol) 
               do j = 1, iring(imol)
                  if (i.eq.j) goto 110
                  do bin = 1, maxbin
                     read(23,*) bdum,probf(i,j,bin)
                     hist(i,j,bin) = 0
                  enddo                 
 110              continue
               enddo
            enddo
         endif
      enddo

c     --- read information for CBMC bond angle growth
      read(4,*) 
      read(4,*) (nchbna(i),i=1,nmolty),(nchbnb(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            do i = 1,nmolty
               write(iou,*) 'nchbna and nchbnb for molecule type',i,':',
     &              nchbna(i),nchbnb(i)
            enddo
         else
            write(iou,*) 'nchbna ',(nchbna(i),i=1,nmolty)
            write(iou,*) 'nchbnb ',(nchbnb(i),i=1,nmolty)
         endif
      endif

      read(4,*)
      read(4,*) (icbdir(i),i=1,nmolty),(icbsta(i),i=1,nmolty)
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            do i = 1,nmolty
               write(iou,*) 'icbdir for molecule type',i,':',icbdir(i)
            enddo
            do i = 1,nmolty
               write(iou,*) 'icbsta for molecule type',i,':',icbsta(i)
            enddo
         else
            write(iou,*) 'icbdir',(icbdir(i),i=1,nmolty)
            write(iou,*) 'icbsta',(icbsta(i),i=1,nmolty)
         endif
      endif

c     ---- Error checking
      do i=1,nmolty
         if ( nchoi1(i) .gt. nchmax ) then
            write(iou,*) 'nchoi1 gt nchmax'
            ldie = .true.
            return
         endif
         if (nchoi(i) .gt. nchmax ) then
            write(iou,*) 'nchoi gt nchmax'
            ldie = .true.
            return
         endif
         if ( nchoih(i) .ne. 1 .and. nunit(i) .eq. nugrow(i) ) 
     &        then
            write(iou,*) ' nchoih must be one if nunit = nugrow'
            ldie = .true.
            return
         endif
         if ( nchtor(i) .gt. nchtor_max ) then
            write(iou,*) 'nchtor gt nchtor_max'
            ldie = .true.
            return
         endif
         if ( nchbna(i) .gt. nchbn_max ) then
            write(iou,*) 'nchbna gt nchbn_max'
            ldie = .true.
            return
         endif
         if ( nchbnb(i) .gt. nchbn_max ) then
            write(iou,*) 'nchbnb gt nchbn_max'
            ldie = .true.
            return
         endif
         if ( icbsta(i) .gt. numax ) then
            write(iou,*) 'icbsta gt numax'
            ldie = .true.
            return
         endif
      enddo

c --- read exclusion table for intermolecular interactions
      read(4,*)
      read(4,*) nexclu
      do i = 1, nmolty
         do ii = 1, numax
            do j = 1, nmolty
               do jj = 1, numax
                  lexclu(i,ii,j,jj) = .false.
               enddo
            enddo
         enddo
      enddo       
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*)'nexclu:',nexclu
         else 
            write(iou,*) nexclu
         endif
      endif
      if (nexclu .ne. 0) then
         do ndum = 1, nexclu
            read(4,*) i, ii, j, jj
            lexclu(i,ii,j,jj) = .true.
            lexclu(j,jj,i,ii) = .true.
            if (lverbose.and.myid.eq.0) then
               write(iou,*) 'excluding interactions between bead',ii,
     &              ' on chain',i,' and bead',jj,' on chain',j
            endif
         enddo
      else
         read(4,*)
      endif

c --- read inclusion list 
      read(4,*)
      read(4,*) inclnum
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'inclnum:',inclnum
         else 
            write(iou,*) inclnum
         endif
      endif
      if (inclnum .ne. 0) then
         do ndum = 1, inclnum
            read(4,*) inclmol(ndum),inclbead(ndum,1),inclbead(ndum,2)
     &           ,inclsign(ndum),ofscale(ndum),ofscale2(ndum)
            if (lecho.and.myid.eq.0 ) then
               if (lverbose) then
                  if (inclsign(ndum) .eq. 1) then
                     write(iou,*) 'including intramolecular '
     &                    ,'interactions for chain type',inclmol(ndum),
     &                    ' between beads',inclbead(ndum,1),' and',
     &                    inclbead(ndum,2)
                  else
                     write(iou,*) 'excluding intramolecular '
     &                    ,'interactions for chain type',inclmol(ndum),
     &                    ' between beads',inclbead(ndum,1),' and',
     &                    inclbead(ndum,2)

                  endif
               else
                  write(iou,*) inclmol(ndum),inclbead(ndum,1)
     &                 ,inclbead(ndum,2),inclsign(ndum),ofscale(ndum)
     &		       ,ofscale2(ndum)	
               endif
            endif
         enddo
      else
         read(4,*)
      endif

c --- read a15 inclusion list 
      read(4,*)
      read(4,*) ainclnum
      if ( lecho.and.myid.eq.0 ) then
         if (lverbose) then
            write(iou,*) 'ainclnum:',ainclnum
         else 
            write(iou,*) ainclnum
         endif
      endif
      if (ainclnum .ne. 0) then
         do ndum = 1, ainclnum
            read(4,*) ainclmol(ndum),ainclbead(ndum,1),ainclbead(ndum,2)
     &           ,a15t(ndum)
            if (lecho.and.myid.eq.0) then
               if (lverbose) then
                  write(iou,*) 'repulsive 1-5 OH interaction for',
     &                 ' chain type',ainclmol(ndum),' between beads',
     &                 ainclbead(ndum,1),' and',ainclbead(ndum,2),
     &                 ' of type:',a15t(ndum)
               else
                  write(iou,*) ainclmol(ndum),ainclbead(ndum,1)
     &                 ,ainclbead(ndum,2),a15t(ndum)
               endif
            endif
         enddo
      else
         read(4,*)
      endif

c - set up the inclusion table
      call inclus( inclnum,inclmol,inclbead,inclsign,ncarbon,
     &     ainclnum,ainclmol,ainclbead,a15t,ofscale,ofscale2)

c -  read in information on the chemical potential checker
      read(4,*)
      read(4,*) lucall
      read(4,*)
      read(4,*) (ucheck(jj),jj=1,nmolty)
      if (lecho.and.myid.eq.0) then
         if (lverbose) then
            write(iou,*) 'lucall:',lucall
            do jj = 1,nmolty
               write(iou,*) '   ucheck for molecule type',
     &              jj,':',ucheck(jj)
            enddo
         else 
            write(iou,*) lucall,(ucheck(jj),jj=1,nmolty)
         endif
      endif

c -   read information for virial coefficient calculation
      read(4,*) 
      read(4,*) nvirial,starvir,stepvir
      if (lecho.and.myid.eq.0) then
         if (lverbose) then
            write(iou,*) 'nvirial:',nvirial
            write(iou,*) 'starvir:',starvir
            write(iou,*) 'stepvir:',stepvir
         else 
            write(iou,*) nvirial,starvir,stepvir
         endif
      endif

      if (lvirial) then
         if ( nvirial .gt. maxvir ) then
            write(iou,*) 'nvirial .gt. maxvir'
            ldie = .true.
            return
         endif

         read(4,*)
         read(4,*) ntemp,(virtemp(jj),jj=1,ntemp)
         if (lecho.and.myid.eq.0) then
            if (lverbose) then
               write(iou,*) 'ntemp:',ntemp
               write(iou,*) 'calculation of virial coefficient ',
     &              'at the following temperatures:',
     &              (virtemp(jj),jj=1,ntemp)
            else 
               write(iou,*) ntemp,
     &              'Calculation of virial coefficient ',
     &              'at the following temperatures:',
     &              (virtemp(jj),jj=1,ntemp)
            endif
         endif
      endif

c --- JLR 11-11-09 
c --- reading in extra variables for RPLC simulations
c --- if lideal=.true. then intermolecular interactions are not computed    
c --- if ltwice=.true. then mimage is applied twice  
c --- if lrplc=.true. there are some special rules in CBMC for how to grow chains  
      read(4,*)
      read(4,*) (lideal(i),i=1,nbox)  
      if (lecho.and.myid.eq.0)  then
         write(iou,*) 'lideal: ', (lideal(i),i=1,nbox)  
      endif
      do i = 1,nbox 
         if (lideal(i) .and. lexpee) then
            write(iou,*) 'cannot have lideal and lexpee both true'
            write(iou,*) 'If you want this you will have change code'
            ldie = .true.
            return
         endif
      enddo
      read(4,*)
      read(4,*) (ltwice(i),i=1,nbox)
      read(4,*)
      read(4,*) (lrplc(i),i=1,nmolty)
      if (lecho.and.myid.eq.0) then
         write(iou,*) 'ltwice: ', (ltwice(i),i=1,nbox) 
         write(iou,*) 'lrplc: ', (lrplc(i),i=1,nmolty)
      endif
c --- END JLR 11-11-09                         

c--- JLR 11-11-09
c--- skip this part if analysis will not be done (if ianalyze .gt. nstep)
c KM 01/10 remove analysis 
c      if (ianalyze.lt.nstep) then
c         nhere = 0
c         do zz=1,nntype
c            if ( lhere(zz) ) then
c               nhere = nhere + 1
c               temphe(nhere) = zz
c               beadtyp(nhere)=zz
c            endif
c         enddo
c         do zz = 1,nhere 
c            atemp = temphe(zz)
c            decode(atemp) = zz
c         enddo
c      endif
c --- END JLR 11-11-09


C -------------------------------------------------------------------
 
c *** restart saver
      if ( nstep .gt. 100 ) then
         irsave = nstep / 20
      else
         irsave = nstep + 1
      endif

c -------------------------------------------------------------------
 
c * check information of .INC-files 
      if (myid.eq.0) then
         write(iou,*)
         write(iou,*) '***** program   =  THE MAGIC BLACK BOX'
         iensem = 0
         if ( lgrand ) then
            iensem = iensem + 1
            write(iou,*) 'grand-canonical ensemble'
         endif
         if ( lvirial ) then
            write(iou,*) 'Computing Second Virial Coefficient'
            if ( nchain .ne. 2) then
               write(iou,*) 'nchain must equal 2'
               ldie = .true.
               return
            endif
         endif
         if ( lgibbs ) then
            iensem = iensem + 1
            if ( lnpt ) then
               write(iou,*) 'NPT Gibbs ensemble'
            else
               write(iou,*) 'NVT Gibbs ensemble'
            endif
         endif
         if ( iensem .eq. 0 ) then
            if ( lnpt ) then
               write(iou,*) 'Isobaric-isothermal ensemble'
               iensem = iensem + 1
            else
               write(iou,*) 'Canonical ensemble'
               iensem = iensem + 1
            endif
         endif
         if ( iensem .gt. 1 ) then
            write(iou,*) 'INCONSISTENT ENSEMBLE SPECIFICATION'
            ldie = .true.
            return
         endif
         
         inpbc = 0
         if ( lpbc ) then
            write(iou,*) 'using periodic boundaries'
            if ( lpbcx ) then
               inpbc = inpbc + 1
               write(iou,*) 'in x-direction'
            endif
            if ( lpbcy ) then
               inpbc = inpbc + 1
               write(iou,*) 'in y-direction'
            endif
            if ( lpbcz ) then
               inpbc = inpbc + 1
               write(iou,*) 'in z-direction'
            endif
            if ( inpbc .eq. 0 ) then
               write(iou,*) 'INCONSISTENT PBC SPECIFICATION'
               ldie = .true.
               return
            endif
            write(iou,*) inpbc,'-dimensional periodic box'
         else
            write(iou,*) 'cluster mode (no pbc)'
            if ( lgibbs .or. lgrand ) then
               write(iou,*) 'INCONSISTENT SPECIFICATION OF LPBC 
     &              AND ENSEMBLE'
               ldie = .true.
               return
            endif
         endif
         
         if ( lfold ) then
            write(iou,*) 'particle coordinates are folded 
     &           into central box'
            if ( .not. lpbc ) then 
               write(iou,*) 'INCONSISTENT SPECIFICATION OF LPBC AND 
     &              LFOLD'
               ldie = .true.
               return
            endif
         endif
         
         if ( lijall ) then
            write(iou,*) 'all i-j-interactions are considered',
     &           '  (no potential cut-off)'
         endif
         
         if ( lchgall ) then
            write(iou,*) 'all the inter- and intramolecular Coulombic',
     &           ' interactions are considered (no group-based cutoff)'
         endif
         if ( lcutcm ) then
            write(iou,*) 'additional (COM) cutoff on computed rcmu'
            if ( lijall ) then
               write(iou,*) 'cannot have lijall with lcutcm'
               ldie = .true.
               return
            endif
c         if ( lchgall ) stop 'cannot have lchgall with lcutcm'
         endif
         if ( ldual ) then
            write(iou,*) 'Dual Cutoff Configurational-bias Monte Carlo'
         endif
         
         write(iou,*) 
     &        'CBMC simultaneously grows all beads conected to the 
     &        same bead'
         write(iou,*) 'with bond angles generated from Gaussian',
     &        ' distribution'
         
         if ( ljoe ) then
            write(iou,*) 'external 12-3 potential for SAM (Joe Hautman)'
         endif
         if (lslit) then
            write(iou,*) '10-4-3 graphite slit pore potential'
         endif
         if ( lsami ) then
            write(iou,*) 'external potential for Langmuir films (Sami)'
            write(iou,*) 'WARNING: LJ potential defined in SUSAMI'
            write(iou,*) 'WARNING: sets potential cut-off to 2.5sigma'
            write(iou,*) 'WARNING: has build-in tail corrections'
            if ( ltailc ) then
               write(iou,*) 'INCONSISTENT SPECIFICATION OF LTAILC AND 
     &              LSAMI'
               ldie = .true.
               return
            endif
         endif
         
         if ( lexzeo ) then
            write(iou,*) 'external potential for zeolites (Berend)'
            write(iou,*) 'WARNING: potential defined in SUZEO'
         endif
         
         write(iou,*) 'Program will call Explct.f if needed'
         
         if ( lexpsix ) then
            write(iou,*) 'Exponential-6 potential for bead-bead'
         elseif ( lmmff ) then
            write(iou,*) 'Buffered 14-7 potential for bead-bead'
         elseif (lninesix) then
            write(iou,*) '9-6 potential for bead-bead'
         elseif (lgenlj) then
            write(iou,*) 'Generalized Lennard Jones potential
     &           for bead-bead'
         elseif (lgaro) then
            write(iou,*) 'Feuston-Garofalini potential for bead-bead'
         else
            write(iou,*) 'Lennard-Jones potential for bead-bead'
         endif
         
         if ( ltailc ) then
            write(iou,*) 'with additional tail corrections'
            if (lshift) then
               write(iou,*) ' lshift.and.ltailc!'
               ldie = .true.
               return
            endif
         endif
         
         if ( lshift) then
            write(iou,*) 'using a shifted potential'
            if (ltailc) then
               write(iou,*) ' lshift.and.ltailc!'
               ldie = .true.
               return
            endif
         endif
         
         write(iou,*) 'Coulombic inter- and intramolecular interactions'
         if ( lewald ) write(iou,*)
     &        'Ewald-sum will be used to calculate Coulombic 
     &        interactions'
         write(iou,*)
         write(iou,*) 'MOLECULAR MASS:', (masst(i),i=1,nmolty)
      endif

C -------------------------------------------------------------------

c *** read/produce initial/starting configuration ***
c *** zeolite external potential
      if ( lexzeo ) then
         call suzeo(rcut(1),newtab)
         boxlx(1)=zeorx
         boxly(1)=zeory
         boxlz(1)=zeorz
         if (myid.eq.0) write(iou,*) ' note zeolite determines 
     &        the box size !'
      endif

C ----------------------------------------------------------------      
      if (lgrand .and. .not.(lslit)) then
c ---     volume ideal gas box is set arbitry large!
         boxlx(2)=1000*boxlx(1)          
         boxly(2)=1000*boxly(1)
         boxlz(2)=1000*boxlz(1)
      endif
      if (linit) then
         do ibox = 1,nbox
            if (lsolid(ibox) .and. .not. lrect(ibox)) then
               write(iou,*) 'Cannot initialize non-rectangular system'
               ldie = .true.
               return
            endif
         enddo
      endif

      if ( linit ) then
         call initia(qelect)
         nnstep = 0
      else
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM needed for long input records, like
c nboxi and moltyp for a lot of molecules
         OPEN(77,FILE='fort.77',RECL=4096)
         read (77,*) nnstep
         read(77,*) Armtrax, Armtray, Armtraz
         do im = 1,nbox
            do imol = 1,nmolty
               read(77,*) rmtrax(imol,im), rmtray(imol,im)
     &              , rmtraz(imol,im)
               read(77,*) rmrotx(imol,im), rmroty(imol,im)
     &              , rmrotz(imol,im)
            enddo
         enddo
         if (myid.eq.0) then 
            write(iou,*) 
     &           'new maximum displacements read from restart-file'
            do im = 1,nbox
               if (myid.eq.0) write(iou,*)'box      #',im
               do imol = 1,nmolty
                  write(iou,*) 'molecule type',imol
                  write(iou,1101) rmtrax(imol,im), rmtray(imol,im)
     &                 , rmtraz(imol,im)
                  write(iou,1102) rmrotx(imol,im), rmroty(imol,im)
     &                 , rmrotz(imol,im)
               enddo
            enddo
         endif

         do im=1,nbox
            read (77,*) (rmflcq(i,im),i=1,nmolty)
         enddo
         if ( lecho.and.myid.eq.0) then
            do im=1,nbox
               write(iou,*) 'maximum fluc q displacements: Box #',im
               write(iou,*) (rmflcq(i,im),i=1,nmolty)
            enddo
         endif
c -- changed so fort.77 the same for all ensembles
c -- 06/08/09 KM
c         if ( lgibbs .or. lgrand .or. lnpt ) then
         read (77,*) (rmvol(ibox), ibox = 1,nbox)
         if (myid.eq.0) then
            write(iou,1103) (rmvol(ibox), ibox = 1,nbox)
            write(iou,*)
            write(iou,*) 
         endif
         do ibox = 1,nbox
            if (lsolid(ibox) .and. .not. lrect(ibox)) then
               read(77,*) (rmhmat(ibox,j),j=1,9)
               if (myid.eq.0) write(iou,*) (rmhmat(ibox,j),j=1,9)
            endif
         enddo
         
         if (lexzeo) then
            do ibox = 1,nbox
               read(77,*) dum,dum,dum
            enddo
         else
            if (myid.eq.0) then
               write(iou,*) 
               write(iou,*) 'new box size read from restart-file'
               write(iou,*)
            endif
            
            do ibox = 1,nbox
               
               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  
                  read(77,*) (hmat(ibox,j),j=1,9)
                  if (myid.eq.0) then
                     write(iou,*)
                     write(iou,*) 'HMAT COORDINATES FOR BOX',ibox
                     write(iou,*) (hmat(ibox,j),j=1,3)
                     write(iou,*) (hmat(ibox,j),j=4,6)
                     write(iou,*) (hmat(ibox,j),j=7,9)
                     write(iou,*)
                  endif
                  
                  call matops(ibox) 
                  
                  if (myid.eq.0) then
                     write(iou,*)    
                     write(iou,1106) ibox,min_width(ibox,1),min_width
     &                    (ibox,2),min_width(ibox,3)
                  endif

                  w(1) = min_width(ibox,1)
                  w(2) = min_width(ibox,2)
                  w(3) = min_width(ibox,3)  
                  
                  if (rcut(ibox)/w(1) .gt. 0.5d0 .or. 
     &                 rcut(ibox)/w(2) .gt. 0.5d0 .or. 
     &                 rcut(ibox)/w(3) .gt. 0.5d0) then
                     write(iou,*) 'rcut > half cell width'
                     ldie = .true.
                     return
                  endif
                  
                  if (myid.eq.0) then
                     write(iou,*)
                     write(iou,*) 'ibox:  ', ibox
                     write(iou,'("cell length |a|:",2x,f12.3)')
     &                    cell_length(ibox,1)
                     write(iou,'("cell length |b|:",2x,f12.3)') 
     &                    cell_length(ibox,2)
                     write(iou,'("cell length |c|:",2x,f12.3)') 
     &                    cell_length(ibox,3) 
                     
                     write(iou,*)
                     write(iou,'("cell angle alpha:",2x,f12.3)') 
     &                    cell_ang(ibox,1)*180.0d0/onepi
                     write(iou,'("cell angle beta: ",2x,f12.3)')
     &                    cell_ang(ibox,2)*180.0d0/onepi
                     write(iou,'("cell angle gamma:",2x,f12.3)')
     &                    cell_ang(ibox,3)*180.0d0/onepi
                   
c                     write(iou,1105) ibox,cell_ang(ibox,1)*180.0d0/onepi,
c     &                cell_ang(ibox,2)*180.0d0/onepi,cell_ang(ibox,3)
c     &                *180/onepi     
                  endif
               else

                  read (77,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
                  if (myid.eq.0) then
                     write(iou,*)
                     write(iou,*) 
                     write(iou,1104) ibox,
     &                    boxlx(ibox),boxly(ibox),boxlz(ibox)
                  endif

                  do i = 1, nbox
                     if( (rcut(i)/boxlx(i) .gt. 0.5d0).or.
     &                    (rcut(i)/boxly(i) .gt. 0.5d0).or.
     &                    (rcut(i)/boxlz(i) .gt. 0.5d0)) then
                        write(iou,*) 'rcut > 0.5*boxlx'
                        ldie = .true.
                        return
                     endif
                  enddo
                  
                  
                  
               endif
               
            enddo
         endif

c         endif ! end if ( lgibbs .or. lgrand .or. lnpt )
         if (myid.eq.0) then
            write(iou,*)
            write(iou,*) 'Finished writing simulation box related info'
         endif

         read (77,*) ncres
         read (77,*) nmtres
c --- check that number of particles in fort.4 & fort.77 agree ---
         if ( ncres .ne. nchain .or. nmtres .ne. nmolty ) then
            write(iou,*)
     +         'conflicting information in restart and control files'
            write(iou,*) 'nchain',nchain,'ncres',ncres
            write(iou,*) 'nmolty',nmolty,'nmtres',nmtres
            ldie = .true.
            return
         endif
         read (77,*) (nures(i),i=1,nmtres)

c         do i = 1, nmtres
c            if ( nures(i) .ne. nunit(i) ) then
c               write(iou,*)
c     +           'conflicting information in restart and control files'
c               write(iou,*) 'unit',i,'nunit',nunit(i),'nures',nures(i)
c               stop
c            endif
c         enddo
c         write(iou,*) 'ncres',ncres,'   nmtres',nmtres

         read (77,*) (moltyp(i),i=1,ncres)
         read (77,*) (nboxi(i),i=1,ncres)
         if ( lee ) then
            do i = 1, nmtres
               if ( lexpand(i) ) read(77,*) eetype(i)
            enddo
            do i = 1, nmtres
               if ( lexpand(i) ) read(77,*) rmexpc(i)
            enddo
         endif

c         write(iou,*) 'start reading coordinates'
c --- check that particles are in correct boxes ---
c --- obtain ncmt values
         do ibox = 1,nbox

            nchbox(ibox) = 0
 
         enddo
         do i = 1, nmolty
            do ibox = 1,nbox
               ncmt(ibox,i) = 0
            enddo
            if ( lexpand(i) ) then

c ??? problem in expand ensemble

               do j = 1, numcoeff(i)
                  do ibox = 1,2
                     ncmt2(ibox,i,j) = 0
                  enddo
               enddo
            endif
         enddo


         do i = 1, ncres

            if( nboxi(i) .le. nbox ) then
               ibox = nboxi(i)
               if ( ibox .ne. 1 .and. .not. lgibbs 
     &              .and. .not.lgrand) then 
                  write(iou,*) 'Particle found outside BOX 1'
                  ldie = .true.
                  return
               endif
               nchbox(ibox) = nchbox(ibox) + 1
               imolty = moltyp(i)
               ncmt(ibox,imolty) = ncmt(ibox,imolty) + 1
               if ( lexpand(imolty) ) then
                  if ( ibox .gt. 2 ) then
                     write(iou,*) 'put in box 1 and 2 for such 
     &                    molecules'
                     ldie = .true.
                     return
                  endif
                  itype = eetype(imolty)
                  ncmt2(ibox,imolty,itype) = 
     &                 ncmt2(ibox,imolty,itype) + 1
                  do j = 1,nunit(imolty)
                     sigma(imolty,j) = sigm(imolty,j,itype)
                     epsilon(imolty,j) = epsil(imolty,j,itype)
                  enddo
               endif 
            else
               write(iou,*) 'i:',i,'nboxi(i)',nboxi(i)
               write(iou,*) 'Particle found in ill-defined box'
               ldie = .true.
               return
            endif

         enddo

c --- check that number of particles of each type is consistent
         do i = 1, nmolty
            tcount = 0
            do ibox = 1, nbox
               tcount = tcount + ncmt(ibox,i)
            enddo
            if ( tcount .ne. temtyp(i) ) then
               write(iou,*) 'Particle type number inconsistency'
               write(iou,*) 'type',i
               write(iou,*) 'ncmt',(ncmt(ibox,i), ibox = 1,nbox)
               write(iou,*) 'temtyp', temtyp(i)
               ldie = .true.
               return
            endif
         enddo

c         write(iou,*) 'particles found in correct box with correct type'

         do i = 1,nbxmax
            qbox(i) = 0.0d0
         enddo

         do 22 i = 1, nchain
c            write(iou,*) 'reading coord of chain i',i
            imolty = moltyp(i)
            do j = 1, nunit(imolty)
               read (77,*) rxu(i,j), ryu(i,j), rzu(i,j),qqu(i,j)
               if (.not. lreadq) then
                  qqu(i,j) = qelect(ntype(imolty,j))
               endif
               qbox(nboxi(i)) = qbox(nboxi(i)) + qqu(i,j)
            enddo
            
 22      continue

         close(77)

          do i = 1,nbxmax
            if ( dabs(qbox(i)) .gt. 1d-6 ) then
               write(iou,*) 'box',i,' has a net charge of',qbox(i)
               ldie = .true.
               return
            endif
         enddo

      endif
 
      if (L_Ewald_Auto) then
        if ( lchgall ) then
           if (lewald) then
               do ibox = 1,nbox
                  if (lsolid(ibox).and.(.not.lrect(ibox))) then
                     min_boxl = min(min_width(ibox,1),min_width(ibox,2),
     &                       min_width(ibox,3))
                  else
                     min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
                  endif
                  kalp(ibox) = 6.40d0
                  calp(ibox) = kalp(ibox)/min_boxl
               enddo 
           else
              write(iou,*) 'lewald should be true when lchgall is true'
              ldie = .true.
              return
           endif
        else
           do ibox = 1, nbox 
              if (lsolid(ibox).and.(.not.lrect(ibox))) then
                  min_boxl = min(min_width(ibox,1),min_width(ibox,2),
     &                       min_width(ibox,3))
              else
                  min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
              endif
!              rcut(ibox) = 0.4d0*min_boxl
              calp(ibox) = 3.2d0/rcut(ibox) 
           enddo 
        endif
      else
        if ( lchgall ) then
c *** if real space term are summed over all possible pairs in the box
c *** kalp(1) & kalp(2) are fixed while calp(1) & calp(2) change
c *** according to the boxlength, kalp(1) should have a value greater than
c *** 5.0
          do ibox = 1, nbox
             if (lsolid(ibox).and.(.not.lrect(ibox))) then
                  min_boxl = min(min_width(ibox,1),min_width(ibox,2),
     &                       min_width(ibox,3))
              else
                  min_boxl = min(boxlx(ibox),boxly(ibox),boxlz(ibox))
              endif
             calp(ibox) = kalp(ibox)/min_boxl
             if ( lewald ) then
                if ( kalp(ibox) .lt. 5.6d0 ) then
                   write(iou,*) 'Warning, kalp is too small'
                   ldie = .true.
                   return
                endif
             else
ckea
                if(.not. lgaro) then
                   write(iou,*) 'lewald should be true when lchgall 
     &                  is true'
                   ldie = .true.
                   return
                endif
             endif
          enddo
        else
c *** if not lchgall, calp(1) & calp(2) are fixed
         do ibox = 1, nbox
            calp(ibox) = kalp(ibox)
            if ( lewald ) then
               if (calp(ibox)*rcut(ibox).lt.3.2d0.and.myid.eq.0) then
                  write(iou,*) 'Warning, kalp too small in box',ibox
                  write(iou,*) ibox,calp(ibox),rcut(ibox)
ccc --- JLR 11-24-09  
ccc --- you may want a smaller kalp, e.g. when comparing to previous work
ccc --- This does not need to be a stop  
c                  stop 'kalp is too small, set to 3.2/rcutchg'
                  write(iou,*) 'kalp is too small, set to 3.2/rcutchg'
ccc --- END JLR 11-24-09
               endif
            endif
         enddo
        endif
      endif

      if (L_Ewald_Auto) then
        do ibox = 1,nbox
           if ( (.not. lsolid(ibox)) .or. lrect(ibox) )  then
              k_max_l(ibox) = dint(boxlx(ibox)*calp(ibox))+1
              k_max_m(ibox) = dint(boxly(ibox)*calp(ibox))+1
              k_max_n(ibox) = dint(boxlz(ibox)*calp(ibox))+1
           else
              k_max_l(ibox) = dint(hmat(ibox,1)*calp(ibox))+2
              k_max_m(ibox) = dint(hmat(ibox,5)*calp(ibox))+2
              k_max_n(ibox) = dint(hmat(ibox,9)*calp(ibox))+2
           endif
        enddo  
      endif

      if (L_Ewald_Auto.and.myid.eq.0) then
         write(iou,*) 
         write(iou,*) '****Ewald Parameters*****'   
         write(iou,*) 'ibox   calp(ibox)  kmaxl(ibox)   kmaxm(ibox)',
     &                       '   kmaxn(ibox)   rcut(ibox)'   
         do ibox = 1,nbox
            write(iou,'(i4,5x,f12.6,3i12,12x,f12.4)') ibox, calp(ibox), 
     &                 k_max_l(ibox),k_max_m(ibox),
     &                 k_max_n(ibox),rcut(ibox) 
         enddo 
         write(iou,*)
      endif

c--- book keeping arrays
      do ibox=1,nbox
         do imolty=1,nmolty
            idummy(imolty)=0
         enddo      
         do i = 1, nchain
            if ( nboxi(i) .eq. ibox ) then
               imolty=moltyp(i)
               idummy(imolty) = idummy(imolty)+1
               parbox(idummy(imolty),ibox,imolty)=i
c               pparbox(i,ibox) = nmcount
            endif
         enddo
         nmcount = 0
         do imolty=1,nmolty
            nmcount = nmcount + idummy(imolty)
         enddo
         if ( nmcount .ne. nchbox(ibox) ) then
            write(iou,*) 'Readdat: nmcount ne nchbox', nmcount, nchbox
            ldie = .true.
            return
         endif
      enddo
c     set idummy counter to 0
      do imolty=1,nmolty
         idummy(imolty)=0
      enddo
c     set up parall
      do i =1,nchain
         imolty = moltyp(i)
         idummy(imolty) = idummy(imolty) + 1
         parall(imolty,idummy(imolty)) = i
      enddo

C -------------------------------------------------------------------
 
c - reordering of numbers for charmm
c      do i = 1, nchain
c         imolty = moltyp(i)
c         temtyp(i) = imolty
c         do j = 1, nunit(imolty)
c            temx(i,j) = rxu(i,j)
c            temy(i,j) = ryu(i,j)
c            temz(i,j) = rzu(i,j)
c         enddo
c      enddo
c      innew = 0
c      do it = 1, nmolty
c         do i = 1, nchain
c            imolty = temtyp(i)
c            if ( imolty .eq. it ) then
c               innew = innew + 1
c               moltyp(innew) = it
c               do j = 1, nunit(imolty)
c                  rxu(innew,j) = temx(i,j)
c                  ryu(innew,j) = temy(i,j)
c                  rzu(innew,j) = temz(i,j)
c               enddo
c            endif
c         enddo
c         write(iou,*) 'it =',it,'   innew =',innew
c      enddo

C -------------------------------------------------------------------

c --- set the centers of mass if LFOLD = .TRUE.

      if ( lfold ) then
         do ibox = 1, nbox
            call ctrmas(.true.,ibox,0,6)
         enddo
      endif

c * check that rintramax is really valid
      if (licell) then
         do i = 1,nchain
            if (2.0d0*rcmu(i) .gt. rintramax) then
               write(iou,*) 'rintramax for the linkcell list too small'
               ldie = .true.
               return
            endif
         enddo
      endif

c * calculate number of frames *
      nnframe = nstep / imv
 
c *** write out movie-header ***
      if (myid.eq.0) then
         open(unit=10, file=file_movie, status='unknown')
         if ( nnframe .ne. 0 ) then
            write(10,*) nnframe,nchain,nmolty,(rcut(ibox),ibox=1,nbox)
            nhere = 0
            do zz=1,nntype
               if ( lhere(zz) ) then
                  nhere = nhere + 1
                  temphe(nhere) = zz
               endif
            enddo

            write(10,*) nhere,(temphe(zz),zz=1,nhere)
         
            do imolty = 1,nmolty
               write(10,*) nunit(imolty)
c     output bond connectivity information
               do ii=1,nunit(imolty)
                  write(10,*) invib(imolty,ii)
     &                 ,(ijvib(imolty,ii,z),z=1,invib(imolty,ii))
               enddo

c     output torsional connectivity information
               do j = 1,nunit(imolty)
                  write(10,*) intor(imolty,j)
                  do ii = 1,intor(imolty,j)
                     write(10,*) ijtor2(imolty,j,ii)
     &                    ,ijtor3(imolty,j,ii),ijtor4(imolty,j,ii)
                  enddo
               enddo
            enddo

         endif
      endif

c     --- write out isolute movie header
      lprint = .false.
      do imol = 1,nmolty
         if (lsolute(imol)) then
            lprint = .true.
         endif
      enddo

      if (lprint.and.myid.eq.0) then
         write(11,*) nmolty
         do imol = 1, nmolty
            write(11,*) imol,nunit(imol),(nstep / isolute(imol))
     &           * temtyp(imol)
         enddo
      endif


c *** write out initial configuration for first movie frame ***
      if (nnstep .eq. 0) then
         dum = 1.0d0
            call monper(dum,dum,dum,dum,dum,dum,dum,dum
c * fixed by adding nbox, why the hell didn't this cause errors before?
c * KM fixed by adding dum for acsolpar
     &       ,dum,dum,dum,nbox,nnstep,dum,.false.,.false.,.false.
     &       ,.true.,.false.,.false.,lratfix,lsolute,dum,0.0d0,0.0d0)
      endif

c *** calculate constants for lmuir external potential ***
      if ( lmuir ) then
         sigpri = 0.715d0 * dsqrt( 3.8d0 * 3.93d0 )
         c9ch2 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*47.0d0) ) * sigpri**9
         c3ch2 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*47.0d0) ) * sigpri**3
         c9ch3 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*114.0d0) ) * sigpri**9
         c3ch3 = 4.0d0 * ( 1.43d0 * dsqrt(80.0d0*114.0d0) ) * sigpri**3
         zprmin = ( 3.0d0**(1/6.0d0) ) * sigpri
         v2prmin = c9ch2 / zprmin**9 - c3ch2 / zprmin**3
         v3prmin = c9ch3 / zprmin**9 - c3ch3 / zprmin**3
         betac2 = beta1 - v2prmin
         betac3 = beta1 - v3prmin
         write(iou,*) 'external potential for Langmuir monolayers used'
         write(iou,*) 'zprmin',zprmin
         write(iou,*) 'v2prmin',v2prmin,'v3prmin',v3prmin
      endif

C -------------------------------------------------------------------

c * write out connectivity and bonded interactions
      if (lverbose.and.myid.eq.0) then
         do imol = 1,nmolty
            write(iou,*) 'molecule type',imol
            if (nunit(imol) .gt. 1) then
               write(iou,*) '   i   j   type_i type_j   bond length',
     &              '        k/2'
            endif
            do i = 1,nunit(imol)
               do j = 1, invib(imol,i)
                  write(iou,1014) i,ijvib(imol,i,j),ntype(imol,i),
     &                 ntype(imol,ijvib(imol,i,j)),
     &                 brvib(itvib(imol,i,j)),brvibk(itvib(imol,i,j))
               enddo
            enddo

            if (nunit(imol) .gt. 2) then
               write(iou,*)
               write(iou,*) '   i   j   k   type_i type_j type_k',
     &              '     angle      k/2'
            endif
            do i = 1,nunit(imol)
               do j = 1,inben(imol,i)
                  write(iou,1015) i,ijben2(imol,i,j),
     &                 ijben3(imol,i,j),ntype(imol,i),
     &                 ntype(imol,ijben2(imol,i,j)),
     &                 ntype(imol,ijben3(imol,i,j)),
     &                 brben(itben(imol,i,j))*180.0d0/onepi,
     &                 brbenk(itben(imol,i,j))
               enddo
            enddo

            if (nunit(imol) .gt. 3) then
               write(iou,*)
               write(iou,*) '   i   j   k   l    type_i type_j type_k',
     &              ' type_l     torsion type'
            endif

            do i = 1,nunit(imol)
               do j = 1, intor(imol,i)
                  write(iou,1016) i,ijtor2(imol,i,j),ijtor3(imol,i,j),
     &                 ijtor4(imol,i,j),ntype(imol,i),
     &                 ntype(imol,ijtor2(imol,i,j)),
     &                 ntype(imol,ijtor3(imol,i,j)),
     &                 ntype(imol,ijtor4(imol,i,j)),
     &                 ittor(imol,i,j)
               enddo
            enddo

         enddo
      endif

c * write out non-bonded interaction table
      if (myid.eq.0) then
         write(iou,*)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
c added a flag for lgaro...need to add for all potentials, probably
         if ((.not.lexpsix).and.(.not.lmmff).and.(.not.lgaro)) then
            if (lninesix) then
               write(iou,*)
     & '     i   j    r_0,ij     epsij         q0(i)          q0(j)'
            elseif (lgenlj) then
               write(iou,*)
     & '     i   j    sig2ij     epsij         q0(i)          q0(j)'
            else
               write(iou,*) 
     & '     i   j    sig2ij     epsij         q0(i)          q0(j)'
            endif
            do i = 1, nntype
               do j = 1,nntype
                  if ( lhere(i) .and. lhere(j) ) then
                     if (lninesix) then
                        ij = (i-1)*nxatom + j
                        write(iou,'(3x,2i4,2f10.5,2f15.6)') i,j
     &                       ,rzero(ij),epsnx(ij),qelect(i),qelect(j)
                     elseif (lgenlj) then
                        ij = (i-1)*nntype + j
                        write(2,'(3x,2i4,2f10.5,2f15.6)') i,j
     &                       ,dsqrt(sig2ij(ij))
     &                    , epsij(ij),qelect(i),qelect(j)
                     else
                        ij = (i-1)*nntype + j
                        write(iou,'(3x,2i4,2f10.5,2f15.6)') i,j
     &                    ,dsqrt(sig2ij(ij))
     &                    , epsij(ij),qelect(i),qelect(j)
                     endif
                  endif
               enddo
            enddo
         endif
      endif

      if ( lneigh ) then
c *** If using neighbour list make sure the rcut & rcutnn is the same
c *** for all the boxes
         do ibox = 2,nbox
            if ((dabs(rcut(1)-rcut(ibox)).gt.1.0d-10).and.
     &           (dabs(rcutnn(1)-rcut(ibox)).gt.1.0d-10)) then
               write(iou,*) 'Keep rcut and rcutnn for all the 
     &              boxes same'
               ldie = .true.
               return
            endif
         enddo
 
c *** calculate squares of nn-radius and max. update displacement ***
         rcnnsq = rcutnn(1) * rcutnn(1)
         upnn = ( rcutnn(1) - rcut(1) ) / 3.0d0
         upnnsq = upnn * upnn
 
c *** calculate max. angular displacement that doesn't violate upnn ***
c *** calculate max. all-trans chain length ( umatch ) ***
         umatch = 0.0d0
         do j = 1, nunit(1) - 1
            umatch = umatch + brvib(1)
         enddo
         upnndg = asin( upnn / umatch )
 
        do im=1,2
           do imol=1,nmolty
              if ( rmtrax(imol,im) .gt. upnn ) then
                 write(iou,*) ' rmtrax greater than upnn',im,imol
                 rmtrax(imol,im) = upnn
              endif
              if ( rmtray(imol,im) .gt. upnn ) then
                 write(iou,*) ' rmtray greater than upnn',im,imol
                 rmtray(imol,im) = upnn
              endif
              if ( rmtraz(imol,im) .gt. upnn ) then
                 write(iou,*) ' rmtraz greater than upnn',im,imol
                 rmtraz(imol,im) = upnn
              endif
              
              if ( rmrotx(imol,im) .gt. upnndg ) then
                 write(iou,*) ' rmrotx greater than upnndg',im,imol
                 rmrotx(imol,im) = upnndg
              endif
              if ( rmroty(imol,im) .gt. upnndg ) then
                 write(iou,*) ' rmroty greater than upnndg',im,imol
                 rmroty(imol,im) = upnndg
              endif
              if ( rmrotz(imol,im) .gt. upnndg ) then
                 write(iou,*) ' rmrotz greater than upnndg',im,imol
                 rmrotz(imol,im) = upnndg
              endif
           enddo
        enddo
      endif

c *** write input data to unit 6 for control ***
      if (myid.eq.0) then
         write(iou,*)
         write(iou,*) 'number of mc cycles:            ', nstep
         write(iou,*) 'number of chains:               ', nchain
         write(iou,*)
         write(iou,*) 'temperature:                    ', temp
         write(iou,*)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MJM
         write(iou,*) 'ex-pressure:                    ', 
     &        (express(ibox),ibox=1,nbox)
         write(iou,*)
      endif


C -------------------------------------------------------------------

      do i = 1,nbox   
         if ( rcut(i) .ge. rcutnn(i) .and. lneigh ) then
            write(iou,*) ' rcut greater equal rcutnn for box',i
            ldie = .true.
            return
         endif
      enddo
 
 1012 format(3(1x,f10.6),2i5)
 1013 format(a20,i3,a13,f9.3,a5,f9.1)
 1014 format(i5,i4,i7,i7,f13.4,f14.1)
 1015 format(i5,i4,i4,i7,i7,i7,f12.2,f12.1)
 1016 format(i5,i4,i4,i4,i8,i7,i7,i7,i14)
 1017 format(1x,a10,i4,a20,a8,i4,a10,i4)
 1018 format(1x,a10,i3,a22,a8,i3,a4,i3,a10,i3,a19,i4)
 1019 format(a41,f8.4,f8.4,f8.4)
 1020 format(a36,a18,i5,i5,i5)
 1021 format(1x,a41,a5,i4,a10,f8.4)
      
 1101 format(' max trans. displacement:        ',3f10.6)
 1102 format(' max rot. displacement:          ',3f10.6)
 1103 format(' max volume displacement:        ',3e12.4)
 1104 format(' dimension box ',i1,'  :','  a:   ',f12.6,'   b:   ',f12.6
     & ,'   c   :  ' ,f12.6)
 1105 format(' angle of  box ',i1,'  :  ','   alpha:   ',f12.6,' 
     &  beta: ', f12.6, '    gamma:   ',f12.6)
 1106 format(' width of  box ',i1,':',3f12.6)


      return
      end





