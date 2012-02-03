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
      include 'gor.inc'

      integer temnc, imol, iutemp, imolty, itype,ipair,bdum,bin,histtot
      integer idummy(ntmax), atemp 

      integer i,j,k,ncres, nmtres, iensem, inpbc, nmcount
      integer im,nures, ibox,  ij, tcount,ucheck,nnframe
      integer nijspecial,ispecial,jspecial,ji,ii,jj,nexclu,ndum,ntii

      integer zz,temphe,z,itemp,zzz
      integer nvirial
      integer inclnum,inclmol,inclbead,inclsign,ncarbon
      dimension inclmol(ntmax*numax*numax),inclsign(ntmax*numax*numax)
      dimension inclbead(ntmax*numax*numax,2)

      integer ainclnum,ainclmol,ainclbead,a15t
      dimension ainclmol(ntmax*numax*numax)
      dimension ainclbead(ntmax*numax*numax,2)
      dimension a15t(ntmax*numax*numax)

      double precision starvir,stepvir,fqtemp,qelect,qbox,vol,v(3),w(3)

      double precision pie2,rcnnsq,umatch,aspecd,bspecd,dum,pm,pcumu
      logical newtab,lnrtab,lucall,lpolar,lqqelect,lee,lratfix,lreadq
      logical  linit, lecho, lmixlb, lmixjo, lhere,lsetup,lsolute
      logical lprint,lverbose,lxyz

      dimension lratfix(ntmax)
      dimension qbox(nbxmax)

c      double precision temx,temy,temz
      
      dimension nures(ntmax)
      dimension ncarbon(ntmax)
      dimension lhere(nntype)
      dimension lsolute(ntmax)
      dimension ucheck(ntmax)
      dimension qelect(nntype)
      dimension temphe(nntype)
c      dimension temx(nmax,numax),temy(nmax,numax),temz(nmax,numax)


c -- reads input data and initializes the positions
c
C --------------------------------------------------------------------
 
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

c *** read echoing and long output flags
      read(4,*)
      read(4,*) lecho,lverbose


c - read run information
      read(4,*)
      read(4,*) nstep, lstop, lpresim, iupdatefix
      if ( lecho) then
         if (lverbose) then
            if (lstop) then
               write(6,*) 'number of steps:',nstep
            else
               write(6,*) 'number of cycles:',nstep
            endif
            write(6,*) 'lstep:',lstop
            write(6,*) 'lpresim:',lpresim
            write(6,*) 'iupdatefix:',iupdatefix
         else
            write(6,*) nstep, lstop, lpresim, iupdatefix
         endif
      endif
      read(4,*)
      read(4,*) temp, express, fqtemp
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'temperature:',temp,' K'
            write(6,*) 'external pressure:',express,' s.u.'
            write(6,*) 'fluctuating charge temperature:',fqtemp,' K'
         else 
            write(6,*) temp, express, fqtemp
         endif
      endif

C$$   read the analysis information
      read(4,*)
      read(4,*)ianalyze,nbin,lrdf,lintra,lstretch,lgvst,lbend,lete,
     &           lrhoz,bin_width
      if (lecho) then
         if(lverbose) then
            write(6,*) 'ianalyze:', ianalyze
	    write(6,*) 'nbin', nbin
            write(6,*) 'lrdf:',lrdf
            write(6,*) 'lintra:',lintra
            write(6,*) 'lstretch:', lstretch
            write(6,*) 'lgvst:' ,lgvst
            write(6,*) 'lbend:' ,lbend
            write(6,*) 'lete:' ,lete
            write(6,*) 'lrhoz:', lrhoz 
            write(6,*) 'bin_width:' ,bin_width
         else
            write(6,*) ianalyze,nbin,lrdf,lintra,lstretch,lgvst,lbend
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
      
      beta = 1.0d0 / temp
      fqbeta = 1.0d0 / fqtemp
 
C -------------------------------------------------------------------

      read(4,*)
      read(4,*) iprint, imv, iratio, iblock, idiele
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'iprint:',iprint
            write(6,*) 'imv:',imv
            write(6,*) 'iratio:',iratio
            write(6,*) 'iblock:',iblock
            write(6,*) 'idiele:',idiele
         else
            write(6,*) iprint, imv, iratio, iblock, idiele
         endif
      endif
c - read information for histogram output (added 8/30/99)
      if (lgrand) then
         read(4,*)
         read(4,*) nequil,ninstf, ninsth, ndumph, run_num, suffix
         if ( lecho ) then
            if (lverbose) then
               write(6,*) 'nequil:',nequil
               write(6,*) 'ninstf:',ninstf
               write(6,*) 'ninsth:',ninsth
               write(6,*) 'ndumph:',ndumph
               write(6,*) 'run_num:',run_num
               write(6,*) 'suffix:',suffix
            else
               write(6,*) nequil,ninstf,ninsth,ndumph,run_num,suffix
            endif
         endif
      end if

      if (dint(dble(nstep)/dble(iblock)) .gt. 100) 
     &     stop 'too many blocks'

c - read system information

      read(4,*) 
      read(4,*) nbox
      if ( lecho ) write(6,*) 'number of boxes in the system:',nbox
      do i = 1,nbox
         read(4,*)
         read(4,*) boxlx(i),boxly(i),boxlz(i),lsolid(i),lrect(i),
c     +        kalp(i),kmax(i),rcutchg(i)
     +        kalp(i),rcutchg(i)
         if ( lecho ) then
            if (lverbose) then
               write(6,*) 'box:',i
               write(6,*) '   boxlx:',boxlx(i),' A'
               write(6,*) '   boxly:',boxly(i),' A'
               write(6,*) '   boxlz:',boxlz(i),' A'
               write(6,*) '   lsolid:',lsolid(i)
               write(6,*) '   lrect:',lrect(i)
               write(6,*) '   kalp:',kalp(i)
               write(6,*) '   rcutchg:',rcutchg(i),' A'
            else
               write(6,*) boxlx(i),boxly(i),boxlz(i),
c     +        lsolid(i),lrect(i),kalp(i),kmax(i),rcutchg(i)
     +        lsolid(i),lrect(i),kalp(i),rcutchg(i)
            endif
         endif
      enddo

      read(4,*) 
      read(4,*) nchain, nmolty
      if ( nmolty .gt. ntmax ) stop 'nmolty gt ntmax' 
      if ( nchain .gt. nmax-2 ) stop 'nchain gt nmax-2'
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'number of chains:',nchain
            write(6,*) 'number of molecule types:',nmolty
         else 
            write(6,*) nchain, nmolty
         endif
      endif
      read(4,*)
      read(4,*) (temtyp(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            do i = 1,nmolty
               write(6,*) 'number of chains of molecule type',i,':',
     &              temtyp(i)
            enddo
         else 
            write(6,*) (temtyp(i),i=1,nmolty)
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
         if (lecho) then
            if (lverbose) then
               do i = 1,nmolty
                  write(6,*) 'chemical potential for molecule type',
     &                 i,':',B(i)
               enddo
            else 
               write(6,*) "B ", (B(i),i=1,nmolty)
            endif
         endif
c - convert chemical potentials to activities
         do i=1,nmolty
            B(i) = exp(B(i)/temp) 
         enddo
         
      end if

c      if ( lecho ) write(6,*) 'moltyp',(moltyp(i),i=1,nchain)

c      if (lgrand) then
c         nchain=nmax
c         write(6,*) ' in Grand Can total number of chains set by NMAX!'
c      endif

      read(4,*)
      read(4,*) lmixlb, lmixjo
      if (lmixlb .and. lmixjo) stop 'cant use both combining rules!'
      if ( lecho ) then
         if (lverbose) then
            if (lmixlb) then
               write(6,*) 'Lorentz-Berthelot combining rules apply'
            else
               write(6,*) 'Jorgensen combining rules apply'
            endif
            write(6,*) '   lmixlb:',lmixlb,' lmixjo:',lmixjo
         else 
            write(6,*) lmixlb,lmixjo
         endif
      endif

c --- read special combining rule information
      read(4,*)
      read(4,*) nijspecial
      if (lecho) then
         if (lverbose) then
            write(6,*) 'number of special combining parameters:',
     &           nijspecial
         else
            write(6,*) nijspecial
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
            if (lecho) then
               if (lverbose) then
                  write(6,*) 'special parameter number',i
                  write(6,*) '   ispecial:',ispecial
                  write(6,*) '   jspecial:',jspecial
                  write(6,*) '   aspecd:',aspecd
                  write(6,*) '   bspecd:',bspecd
               else 
                  write(6,*) ispecial,jspecial,aspecd,bspecd
               endif
            endif
         enddo
      endif

      read(4,*) 
      read(4,*) rmin, rcut, rcutnn,softcut,rcutin,
     &     rbsmax,rbsmin
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'minimum cutoff (rmin):',rmin,' A'
            write(6,*) 'potential cutoff (rcut):',rcut,' A'
            write(6,*) 'neighbor list cutoff (rcutnn):',rcutnn,' A'
            write(6,*) 'softcut:',softcut
            write(6,*) 'CBMC inner cutoff (rcutin):',rcutin,' A'
            write(6,*) 'AVBMC outer cutoff (rbsmax):',rbsmax,' A'
            write(6,*) 'AVBMC inner cutoff (rbsmin):',rbsmin,' A'
         else 
            write(6,*) rmin, rcut, rcutnn, softcut
     &           ,rcutin,rbsmax,rbsmin 
         endif
      endif
      do i = 1, nbox
         if( rcut/boxlx(i) .gt. 0.5d0) then
            write(6,*) 'rcut > 0.5*boxlx'
            stop
         endif
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
         

         if ( nunit(imol) .gt. numax ) stop 'nunit gt numax'
         if ( lflucq(imol)
     &        .and. (.not. lelect(imol) ) )
     &        stop 'lelect must be true if flucq is true' 

         if ( lqtrans(imol) ) then
            if (.not. lflucq(imol) )
     &           stop 'lflucq must be true if interm. CT is allowed'
            write(6,*) 'Intermolecular Charge Transfer is allowed'
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
         if ( lavbmc2(imol) .or. lavbmc3(imol) ) lneighbor = .true.

         if ( lecho ) then
            if (lverbose) then
               write(6,*) 'molecule type:',imol
               write(6,*) '   number of units:',nunit(imol)
               write(6,*) '   number of units for CBMC growth:',
     &              nugrow(imol)
               write(6,*) '   number of carbons for EH alkane:',
     &              ncarbon(imol)
               write(6,*) '   maximum number of units for CBMC:',
     &              nmaxcbmc(imol)
               write(6,*) '   iurot:',iurot(imol)
               write(6,*) '   lelect:',lelect(imol)
               write(6,*) '   lflucq:',lflucq(imol)
               write(6,*) '   lqtrans:',lqtrans(imol)
               write(6,*) '   lexpand:',lexpand(imol)
               write(6,*) '   lavbmc1:',lavbmc1(imol)
               write(6,*) '   lavbmc2:',lavbmc2(imol)
               write(6,*) '   lavbmc3:',lavbmc3(imol)
               write(6,*) '   fqegp:',fqegp(imol)
               write(6,*) '   lsetup:',lsetup
               do i = 1,nbox
                  write(6,*) '   energy offset for box',
     &                 i,':',eta2(i,imol),' K'
               enddo
            else 
               write(6,*) nunit(imol),nugrow(imol),ncarbon(imol)
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
               if ( lecho ) 
     &              write(6,*) '   bead ',j,' beadtype ',ntype(imol,i),
     &              ' charge leader ',leaderq(imol,i)
               if ( leaderq(imol,i) .gt. j .and. .not. lchgall)
     &              stop 'group-based cut-off screwed for qq'
            else
               read(4,*) j, ntype(imol,i)
               if ( lecho ) then
                  if (lverbose) then
                     write(6,*) '   bead ',j,' beadtype ',ntype(imol,i),
     &                    chname(ntype(imol,i))
                  else 
                     write(6,*) '   bead ',j,' beadtype ',ntype(imol,i)
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
               write(6,*) 'imol',imol,'   i',i,'   invib',invib(imol,i)
               stop 'too many vibrations'
            endif
            do j = 1, invib(imol,i)
               read(4,*) ijvib(imol,i,j),itvib(imol,i,j)
               if (lverbose) then
c                  write(6,*) '      bead',i,' bonded to bead',
c     &                 ijvib(imol,i,j),' with bond type:',
c     &                 itvib(imol,i,j)
                  write(6,*) '      bead',i,' bonded to bead',
     &                 ijvib(imol,i,j)
                  write(6,1013) '          bond type:',
     &                 itvib(imol,i,j),' bond length:',
     &                 brvib(itvib(imol,i,j)),' k/2:',
     &                 brvibk(itvib(imol,i,j))
               endif
            enddo
c - bond bending -
            read(4,*)
            read(4,*) inben(imol,i)
            if ( inben(imol,i) .gt. 12 ) stop 'too many bends'
            do j = 1, inben(imol,i)
               read(4,*) ijben2(imol,i,j),ijben3(imol,i,j)
     &              ,itben(imol,i,j)
               if (lverbose) then
c                  write(6,*) '      bead',i, ' bending interaction',
c     &                 ' through',ijben2(imol,i,j),' with bead',
c     &                 ijben3(imol,i,j),' of bend type:',itben(imol,i,j)
                  write(6,1017) '      bead',i, ' bending interaction',
     &                 ' through',ijben2(imol,i,j),' with bead',
     &                 ijben3(imol,i,j)
                  write(6,1013) '          bend type:',itben(imol,i,j),
     &                 ' bend angle :',
     &                 brben(itben(imol,i,j))*180.0d0/onepi,
     &                 ' k/2:',brbenk(itben(imol,i,j))
               endif
            enddo
c - bond torsion -
            read(4,*)
            read(4,*) intor(imol,i)
            if ( intor(imol,i) .gt. 10 ) stop 'too many torsions'
            do j = 1, intor(imol,i)
               read(4,*) ijtor2(imol,i,j),ijtor3(imol,i,j),
     &              ijtor4(imol,i,j),ittor(imol,i,j)
               if (lverbose) then
                  write(6,1018) '      bead',i, ' torsional interaction'
     &                 ,' through',ijtor2(imol,i,j),' and',
     &                 ijtor3(imol,i,j),' with bead',ijtor4(imol,i,j),
     &                 ' of torsional type:',ittor(imol,i,j)
               endif
            enddo
         enddo

 112     continue

         if ( lexpand(imol) ) then
            if ( temtyp(imol) .gt. 1 ) then
               write(6,*) 'Only one molecule of this type is allowed!'
               stop
            endif
            lee = .true.
            read(7,*) 
            read(7,*) numcoeff(imol)
            do j = 1,numcoeff(imol)
               read(7,*)
               read(7,*) (epsil(imol,ii,j),ii=1,nunit(imol))
               write(6,*) 'itype:',j
               write(6,*) (epsil(imol,ii,j),ii=1,nunit(imol))
               read(7,*) (sigm(imol,ii,j),ii=1,nunit(imol))
               write(6,*) (sigm(imol,ii,j),ii=1,nunit(imol))
               read(7,*) (qcharge(imol,ii,j),ii=1,nunit(imol))
               write(6,*) (qcharge(imol,ii,j),ii=1,nunit(imol))
               read(7,*)
               read(7,*) (eta(ii,imol,j),ii=1,2)
               write(6,*) 'eta:',(eta(ii,imol,j),ii=1,2)
            enddo
         endif
         if ( lbias(imol) ) then
            read(4,*)
            read(4,*) pmbias(imol),(pmbsmt(ii),ii=1,nmolty)
     &           ,pmbias2(imol)
            if ( lecho ) then
               if (lverbose) then
                  write(6,*) '   AVBMC pmbias',pmbias(imol)
                  do ii = 1,nmolty
                     write(6,*) '   AVBMC2 and 3 probability for',
     &                    ' molecule',' type',ii,':',pmbsmt(ii)
                  enddo
                  write(6,*) '   AVBMC3 pmbias2:',pmbias2(imol)
               else
                  write(6,*) '   AVBMC bias for cluster formation and',
     &                 ' destruction'
                  write(6,*) pmbias(imol),(pmbsmt(ii),ii=1,nmolty)
     &                 ,pmbias2(imol)
               endif
            endif
            if (rbsmax .lt. rbsmin)
     &           stop 'rbsmax should be greater than rbsmin'
         endif
      enddo

c -- check whether there is a polarizable molecule

      lpolar = .false.
      lqqelect = .false.
      do imol = 1, nmolty
         if (lflucq(imol)) lpolar = .true.
         if (lelect(imol)) lqqelect = .true.
      enddo
      if ( .not. lqqelect ) then
         if ( lewald .or. lchgall ) 
     &        stop 'no charges in the system and turn off lewald'
         do ibox = 1, nbox
            rcutchg(ibox) = 0.0d0
         enddo
      endif

      if ( .not. lpolar ) then
         if ( lanes  ) 
     &        stop 'lanes should be false for nonpolarizable systems!'
      if ( lfepsi ) 
     &        stop 'lfepsi should be false for nonpolarizable systems!'
      endif

c - read linkcell information
      read(4,*)
      read(4,*) licell,rintramax,boxlink

      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'licell:',licell
            write(6,*) 'rintramax:',rintramax,' A'
            write(6,*) 'boxlink:',boxlink
         else 
            write(6,*) licell,rintramax,boxlink
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
      if ( lecho ) then
         if (lverbose) then
            write(6,1019) ' initial maximum x, y and z displacement:',
     &           rmtrax(1,1),rmtray(1,1),rmtraz(1,1)
         else 
            write(6,*) rmtrax(1,1), rmtray(1,1), rmtraz(1,1)
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
      if ( lecho ) then
         if (lverbose) then
            write(6,1019) ' initial maximum x, y and z rotation:    ',
     &           rmrotx(1,1),rmroty(1,1),rmrotz(1,1)
         else 
            write(6,*) rmrotx(1,1), rmroty(1,1), rmrotz(1,1)
         endif
      endif
      read(4,*) 
      read(4,*) tatra,tarot
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'target translational acceptance ratio:',tatra
            write(6,*) 'target rotational acceptance ratio:',tarot
         else 
            write(6,*) tatra, tarot
         endif
      endif

c - read initial setup information
      read(4,*)
      read(4,*) linit,newtab,lreadq,qscale
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'linit:',linit
            write(6,*) 'lnewtab:',newtab
            write(6,*) 'lreadq:',lreadq
            write(6,*) 'qscale:',qscale
         else 
            write(6,*) linit, newtab, lreadq,qscale
         endif
      endif
      read(4,*)
      read(4,*) (lbranch(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            do i = 1,nmolty
               write(6,*) 'lbranch for molecule type',i,':',lbranch(i)
            enddo
         else 
            write(6,*) (lbranch(i),i=1,nmolty)
         endif
      endif
      do i = 1, nbox
         read(4,*)
         read(4,*) (ininch(j,i),j=1,nmolty)
         if ( lecho ) then
            if (lverbose) then
               write(6,*) 'box:',i
               do j = 1,nmolty
                  write(6,*) '   initial number of chains of type',
     &                 j,':',ininch(j,i)
               enddo
            else 
               write(6,*) 'box:',i,(ininch(j,i),j=1,nmolty)
            endif
         endif
         read(4,*)
         read(4,*) inix(i),iniy(i),iniz(i),inirot(i),inimix(i),
     &        zshift(i),dshift(i),nchoiq(i)
         if ( lecho ) then
            if (lverbose) then
               write(6,1020) '    initial number of chains in x, y',
     &              ' and z directions:',inix(i),iniy(i),iniz(i)
               write(6,*) '   initial rotational displacement:',
     &              inirot(i)
               write(6,*) '   inimix:',inimix(i)
               write(6,*) '   zshift:',zshift(i)
               write(6,*) '   dshift:',dshift(i)
               write(6,*) '   nchoiq:',nchoiq(i)
            else 
               write(6,*) inix(i),iniy(i),iniz(i),inirot(i),
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
      if ( lecho ) then
         if (lverbose) then
            do zz = 1,nbox
               write(6,*) 'initial maximum volume displacement (rmvol)',
     &              ' in box',zz,':',rmvol(zz)
            enddo
            write(6,*) 'target volume acceptance ratio (tavol):',tavol
            write(6,*) 'iratv:',iratv
            write(6,*) 'iratp:',iratp
            do zz = 1,nmolty
               do zzz = 1,nbox
                  write(6,*) 'initial maximum fluct. charge',
     &                 ' displ. for chain type',zz,' box',
     &                 zzz,':',rmflcq(zz,zzz)
               enddo
            enddo
            write(6,*) 'target fluctuating charge acceptance ratio',
     &           ' (taflcq):',taflcq
         else
            write(6,*) rmvol, tavol, iratv, iratp
            write(6,*) 'rmflcq',
     &           ((rmflcq(zz,zzz),zzz=1,nbox),zz=1,nmolty),taflcq
         endif
      endif

      read(4,*)
      read(4,*) pmvol,(pmvlmt(j),j=1,nbox)
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmvol:',pmvol
            do j = 1,nbox
               write(6,*) '   pmvlmt for box',j,':',pmvlmt(j)
            enddo
         else 
            write(6,*) 'pmvol',pmvol,(pmvlmt(j),j=1,nbox)
         endif
      endif
      read(4,*)
      read(4,*) nvolb,(pmvolb(j),j=1,nvolb)
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'nvolb:',nvolb
            do j = 1,nvolb
               write(6,*) '   pmvolb:',pmvolb(j)
            enddo
         else 
            write(6,*) '   nvolb',nvolb,(pmvolb(j),j=1,nvolb)
         endif
      endif
      read(4,*)
      do j = 1,nvolb
         read(4,*) box5(j),box6(j)
         if ( lecho ) then
            if (lverbose) then
               write(6,*) '   box pair for volume move number',j,':',
     &              box5(j),box6(j)
            else 
               write(6,*) box5(j),box6(j)
            endif
         endif
      enddo
      
      lxyz = .false.
      do j = 1,nbox
         if (lsolid(j) .and. .not. lxyz) then
            lxyz = .true.
            read(4,*)
            read(4,*) pmvolx,pmvoly
            if (lecho) then
               if (lverbose) then
                  write(6,*) 'pmvolx:',pmvolx
                  write(6,*) 'pmvoly:',pmvoly
               else
                  write(6,*) pmvolx,pmvoly
               endif
            endif
         endif
      enddo

c     --- read swatch information
      read(4,*)
      read(4,*) pmswat,nswaty
      if ( nswaty .gt. npamax ) then
         write(6,*) 'nswaty gt npamax'
         stop
      endif

      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmswat:',pmswat
            write(6,*) '   number of swatch pairs (nswaty):',nswaty
         else 
            write(6,*) 'pmswat, nswaty',pmswat,nswaty
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
               write(6,*) 'nswaty ',i,' has identical moltyp'
               stop 'cannot swatch identical moltyp'
            endif
         enddo

         if( lecho ) then
            if (lverbose) then
               do i = 1,nswaty
                  write(6,*) '   swatch molecule type pairs:',
     &                 (nswatb(i,j),j=1,2)
               enddo
               do i = 1,nswaty
                  write(6,*) '   probability of each swatch pair:',
     &                 pmsatc(i)
               enddo
            else
               do i=1,nswaty
                  write(6,*) 'nswatb', (nswatb(i,j),j=1,2)
                  write(6,*) 'pmsatc',pmsatc(i)
               enddo
            endif
         endif
         
         do i=1,nswaty
c           --- number of beads that remain in the same position
            read(4,*)
            read(4,*) nsampos(i),(ncut(i,j),j=1,2)
            if (lecho) then
               if (lverbose) then
                  write(6,*) '   nsampos:',nsampos(i)
                  write(6,*) '   ncut:',(ncut(i,j),j=1,2)
               else 
                  write(6,*) 'nsampos',nsampos(i),' ncut',
     &                 (ncut(i,j),j=1,2)
               endif
            endif
c           --- bead number
            read(4,*)
            do j = 1,nsampos(i)
               read(4,*) (splist(i,j,k),k=1,2)
               if (lecho) then
                  if (lverbose) then
                     write(6,*) '   splist:',
     &                    (splist(i,j,k),k=1,2)
                  else 
                     write(6,*) 'splist',(splist(i,j,k),k=1,2)
                  endif
               endif
            enddo
            
            read(4,*)
            read(4,*) (( gswatc(i,j,k), k=1,2*ncut(i,j) ), j=1,2 )
!!            if (lecho) then
!!               if (lverbose) then
!!                  do zz = 1,ncut(i,j)
!!                     write(6,*) '   grow from and prev for ncut',zz,':',
!!     &                    (gswatc(i,j,zz),j=1,2)
!!                  enddo
!!               else 
!!                  write(6,*) 'gswatc',(( gswatc(i,j,k), 
!!     &                 k=1,2*ncut(i,j) ), j=1,2 )
!!               endif
!!            endif

            read(4,*) 
            read(4,*) nswtcb(i), (pmswtcb(i,ipair), ipair=1,nswtcb(i))
            read(4,*)
            do ipair = 1,nswtcb(i)
               read(4,*) box3(i,ipair),box4(i,ipair)
               if (lecho) then
                  if (lverbose) then
                     write(6,*) '   box pair:',
     &                    box3(i,ipair),box4(i,ipair)
                  else
                     write(6,*) box3(i,ipair),box4(i,ipair)
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
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmswap:',pmswap
            do i = 1,nmolty
               write(6,1021) '   swap probability for molecule type    '
     &              ,'     ',i,' (pmswmt):',pmswmt(i)
            enddo
         else 
            write(6,*) 'pmswap',pmswap,(pmswmt(i),i=1,nmolty)
         endif
      endif
      do i = 1, nmolty
         read(4,*)
         read(4,*) nswapb(i), (pmswapb(i,ipair),ipair=1,nswapb(i))
         if ( lecho ) then
            if (lverbose) then
               write(6,*) '   number of swap moves for molecule type',i,
     &              ':',nswapb(i)
               do ipair = 1,nswapb(i)
                  write(6,*) '      pmswapb:',pmswapb(i,ipair)
               enddo
            else 
               write(6,*) nswapb(i), 
     &              (pmswapb(i,ipair),ipair=1,nswapb(i))
            endif
         endif
         read(4,*)
         do ipair = 1, nswapb(i)
            read(4,*) box1(i,ipair), box2(i,ipair)
            if ( lecho ) then
               if (lverbose) then
                  write(6,*) '      box pair:',
     &                 box1(i,ipair), box2(i,ipair)
               else 
                  write(6,*) 'box pair:',
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
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmcb:',pmcb
            do i = 1,nmolty
               write(6,1021) '   CBMC probability for molecule type    '
     &              ,'     ',i,' (pmcbmt):',pmcbmt(i)
            enddo
            do i = 1,nmolty
               write(6,*) '   pmall for molecule type',i,':',pmall(i)
            enddo
         else 
            write(6,*) 'pmcb',pmcb,(pmcbmt(i),i=1,nmolty),
     &           'pmall',(pmall(i),i=1,nmolty)
         endif
      endif
      read(4,*)
      read(4,*) (pmfix(i),i=1,nmolty)
      if (lecho) then
         if (lverbose) then
            do i = 1,nmolty
               write(6,1021) '   probability for SAFE-CBMC for molecule'
     &              ,' type',i,'  (pmfix):',pmfix(i)
            enddo
         else
            write(6,*) (pmfix(i),i=1,nmolty)
         endif
      endif
      do i = 1, nmolty
         if (lring(i).and.pmfix(i).lt.0.999.and.
     &        .not.lrig(i)) then
            stop 'a ring can only be used with safe-cbmc'
         endif
      enddo
c --- read fluctuating charge info
      read(4,*)
      read(4,*) pmflcq, (pmfqmt(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmflcq:',pmflcq
            do i = 1,nmolty
               write(6,1021) '   flcq probability for molecule type    '
     &              ,'     ',i,' (pmfqmt):',pmfqmt(i)
            enddo
         else 
            write(6,*) 'pmflcq',pmflcq,(pmfqmt(i),i=1,nmolty)
         endif
      endif
c --- read expanded-coefficient move info
      read(4,*)
      read(4,*) pmexpc, (pmeemt(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmexpc:',pmexpc
            do i = 1,nmolty
               write(6,1021) '   expanded ens. prob. for molecule type '
     &              ,'     ',i,' (pmeemt):',pmeemt(i)
            enddo
         else 
            write(6,*) 'pmexpc',pmexpc,(pmeemt(i),i=1,nmolty)
         endif
      endif
c --- read translation info
      read(4,*)
      read(4,*) pmtra,(pmtrmt(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmtra:',pmtra
            do i = 1,nmolty
               write(6,1021) '   translation probability for molecule  '
     &              ,' type',i,' (pmtrmt):',pmtrmt(i)
            enddo
         else 
            write(6,*) 'pmtra',pmtra,(pmtrmt(i),i=1,nmolty)
         endif
      endif
c --- read rotation info
      read(4,*)
      read(4,*) (pmromt(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'pmrot:',1.0d0
            do i = 1,nmolty
               write(6,1021) '   rotational probability for molecule   '
     &              ,' type',i,' (pmromt):',pmromt(i)
            enddo
         else 
            write(6,*) (pmromt(i),i=1,nmolty)
         endif
      endif

c *** writeout probabilities, in percent
      if (lverbose) then
         write(6,*)
         write(6,*) 'percentage move probabilities:'
         write(6,'(1x,a19,f8.2,a2)') 'volume move       :',
     &      100.0d0*pmvol,' %'
         pcumu = pmvol
         if (pmswat .gt. pmvol) then
            pm = pmswat - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(6,'(1x,a19,f8.2,a2)') 'swatch move       :',
     &      100.0d0*pm,' %'
         if (pmswap .gt. pmswat) then
            pm = pmswap - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(6,'(1x,a19,f8.2,a2)') 'swap move         :',
     &      100.0d0*pm,' %'
         if (pmcb .gt. pmswap) then
            pm = pmcb - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(6,'(1x,a19,f8.2,a2)') 'CBMC move         :',
     &      100.0d0*pm,' %'
         if (pmflcq .gt. pmcb) then
            pm = pmflcq - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(6,'(1x,a19,f8.2,a2)') 'fluct charge move :',
     &      100.0d0*pm,' %'
         if (pmexpc .gt. pmflcq) then
            pm = pmexpc - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(6,'(1x,a19,f8.2,a2)') 'expanded ens move :',
     &      100.0d0*pm,' %'
         if (pmtra .gt. pmexpc) then
            pm = pmtra - pcumu
            pcumu = pcumu + pm
         else
            pm = 0.0d0
         endif
         write(6,'(1x,a19,f8.2,a2)') 'translation move  :',
     &      100.0d0*pm,' %'
         pm = 1.0d0 - pmtra
         write(6,'(1x,a19,f8.2,a2)') 'rotation move     :',
     &      100.0d0*pm,' %'

         write(6,*) 
      endif

c --- read growth details
      read(4,*)
      read(4,*) (nchoi1(i),i=1,nmolty),(nchoi(i),i=1,nmolty)
     &     ,(nchoir(i),i=1,nmolty)
     &     ,(nchoih(i),i=1,nmolty),(nchtor(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'molecule type :',(i,'  ',i=1,nmolty)
            write(6,*) '     nchoi1   :',(nchoi1(i),' ',i=1,nmolty)
            write(6,*) '     nchoi    :',(nchoi(i),' ',i=1,nmolty)
            write(6,*) '     nchoir   :',(nchoir(i),' ',i=1,nmolty)
            write(6,*) '     nchoih   :',(nchoih(i),' ',i=1,nmolty)
            write(6,*) '     nchtor   :',(nchtor(i),i=1,nmolty)
         else
            write(6,*) 'nchoi1',(nchoi1(i),i=1,nmolty)
            write(6,*) 'nchoi',(nchoi(i),i=1,nmolty)
            write(6,*) 'nchoir',(nchoir(i),i=1,nmolty)
            write(6,*) 'nchoih',(nchoih(i),i=1,nmolty)
            write(6,*) 'nchtor',(nchtor(i),i=1,nmolty)
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
      if ( lecho ) then
         if (lverbose) then
            do i = 1,nmolty
               write(6,*) 'nchbna and nchbnb for molecule type',i,':',
     &              nchbna(i),nchbnb(i)
            enddo
         else
            write(6,*) 'nchbna ',(nchbna(i),i=1,nmolty)
            write(6,*) 'nchbnb ',(nchbnb(i),i=1,nmolty)
         endif
      endif

      read(4,*)
      read(4,*) (icbdir(i),i=1,nmolty),(icbsta(i),i=1,nmolty)
      if ( lecho ) then
         if (lverbose) then
            do i = 1,nmolty
               write(6,*) 'icbdir for molecule type',i,':',icbdir(i)
            enddo
            do i = 1,nmolty
               write(6,*) 'icbsta for molecule type',i,':',icbsta(i)
            enddo
         else
            write(6,*) 'icbdir',(icbdir(i),i=1,nmolty)
            write(6,*) 'icbsta',(icbsta(i),i=1,nmolty)
         endif
      endif

c     ---- Error checking
      do i=1,nmolty
         if ( nchoi1(i) .gt. nchmax ) stop 'nchoi1 gt nchmax'
         if (nchoi(i) .gt. nchmax ) stop 'nchoi gt nchmax'
         if ( nchoih(i) .ne. 1 .and. nunit(i) .eq. nugrow(i) ) 
     &        stop ' nchoih must be one if nunit = nugrow'
         if ( nchtor(i) .gt. nchtor_max ) 
     &        stop 'nchtor gt nchtor_max'
         if ( nchbna(i) .gt. nchbn_max ) 
     &        stop 'nchbna gt nchbn_max'
         if ( nchbnb(i) .gt. nchbn_max ) 
     &        stop 'nchbnb gt nchbn_max'
         if ( icbsta(i) .gt. numax ) stop 'icbsta gt numax'
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
      if ( lecho ) then
         if (lverbose) then
            write(6,*)'nexclu:',nexclu
         else 
            write(6,*) nexclu
         endif
      endif
      if (nexclu .ne. 0) then
         do ndum = 1, nexclu
            read(4,*) i, ii, j, jj
            lexclu(i,ii,j,jj) = .true.
            lexclu(j,jj,i,ii) = .true.
            if (lverbose) then
               write(6,*) 'excluding interactions between bead',ii,
     &              ' on chain',i,' and bead',jj,' on chain',j
            endif
         enddo
      else
         read(4,*)
      endif

c --- read inclusion list 
      read(4,*)
      read(4,*) inclnum
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'inclnum:',inclnum
         else 
            write(6,*) inclnum
         endif
      endif
      if (inclnum .ne. 0) then
         do ndum = 1, inclnum
            read(4,*) inclmol(ndum),inclbead(ndum,1),inclbead(ndum,2)
     &           ,inclsign(ndum)
            if (lecho ) then
               if (lverbose) then
                  if (inclsign(ndum) .eq. 1) then
                     write(6,*) 'including intramolecular interactions',
     &                    'for chain type',inclmol(ndum),
     &                    ' between beads',inclbead(ndum,1),' and',
     &                    inclbead(ndum,2)
                  else
                     write(6,*) 'excluding intramolecular interactions',
     &                    'for chain type',inclmol(ndum),
     &                    ' between beads',inclbead(ndum,1),' and',
     &                    inclbead(ndum,2)

                  endif
               else
                  write(6,*) inclmol(ndum),inclbead(ndum,1)
     &                 ,inclbead(ndum,2),inclsign(ndum)
               endif
            endif
         enddo
      else
         read(4,*)
      endif

c --- read a15 inclusion list 
      read(4,*)
      read(4,*) ainclnum
      if ( lecho ) then
         if (lverbose) then
            write(6,*) 'ainclnum:',ainclnum
         else 
            write(6,*) ainclnum
         endif
      endif
      if (ainclnum .ne. 0) then
         do ndum = 1, ainclnum
            read(4,*) ainclmol(ndum),ainclbead(ndum,1),ainclbead(ndum,2)
     &           ,a15t(ndum)
            if (lecho) then
               if (lverbose) then
                  write(6,*) 'repulsive 1-5 OH interaction for',
     &                 ' chain type',ainclmol(ndum),' between beads',
     &                 ainclbead(ndum,1),' and',ainclbead(ndum,2),
     &                 ' of type:',a15t(ndum)
               else
                  write(6,*) ainclmol(ndum),ainclbead(ndum,1)
     &                 ,ainclbead(ndum,2),a15t(ndum)
               endif
            endif
         enddo
      else
         read(4,*)
      endif

c - set up the inclusion table
      call inclus( inclnum,inclmol,inclbead,inclsign,ncarbon,
     &     ainclnum,ainclmol,ainclbead,a15t)

c -  read in information on the chemical potential checker
      read(4,*)
      read(4,*) lucall
      read(4,*)
      read(4,*) (ucheck(jj),jj=1,nmolty)
      if (lecho) then
         if (lverbose) then
            write(6,*) 'lucall:',lucall
            do jj = 1,nmolty
               write(6,*) '   ucheck for molecule type',
     &              jj,':',ucheck(jj)
            enddo
         else 
            write(6,*) lucall,(ucheck(jj),jj=1,nmolty)
         endif
      endif

c -   read information for virial coefficient calculation
      read(4,*) 
      read(4,*) nvirial,starvir,stepvir
      if (lecho) then
         if (lverbose) then
            write(6,*) 'nvirial:',nvirial
            write(6,*) 'starvir:',starvir
            write(6,*) 'stepvir:',stepvir
         else 
            write(6,*) nvirial,starvir,stepvir
         endif
      endif

      if (lvirial) then
         if ( nvirial .gt. maxvir ) stop 'nvirial .gt. maxvir'
      
         read(4,*)
         read(4,*) ntemp,(virtemp(jj),jj=1,ntemp)
         if (lecho) then
            if (lverbose) then
               write(6,*) 'ntemp:',ntemp
               write(6,*) 'calculation of virial coefficient ',
     &              'at the following temperatures:',
     &              (virtemp(jj),jj=1,ntemp)
            else 
               write(6,*) ntemp,
     &              'Calculation of virial coefficient ',
     &              'at the following temperatures:',
     &              (virtemp(jj),jj=1,ntemp)
            endif
         endif
      endif

      nhere = 0
      do zz=1,nntype
         if ( lhere(zz) ) then
             nhere = nhere + 1
             temphe(nhere) = zz
             beadtyp(nhere)=zz
         endif
      enddo


      do zz = 1,nhere 
          atemp = temphe(zz)
          decode(atemp) = zz
c	  print*, 'readdat'
c	  print*, 'atemp=',atemp, 'decode=',decode(atemp)

      enddo



C -------------------------------------------------------------------
 
c *** restart saver
      if ( nstep .gt. 100 ) then
         irsave = nstep / 10
      else
         irsave = nstep + 1
      endif

c -------------------------------------------------------------------
 
c * check information of .INC-files 
      write(6,*)
      write(6,*) '***** program   =  THE MAGIC BLACK BOX'
      iensem = 0
      if ( lgrand ) then
         iensem = iensem + 1
         write(6,*) 'grand-canonical ensemble'
      endif
      if ( lvirial ) then
         write(6,*) 'Computing Second Virial Coefficient'
         if ( nchain .ne. 2) stop 'nchain must equal 2'
      endif
      if ( lgibbs ) then
         iensem = iensem + 1
         if ( lnpt ) then
            write(6,*) 'NPT Gibbs ensemble'
         else
            write(6,*) 'NVT Gibbs ensemble'
         endif
      endif
      if ( iensem .eq. 0 ) then
         if ( lnpt ) then
            write(6,*) 'Isobaric-isothermal ensemble'
            iensem = iensem + 1
         else
            write(6,*) 'Canonical ensemble'
            iensem = iensem + 1
         endif
      endif
      if ( iensem .gt. 1 ) stop 'INCONSISTENT ENSEMBLE SPECIFICATION'

      inpbc = 0
      if ( lpbc ) then
         write(6,*) 'using periodic boundaries'
         if ( lpbcx ) then
            inpbc = inpbc + 1
            write(6,*) 'in x-direction'
         endif
         if ( lpbcy ) then
            inpbc = inpbc + 1
            write(6,*) 'in y-direction'
         endif
         if ( lpbcz ) then
            inpbc = inpbc + 1
            write(6,*) 'in z-direction'
         endif
         if ( inpbc .eq. 0 ) stop 'INCONSISTENT PBC SPECIFICATION'
         write(6,*) inpbc,'-dimensional periodic box'
      else
         write(6,*) 'cluster mode (no pbc)'
         if ( lgibbs .or. lgrand ) stop
     &        'INCONSISTENT SPECIFICATION OF LPBC AND ENSEMBLE'
      endif

      if ( lfold ) then
         write(6,*) 'particle coordinates are folded into central box'
         if ( .not. lpbc ) stop 
     &        'INCONSISTENT SPECIFICATION OF LPBC AND LFOLD'
      endif

      if ( lijall ) then
         write(6,*) 'all i-j-interactions are considered',
     &        '  (no potential cut-off)'
      endif

      if ( lchgall ) then
         write(6,*) 'all the inter- and intramolecular Coulombic',
     &        ' interactions are considered (no group-based cutoff)'
      endif
      if ( lcutcm ) then
         write(6,*) 'additional (COM) cutoff on computed rcmu'
         if ( lijall ) stop 'cannot have lijall with lcutcm'
c         if ( lchgall ) stop 'cannot have lchgall with lcutcm'
      endif
      if ( ldual ) then
         write(6,*) 'Dual Cutoff Configurational-bias Monte Carlo'
      endif

      write(6,*) 
     & 'CBMC simultaneously grows all beads conected to the same bead'
      write(6,*) 'with bond angles generated from Gaussian distribution'

      if ( ljoe ) then
         write(6,*) 'external 12-3 potential for SAM (Joe Hautman)'
      endif

      if ( lsami ) then
         write(6,*) 'external potential for Langmuir films (Sami)'
         write(6,*) 'WARNING: LJ potential defined in SUSAMI'
         write(6,*) 'WARNING: sets potential cut-off to 2.5sigma'
         write(6,*) 'WARNING: has build-in tail corrections'
         if ( ltailc ) stop
     &        'INCONSISTENT SPECIFICATION OF LTAILC AND LSAMI'
      endif

      if ( lexzeo ) then
         write(6,*) 'external potential for zeolites (Berend)'
         write(6,*) 'WARNING: potential defined in SUZEO'
      endif

      write(6,*) 'Program will call Explct.f if needed'

      if ( lexpsix ) then
         write(6,*) 'Exponential-6 potential for bead-bead'
      elseif ( lmmff ) then
         write(6,*) 'Buffered 14-7 potential for bead-bead'
      elseif (lninesix) then
         write(6,*) '9-6 potential for bead-bead'
      else
         write(6,*) 'Lennard-Jones potential for bead-bead'
      endif

      if ( ltailc ) then
         write(6,*) 'with additional tail corrections'
         if (lshift) stop ' lshift.and.ltailc!'
      endif

      if ( lshift) then
         write(6,*) 'using a shifted potential'
         if (ltailc) stop ' lshift.and.ltailc!'
      endif
         
      write(6,*) 'Coulombic inter- and intramolecular interactions'
      if ( lewald ) write(6,*)
     &     'Ewald-sum will be used to calculate Coulombic interactions'
      write(6,*)
      write(6,*) 'MOLECULAR MASS:', (masst(i),i=1,nmolty)

C -------------------------------------------------------------------

c *** read/produce initial/starting configuration ***
c *** zeolite external potential
      if ( lexzeo ) then
         call suzeo(rcut,newtab)
         boxlx(1)=zeorx
         boxly(1)=zeory
         boxlz(1)=zeorz
         write(6,*) ' note zeolite determines the box size !'
      endif

C ----------------------------------------------------------------
      if (lgrand) then
c ---     volume ideal gas box is set arbitry large!
         boxlx(2)=1000*boxlx(1)          
         boxly(2)=1000*boxly(1)
         boxlz(2)=1000*boxlz(1)
      endif
      if (linit) then
         do ibox = 1,nbox
            if (lsolid(ibox) .and. .not. lrect(ibox)) then
               write(6,*) 'Cannot initialize non-rectangular system'
               stop
            endif
         enddo
      endif
      if ( linit ) then
         call initia(qelect)
         nnstep = 0
      else
         read (77,*) nnstep
         do im = 1,nbox
            do imol = 1,nmolty
               read(77,*) rmtrax(imol,im), rmtray(imol,im)
     &              , rmtraz(imol,im)
               read(77,*) rmrotx(imol,im), rmroty(imol,im)
     &              , rmrotz(imol,im)
            enddo
         enddo
         write(6,*) 
     +        'new maximum displacements read from restart-file'
         do im = 1,nbox
            write(6,*)'box      #',im
            do imol = 1,nmolty
               write(6,*) 'molecule type',imol
               write(6,1101) rmtrax(imol,im), rmtray(imol,im)
     &              , rmtraz(imol,im)
               write(6,1102) rmrotx(imol,im), rmroty(imol,im)
     &              , rmrotz(imol,im)
            enddo
         enddo

         do im=1,nbox
            read (77,*) (rmflcq(i,im),i=1,nmolty)
         enddo
         if ( lecho ) then
            do im=1,nbox
               write(6,*) 'maximum fluc q displacements: Box #',im
               write(6,*) (rmflcq(i,im),i=1,nmolty)
            enddo
         endif

         if ( lgibbs .or. lgrand .or. lnpt ) then
            read (77,*) (rmvol(ibox), ibox = 1,nbox)
            write(6,1103) (rmvol(ibox), ibox = 1,nbox)
            write(6,*)

            do ibox = 1,nbox
               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  read(77,*) (rmhmat(ibox,j),j=1,9)
                  write(6,*) (rmhmat(ibox,j),j=1,9)
               endif
            enddo

            if (lexzeo) then
               do ibox = 1,nbox
                  read(77,*) dum,dum,dum
               enddo
            else
               write(6,*) 'new box size read from restart-file'

               do ibox = 1,nbox

                  if (lsolid(ibox) .and. .not. lrect(ibox)) then

                     read(77,*) (hmat(ibox,j),j=1,9)

                     write(6,*) 'HMAT COORDINATES FOR BOX',ibox
                     write(6,*) (hmat(ibox,j),j=1,3)
                     write(6,*) (hmat(ibox,j),j=4,6)
                     write(6,*) (hmat(ibox,j),j=7,9)

                     vol = (hmat(ibox,1) * (hmat(ibox,5)
     +                    * hmat(ibox,9) - hmat(ibox,8)
     +                    * hmat(ibox,6)) + hmat(ibox,4)
     +                    * (hmat(ibox,8) * hmat(ibox,3)
     +                    - hmat(ibox,2) * hmat(ibox,9))
     +                    + hmat(ibox,7) * (hmat(ibox,2)
     +                    * hmat(ibox,6) - hmat(ibox,5)
     +                    *hmat(ibox,3)))

                     hmati(ibox,1) = (hmat(ibox,5)*hmat(ibox,9)
     &                    -hmat(ibox,8)*hmat(ibox,6))/vol
                     hmati(ibox,5) = (hmat(ibox,1)*hmat(ibox,9)
     &                    -hmat(ibox,7)*hmat(ibox,3))/vol
                     hmati(ibox,9) = (hmat(ibox,1)*hmat(ibox,5)
     &                    -hmat(ibox,4)*hmat(ibox,2))/vol
                     hmati(ibox,4) = (hmat(ibox,7)*hmat(ibox,6)
     &                    -hmat(ibox,4)*hmat(ibox,9))/vol
                     hmati(ibox,2) = (hmat(ibox,3)*hmat(ibox,8)
     &                    -hmat(ibox,2)*hmat(ibox,9))/vol
                     hmati(ibox,7) = (hmat(ibox,4)*hmat(ibox,8)
     &                    -hmat(ibox,7)*hmat(ibox,5))/vol
                     hmati(ibox,3) = (hmat(ibox,2)*hmat(ibox,6)
     &                    -hmat(ibox,3)*hmat(ibox,5))/vol
                     hmati(ibox,8) = (hmat(ibox,7)*hmat(ibox,2)
     &                    -hmat(ibox,8)*hmat(ibox,1))/vol
                     hmati(ibox,6) = (hmat(ibox,3)*hmat(ibox,4)
     &                    -hmat(ibox,6)*hmat(ibox,1))/vol

                     if ( vol .eq. 0 )
     +                    stop 'The volume of this box is 0!'
c * compare boxwidths to cutoff
c * v = u2 x u3                     
                     v(1) = hmat(ibox,5)*hmat(ibox,9) - 
     &                    hmat(ibox,6)*hmat(ibox,8)
                     v(2) = hmat(ibox,6)*hmat(ibox,7) - 
     &                    hmat(ibox,4)*hmat(ibox,9)
                     v(3) = hmat(ibox,4)*hmat(ibox,8) - 
     &                    hmat(ibox,5)*hmat(ibox,7)

                     w(1) = vol / dsqrt(v(1)**2 + v(2)**2 + v(3)**2)

c * v = u3 x u1                     
                     v(1) = hmat(ibox,8)*hmat(ibox,3) - 
     &                    hmat(ibox,9)*hmat(ibox,2)
                     v(2) = hmat(ibox,9)*hmat(ibox,1) - 
     &                    hmat(ibox,7)*hmat(ibox,3)
                     v(3) = hmat(ibox,7)*hmat(ibox,2) - 
     &                    hmat(ibox,8)*hmat(ibox,1)

                     w(2) = vol / dsqrt(v(1)**2 + v(2)**2 + v(3)**2)

c * v = u1 x u2
                     v(1) = hmat(ibox,2)*hmat(ibox,6) - 
     &                    hmat(ibox,3)*hmat(ibox,5)
                     v(2) = hmat(ibox,3)*hmat(ibox,4) - 
     &                    hmat(ibox,1)*hmat(ibox,6)
                     v(3) = hmat(ibox,1)*hmat(ibox,5) - 
     &                    hmat(ibox,2)*hmat(ibox,4)

                     w(3) = vol / dsqrt(v(1)**2 + v(2)**2 + v(3)**2)

                     write(6,1106) ibox,w(1),w(2),w(3)

                     if (rcut/w(1) .gt. 0.5d0 .or. 
     &                    rcut/w(2) .gt. 0.5d0 .or. 
     &                    rcut/w(3) .gt. 0.5d0) then
                        write(6,*) 'rcut > half cell width'
                        stop
                     endif
                     
                     boxlx(ibox) = hmat(ibox,1)
                     boxly(ibox) = dsqrt(hmat(ibox,4)*hmat(ibox,4)
     &                    +hmat(ibox,5)*hmat(ibox,5))
                     boxlz(ibox) = dsqrt(hmat(ibox,7)*hmat(ibox,7)
     &                    +hmat(ibox,8)*hmat(ibox,8)+
     &                    hmat(ibox,9)*hmat(ibox,9))
                     write(6,1104) ibox,boxlx(ibox),boxly(ibox),
     &                    boxlz(ibox)
                     write(6,1105) ibox,dacos(hmat(ibox,4)/boxly(ibox))
     &                    *180.0d0/onepi,180.0d0/onepi*
     &                    dacos(hmat(ibox,7)/boxlz(ibox)),180.0d0/onepi*
     &                    dacos((hmat(ibox,4)*hmat(ibox,7)
     &                    +hmat(ibox,5)*hmat(ibox,8))/
     &                    boxly(ibox)/boxlz(ibox))
                     
                  else

                     read (77,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
                     write(6,1104) ibox,
     &                    boxlx(ibox),boxly(ibox),boxlz(ibox)

                     do i = 1, nbox
                          if( (rcut/boxlx(i) .gt. 0.5d0).or.
     &                      (rcut/boxly(i) .gt. 0.5d0).or.
     &                      (rcut/boxlz(i) .gt. 0.5d0)) then
                             write(6,*) 'rcut > 0.5*boxlx'
                             stop
                          endif
                     enddo
   


                  endif

               enddo
            endif

         endif

         write(6,*)

         read (77,*) ncres
         read (77,*) nmtres
c --- check that number of particles in fort.4 & fort.77 agree ---
         if ( ncres .ne. nchain .or. nmtres .ne. nmolty ) then
            write(6,*)
     +         'conflicting information in restart and control files'
            write(6,*) 'nchain',nchain,'ncres',ncres
            write(6,*) 'nmolty',nmolty,'nmtres',nmtres
            stop
         endif
         read (77,*) (nures(i),i=1,nmtres)

c         do i = 1, nmtres
c            if ( nures(i) .ne. nunit(i) ) then
c               write(6,*)
c     +           'conflicting information in restart and control files'
c               write(6,*) 'unit',i,'nunit',nunit(i),'nures',nures(i)
c               stop
c            endif
c         enddo
c         write(6,*) 'ncres',ncres,'   nmtres',nmtres

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

c         write(6,*) 'start reading coordinates'
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
     &              .and. .not.lgrand) stop 
     &              'Particle found outside BOX 1'
               
               nchbox(ibox) = nchbox(ibox) + 1
               imolty = moltyp(i)
               ncmt(ibox,imolty) = ncmt(ibox,imolty) + 1
               if ( lexpand(imolty) ) then
                  if ( ibox .gt. 2 ) 
     &                 stop 'put in box 1 and 2 for such molecules'
                  itype = eetype(imolty)
                  ncmt2(ibox,imolty,itype) = 
     &                 ncmt2(ibox,imolty,itype) + 1
                  do j = 1,nunit(imolty)
                     sigma(imolty,j) = sigm(imolty,j,itype)
                     epsilon(imolty,j) = epsil(imolty,j,itype)
                  enddo
               endif 
            else
               write(6,*) 'i:',i,'nboxi(i)',nboxi(i)
               stop 'Particle found in ill-defined box'
            endif

         enddo

c --- check that number of particles of each type is consistent
         do i = 1, nmolty
            tcount = 0
            do ibox = 1, nbox
               tcount = tcount + ncmt(ibox,i)
            enddo
            if ( tcount .ne. temtyp(i) ) then
               write(6,*) 'Particle type number inconsistency'
               write(6,*) 'type',i
               write(6,*) 'ncmt',(ncmt(ibox,i), ibox = 1,nbox)
               write(6,*) 'temtyp', temtyp(i)
               stop
            endif
         enddo

c         write(6,*) 'particles found in correct box with correct type'

         do i = 1,nbxmax
            qbox(i) = 0.0d0
         enddo

         do 22 i = 1, nchain
c            write(6,*) 'reading coord of chain i',i
            imolty = moltyp(i)
            do j = 1, nunit(imolty)
               read (77,*) rxu(i,j), ryu(i,j), rzu(i,j),qqu(i,j)
               if (.not. lreadq) then
                  qqu(i,j) = qelect(ntype(imolty,j))
               endif
               qbox(nboxi(i)) = qbox(nboxi(i)) + qqu(i,j)
            enddo
            
 22      continue

         do i = 1,nbxmax
            if ( dabs(qbox(i)) .gt. 1d-6 ) then
               write(6,*) 'box',i,' has a net charge of',qbox(i)
               stop
            endif
         enddo

      endif

      if ( lchgall ) then
c *** if real space term are summed over all possible pairs in the box
c *** kalp(1) & kalp(2) are fixed while calp(1) & calp(2) change
c *** according to the boxlength, kalp(1) should have a value greater than
c *** 5.0
         do ibox = 1, nbox 
            calp(ibox) = kalp(ibox)/boxlx(ibox)
            if ( lewald ) then
               if ( kalp(ibox) .lt. 5.0d0 ) then
                  write(6,*) 'Warning, kalp is too small'
                  stop
               endif
c * removing kmax from the code
c               if (kmax(ibox) .lt. kalp(ibox) ) then
c                  write(6,*) 'Warning, kmax should not be',
c     &                 ' smaller than kalp'
c                  stop
c               endif
c * Ewald sum seems to work fine now
c               if ( dabs(boxlx(ibox)-boxly(ibox)) .gt. 1d-6 .or.
c     &              dabs(boxly(ibox)-boxlz(ibox)) .gt. 1d-6 ) then
c                  write(6,*) 
c     &            'Ewald-sum in this code for non-cubic still testing!'
c                  write(6,*) 'box dimension',boxlx(ibox),boxly(ibox)
c     &                 ,boxlz(ibox)
c                  stop
c               endif
            else
               stop 'lewald should be true when lchgall is true'
            endif
         enddo
      else
c *** if not lchgall, calp(1) & calp(2) are fixed
         do ibox = 1, nbox
            calp(ibox) = kalp(ibox)
            if ( lewald ) then

               write(6,*) '##################'
               write(6,*) 'Warning, the case',
     &              ' lchgall=false lewald=true'
               write(6,*) 'has not recently been checked for accuracy',
     &              ' in this code.'
               write(6,*) 'Do you really need it?'
               write(6,*) '##################'

               if ( calp(ibox)*rcutchg(ibox) .lt. 2.75d0 ) then
                  write(6,*) 'Warning, kalp too small in box',ibox
                  write(6,*) ibox,calp(ibox),rcutchg(ibox)
c                  stop 'kalp is too small'
               endif
               if ( kmax(ibox) .lt. calp(ibox)*boxlx(ibox) ) then
                  write(6,*) 'Warning, kmax is smaller than kalp*boxlx'
                  write(6,*) ibox,kmax(ibox),calp(ibox)*boxlx(ibox)
c                  stop 'kmax should not be smaller than kalp*boxlx'
               endif
            endif
         enddo
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
            write(6,*) 'Readdat: nmcount ne nchbox', nmcount, nchbox
            stop
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
c         write(6,*) 'it =',it,'   innew =',innew
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
               write(6,*) 'rintramax for the linkcell list too small'
               stop
            endif
         enddo
      endif

c * calculate number of frames *
      nnframe = nstep / imv
 
c *** write out movie-header ***
      if ( nnframe .ne. 0 ) then
         write(10,*) nnframe,nchain,nmolty,(rcut,ibox=1,nbox)
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
     &              ,(ijvib(imolty,ii,z),z=1,invib(imolty,ii))
            enddo

c     output torsional connectivity information
            do j = 1,nunit(imolty)
               write(10,*) intor(imolty,j)
               do ii = 1,intor(imolty,j)
                  write(10,*) ijtor2(imolty,j,ii)
     &                 ,ijtor3(imolty,j,ii),ijtor4(imolty,j,ii)
               enddo
            enddo
         enddo

      endif

c     --- write out isolute movie header
      lprint = .false.
      do imol = 1,nmolty
         if (lsolute(imol)) then
            lprint = .true.
         endif
      enddo

      if (lprint) then
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
     &       ,dum,dum,dum,nbox,nnstep,dum,.false.,.false.,.false.
     &       ,.true.,.false.,.false.,lratfix,lsolute)
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
         write(6,*) 'external potential for Langmuir monolayers used'
         write(6,*) 'zprmin',zprmin
         write(6,*) 'v2prmin',v2prmin,'v3prmin',v3prmin
      endif

C -------------------------------------------------------------------

c * write out connectivity and bonded interactions
      if (lverbose) then
         do imol = 1,nmolty
            write(6,*) 'molecule type',imol
            if (nunit(imol) .gt. 1) then
               write(6,*) '   i   j   type_i type_j   bond length',
     &              '        k/2'
            endif
            do i = 1,nunit(imol)
               do j = 1, invib(imol,i)
                  write(6,1014) i,ijvib(imol,i,j),ntype(imol,i),
     &                 ntype(imol,ijvib(imol,i,j)),
     &                 brvib(itvib(imol,i,j)),brvibk(itvib(imol,i,j))
               enddo
            enddo

            if (nunit(imol) .gt. 2) then
               write(6,*)
               write(6,*) '   i   j   k   type_i type_j type_k',
     &              '     angle      k/2'
            endif
            do i = 1,nunit(imol)
               do j = 1,inben(imol,i)
                  write(6,1015) i,ijben2(imol,i,j),
     &                 ijben3(imol,i,j),ntype(imol,i),
     &                 ntype(imol,ijben2(imol,i,j)),
     &                 ntype(imol,ijben3(imol,i,j)),
     &                 brben(itben(imol,i,j))*180.0d0/onepi,
     &                 brbenk(itben(imol,i,j))
               enddo
            enddo

            if (nunit(imol) .gt. 3) then
               write(6,*)
               write(6,*) '   i   j   k   l    type_i type_j type_k',
     &              ' type_l     torsion type'
            endif

            do i = 1,nunit(imol)
               do j = 1, intor(imol,i)
                  write(6,1016) i,ijtor2(imol,i,j),ijtor3(imol,i,j),
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
      write(6,*)
      if ((.not. lexpsix) .and. (.not. lmmff)) then
         if (lninesix) then
            write(6,*)
     & '     i   j    r_0,ij     epsij         q0(i)          q0(j)'
         else
            write(6,*) 
     & '     i   j    sig2ij     epsij         q0(i)          q0(j)'
         endif
         do i = 1, nntype
            do j = 1,nntype
               if ( lhere(i) .and. lhere(j) ) then
                  if (lninesix) then
                     ij = (i-1)*nxatom + j
                     write(6,'(3x,2i4,2f10.5,2f15.6)') i,j
     &                    ,rzero(ij),epsnx(ij),qelect(i),qelect(j)
                  else
                     ij = (i-1)*nntype + j
                     write(6,'(3x,2i4,2f10.5,2f15.6)') i,j
     &                    ,dsqrt(sig2ij(ij))
     &                    , epsij(ij),qelect(i),qelect(j)
                  endif
               endif
            enddo
         enddo
      endif

      if ( lneigh ) then
c *** calculate squares of nn-radius and max. update displacement ***
         rcnnsq = rcutnn * rcutnn
         upnn = ( rcutnn - rcut ) / 3.0d0
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
                 write(6,*) ' rmtrax greater than upnn',im,imol
                 rmtrax(imol,im) = upnn
              endif
              if ( rmtray(imol,im) .gt. upnn ) then
                 write(6,*) ' rmtray greater than upnn',im,imol
                 rmtray(imol,im) = upnn
              endif
              if ( rmtraz(imol,im) .gt. upnn ) then
                 write(6,*) ' rmtraz greater than upnn',im,imol
                 rmtraz(imol,im) = upnn
              endif
              
              if ( rmrotx(imol,im) .gt. upnndg ) then
                 write(6,*) ' rmrotx greater than upnndg',im,imol
                 rmrotx(imol,im) = upnndg
              endif
              if ( rmroty(imol,im) .gt. upnndg ) then
                 write(6,*) ' rmroty greater than upnndg',im,imol
                 rmroty(imol,im) = upnndg
              endif
              if ( rmrotz(imol,im) .gt. upnndg ) then
                 write(6,*) ' rmrotz greater than upnndg',im,imol
                 rmrotz(imol,im) = upnndg
              endif
           enddo
        enddo
      endif

c *** write input data to unit 6 for control ***
      write(6,*)
      write(6,*) 'number of mc cycles:            ', nstep
      write(6,*) 'number of chains:               ', nchain
      write(6,*)
      write(6,*) 'temperature:                    ', temp
      write(6,*)
      write(6,*) 'ex-pressure:                    ', express
      write(6,*)


C -------------------------------------------------------------------
 
      if ( rcut .ge. rcutnn .and. lneigh ) then
         write(6,*) ' rcut greater equal rcutnn '
         stop
      endif
 
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
 1104 format(' dimension box ',i1,':                ',3f12.6)
 1105 format(' angle of  box ',i1,':                ',3f12.6)
 1106 format(' width of  box ',i1,':                ',3f12.6)

      return
      end





