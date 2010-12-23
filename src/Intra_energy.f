      subroutine Intra_energy ( i,imolty, v, vintra, vinter,vext
     &     ,velect,vewald,flagon,ibox, istart, iend,lljii,ovrlap
     &     ,ltors,vtors,lcharge_table,lfavor,vvib,vbend,vtg)

 
!    *******************************************************************
!    ** calculates the total potential energy for a configuration.    **
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

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'system.inc'
!$$$      include 'neigh.inc'
!$$$      include 'poten.inc'
!$$$      include 'coord2.inc' 
!$$$      include 'external.inc'
!$$$      include 'connect.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'fepsi.inc'
!$$$      include 'qqlist.inc'
!$$$      include 'clusterbias.inc'
!$$$      include 'nsix.inc'
!$$$      include 'peboco.inc'
!$$$      include 'cell.inc'
!$$$      include 'tabulated.inc'

      logical::lqimol,lqjmol,lexplt,lcoulo,lfavor,lij2,liji,lqchgi
      logical::lljii,ovrlap,ltors,lcharge_table,lt,lfound

      integer(KIND=int)::growii,growjj,k,cellinc,jcell,ic,nmole
      integer(KIND=int)::i,ibox, istart, iend,ii,ntii,flagon,jjj,iii
     &       ,j,jj,ntjj,ntij,ntj,imolty,jmolty,ncell
      integer(KIND=int)::iivib,jjtor,ip1,ip2,ip3,it,nchp2,acellinc

      integer(KIND=int)::jjvib,jjben  

      real(KIND=double_precision)::vvib,vbend,vtg,theta

      real(KIND=double_precision)::ljsami,ljpsur,ljmuir,v,vintra, vinter,vext 
     &                ,rcutsq,rminsq,rxui,rzui,ryui,rxuij,rcinsq
     &                ,ryuij,rzuij,sr2,sr6,rij,rijsq,dzui,dz3,dz12
     &                ,exgrph,exsami,exmuir,exzeo,vtors,exsix,velect
     &                ,vewald,mmff,rbcut,rvdwsq,rchgsq,ninesix, genlj
      real(KIND=double_precision)::erfunc,qave
      real(KIND=double_precision)::xvec,yvec,zvec,xaa1,yaa1,zaa1,xa1a2,ya1a2,za1a2
     &     ,daa1,da1a2,dot,thetac,vtorso
      real(KIND=double_precision)::xcmi,ycmi,zcmi,rcmi,rcm,rcmsq,epsilon2,sigma2
      real(KIND=double_precision)::sx,sy,sz
      real(KIND=double_precision)::slitpore	
      real(KIND=double_precision)::distij(numax,numax)
      real(KIND=double_precision)::xcc,ycc,zcc,tcc,spltor, tabulated_vib

      dimension xvec(numax,numax),yvec(numax,numax),zvec(numax,numax)
      dimension lcoulo(numax,numax),cellinc(cmax),jcell(nmax)
      dimension acellinc(numax,27)

! --------------------------------------------------------------------

!      write(iou,*) 'start ENERGY'
      if ( lpbc ) call setpbc (ibox)

      rcutsq = rcut(ibox) * rcut(ibox)
      rbcut = rcut(ibox)
      if (ldual) rcinsq = rcutin*rcutin
      rminsq = rmin * rmin

      v = 0.0d0
      vvib = 0.0d0
      vbend = 0.0d0
      vtg = 0.0d0
 
! - branched and linear molecules with connectivity table -
! - go through entire chain -
! - calculate all bonds vectors and lengths
! - calculate all stretching, bending, and torsional potentials
! - that have an end-bead with an index smaller than the current bead
             do ii = 1, nunit(imolty)
                rxui=rxu(i,ii)
                ryui=ryu(i,ii)
                rzui=rzu(i,ii)
                do iivib = 1, invib(imolty,ii)
                   jj = ijvib(imolty,ii,iivib)
!                   xvec(ii,jj) = rxu(i,jj) - rxui
!                   yvec(ii,jj) = ryu(i,jj) - ryui
!                   zvec(ii,jj) = rzu(i,jj) - rzui
                   xvec(ii,jj) = rxu(i,jj) - rxui
                   yvec(ii,jj) = ryu(i,jj) - ryui
                   zvec(ii,jj) = rzu(i,jj) - rzui 
                   distij(ii,jj) = dsqrt( xvec(ii,jj)**2
     &                 + yvec(ii,jj)**2 + zvec(ii,jj)**2 )

                   if ( nunit(imolty) .ne. nugrow(imolty) )then
!                  --- account for explct atoms in opposite direction
                      xvec(jj,ii)   = -xvec(ii,jj)
                      yvec(jj,ii)   = -yvec(ii,jj)
                      zvec(jj,ii)   = -zvec(ii,jj)
                      distij(jj,ii) = distij(ii,jj)
                   end if
                end do
             end do

! - stretching -
!             if ( brvibk(1) .gt. 0.01d0 .or. lninesix  ) then
                do j = 2, nunit(imolty)
                   do jjvib = 1, invib(imolty,j)
                      ip1 = ijvib(imolty,j,jjvib)
                      it  = itvib(imolty,j,jjvib)
                      if (L_vib_table) then
                         call lininter_vib(distij(ip1,j), 
     &                        tabulated_vib, it)
                         vvib = vvib + tabulated_vib
!                         write(2,*) 'INTRA_ENERGY VVIB: ', 
!     &                        tabulated_vib
                      end if
                      if ( ip1 .lt. j .and..not.L_vib_table) vvib = vvib
     &                 + brvibk(it) * ( distij(ip1,j) - brvib(it) )**2
                   end do
                end do
!             end if


! - bending -
! ### molecule with bond bending
             do j = 2, nunit(imolty)
                do jjben = 1, inben(imolty,j)
                   ip2 = ijben3(imolty,j,jjben)
                   if ( ip2 .lt. j ) then
                      ip1 = ijben2(imolty,j,jjben)
                      it  = itben(imolty,j,jjben)
                      thetac = ( xvec(ip1,j)*xvec(ip1,ip2) +
     &                     yvec(ip1,j)*yvec(ip1,ip2) +
     &                     zvec(ip1,j)*zvec(ip1,ip2) ) /
     &                     ( distij(ip1,j)*distij(ip1,ip2) )
                      if ( thetac .ge. 1.0d0 ) thetac = 1.0d0
                      if ( thetac .le. -1.0d0 ) thetac = -1.0d0

                      theta = dacos(thetac)
                      vbend = vbend +
     &                     brbenk(it) * (theta-brben(it))**2

!                      write(iou,*) 'ip2,ip1,j',ip2,ip1,j
!                      write(iou,*) 'bend energy, theta '
!     &                     ,brbenk(it) * (theta-brben(it))**2,theta
                   end if
                end do
             end do

! - torsions -
! ### molecule with dihedral potenials ###
             do j = 2, nunit(imolty)
                do jjtor = 1, intor(imolty,j)
                   ip3 = ijtor4(imolty,j,jjtor)
                   if ( ip3 .lt. j ) then
                      ip1 = ijtor2(imolty,j,jjtor)
                      ip2 = ijtor3(imolty,j,jjtor)
                      it  = ittor(imolty,j,jjtor)
!*** calculate cross products d_a x d_a-1 and d_a-1 x d_a-2 ***
                      xaa1 = yvec(ip1,j) * zvec(ip2,ip1) +
     &                     zvec(ip1,j) * yvec(ip1,ip2)
                      yaa1 = zvec(ip1,j) * xvec(ip2,ip1) +
     &                     xvec(ip1,j) * zvec(ip1,ip2)
                      zaa1 = xvec(ip1,j) * yvec(ip2,ip1) +
     &                     yvec(ip1,j) * xvec(ip1,ip2)
                      xa1a2 = yvec(ip1,ip2) * zvec(ip2,ip3) +
     &                     zvec(ip1,ip2) * yvec(ip3,ip2)
                      ya1a2 = zvec(ip1,ip2) * xvec(ip2,ip3) +
     &                     xvec(ip1,ip2) * zvec(ip3,ip2)
                      za1a2 = xvec(ip1,ip2) * yvec(ip2,ip3) +
     &                     yvec(ip1,ip2) * xvec(ip3,ip2)
! *** calculate lengths of cross products ***
                      daa1 = dsqrt(xaa1**2+yaa1**2+zaa1**2)
                      da1a2 = dsqrt(xa1a2**2+ya1a2**2+za1a2**2)
! *** calculate dot product of cross products ***
                      dot = xaa1*xa1a2 + yaa1*ya1a2 + zaa1*za1a2
                      thetac = - dot / ( daa1 * da1a2 )
!     KEA -- added for extending range to +/- 180
!     additional definitions for torsions
                     if (L_tor_table) then
!     *** calculate cross product of cross products ***
                       xcc = yaa1*za1a2 - zaa1*ya1a2
                       ycc = zaa1*xa1a2 - xaa1*za1a2
                       zcc = xaa1*ya1a2 - yaa1*xa1a2
!     *** calculate scalar triple product ***
                       tcc = xcc*xvec(ip1,ip2) + ycc*yvec(ip1,ip2)
     &                    + zcc*zvec(ip1,ip2)
                       theta = dacos(thetac)
                       if (tcc .lt. 0.0d0) theta = -theta
                       if (L_spline) then
                          call splint(theta,spltor,it)
                       elseif(L_linear) then
                          call lininter(theta,spltor,it)
                       end if
                       vtg = vtg + spltor
                     else
                      vtg = vtg + vtorso( thetac, it )
!                       write(17,*) j,it,vtg,
!     &                           vtorso( thetac, it )
                     end if
                   end if
                end do
             end do

!----------------------------------------
     
            velect = velect*qqfact
            vewald = vewald*qqfact

 
!      velect = velect * qqfact
!      vewald = vewald * qqfact

!     note that vintra is only computed when the flag lljii is true
      v = vinter + vext + vintra + velect + vewald + vvib + vbend + vtg 
!  NEERAJ: Debugging start

!      write(iou,*) 'vinter:',vinter,'vext:',vext,'vintra:',vintra,'velect'
!     & ,velect,'vewald:'vewald

!  NEERAJ: DEbugging end

!      write(iou,*) 'end ENERGY'

      return
      end





