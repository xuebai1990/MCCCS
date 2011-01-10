      module poten
      implicit none
      save
      logical::lq14scale
      logical::L_Coul_CBMC
      logical::L_tor_table, L_spline, L_linear
      integer::a15type
      logical,allocatable::lspecial(:),lpl(:),lij(:),lqchg(:)
      real(8),allocatable::sig2ij(:),epsij(:),epsi(:),sigi(:),qelect(:)
     &     ,xiq(:),jayself(:),extc12(:),extc3(:),extz0(:),aspecial(:)
     &     ,bspecial(:),q1(:),ecut(:),jayq(:)
      real(8)::qqfact,fqegp,a15,qscale,epsilon,sigma
      real(8)::ljscale,qscale2
      
      dimension fqegp(ntmax),a15(2),a15type(ntmax,numax,numax)
      dimension epsilon(2,numax),sigma(2,numax)
      dimension ljscale(ntmax,numax,numax),qscale2(ntmax,numax,numax)
      dimension lq14scale(ntmax),qscale(ntmax)
      real(8)::n0,n1

! ***********************************************************
! *** parameters for Generalized Lennard Jones Potential  ***
! repulsive part
      parameter(n0=12.0d0)
!   attractive part
      parameter(n1=6.0d0)  
! ****  Ref:  J. Chem. Phys. 120, 4994 (2004)         ****
! ***********************************************************
 
! *** common blocks *** 
!$$$      common /ljpot/ sig2ij,epsij,ecut,epsilon,sigma,epsi,sigi
!$$$     &   ,q1,qscale,a15,a15type
!$$$      common /ljepot/ ljscale,qscale2,lq14scale,qelect
!$$$      common /extpot/ extc12, extc3, extz0
!$$$      common /specail/ aspecial,bspecial,lspecial,lpl
!$$$      common /electr/ qqfact,xiq,jayq,jayself,fqegp,lqchg,lij
!$$$      common /etects/ L_Coul_CBMC,L_tor_table, L_spline, L_linear

      end module poten
