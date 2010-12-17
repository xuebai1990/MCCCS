! SYSTEM.INC

      logical::lsolid,lrect,ltwice,lideal
      integer::icell,nicell,iucell,solcount
      real(8)::boxlx, boxly, boxlz, rcut, rmin
     &   ,rcutnn, softcut,softlog,rcutin
     &   ,avsolinter,avsolintra,avsolbend,avsoltor,avsolelc
      dimension icell(nmax),nicell(cmax),iucell(cmax,cmaxa)
     &   ,solcount(nbxmax,ntmax)
      dimension avsolinter(nbxmax,ntmax),avsolintra(nbxmax,ntmax)
     &   ,avsolbend(nbxmax,ntmax),avsoltor(nbxmax,ntmax)
     &   ,avsolelc(nbxmax,ntmax)

      common /boxpha/ lsolid(nbxmax),lrect(nbxmax),ltwice(nbxmax)
      common /boxdim/ boxlx(nbxmax),boxly(nbxmax),boxlz(nbxmax),
     &                 rcut(nbxmax),rcutnn(nbxmax),lideal(nbxmax)

      common /solav/ avsolinter,avsolintra,avsolbend,avsoltor,avsolelc
      common /input1/ rmin, softcut,softlog
     &  ,rcutin,solcount

      common /linkc/ icell,nicell,iucell
