! INPAR.INC

      logical::lbranch
      integer::ininch,inix, iniy, iniz, inirot, inimix
      real(8)::zshift, dshift
      dimension inix(nbxmax),iniy(nbxmax),iniz(nbxmax),
     &          inirot(nbxmax),inimix(nbxmax)
      dimension zshift(nbxmax),dshift(nbxmax)
      dimension ininch(ntmax,nbxmax)
      dimension lbranch(ntmax)
      common /input2/ zshift,dshift,inix,iniy,iniz,inirot,inimix
     &               ,ininch,lbranch
