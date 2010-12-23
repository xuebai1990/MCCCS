! FIX.INC

      integer::  iend,ipast,endnum,pastnum,fclose,fcount,iwbef
     &   ,ibef,befnum,wbefnum,nextnum,inext

      logical::lcrank

      real(8)::xx,yy,zz,xvec,yvec,zvec,distij,vtgtr,vtbend
     &   ,bsum_tor,vtvib

      dimension iend(numax),ipast(numax,numax)
     &   ,pastnum(numax),fclose(numax,numax)
     &   ,fcount(numax),iwbef(numax),ibef(numax,numax)
     &   ,befnum(numax),xx(numax),yy(numax),zz(numax)
     &   ,xvec(numax,numax),yvec(numax,numax)
     &   ,zvec(numax,numax),distij(numax,numax)
     &   ,nextnum(numax),inext(numax,numax),vtvib(nchmax)
     &   ,vtgtr(nchmax),vtbend(nchmax),bsum_tor(nchmax)

      common /fecba/ iend,ipast,endnum,pastnum,fclose,fcount
     &   ,befnum,wbefnum,iwbef,ibef,nextnum
     &   ,inext,lcrank
   
      common /fecbb/ xx,yy,zz,xvec,yvec,zvec,distij,vtgtr,vtbend
     &   ,bsum_tor,vtvib

