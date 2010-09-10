      integer::nw,nwell,iratipsw
      parameter (nw = 4000)

      real(8)::dvdl,acdvdl,acipsw,vipsw,pipsw,vwellipsw
     &   ,pwellipsw,vwellipswb,vipswb,etais,lambdais,rxwell,rywell
     &   ,rzwell,awell,bwell,vipswo,vipswn,vwellipswo,vwellipswn
     &   ,lena,lenc,pwellips,pips,dhmat,sxwell,sywell,szwell
     &   ,vwellipswot,vwellipswnt,vipswnt,vipswot

      logical::lwell,lstagea,lstageb,lstagec

      dimension vipswb(nbxmax),vwellipswb(nbxmax),rxwell(nw,ntmax)
     &         ,rywell(nw,ntmax)
     &         ,rzwell(nw,ntmax),pwellips(3,3),pips(3,3),nwell(ntmax)
     &         ,lwell(ntmax),awell(numax,numax,ntmax),dhmat(3,3)
     &         ,sxwell(nw,ntmax),sywell(nw,ntmax),szwell(nw,ntmax)
     &         ,vwellipswot(nchmax),vwellipswnt(nchmax)
     &         ,vipswot(nchmax),vipswnt(nchmax)

      common /paripsw1/ dvdl,acdvdl,acipsw,vipsw,pipsw,vwellipsw
     &   ,pwellips,pips,dhmat,sxwell,sywell,szwell
     &   ,vwellipswot,vwellipswnt,vipswnt,vipswot
     &   ,pwellipsw,vwellipswb,vipswb,etais,lambdais,rxwell,rywell

      common /paripsw2/ rzwell,nwell

      common /paripsw5/ awell,bwell
          
      common /paripsw3/ vipswo,vipswn,vwellipswo,vwellipswn,
     &                lena,lenc,iratipsw

      common /paripsw4/ lwell,lstagea,lstageb,lstagec
