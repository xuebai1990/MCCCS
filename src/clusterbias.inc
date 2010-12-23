! CLUSTERBIAS.INC

      logical::lbias(ntmax),lavbmc1(ntmax),lavbmc2(ntmax),lavbmc3(ntmax)
      real(8)::rbsmax,rbsmin,vol_eff,pmbias(ntmax),pmbsmt(ntmax)
     &                 ,pmbias2(ntmax)
      common /cbiasa/ lbias,lavbmc1,lavbmc2,lavbmc3
      common /cbiasb/ rbsmax,rbsmin,vol_eff,pmbias,pmbias2,pmbsmt
