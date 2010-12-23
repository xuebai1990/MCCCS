! MERCK.INC

      integer::natomtyp
      parameter (natomtyp = 3)
      real(8)::epsimmff,sigimmff,smmff,sigisq,ammff,nmmff
      real(8)::gmmff,alphammff,coru_cons,corp_cons 
      dimension epsimmff(natomtyp),sigimmff(natomtyp),smmff(natomtyp)
      dimension ammff(natomtyp),nmmff(natomtyp),gmmff(natomtyp)
      dimension sigisq(natomtyp),alphammff(natomtyp)
      dimension coru_cons(natomtyp),corp_cons(natomtyp) 
      common /merk/ epsimmff,sigimmff,smmff,sigisq,ammff,nmmff
     &   ,alphammff,gmmff,coru_cons,corp_cons 

 
