! BNBSMA.INC

      real(8)::bsvol,bnvol,bntrax,bntray,bntraz,bshmat
     &  ,bstrax,bstray,bstraz,bnrotx,bnroty,bnrotz,bnhmat
     &  ,bsrotx,bsroty,bsrotz,bnflcq,bsflcq,bnexpc,bsexpc
     &  ,bnflcq2,bsflcq2

       real(8)::Abntrax,Abntray,Abntraz,   
     &                  Abstrax,Abstray,Abstraz

      dimension bsvol(nbxmax), bnvol(nbxmax),bshmat(nbxmax,9),
     &	        bnhmat(nbxmax,9)
      dimension bntrax(ntmax,nbxmax),bntray(ntmax,nbxmax),
     &          bntraz(ntmax,nbxmax)
      dimension bstrax(ntmax,nbxmax),bstray(ntmax,nbxmax),
     &          bstraz(ntmax,nbxmax)
      dimension bnrotx(ntmax,nbxmax),bnroty(ntmax,nbxmax),
     &          bnrotz(ntmax,nbxmax)
      dimension bsrotx(ntmax,nbxmax),bsroty(ntmax,nbxmax),
     &          bsrotz(ntmax,nbxmax)
      dimension bnexpc(ntmax,nbxmax),bsexpc(ntmax,nbxmax)
      dimension bnflcq(ntmax,nbxmax),bsflcq(ntmax,nbxmax),
     &          bnflcq2(ntmax,nbxmax),bsflcq2(ntmax,nbxmax)


      common /bvolum/ bnvol,bsvol,bshmat,bnhmat
      common /btrans/ bntrax,bntray,bntraz,bstrax,bstray,bstraz
      common /brotat/ bnrotx,bnroty,bnrotz,bsrotx,bsroty,bsrotz
      common /bflucq/ bnflcq,bsflcq,bnexpc,bsexpc,bnflcq2,bsflcq2
      common /atomtra/ Abntrax,Abntray,Abntraz,Abstrax,
     &          Abstray,Abstraz  
