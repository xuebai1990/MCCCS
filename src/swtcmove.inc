! SWTCMOVE.INC

      integer::nswaty,nswatb,gswatc,nswach
      integer::nsampos,splist,nswtcb,box3,box4,ncut

      real(8)::bnswat,bniswat,pmswtcb
      real(8)::bnswat_empty	

      dimension nswach(npamax),nsampos(npamax)
      dimension nswatb(npamax,npabmax)
      dimension bnswat(npamax,npabmax),ncut(npamax,2)
      dimension gswatc(npamax,2,2*npamax),splist(npamax,numax,2)
      dimension nswtcb(npamax)
      dimension pmswtcb(npamax,npabmax)
      dimension box3(npamax,npabmax)
      dimension box4(npamax,npabmax)
      dimension bnswat_empty(npamax,npabmax)

      common /swtcma/ nswaty,nswtcb
      common /swtcmb/ nswach,nsampos
      common /swtcmc/ box3,box4,bnswat,bnswat_empty,nswatb
      common /swtcmd/ gswatc,splist,ncut
      common /swtcme/ pmswtcb
