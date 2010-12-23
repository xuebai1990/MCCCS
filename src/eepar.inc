! EEPAR.INC

      logical::leemove,lmstate, leeacc
      integer::fmstate,sstate1,sstate2,mstate,nstate,ee_moltyp
     &       ,box_state,eepointp,eeirem,boxrem1,boxins1,ee_prob,nmolty1
      real(8)::wee_ratio,psi,ee_qqu,um_markov,eeratio,rminee

      dimension psi(smax),ee_qqu(numax,smax),box_state(smax)
     &         ,ee_moltyp(smax),um_markov(smax,smax),ee_prob(smax)
     &         ,rminee(nntype*nntype)

      common /expepar/ rminee,ee_qqu,psi,wee_ratio,lmstate,leemove
     &      ,fmstate,sstate1,sstate2,mstate,nstate,ee_moltyp,box_state
     &      ,eepointp,eeirem,boxrem1,boxins1,leeacc,um_markov,eeratio
     &      ,ee_prob,nmolty1
