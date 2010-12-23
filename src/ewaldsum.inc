! EWALDSUM.INC

      logical::L_Ewald_Auto
      integer::kmax,numvect
! - kalp, a parameter to control the real::space sum
! - kmax, a parameter to control the total number of reciprocal vectors
! - numvect, the total number of reciprocal vectors
      dimension numvect(nbxmax)
      real(8)::kx,ky,kz,prefact,ssumr,ssumi,ssumrn,ssumin,
     &                 ssumro,ssumio
      real(8)::sself,correct,kalp,calp,dipolex,dipoley,dipolez,
     &                 dipolexo,dipoleyo,dipolezo
! - calp = kalp / boxlen
      dimension kx(vectormax,nbxmax),ky(vectormax,nbxmax),
     &          kz(vectormax,nbxmax),prefact(vectormax,nbxmax),
     &          ssumr(vectormax,nbxmax),ssumi(vectormax,nbxmax),
     &          ssumrn(vectormax,nbxmax),ssumin(vectormax,nbxmax),
     &          calp(nbxmax),kmax(nbxmax),kalp(nbxmax),
     &          ssumro(vectormax,nbxmax),ssumio(vectormax,nbxmax)
      dimension dipolex(nbxmax),dipoley(nbxmax),dipolez(nbxmax)

      common /ewaldl/ L_Ewald_Auto 
      common /numvect/ numvect,kmax 
      common /vector/ kx,ky,kz
      common /ssum/ calp,kalp,ssumr,ssumi,prefact,ssumrn,
     &                  ssumin,ssumro,ssumio
      common /selfterm/ sself,correct
      common /dipoterm/ dipolex,dipoley,dipolez,dipolexo,dipoleyo,
     &                  dipolezo
