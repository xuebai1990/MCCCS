      program dielectric

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   calculates the dielectric constant
ccc   input file is dielectinput, contains T and V
ccc   gets dipole data from fort.27 - dipx, dipy, dipz
ccc   modified 10/27/08 to include extra term for when there is
ccc   an applied electric field
ccc   (see P. Kusalik, Mol. Phys. 81, 1994, 199-216)
ccc   revised and checked 3/19/09
ccc   included the forgotten cross terms in <mu^4> and <mu^2>^2!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer i, j, k, nmax, nnn, istep, nread, count
      parameter (nmax=8000000)
      logical lfield
      double precision field
      double precision dielect, eps0, V, qqfact, pi
      double precision const, const2, T, avgM, avgM2, avgM4, boxlx
      double precision avgMx, avgMy, avgMz
      double precision avgMx2, avgMy2, avgMz2
      double precision avgMx4, avgMy4, avgMz4
       
      double precision dipole(nmax), dipole2(nmax), dipole4(nmax)
      double precision avgdip(nmax), avgdip2(nmax), avgdip4(nmax)
      double precision dipx(nmax),dipy(nmax), dipz(nmax)
      double precision dipx2(nmax), dipy2(nmax), dipz2(nmax)
      double precision dipx4(nmax), dipy4(nmax), dipz4(nmax)
      double precision dip(nmax), dip2(nmax), dip4(nmax), tmp
      double precision avgdipx(nmax), avgdipy(nmax), avgdipz(nmax)
      double precision avgdipx2(nmax), avgdipy2(nmax), avgdipz2(nmax)
      double precision avgdipx4(nmax), avgdipy4(nmax), avgdipz4(nmax)
      double precision avgdipxy2(nmax), avgdipxz2(nmax), avgdipyz2(nmax)
      double precision avgdip22(nmax), avgMxy2, avgMxz2, avgMyz2
      double precision fluc(nmax), fluc2(nmax), dielect2
      double precision temp1, temp2, temp3

      pi = 4.0d0*datan(1.0d0)

      open(4,file='dielectinput')
c     read temperature
      read(4,*)
      read(4,*) T
c     read boxlength (calculation of dielect const only for NVT simulation)
      read(4,*)
      read(4,*) boxlx
      read(4,*)
c     read istep
      read(4,*) istep
      read(4,*)
c     read logical
      read(4,*) lfield
      if (lfield) then
c     read field strength (in V/A)
         read(4,*)
         read(4,*) field
c  convert V to K/e
c  1 V = J/C
c  1.602*10^-19 C = 1 e
c  1 J *(1 kJ/1000 J) * (6.022*10^23/mol) * (120.27 K/1 kJ/mol) = 7.2426*10^22 K
c  1 V = 11602.74036 K
         field = field*11602.74036d0
      endif
      
      V = boxlx**3

      do i=1,nmax
         avgdip(i) = 0.0d0
         avgdip2(i) = 0.0d0
         avgdip4(i) = 0.0d0
         avgdipx(i) = 0.0d0
         avgdipy(i) = 0.0d0
         avgdipz(i) = 0.0d0
         avgdipx2(i) = 0.0d0
         avgdipy2(i) = 0.0d0
         avgdipz2(i) = 0.0d0
         avgdipx4(i) = 0.0d0
         avgdipy4(i) = 0.0d0
         avgdipz4(i) = 0.0d0
         avgdipxy2(i) = 0.0d0
         avgdipxz2(i) = 0.0d0
         avgdipyz2(i) = 0.0d0
      enddo
      avgM = 0.0d0
      avgM2 = 0.0d0
      avgM4 = 0.0d0

      count=0

      avgMx = 0.0d0
      avgMy = 0.0d0
      avgMz = 0.0d0
      
      avgMx2 = 0.0d0
      avgMy2 = 0.0d0
      avgMz2 = 0.0d0

      avgMx4 = 0.0d0
      avgMy4 = 0.0d0
      avgMz4 = 0.0d0

      avgMxy2 = 0.0d0
      avgMxz2 = 0.0d0
      avgMyz2 = 0.0d0

c      write(6,*) 'pi ', pi
      write(6,*) 'T ', T, ' K'
      write(6,*) 'V ', V, ' A^3'
      write(6,*) 'lfield ?', lfield
      if (lfield) then
         write(6,*) 'field strength ', field/11602.74036d0
      endif

c --- eps0 in e^2/K*Angstrom
      eps0 = 4.762445019d-7
      
      qqfact = 1.0d0/(4.0d0*pi*eps0)
c      write(6,*) '1/(4*pi*eps0)', qqfact

c ---  kB = 1, V in A^3, T in K
      const = (4.0d0*pi)/(3*V*T)
      const2 = (4.0d0*pi)/(90*V*T**3)
      
c      write(6,*) '(4*pi)/(3*V*T) ', const
c      write(6,*) '(4*pi)/(90*V*T**3) ', const2

c      open(25)

      do i=1, nmax
         read(27,*,end=10) dipx(i), dipy(i), dipz(i)
c         write(6,*) i, dipx(i), dipy(i), dipz(i)
      enddo

 10   nread = i-1
      write(6,*) 'read ', nread, ' values of the dipole moment'      
       
      
      open(unit=23,file='dielect.dat')

      do j=1, nread

c ---    calculate <M> and <M^2> in order to calculate standard deviation
c         dipole(j) = dsqrt(dipx(j)**2+dipy(j)**2+dipz(j)**2)
c         dipole(j) = dipx(j)+dipy(j)+dipz(j) -- NO
         dipole2(j) = dipx(j)**2+dipy(j)**2+dipz(j)**2

         avgM = avgM+dipole(j)
         avgdip(j) = avgM/dble(j)

         avgM2 =  avgM2 + dipole2(j)
         avgdip2(j) = avgM2/dble(j)


c ---    calculate <Mx>, <My>, <Mz> and <Mx^2>, <My^2>, <Mz^2>
c ---    in order to calculate the dielectric constant

         avgMx = avgMx + dipx(j)
         avgMy = avgMy + dipy(j)
         avgMz = avgMz + dipz(j)

         avgMx2 = avgMx2 + dipx(j)**2
         avgMy2 = avgMy2 + dipy(j)**2
         avgMz2 = avgMz2 + dipz(j)**2

         avgdipx(j) = avgMx/dble(j)
         avgdipy(j) = avgMy/dble(j)
         avgdipz(j) = avgMz/dble(j)

         avgdipx2(j) = avgMx2/dble(j)
         avgdipy2(j) = avgMy2/dble(j)
         avgdipz2(j) = avgMz2/dble(j)


c --- <mu^2> - <mu>^2
         fluc(j) = avgdipx2(j) + avgdipy2(j) + avgdipz2(j) 
     &        - avgdipx(j)**2 - avgdipy(j)**2 - avgdipz(j)**2
c         fluc(j) = avgdip2(j) - avgdip(j)**2
c --- avgdip2 = avgdipx2 + avgdipy2 + avgdipz2
c --- avgdip**2 != avgdipx**2 + avgdipy**2 + avgdipz**2        


         if (lfield) then
c            dipole4(j) = dipx(j)**4 + dipy(j)**4 + dipz(j)**4
            dipole4(j) = dipole(j)**4

c            avgM4 = avgM4 + dipole4(j)
c            avgdip4(j) = avgM4/dble(j)

            avgMx4 = avgMx4 + dipx(j)**4
            avgMy4 = avgMy4 + dipy(j)**4
            avgMz4 = avgMz4 + dipz(j)**4

            avgdipx4(j) = avgMx4/dble(j)
            avgdipy4(j) = avgMy4/dble(j)
            avgdipz4(j) = avgMz4/dble(j)
            
            avgMxy2 = avgMxy2 + (dipx(j)**2)*(dipy(j)**2)
            avgMxz2 = avgMxz2 + (dipx(j)**2)*(dipz(j)**2)
            avgMyz2 = avgMyz2 + (dipy(j)**2)*(dipz(j)**2)

            avgdipxy2(j) = avgMxy2/dble(j)
            avgdipxz2(j) = avgMxz2/dble(j)
            avgdipyz2(j) = avgMyz2/dble(j)


ccc --- <mu^4> = <mux^4> + <muy^4> + <muz^4> +
ccc ---   2<mux^2muy^2> + 2<mux^2muz^2> + 2<muy^2muz^2>
            avgdip4(j) = avgdipx4(j) + avgdipy4(j) + avgdipz4(j) +
     &           2*avgdipxy2(j) + 2*avgdipxz2(j) + 2*avgdipyz2(j)

ccc --- <mu^2>^2 = <mux^2>^2 + <muy^2>^2 + <muz^2>^2 +
ccc ---   2<mux^2><muy^2> + 2<mux^2><muz^2> + 2<muy^2><muz^2>
            avgdip22(j) = avgdipx2(j)**2 + avgdipy2(j)**2 + 
     &           avgdipz2(j)**2 + 2*avgdipx2(j)*avgdipy2(j) + 
     &           2*avgdipx2(j)*avgdipz2(j) + 
     &           2*avgdipy2(j)*avgdipz2(j)
            

            fluc2(j) = 3*avgdip4(j) - 5*avgdip22(j)

            dielect = 1.0d0 + qqfact*const*fluc(j) + 
     &           qqfact*const2*fluc2(j)*field**2


        
         else

            dielect = 1.0d0 + qqfact*const*fluc(j)
            
         endif  ! end if(lfield)

         if (mod(j,istep).eq.0) then
            write(23,*) j, dielect
         endif

      enddo  ! end do j=1,nread

c --  calculate standard deviation
c --  sigma = dsqrt(<M^2> - <M>^2)

c      open(unit=17, file='dielect.out')

      write(6,*) 'dielectric constant: ', dielect
c         write(6,*) 'standard deviation: ', dsqrt(avgdip2(nread) 
c     &        - avgdip(nread)**2)
      write(6,*) 'average dipole moment x, y, z: ', 
     &     avgdipx(nread), avgdipy(nread), avgdipz(nread)
      write(6,*) 'first term: ', avgdipx2(nread)+avgdipy2(nread)+
     &     avgdipz2(nread)
      write(6,*) 'second term: ', avgdipx(nread)**2+avgdipy(nread)**2
     &     +avgdipz(nread)**2
      write(6,*) 'fluc ', fluc(nread)
         
      write(23,*)

      close(4)
      close(27)
      close(23)


      end
