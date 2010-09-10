      function v_elect_field(i, j, rzfield,E)
c **********************************************************************
c **  calculates interaction of molecule i with an external field E  ***
c **  added 06/24/07 by KM                                           ***
c **********************************************************************
      implicit none

c *** common blocks ***
      include 'control.inc'
      include 'external.inc'
      include 'coord.inc'

      real(8)::v_elect_field, convert, rzfield, E
      integer::i, j,ibox


c ********************************************
c **  units
c **  E in V/A, q in e, rz in A
c **  E*q*rz = V*e
c **  1 V*e = 11600 K
c ********************************************

      v_elect_field = -E*rzfield*qqu(i,j)  
      
c      write(6,*) 'E ', E, ' exfield ', exfield

      return
      end
