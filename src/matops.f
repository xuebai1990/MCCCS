      subroutine matops(ibox)
     

!     calculates boxlengths: cell_length ,minimum boxwidths: min_length
!     ,boxvolume: cell_vol,inverse h matrix for 
!     non cubic simulation cell
!     Written by Neeraj Rai (in Merck Apr 2005)

      implicit none
      include 'control.inc'
      include 'cell.inc'
      include 'conver.inc'       
      include 'system.inc' 
 
      real(8)::abx,aby,abz,bcx,bcy,bcz,cax,cay,caz
      real(8)::elem(9),inv_vol,adj(9),cosa,cosb,cosg
      integer::ibox,i,j

      do i=1,9
         elem(i)=hmat(ibox,i)
      enddo
   
!         -- calculating the length of cell vectors 
           
      cell_length(ibox,1)=sqrt(elem(1)*elem(1)+elem(2)*elem(2)+
     +                          elem(3)*elem(3))
      cell_length(ibox,2)=sqrt(elem(4)*elem(4)+elem(5)*elem(5)+
     +                          elem(6)*elem(6))
      cell_length(ibox,3)=sqrt(elem(7)*elem(7)+elem(8)*elem(8)+
     +                          elem(9)*elem(9))

      boxlx(ibox) = cell_length(ibox,1)
      boxly(ibox) = cell_length(ibox,2)
      boxlz(ibox) = cell_length(ibox,3) 

!    -- calculating cross product of cell vectors

      abx=elem(2)*elem(6)-elem(3)*elem(5)
      aby=elem(3)*elem(4)-elem(1)*elem(6)
      abz=elem(1)*elem(5)-elem(2)*elem(4)
      bcx=elem(5)*elem(9)-elem(6)*elem(8)
      bcy=elem(6)*elem(7)-elem(4)*elem(9)
      bcz=elem(4)*elem(8)-elem(5)*elem(7)
      cax=elem(8)*elem(3)-elem(2)*elem(9)
      cay=elem(1)*elem(9)-elem(3)*elem(7)
      caz=elem(2)*elem(7)-elem(1)*elem(8)

!    -- calculating cell volume
   
      cell_vol(ibox) = elem(1)*bcx+elem(2)*bcy+elem(3)*bcz
      
      if(abs(cell_vol(ibox)).lt.1D-16) then 
         call cleanup('Volume of cell negligible, check input H matrix')
      endif
              
      inv_vol = 1.0d0/cell_vol(ibox) 

!    -- calculating minimum cell widths

      min_width(ibox,1) = cell_vol(ibox)/sqrt(bcx*bcx+bcy*bcy+
     +                                                 bcz*bcz)
      min_width(ibox,2) = cell_vol(ibox)/sqrt(cax*cax+cay*cay+
     +                                                 caz*caz)
      min_width(ibox,3) = cell_vol(ibox)/sqrt(abx*abx+aby*aby+
     +                                                 abz*abz)

!    -- calculating adjoint for inverting the h-matrix
               
      adj(1)=elem(5)*elem(9)-elem(3)*elem(8)
      adj(2)=elem(3)*elem(8)-elem(2)*elem(9)
      adj(3)=elem(2)*elem(6)-elem(3)*elem(5)
      adj(4)=elem(6)*elem(7)-elem(4)*elem(9)
      adj(5)=elem(1)*elem(9)-elem(3)*elem(7)
      adj(6)=elem(3)*elem(4)-elem(1)*elem(6)
      adj(7)=elem(4)*elem(8)-elem(5)*elem(7)
      adj(8)=elem(2)*elem(7)-elem(1)*elem(8)
      adj(9)=elem(1)*elem(5)-elem(2)*elem(4)
   
!     --  inverting the matrix


      hmati(ibox,1) = inv_vol * adj(1)
      hmati(ibox,2) = inv_vol * adj(2)
      hmati(ibox,3) = inv_vol * adj(3)
      hmati(ibox,4) = inv_vol * adj(4)
      hmati(ibox,5) = inv_vol * adj(5)
      hmati(ibox,6) = inv_vol * adj(6)
      hmati(ibox,7) = inv_vol * adj(7)
      hmati(ibox,8) = inv_vol * adj(8)
      hmati(ibox,9) = inv_vol * adj(9)
           
!   ---  calculating alpha, beta and gamma using the dot product rule

!   ---  cell_ang(ibox,1)=alpha;2= beta; 3= gamma
!   ---  In crystallography Literature angle between b & c is alpha (1),
!        angle between c & a is Beta (2) and angle between a & b is gamma (3)


      cosa = (elem(4)*elem(7)+elem(5)*elem(8)+elem(6)*elem(9))/
     +           (cell_length(ibox,2)*cell_length(ibox,3))

      cosb = (elem(1)*elem(7)+elem(2)*elem(8)+elem(3)*elem(9))/
     +           (cell_length(ibox,1)*cell_length(ibox,3))

      cosg = (elem(4)*elem(1)+elem(5)*elem(2)+elem(6)*elem(3))/
     +           (cell_length(ibox,2)*cell_length(ibox,1))

      cell_ang(ibox,1) = dacos(cosa)
      cell_ang(ibox,2) = dacos(cosb)
      cell_ang(ibox,3) = dacos(cosg)
          
      return
      end subroutine matops

