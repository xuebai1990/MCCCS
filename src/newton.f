	subroutine newton(x_in,x_out)

C  --  This subroutine calcultes the roots of non-linear equation using 
C  --  Newton-Raphson method.

       implicit none
  
       integer::i,j,maxiter
       real(8)::x_in,x_out,erfunc
       real(8)::root(11),sqrtpi,deno

       parameter (maxiter = 10)
 
       root(1) = x_in
 
       sqrtpi = sqrt(4.0d0*datan(1.0d0)) 

       do i = 1,maxiter
          deno = (-2.0d0/sqrtpi)*exp(-(root(i))**2)
         root(i+1) = root(i) * erfunc(root(i))/deno
          if(abs(root(i+1)-root(i)).lt.1.0d-2) then
          x_out = root(i+1)
          return
          endif
       enddo
       x_out = root(maxiter+1)
       return
       end subroutine newton   
        
        
           
      
 
