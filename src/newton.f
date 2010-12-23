	subroutine newton(x_in,x_out)

!  --  This subroutine calcultes the roots of non-linear equation using 
!  --  Newton-Raphson method.

       use global_data
       use var_type
       use const_phys
       use const_math
       use util_math
       use util_string
       use util_files
       use util_timings
       implicit none
  
       integer(KIND=int)::i,j,maxiter
       real(KIND=double_precision)::x_in,x_out,erfunc
       real(KIND=double_precision)::root(11),sqrtpi,deno

       parameter (maxiter = 10)
 
       root(1) = x_in
 
       sqrtpi = sqrt(4.0d0*datan(1.0d0)) 

       do i = 1,maxiter
          deno = (-2.0d0/sqrtpi)*exp(-(root(i))**2)
         root(i+1) = root(i) * erfunc(root(i))/deno
          if(abs(root(i+1)-root(i)).lt.1.0d-2) then
          x_out = root(i+1)
          return
          end if
       end do
       x_out = root(maxiter+1)
       return
       end subroutine newton   
        
        
           
      
 
