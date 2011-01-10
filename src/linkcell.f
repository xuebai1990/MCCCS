      subroutine linkcell(iinit,imol,xcmi,ycmi,zcmi,cellinc)

      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_runtime,only:err_exit
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'system.inc'      

      integer(KIND=normal_int)::i,j,k,n,ncellx,ncelly,ncellz,iinit,ibox ,linkdecode,imol,ic,cellinc(27),ia,ja,ka,ib,jb,kb,count,ncell ,ncello

      real(KIND=double_precision)::dcellx,dcelly,dcellz,rx,ry,rz,xcmi ,ycmi,zcmi

!     *** rintramax is the maximum distance between endpoints 
!     *** in a molecule

      save dcellx,dcelly,dcellz,ncellx,ncelly,ncellz,ncello
      
!      write(io_output,*) 'START LINKCELL IINIT=',iinit

       

      if (iinit.eq.1) then
!     --- we will set up or update our cell sizes

         ibox = boxlink

         if (imol .eq. 0) then
! * called from monola
            ncello = 0
         end if

!     --- determine dcell
!         write(io_output,*) 'linkcell used',rcut,rintramax

         dcellx = rcut(ibox) + rintramax
         dcelly = dcellx
         dcellz = dcellx

!     --- find hypothetical ncell
         ncellx = int( boxlx(ibox) / dcellx ) 
         ncelly = int( boxly(ibox) / dcelly ) 
         ncellz = int( boxlz(ibox) / dcellz )

         
!     --- make dcells larger so each each cell is the same size
         dcellx = boxlx(ibox) / dble(ncellx)
         dcelly = boxly(ibox) / dble(ncelly)
         dcellz = boxlz(ibox) / dble(ncellz)

!     --- now reweight ncell one more time

         ncellx = anint( boxlx(ibox) / dcellx ) 
         ncelly = anint( boxly(ibox) / dcelly ) 
         ncellz = anint( boxlz(ibox) / dcellz ) 

         ncell = ncellx * ncelly * ncellz
 
         if (ncell .ne. ncello) then
            if (imol .eq. 0) then
               write(io_output,*) 'number of linkcells set to',ncell
            else
               write(io_output,*) 'number of linkcells changed to',ncell
            end if
         end if

         ncello = ncell

         if (ncell.gt.cmax) then
            write(io_output,*) 'ncell,cmax',ncell,cmax
            call err_exit('ncell greater than cmax in linkcell')
         end if
         
         do n = 1, ncell
            nicell(n) = 0
         end do

!     --- assign molecules to cells
         do n = 1, nchain

            if (nboxi(n).eq.boxlink) then
               rx = xcm(n) / dcellx
               ry = ycm(n) / dcelly
               rz = zcm(n) / dcellz

               i = int(rx) + 1
               j = int(ry) + 1
               k = int(rz) + 1
            
               ic = linkdecode(i,j,k,ncellx,ncelly,ncellz)


               if (ic.gt.cmax) then
                  write(io_output,*) 'ic,cmax',ic,cmax
                  call err_exit('ic gt cmax')
               end if

               icell(n) = ic
               nicell(ic) = nicell(ic) + 1

               if (nicell(ic).gt.cmaxa) then
                  write(io_output,*) 'nicell,cmaxa',nicell(ic) ,cmaxa
                  call err_exit('nicell gt cmaxa')
               end if

               iucell(ic,nicell(ic)) = n
            else
               icell(n) = 0
            end if
         end do

      else if (iinit.eq.2) then
!     --- we will update our cell's occupants

         ic = icell(imol)

         if (ic.gt.0) then
!     --- first we will remove our molecule
            do n = 1, nicell(ic) 
               if (iucell(ic,n).eq.imol) then
!     --- replace removed occupant with last occupant and erase last spot
                  iucell(ic,n) = iucell(ic,nicell(ic))
                  iucell(ic,nicell(ic)) = 0
                  nicell(ic) = nicell(ic) - 1
                  goto 100
               end if               
            end do
         
            call err_exit('screwup for iinit = 2 for linkcell')
 100        continue
            icell(imol) = 0
         end if

!     --- now we will add the molecule 
         if (nboxi(imol).eq.boxlink) then
            rx = xcm(imol) / dcellx
            ry = ycm(imol) / dcelly
            rz = zcm(imol) / dcellz

            i = int(rx) + 1
            j = int(ry) + 1
            k = int(rz) + 1

            ic = linkdecode(i,j,k,ncellx,ncelly,ncellz)
            icell(imol) = ic

            nicell(ic) = nicell(ic) + 1
         
            if (nicell(ic).gt.cmaxa) call err_exit('nicell too big')

            iucell(ic,nicell(ic)) = imol

         end if
      else
!     --- we will determine the cell neighbors
         
         rx = xcmi / dcellx
         ry = ycmi / dcelly
         rz = zcmi / dcellz

         i = int(rx) + 1
         j = int(ry) + 1
         k = int(rz) + 1
   
         count = 0
         
         do ia = i-1, i+1
            do ja = j-1, j+1
               do ka = k-1, k+1
                  if (ia.gt.ncellx) then
                     ib = ia - ncellx
                  else if (ia.lt.1) then
                     ib = ia + ncellx
                  else
                     ib = ia
                  end if

                  if (ja.gt.ncelly) then
                     jb = ja - ncelly
                  else if (ja.lt.1) then
                     jb = ja + ncelly
                  else
                     jb = ja
                  end if

                  if (ka.gt.ncellz) then
                     kb = ka - ncellz
                  else if (ka.lt.1) then
                     kb = ka + ncellz
                  else
                     kb = ka
                  end if
                  
                  count = count + 1

                  ic = linkdecode(ib,jb,kb,ncellx,ncelly,ncellz)

                  cellinc(count) = ic          
                          
               end do
            end do
         end do
      end if

!      write(io_output,*) 'END LINKCELL IINIT=',iinit

      return
      end





      function linkdecode(i,j,k,ncellx,ncelly,ncellz)
            
      use sim_system
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
      
      integer(KIND=normal_int)::i,j,k,ncellx,ncelly,ncellz,linkdecode
      
!     *** decodes x,y,z to a single number

      linkdecode = (i-1)*ncelly*ncellz + (j-1)*ncellz + k

      return
      
      end
      







