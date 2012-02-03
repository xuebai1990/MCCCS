      subroutine linkcell(iinit,imol,xcmi,ycmi,zcmi,cellinc)

c linkcell
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Copyright (C) 1999-2004 Bin Chen, Marcus Martin, Jeff Potoff, 
c John Stubbs, and Collin Wick and Ilja Siepmann  
c                     
c This program is free software; you can redistribute it and/or
c modify it under the terms of the GNU General Public License
c as published by the Free Software Foundation; either version 2
c of the License, or (at your option) any later version.
c
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to 
c
c Free Software Foundation, Inc. 
c 59 Temple Place - Suite 330
c Boston, MA  02111-1307, USA.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none


      include 'control.inc'
      include 'coord.inc'
      include 'system.inc'      

      integer i,j,k,n,ncellx,ncelly,ncellz,iinit,ibox,linkdecode
     &     ,imolty,imol,ic,cellinc,ia,ja,ka,ib,jb,kb,count,ncell,ncello

      double precision dcellx,dcelly,dcellz,rx,ry,rz,xcmi,ycmi,zcmi

c     *** rintramax is the maximum distance between endpoints 
c     *** in a molecule

      dimension cellinc(27)

      save dcellx,dcelly,dcellz,ncellx,ncelly,ncellz,ncello
      
c      write(2,*) 'START LINKCELL IINIT=',iinit

       

      if (iinit.eq.1) then
c     --- we will set up or update our cell sizes

         ibox = boxlink

         if (imol .eq. 0) then
c * called from monola
            ncello = 0
         endif

c     --- determine dcell
c         write(2,*) 'linkcell used',rcut,rintramax

         dcellx = rcut(ibox) + rintramax
         dcelly = dcellx
         dcellz = dcellx

c     --- find hypothetical ncell
         ncellx = int( boxlx(ibox) / dcellx ) 
         ncelly = int( boxly(ibox) / dcelly ) 
         ncellz = int( boxlz(ibox) / dcellz )

         
c     --- make dcells larger so each each cell is the same size
         dcellx = boxlx(ibox) / dble(ncellx)
         dcelly = boxly(ibox) / dble(ncelly)
         dcellz = boxlz(ibox) / dble(ncellz)

c     --- now reweight ncell one more time

         ncellx = anint( boxlx(ibox) / dcellx ) 
         ncelly = anint( boxly(ibox) / dcelly ) 
         ncellz = anint( boxlz(ibox) / dcellz ) 

         ncell = ncellx * ncelly * ncellz
 
         if (ncell .ne. ncello) then
            if (imol .eq. 0) then
               write(2,*) 'number of linkcells set to',ncell
            else
               write(2,*) 'number of linkcells changed to',ncell
            endif
         endif

         ncello = ncell

         if (ncell.gt.cmax) then
            write(2,*) 'ncell,cmax',ncell,cmax
            stop 'ncell greater than cmax in linkcell'
         endif
         
         do n = 1, ncell
            nicell(n) = 0
         enddo

c     --- assign molecules to cells
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
                  write(2,*) 'ic,cmax',ic,cmax
                  stop 'ic gt cmax'
               endif

               icell(n) = ic
               nicell(ic) = nicell(ic) + 1

               if (nicell(ic).gt.cmaxa) then
                  write(2,*) 'nicell,cmaxa',nicell(ic)
     &                 ,cmaxa
                  stop 'nicell gt cmaxa'
               endif

               iucell(ic,nicell(ic)) = n
            else
               icell(n) = 0
            endif
         enddo

      elseif (iinit.eq.2) then
c     --- we will update our cell's occupants

         ic = icell(imol)

         if (ic.gt.0) then
c     --- first we will remove our molecule
            do n = 1, nicell(ic) 
               if (iucell(ic,n).eq.imol) then
c     --- replace removed occupant with last occupant and erase last spot
                  iucell(ic,n) = iucell(ic,nicell(ic))
                  iucell(ic,nicell(ic)) = 0
                  nicell(ic) = nicell(ic) - 1
                  goto 100
               endif               
            enddo
         
            stop 'screwup for iinit = 2 for linkcell'
 100        continue
            icell(imol) = 0
         endif

c     --- now we will add the molecule 
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
         
            if (nicell(ic).gt.cmaxa) stop 'nicell too big'

            iucell(ic,nicell(ic)) = imol

         endif
      else
c     --- we will determine the cell neighbors
         
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
                  elseif (ia.lt.1) then
                     ib = ia + ncellx
                  else
                     ib = ia
                  endif

                  if (ja.gt.ncelly) then
                     jb = ja - ncelly
                  elseif (ja.lt.1) then
                     jb = ja + ncelly
                  else
                     jb = ja
                  endif

                  if (ka.gt.ncellz) then
                     kb = ka - ncellz
                  elseif (ka.lt.1) then
                     kb = ka + ncellz
                  else
                     kb = ka
                  endif
                  
                  count = count + 1

                  ic = linkdecode(ib,jb,kb,ncellx,ncelly,ncellz)

                  cellinc(count) = ic          
                          
               enddo
            enddo
         enddo
      endif

c      write(2,*) 'END LINKCELL IINIT=',iinit

      return
      end





      function linkdecode(i,j,k,ncellx,ncelly,ncellz)
            
      implicit none
      
      integer i,j,k,ncellx,ncelly,ncellz,linkdecode
      
c     *** decodes x,y,z to a single number

      linkdecode = (i-1)*ncelly*ncellz + (j-1)*ncellz + k

      return
      
      end
      







