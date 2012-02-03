      subroutine rwztab(iswitc)

c rwztab
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

c   read/write zeolite table from/to disk
      implicit none
      integer iswitc
      include 'control.inc'
      include 'grid.inc'
      include 'zeolite.inc'
      if (iswitc.eq.0) then
c ---    read zeolite table from disk
         read(91) ngrx,ngry,ngrz,dgrx,dgry,dgrz,zunitx,zunity,zunitz
     +         ,zunitxi,zunityi,zunitzi,factx,facty,factz
     +         ,xzz,yzz,zzz
     +         ,egrid 
         write(iou,1000) dgrx,ngrx,dgry,ngry,dgrz,ngrz,zunitx,zunity
     +                ,zunitz
         rewind(91)         
      else
        print*,' write table zeolite to 91 '
        write(91) ngrx,ngry,ngrz,dgrx,dgry,dgrz,zunitx,zunity,zunitz
     +         ,zunitxi,zunityi,zunitzi,factx,facty,factz
     +         ,xzz,yzz,zzz
     +         ,egrid  
      endif
      return
 1000 format(' Zeolite table from disk: ',/,
     +       '   x-dir: grid size : ',f7.4,' nr. point ',i6,/,
     +       '   y-dir: grid size : ',f7.4,' nr. point ',i6,/,
     +       '   z-dir: grid size : ',f7.4,' nr. point ',i6,/,/,
     +       ' Size unit-cell zeolite: ',f7.4,' x ',f7.4,' x ',f7.4,/)
      end

   

