      subroutine updnn( i, ibox )

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'system.inc'
!$$$      include 'neigh.inc'
!$$$      include 'neigh2.inc'
      integer(KIND=int)::i,j,ii,jj,ibox
      real(KIND=double_precision)::rxuij,ryuij,rzuij,rijsq,rcnnsq
! -----------------------------------------------------------------------------
 
      rcnnsq = rcutnn(ibox)**2

      if ( lpbc ) call setpbc (ibox)

! -----------------------------------------------------------------
 
! --- set i-part of logical::map to .false. ---

      do 9 j = 1, nchain
         lnn(i,j) = .false.
 9    continue

! -----------------------------------------------------------------
 
! --- loop over all beads ii of chain i 
      do 99 ii = 1, nunit(i)

         rxui = rxu(i,ii)
         ryui = ryu(i,ii)
         rzui = rzu(i,ii)

! --- loop over all chains j
         do 98 j = 1, nchain
               
            if ( i .eq. j ) go to 98

            if ( lnn(i,j) ) go to 98

! --- loop over all beads jj of chain j 
            do 97 jj = 1, nunit(j)

               rxuij = rxui - rxu(j,jj)
               ryuij = ryui - ryu(j,jj)
               rzuij = rzui - rzu(j,jj)

! *** minimum image the pair separations ***
               if ( lpbc ) call mimage ( rxuij,ryuij,rzuij,ibox )

               rijsq = rxuij*rxuij + ryuij*ryuij + rzuij*rzuij

               if ( rijsq .le. rcnnsq ) then
                  lnn(i,j) = .true.
                  go to 98
               end if
 
 97         continue
 98      continue
 99   continue
 
! -----------------------------------------------------------------------------
 
! *** set displacementvectors to zero ***
 
      do 1099 j = 1, 3
         disvec(1,i,j) = 0.0d0
         disvec(2,i,j) = 0.0d0
 1099 continue
 
!      write(iou,*) '@@@ control updnn @@@'
!      write(iou,*) 'lnn', i, '   ', (lnn(i,j),j=1,nchain)

      return
      end



