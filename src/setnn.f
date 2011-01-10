!    *******************************************************************
!    ** generates the complete near neighbour bitmap.                 **
!    **                                                               **
!    ** principal variables:                                          **
!    **                                                               **
!    ** logic   lnn               the bitmap                          **
!    ** dble    disvec            the displacementvector              **
!    **                                                               **
!    ** usage:                                                        **
!    **                                                               **
!    ** the subroutine generates the complete near neighbour bitmap   **
!    ** and sets all displacementvectors to zero at the beginning of  **
!    ** the run.                                                      **
!    *******************************************************************
      subroutine setnn ( ibox )
      use sim_system
      use sim_cell,only:setpbc,mimage
      implicit none

!$$$      include 'control.inc'
!$$$      include 'coord.inc'
!$$$      include 'neigh2.inc'
!$$$      include 'neigh.inc'
!$$$      include 'system.inc'
      integer::i,j,ii,jj,ibox
      real::rxui,ryui,rzui,rxuij,ryuij,rzuij ,rijsq,rcnnsq
! -----------------------------------------------------------------------------
 
      rcnnsq = rcutnn(ibox)**2

      if ( lpbc ) call setpbc (ibox)

! --------------------------------------------------------------------
 
! --- set logical::map to .false. ---
      do 10 i = 1, nchain
         do 9 j = 1, nchain
            lnn(i,j) = .false.
 9       continue
 10   continue

! --------------------------------------------------------------------
 
! --- loop over all chains i 
      do 100 i = 1, nchain
 
         if ( nboxi(i) .ne. ibox ) go to 100

! --- loop over all beads ii of chain i 
         do 99 ii = 1, nunit(i)

            rxui = rxu(i,ii)
            ryui = ryu(i,ii)
            rzui = rzu(i,ii)

! --- loop over all chains j
            do 98 j = 1, nchain
               
               if ( ( i .eq. j ) .or.  ( nboxi(i) .ne. nboxi(j) ) ) go to 98

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
 
 97            continue
 98         continue
 99      continue
 100  continue
 
! -----------------------------------------------------------------------------
 
! *** set displacementvectors to zero ***
 
      do 1100 i = 1, nchain
 
         do 1099 j = 1, 3
            disvec(1,i,j) = 0.0d0
            disvec(2,i,j) = 0.0d0
1099     continue
1100  continue
 
!      write(io_output,*) '@@@ control setnn @@@'
!      do 2000 i = 1, nchain
!         write(io_output,*) 'lnn', i, '   ', (lnn(i,j),j=1,nchain)
!2000  continue
 
      write(io_output,*)
      write(io_output,*) 'neighbour table created by SETNN'

      return
      end

