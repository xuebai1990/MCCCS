      subroutine dump  

! -----------------------------------------------------------------
! subroutine dump
! dumps the final configuration before stopping the program (Neeraj).
! -----------------------------------------------------------------

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

!$$$      include 'control.inc'
!$$$      include 'coord.inc'      
!$$$      include 'poten.inc'
!$$$      include 'system.inc'
!$$$      include 'ensemble.inc'
!$$$      include 'cbmc.inc'
!$$$      include 'conver.inc'
!$$$      include 'inpar.inc' 
!$$$      include 'external.inc'
!$$$      include 'zeolite.inc'
!$$$      include 'inputdata.inc'
!$$$      include 'blkavg.inc'
!$$$      include 'bnbsma.inc'
!$$$      include 'swtcmove.inc'
!$$$      include 'ewaldsum.inc'
!$$$      include 'neigh.inc'
!$$$      include 'clusterbias.inc'
!$$$      include 'cell.inc'

      integer(KIND=normal_int)::i,j,im,imolty,ibox

      open (8, file="final-config")
      write(8,*) tmcc
      if ( tmcc .gt. 0 ) then
         write(8,*) Armtrax, Armtray, Armtraz 
         do im=1,nbox
            do imolty=1,nmolty
               write(8,*) rmtrax(imolty,im), rmtray(imolty,im) , rmtraz(imolty,im)
               write(8,*) rmrotx(imolty,im), rmroty(imolty,im) , rmrotz(imolty,im)
            end do
         end do
         do im=1, nbox
            write (8,*) (rmflcq(i,im),i=1,nmolty)
         end do
         if ( lgibbs .or. lgrand .or. lnpt ) then
            write(8,*) (rmvol(ibox),ibox=1,nbox)
            do ibox = 1,nbox
               if (lsolid(ibox) .and. .not. lrect(ibox)) then
                  write(8,*) (rmhmat(ibox,i),i=1,9)
                  write(8,*) (hmat(ibox,i),i=1,9)
               else
                  write(8,*) boxlx(ibox),boxly(ibox),boxlz(ibox)
               end if
            end do
         end if
      end if
      write(8,*) nchain
      write(8,*) nmolty
      write(8,*) (nunit(i),i=1,nmolty)
      write(8,*) (moltyp(i),i=1,nchain)
      write(8,*) (nboxi(i),i=1,nchain)
      do i = 1, nmolty
         if ( lexpand(i) ) write(8,*) eetype(i)
      end do
      do i = 1, nmolty
         if ( lexpand(i) ) write(8,*) rmexpc(i)
      end do
      do  i = 1, nchain
         imolty = moltyp(i)
         do j = 1, nunit(imolty)
            write(8,*) rxu(i,j), ryu(i,j), rzu(i,j), qqu(i,j)
         end do
      end do
      close (8)
      return
      end
