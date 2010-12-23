      subroutine writepdb(nmol,natom,ibox)

      use global_data
      use var_type
      use const_phys
      use const_math
      use util_math
      use util_string
      use util_files
      use util_timings
      implicit none
      include 'common.inc'

!$$$      include 'zeolite.inc'
!$$$      include 'control.inc'
!$$$      include 'coord.inc'

      integer(KIND=normal_int)::im,ia,nmol,natom,iguest,imm,ibox
      character(LEN=4)::atom

      open(unit=48, file='system.pdb', form='formatted')
      rewind(48)
      write(iou,*) 'most likely does not work as before 9-3-96'
      write(iou,100)

      atom = 'C   '
      iguest=0 
      do im=1,nmol
         do ia=1,natom
            imm=parbox(im,ibox,1)
            iguest=iguest+1
            write(48,99) "ATOM  ",iguest,atom,
     &                 rxu(imm,ia),ryu(imm,ia),rzu(imm,ia)
         end do
      end do

   99 format(a6,i5,1x,a4,14x,3f8.3)
  100 format(/,' WRITING GUEST ATOMS TO FILE system.pdb')

      return
      end

