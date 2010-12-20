      subroutine writepdb(nmol,natom,ibox)

      implicit none

      include 'zeolite.inc'
      include 'control.inc'
      include 'coord.inc'

      integer::im,ia,nmol,natom,iguest,imm,ibox
      character::atom*4

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

