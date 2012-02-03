ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Ordering the .res file to input coord for fort.77
c     Readin files: fort.33,source_filename(.res)
c     Format of fort.33, which is from source_filename.inp:
c     source_filename
c     nbead
c     C1 H2 C5 O4 ...
c     Output file: source_filename.frac
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program getfrac

      implicit none

      character (len=30) infile,outfile
      character (len=5) sname(500),celltype
      integer nbead, i, j
      real cell(6),bead(3,500),dum
      
      read(4,*) infile
      read(4,*) nbead

      outfile=trim(infile) // '.frac'

      OPEN(22,FILE=infile,STATUS='old')
      read(22,*)
      read(22,*) celltype,dum,cell
      read(22,*)
      read(22,*)
      do i=1,nbead
         read(22,*) sname(i),dum,bead(1:3,i)
      enddo
      close(22)
      
      print*, nbead
      OPEN(23,FILE=outfile,STATUS='unknown')
      write(23,*) cell
      write(23,*) nbead
      do i=1,nbead
         write(23,*) sname(i), bead(1:3,i)
      enddo
      close(23)
      print*, 'coordinates are in : ', outfile
      end program
