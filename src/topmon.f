      program topmon

      use global_data,only:ierr,myid,numprocs
      use var_type
      use util_files
      use util_timings
      implicit none
      include 'common.inc'
      integer::narg,iarg
      character(LEN=default_path_length)::sarg,infile='topmon.inp'
      logical::lrun=.true.,lsetinput=.false.,lversion=.false.,lusage =.false.
! ----------------------------------------------------------------
! Parse the command line arguments
      narg=command_argument_count()
      print*,narg
      iarg=1
      do while (iarg.le.narg)
         call get_command_argument(iarg,sarg)
         print *,iarg,trim(sarg)
         select case(sarg)
         case('--version','-v')
            lversion=.true.
            lrun=.false.
         case('--help','-h')
            lusage=.true.
            lrun=.false.
         case('--input','-i')
            iarg=iarg+1
            if (iarg.gt.narg) then
               lusage=.true.
               lrun=.false.
               exit
            end if
            lsetinput=.true.
            infile=sarg
         case default
            if (.not.lsetinput.and.iarg.eq.narg) then
               lsetinput=.true.
               infile=sarg
            end if
         end select
         iarg=iarg+1
      end do
! ----------------------------------------------------------------
      if (lversion) write(*,'(A,/,A,/)') 'MCCCS topmon 2011-01-01' ,'branch: zeolite'

      if (lusage) then
         call get_command_argument(0,sarg)
         write(*,'(A,/,T4,A,/)') 'Usage: '//trim(sarg)// ' [--version|-v] [--help|-h] [[--input|-i] /path/to/input/file]' ,'If the input file is the last argument, --input or -i can  be omitted'
      end if

      if (lrun) then
         call MPI_INIT( ierr )
         call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
         call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

! --- call main program
         call monola(infile)

         call MPI_FINALIZE(ierr)
      end if
! ----------------------------------------------------------------
      end program topmon
