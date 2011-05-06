      program topmon

      use global_data,only:ierr,myid,numprocs
      use global_data,only:ierr,myid,numprocs,thread_id,thread_num,thread_num_max,thread_num_proc
      use var_type,only:normal_int,default_path_length
      use util_timings,only:time_init
      implicit none
      include 'common.inc'
      integer(KIND=normal_int)::narg,iarg,jerr
      character(LEN=default_path_length)::sarg,file_in='topmon.inp'
      logical::lrun=.true.,lsetinput=.false.,lversion=.false.,lusage =.false.
! ----------------------------------------------------------------
! Initialize the timer
      call time_init()

! Parse the command line arguments
      narg=command_argument_count()
      iarg=1
      do while (iarg.le.narg)
         call get_command_argument(iarg,sarg)
         select case(sarg)
         case('--version','-v')
            lversion=.true.
            lrun=.false.
         case('--help','-h')
            lusage=.true.
            lrun=.false.
         case('--thread','-t')
            iarg=iarg+1
            if (iarg.gt.narg) then
               lusage=.true.
               lrun=.false.
               exit
            end if
            call get_command_argument(iarg,sarg)
            read(sarg,*,iostat=jerr) thread_num
            if (jerr.ne.0) then
               lusage=.true.
               lrun=.false.
               exit
            end if
!$          thread_num_max=omp_get_max_threads()
!$          if (thread_num.gt.thread_num_max) thread_num=thread_num_max
!$          call omp_set_num_threads(thread_num)
         case('--input','-i')
            iarg=iarg+1
            if (iarg.gt.narg) then
               lusage=.true.
               lrun=.false.
               exit
            end if
            call get_command_argument(iarg,sarg)
            lsetinput=.true.
            file_in=sarg
         case default
            if (.not.lsetinput.and.iarg.eq.narg) then
               lsetinput=.true.
               file_in=sarg
            end if
         end select
         iarg=iarg+1
      end do
! ----------------------------------------------------------------
      if (lversion) write(*,'(A,/,A,/)') 'MCCCS topmon 2011-01-01','branch: zeolite'

      if (lusage) then
         call get_command_argument(0,sarg)
         write(*,'(A,/,T4,A,/,T4,A,/)') 'Usage: '//trim(sarg)// ' [--version|-v] [--help|-h] [(--threads|-t) number_of_threads_per_processor]','[(--input|-i) /path/to/input/file]','If the input file is the last argument, --input or -i can  be omitted'
      end if

      if (lrun) then
         call MPI_INIT( ierr )
         call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
         call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

! --- call main program
         call monola(file_in)

         call MPI_FINALIZE(ierr)
      end if
! ----------------------------------------------------------------
      end program topmon
