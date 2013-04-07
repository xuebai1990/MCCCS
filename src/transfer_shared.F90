module transfer_shared
  use util_runtime,only:err_exit
  use sim_system
  implicit none
  private
  save
  public::read_transfer,update_bias,opt_bias,lopt_bias,freq_opt_bias,read_checkpoint_transfer_shared,write_checkpoint_transfer_shared

  real,allocatable::u_bias_diff(:,:) !<u_bias_diff(i,j) is the bias potential difference that should be applied to box i to achieve equal distribution of particles of type j
  integer,allocatable::num_update_bias(:,:)
  logical,allocatable::lopt_bias(:)
  integer::freq_opt_bias=500
  namelist /transfer/ lopt_bias,freq_opt_bias
contains
  subroutine read_transfer(file_in,lprint)
    character(LEN=*),INTENT(IN)::file_in
    LOGICAL,INTENT(IN)::lprint
    integer::io_input,jerr

    allocate(u_bias_diff(nbox,nmolty),num_update_bias(nbox,nmolty),lopt_bias(nmolty),stat=jerr)
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'read_transfer: memory allocation',jerr)
    end if

    u_bias_diff=0.0_dp
    num_update_bias=0
    lopt_bias=.false.

    io_input=get_iounit()
    open(unit=io_input,access='sequential',action='read',file=file_in,form='formatted',iostat=jerr,status='old')
    if (jerr.ne.0) then
       call err_exit(__FILE__,__LINE__,'cannot open transfer input file',myid+1)
    end if

    read(UNIT=io_input,NML=transfer,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) then
       call err_exit(__FILE__,__LINE__,'reading namelist: transfer',jerr)
    end if
    close(io_input)

    if (lprint) then
       write(io_output,*) 'lopt_bias: ',lopt_bias
       write(io_output,*) 'freq_opt_bias: ',freq_opt_bias
    end if

  end subroutine read_transfer

  subroutine update_bias(u_diff,boxrem,boxins,imolty)
    use util_math,only:update_average
    real,intent(in)::u_diff
    integer,intent(in)::boxrem,boxins,imolty

    num_update_bias(boxins,imolty)=num_update_bias(boxins,imolty)+1
    num_update_bias(boxrem,imolty)=num_update_bias(boxrem,imolty)+1
    call update_average(u_bias_diff(boxins,imolty),u_diff/2.0_dp,num_update_bias(boxins,imolty))
    call update_average(u_bias_diff(boxrem,imolty),-u_diff/2.0_dp,num_update_bias(boxrem,imolty))
  end subroutine update_bias

  subroutine opt_bias
    integer::imolty,ibox

    do imolty=1,nmolty
       if (.not.lopt_bias(imolty)) cycle
       do ibox=1,nbox
          eta2(ibox,imolty)=eta2(ibox,imolty)+u_bias_diff(ibox,imolty)
       end do
    end do
    u_bias_diff=0.0_dp
    num_update_bias=0
  end subroutine opt_bias

  subroutine read_checkpoint_transfer_shared(io_chkpt)
    use util_mp,only:mp_bcast
    integer,intent(in)::io_chkpt
    if (myid.eq.rootid) read(io_chkpt) num_update_bias,u_bias_diff
    call mp_bcast(num_update_bias,nbox*nmolty,rootid,groupid)
    call mp_bcast(u_bias_diff,nbox*nmolty,rootid,groupid)
  end subroutine read_checkpoint_transfer_shared

  subroutine write_checkpoint_transfer_shared(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) num_update_bias,u_bias_diff
  end subroutine write_checkpoint_transfer_shared
end module transfer_shared
