module transfer_shared
  use util_runtime,only:err_exit
  use sim_system
  implicit none
  private
  save
  public::read_transfer,update_bias,opt_bias,lopt_bias,freq_opt_bias

  real,allocatable::u_bias_diff(:,:) !<u_bias_diff(i,j) is the bias potential difference that should be applied to box i to achieve equal distribution of particles of type j
  integer,allocatable::num_update_bias(:,:)
  logical,allocatable::lopt_bias(:)
  integer::freq_opt_bias=500
  namelist /transfer/ lopt_bias,freq_opt_bias
contains
  subroutine read_transfer(io_input)
    integer,intent(in)::io_input
    integer::jerr

    allocate(u_bias_diff(nbox,nmolty),num_update_bias(nbox,nmolty),lopt_bias(nmolty),stat=jerr)
    if (jerr.ne.0) then
       write(io_output,*) 'ERROR ',jerr,' in ',TRIM(__FILE__),':',__LINE__
       call err_exit('read_transfer: memory allocation')
    end if

    u_bias_diff=0.0_double_precision
    num_update_bias=0
    lopt_bias=.true.

    rewind(io_input)
    read(UNIT=io_input,NML=transfer,iostat=jerr)
    if (jerr.ne.0.and.jerr.ne.-1) then
       write(io_output,*) 'ERROR ',jerr,' in ',TRIM(__FILE__),':',__LINE__
       call err_exit('reading namelist: transfer')
    end if

    write(io_output,*) 'lopt_bias: ',lopt_bias
    write(io_output,*) 'freq_opt_bias: ',freq_opt_bias

  end subroutine read_transfer

  subroutine update_bias(u_diff,boxrem,boxins,imolty)
    real,intent(in)::u_diff
    integer,intent(in)::boxrem,boxins,imolty

    u_bias_diff(boxins,imolty)=(u_bias_diff(boxins,imolty)*num_update_bias(boxins,imolty)+u_diff/2.0)/real(num_update_bias(boxins,imolty)+1.0)
    u_bias_diff(boxrem,imolty)=(u_bias_diff(boxrem,imolty)*num_update_bias(boxrem,imolty)-u_diff/2.0)/real(num_update_bias(boxrem,imolty)+1.0)
    num_update_bias(boxins,imolty)=num_update_bias(boxins,imolty)+1
    num_update_bias(boxrem,imolty)=num_update_bias(boxrem,imolty)+1
  end subroutine update_bias

  subroutine opt_bias
    integer::imolty,ibox

    do imolty=1,nmolty
       if (.not.lopt_bias(imolty)) cycle
       do ibox=1,nbox
          eta2(ibox,imolty)=eta2(ibox,imolty)+u_bias_diff(ibox,imolty)
       end do
    end do
    u_bias_diff=0.0_double_precision
    num_update_bias=0
  end subroutine opt_bias

end module transfer_shared
