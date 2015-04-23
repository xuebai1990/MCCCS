!> \file Calculate properties using the Widom test particle insertion method.
!>
!> \warning Widom insertion is known to fail at high densities.

module prop_widom
  implicit none
  private
  save
  public::read_prop_widom,calc_prop_widom,blk_avg_prop_widom,write_prop_widom,write_prop_widom_with_stats,read_checkpoint_prop_widom,write_checkpoint_prop_widom

  real,allocatable::blk_avg_setedist(:,:,:),setedist(:,:),blk_avg_setedist_prev(:,:),blk_avg_Uads(:,:),Uads(:),blk_avg_Uads_prev(:),Wrosen(:),blk_avg_Wrosen_prev(:)
contains
  subroutine read_prop_widom(io_input,lprint,blockm)
    use var_type,only:dp
    use util_runtime,only:err_exit
    use sim_system,only:lgrand,ntmax
    INTEGER,INTENT(IN)::io_input,blockm
    LOGICAL,INTENT(IN)::lprint
    integer::jerr

    lgrand=.true.

    if (allocated(blk_avg_setedist)) deallocate(Wrosen,blk_avg_Wrosen_prev,Uads,blk_avg_Uads_prev,blk_avg_Uads,setedist,blk_avg_setedist_prev,blk_avg_setedist,stat=jerr)
    allocate(Wrosen(ntmax),blk_avg_Wrosen_prev(ntmax),Uads(ntmax),blk_avg_Uads_prev(ntmax),blk_avg_Uads(ntmax,blockm),setedist(4,ntmax),blk_avg_setedist_prev(4,ntmax),blk_avg_setedist(4,ntmax,blockm),stat=jerr)
    if (jerr.ne.0) call err_exit(__FILE__,__LINE__,'read_prop_widom: memory allocation',jerr)

    setedist=0._dp
    Wrosen=0._dp
    Uads=0._dp
  end subroutine read_prop_widom

  subroutine calc_prop_widom(imolty,weight)
    use var_type,only:dp
    use sim_system,only:vnew,rxnew,rynew,rznew,nugrow,ivTot
    integer,intent(in)::imolty
    real,intent(in)::weight
    real::s,r1,r2,arg
    integer::ip

    Wrosen(imolty) = Wrosen(imolty) + weight
    Uads(imolty) = Uads(imolty) + vnew(ivTot)*weight
    s=0.0_dp
    do ip=1,3
       if (ip.eq.1) then
          r1=rxnew(1)
          r2=rxnew(nugrow(imolty))
       else if (ip.eq.2) then
          r1=rynew(1)
          r2=rynew(nugrow(imolty))
       else if (ip.eq.3) then
          r1=rznew(1)
          r2=rznew(nugrow(imolty))
       end if
       arg = (r1-r2)**2
       s = s + arg
       setedist(ip,imolty) = setedist(ip,imolty) + arg*weight
    end do
    setedist(0,imolty) = setedist(0,imolty) + s*weight

  end subroutine calc_prop_widom

  subroutine blk_avg_prop_widom(nblock)
    use util_math,only:store_block_average
    use sim_system,only:nmolty
    integer,intent(in)::nblock
    real::Wrosen_tmp
    integer::itype,j

    do itype = 1, nmolty
       do j = 0, 3
          Wrosen_tmp = blk_avg_Wrosen_prev(itype)
          call store_block_average(blk_avg_setedist(j,itype,nblock),setedist(j,itype)/Wrosen(itype),Wrosen(itype),blk_avg_setedist_prev(j,itype),Wrosen_tmp)
       end do
       call store_block_average(blk_avg_Uads(itype,nblock),Uads(itype)/Wrosen(itype),Wrosen(itype),blk_avg_Uads_prev(itype),blk_avg_Wrosen_prev(itype))
    end do
  end subroutine blk_avg_prop_widom

  subroutine write_prop_widom(io_output)
    use sim_system,only:nmolty
    integer,intent(in)::io_output
    integer::i,j

    Uads(1:nmolty)=Uads(1:nmolty)/Wrosen(1:nmolty)
    do i=0,3
       setedist(i,1:nmolty)=setedist(i,1:nmolty)/Wrosen(1:nmolty)
    end do

    write(io_output,"(/,A)") '  ** Properties from Rosenbluth Sampling **'
    do i=1,nmolty
       write(io_output,"(A,I2,A,1X,F12.3)") ' potential energy of type ',i,'          [K] =',Uads(i)
    end do
    do i=1,nmolty
       write(io_output,"(A,I2,A,4(1X,F12.3))") ' sete & x,y,z len of type ',i,'        [A^2] =',(setedist(j,i),j=0,3)
    end do
  end subroutine write_prop_widom

  subroutine write_prop_widom_with_stats(io_output,nblock)
    use util_math,only:calculate_statistics
    use sim_system,only:nmolty
    integer,intent(in)::io_output,nblock
    real::avg(4),sd(4),se(4)
    integer::itype,j

    write(io_output,"(/,A)") '  ** Properties with Uncertainties from Rosenbluth Sampling **'
    do itype = 1, nmolty
       call calculate_statistics(blk_avg_Uads(itype,1:nblock),avg(0),sd(0),se(0))
       write(io_output,"(2(A,I2),A,3(1X,F12.3))") ' potential energy    itype ',itype,' box ',1,' = ',avg(0),sd(0),se(0)
    end do
    do itype = 1, nmolty
       do j = 0, 3
          call calculate_statistics(blk_avg_setedist(j,itype,1:nblock),avg(j),sd(j),se(j))
       end do
       write(io_output,"(2(A,I2),A,12(1X,F7.3))") ' sete & x,y,z length itype ',itype,' box ',1,' = ',(avg(j),sd(j),se(j),j=0,3)
    end do

  end subroutine write_prop_widom_with_stats

  subroutine read_checkpoint_prop_widom(io_chkpt,blockm)
    use util_mp,only:mp_bcast
    use sim_system,only:myid,rootid,groupid,ntmax
    integer,intent(in)::io_chkpt,blockm
    if (myid.eq.rootid) read(io_chkpt) Wrosen,blk_avg_Wrosen_prev,Uads,blk_avg_Uads_prev,blk_avg_Uads,setedist,blk_avg_setedist_prev,blk_avg_setedist
    call mp_bcast(Wrosen,ntmax,rootid,groupid)
    call mp_bcast(blk_avg_Wrosen_prev,ntmax,rootid,groupid)
    call mp_bcast(Uads,ntmax,rootid,groupid)
    call mp_bcast(blk_avg_Uads_prev,ntmax,rootid,groupid)
    call mp_bcast(blk_avg_Uads,ntmax*blockm,rootid,groupid)
    call mp_bcast(setedist,4*ntmax,rootid,groupid)
    call mp_bcast(blk_avg_setedist_prev,4*ntmax,rootid,groupid)
    call mp_bcast(blk_avg_setedist,4*ntmax*blockm,rootid,groupid)
  end subroutine read_checkpoint_prop_widom

  subroutine write_checkpoint_prop_widom(io_chkpt)
    integer,intent(in)::io_chkpt
    write(io_chkpt) Wrosen,blk_avg_Wrosen_prev,Uads,blk_avg_Uads_prev,blk_avg_Uads,setedist,blk_avg_setedist_prev,blk_avg_setedist
  end subroutine write_checkpoint_prop_widom
end module prop_widom
