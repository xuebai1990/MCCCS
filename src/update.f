!
! this subroutine updates the block averages
!
subroutine update(nblock,ipos,ibox,value,acmove)
  use sim_system
  implicit none

  integer::nblock,ipos,ibox
  real::acmove,dpr,dn,value      
      
  if (nblock.eq.1) then
! first block
! write(io_output,*) 'acmove', acmove       
     dn = acmove
     dpr = value
     if ( dn .lt. 0.5E0_dp ) then
        baver(ipos,ibox,nblock) = 0.0E0_dp
     else
        baver(ipos,ibox,nblock) = dpr / dn
     end if
     nccold(ipos,ibox) = dn
     bccold(ipos,ibox) = dpr
     naccu(ipos,ibox) = naccu(ipos,ibox) + dn
     accum(ipos,ibox) = accum(ipos,ibox) + dpr
  else      
! other blocks
! write(io_output,*) 'acmove', acmove
     dn = acmove - nccold(ipos,ibox)
     dpr = value - bccold(ipos,ibox)
     if ( dn .lt. 0.5E0_dp ) then
        baver(ipos,ibox,nblock) = 0.0E0_dp
     else
        baver(ipos,ibox,nblock) = dpr / dn
     end if
     nccold(ipos,ibox) = acmove
     bccold(ipos,ibox) = value
     naccu(ipos,ibox) = naccu(ipos,ibox) + dn
     accum(ipos,ibox) = accum(ipos,ibox) + dpr
  end if
  return
end subroutine update
