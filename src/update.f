!
! this subroutine updates the block averages
!
subroutine update(nblock,ipos,ibox,value,acmove)
  use sim_system
  implicit none

  integer::nblock,ipos,ibox
  real::acmove,dp,dn,value      
      
  if (nblock.eq.1) then
! first block
! write(io_output,*) 'acmove', acmove       
     dn = acmove
     dp = value
     if ( dn .lt. 0.5d0 ) then
        baver(ipos,ibox,nblock) = 0.0d0
     else
        baver(ipos,ibox,nblock) = dp / dn
     end if
     nccold(ipos,ibox) = dn
     bccold(ipos,ibox) = dp
     naccu(ipos,ibox) = naccu(ipos,ibox) + dn
     accum(ipos,ibox) = accum(ipos,ibox) + dp
  else      
! other blocks
! write(io_output,*) 'acmove', acmove
     dn = acmove - nccold(ipos,ibox)
     dp = value - bccold(ipos,ibox)
     if ( dn .lt. 0.5d0 ) then
        baver(ipos,ibox,nblock) = 0.0d0
     else
        baver(ipos,ibox,nblock) = dp / dn
     end if
     nccold(ipos,ibox) = acmove
     bccold(ipos,ibox) = value
     naccu(ipos,ibox) = naccu(ipos,ibox) + dn
     accum(ipos,ibox) = accum(ipos,ibox) + dp
  end if
  return
end subroutine update
