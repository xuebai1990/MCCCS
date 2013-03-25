!
! this subroutine updates the block averages
!
subroutine update1(nblock,ipos,value,acmove,ibox,jbox)
  use sim_system
  implicit none

  integer::nblock,ipos,ibox,jbox
  real::acmove,dp,dn,value
      
  if (nblock.eq.1) then
! first block
     dn = acmove
     dp = value
     if ( dn .lt. 0.5d0 ) then
        baver1(ipos,ibox,jbox,nblock) = 0.0d0
     else
        baver1(ipos,ibox,jbox,nblock) = dp / dn
     end if
     nccold1(ipos,ibox,jbox) = dn
     bccold1(ipos,ibox,jbox) = dp
     naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
     accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dp
  else      
! other blocks
     dn = acmove - nccold1(ipos,ibox,jbox)
     dp = value - bccold1(ipos,ibox,jbox)
     if ( dn .lt. 0.5d0 ) then
        baver1(ipos,ibox,jbox,nblock) = 0.0d0
     else
        baver1(ipos,ibox,jbox,nblock) = dp / dn
     end if
     nccold1(ipos,ibox,jbox) = acmove
     bccold1(ipos,ibox,jbox) = value
     naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
     accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dp
  end if
  return
end subroutine update1



















