!
! this subroutine updates the block averages
!
subroutine update1(nblock,ipos,value,acmove,ibox,jbox)
  use sim_system
  implicit none

  integer::nblock,ipos,ibox,jbox
  real::acmove,dpr,dn,value
      
  if (nblock.eq.1) then
! first block
     dn = acmove
     dpr = value
     if ( dn .lt. 0.5E0_dp ) then
        baver1(ipos,ibox,jbox,nblock) = 0.0E0_dp
     else
        baver1(ipos,ibox,jbox,nblock) = dpr / dn
     end if
     nccold1(ipos,ibox,jbox) = dn
     bccold1(ipos,ibox,jbox) = dpr
     naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
     accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dpr
  else      
! other blocks
     dn = acmove - nccold1(ipos,ibox,jbox)
     dpr = value - bccold1(ipos,ibox,jbox)
     if ( dn .lt. 0.5E0_dp ) then
        baver1(ipos,ibox,jbox,nblock) = 0.0E0_dp
     else
        baver1(ipos,ibox,jbox,nblock) = dpr / dn
     end if
     nccold1(ipos,ibox,jbox) = acmove
     bccold1(ipos,ibox,jbox) = value
     naccu1(ipos,ibox,jbox) = naccu1(ipos,ibox,jbox) + dn
     accum1(ipos,ibox,jbox) = accum1(ipos,ibox,jbox) + dpr
  end if
  return
end subroutine update1



















