  if (.not.allocated(p)) then
     allocate(p(pos:pos))
  else if (pos.gt.ubound(p,1)) then
     call reallocate(p,lbound(p,1),pos)
  else if (pos.lt.lbound(p,1)) then
     call reallocate(p,pos,ubound(p,1))
  end if
  p(pos)=val
