! ROSEN.INC

      integer::growfrom,growprev,growlist,grownum
      dimension growfrom(numax),growprev(numax),grownum(numax)
      dimension growlist(numax,numax)

      real(8)::weight,  weiold, voldt,voldbb,voldtg,voldext
     &     ,voldintra,voldinter,voldbvib,voldelect,voldewald
     &     ,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib
     &     ,vnewinter,vnewelect,vnewewald,vneworient
     &     ,vnewinterr,vnewextr,vnewelectr,voldinterr
     &     ,voldextr,voldelectr,voldorient
      common /rosena/ growfrom,growprev,growlist,grownum

      common /rosenb/ weight,  weiold, voldt,voldbb,voldtg,voldext
     &     ,voldintra,voldinter,voldbvib,voldelect,voldewald
     &     ,vnewt,vnewbb,vnewtg,vnewext,vnewintra,vnewbvib
     &     ,vnewinter,vnewelect,vnewewald,voldelectr
     &     ,vnewinterr,vnewextr,vnewelectr,voldinterr
     &     ,voldextr,vneworient,voldorient
    

