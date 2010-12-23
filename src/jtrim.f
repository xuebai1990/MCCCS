      character*(*) function JTRIM(str)
      character*(*) str
      JTRIM = str(1:len_trim(str))
      end function
