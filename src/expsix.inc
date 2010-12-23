! EXPSIX.INC

      integer::natom
      parameter (natom = 3)
      real(8)::aexsix,bexsix,cexsix,sexsix,consp,consu
      dimension aexsix(natom),bexsix(natom),cexsix(natom)
     &   ,sexsix(natom),consu(natom),consp(natom)

      common /stuff/ aexsix,bexsix,cexsix,sexsix,consu,consp
