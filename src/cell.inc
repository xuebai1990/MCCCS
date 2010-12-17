! CELL.INC

      real(8)::hmat(nbxmax,9),hmati(nbxmax,9)
      real(8)::cell_length(nbxmax,3),min_width(nbxmax,3)
      real(8)::cell_vol(nbxmax) 
      real(8)::cell_ang(nbxmax,3)
      common /cella/ hmat,hmati
      common /cellb/ cell_length,min_width,cell_vol,cell_ang
