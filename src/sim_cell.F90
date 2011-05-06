MODULE sim_cell
  use var_type,only:double_precision,RealPtr
  implicit none
  private
  public::CellType,CellMaskType

  type CellType
     logical::ortho,solid ! ortho=.true. if the simulation box is orthorhombic
                          ! solid=.true. if the simulation box can change shape
     logical,dimension(3)::pbc
     real(KIND=double_precision)::cut,vol,calp,boxl(3),ang(3),height(3),hmat(3,3),hmati(3,3)
  end type CellType

  type CellMaskType
     logical,pointer::ortho,solid ! ortho=.true. if the simulation box is orthorhombic
                          ! solid=.true. if the simulation box can change shape
     logical,dimension(3)::pbc
     real(KIND=double_precision),pointer::cut,vol,calp
     type(RealPtr)::boxl(3),ang(3),height(3),hmat(3,3),hmati(3,3)
  end type CellMaskType
end MODULE sim_cell
