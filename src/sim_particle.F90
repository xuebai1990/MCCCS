MODULE sim_particle
  use var_type,only:double_precision,default_string_length
  implicit none
  private
  public::BeadType,AtomType,MoleculeType

  type BeadType
     integer::type
     real(KIND=double_precision)::coord(3)
  end type BeadType

  type AtomType
     logical::lqchg,llj
     integer::atomID
     real(KIND=double_precision)::mass,charge,eps,sig
     character(LEN=default_string_length)::name,desc !name is the atomic symobl, desc is the full description of the atom
  end type AtomType

  type MoleculeType
     integer::nbead
     type(BeadType),allocatable::bead(:)
  end type MoleculeType
end MODULE sim_particle
