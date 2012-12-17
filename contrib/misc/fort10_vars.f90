!******************!
module fort10_vars !
!*****************************************************************************
! stores all variables in fort.10 movie file so analysis subroutine 
! calls do not have pass :)
!*****************************************************************************

  implicit none

  integer::nbox,nchain,nmolty,nhere,nframe,most_units,most_vibs,most_tors
  integer,dimension(:),allocatable::nunit,temphere,moltyp,ncycles
  integer,dimension(:,:),allocatable::nvib,ntor,molbox,beadtype
  integer,dimension(:,:,:),allocatable::ncmt,ijvib,itor2,itor3,itor4
  real,dimension(:),allocatable::rcut,maxl
  real,dimension(:,:),allocatable::boxlx,boxly,boxlz,xcom,ycom,zcom
  real,dimension(:,:,:),allocatable::xbead,ybead,zbead,qbead
      
!***VARIABLE DESCRIPTIONS***!
!   nbox            : number of boxes in the simulation
!   nchain          : number chains (molecules) in the simulation
!   nmolty          : number of molecule types 
!   nhere           : number of bead types
!   nframe          : number frames or snapshots from the simulation
!   nunit(i)        : number of beads in the molecule type i
!   temphere(i)     : type of bead i
!   moltyp(i)       : molecule type of chain i
!   ncycles(i)      : number of cycles at frame was taken
!   nvib(i,j)       : number of vibration for chain i bead j
!   ntor(i,j)       : number of torsions for chain i bead j
!   molbox(i,j)     : box of chain j in frame i
!   beadtype(i,j)   : type of bead for molecule type i, bead j
!   ncmt(i,j,k)     : number of chains k in box j for frame i
!   ijvib(i,j,k)    : bead involved in vibration with chain i, bead j, vib k
!   itor2(i,j,k)    : 1st bead involved in torsion with chain i, bead j, tors k
!   itor3(i,j,k)    : 2nd bead involved in torsion with chain i, bead j, tors k
!   itor4(i,j,k)    : 3rd bead involved in torsion with chain i, bead j, tors k
!   rcut(i)         : potential cutoff in box i
!   boxlx(i,j)      : x egdelength of box j in frame i
!   boxly(i,j)      : y egdelength of box j in frame i
!   boxlz(i,j)      : z egdelength of box j in frame i
!   xcom(i,j)       : x center of mass position for chain j in frame i
!   ycom(i,j)       : y center of mass position for chain j in frame i
!   zcom(i,j)       : z center of mass position for chain j in frame i
!   xbead(i,j,k)    : x bead position for bead k of chain j in frame i
!   ybead(i,j,k)    : y bead position for bead k of chain j in frame i
!   zbead(i,j,k)    : z bead position for bead k of chain j in frame i
!   qbead(i,j,k)    : charge for bead k of chain j in frame i
!   most_units      : most beads in any molecule
!   most_tors       : most torsions in any molecule
!   most_vibs       : most vibrations in any molecule
!   maxl            : maximum box lengths in any direction (x,y,z)
end module fort10_vars
