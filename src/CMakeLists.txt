# CMake project file for MCCCS Topmon

cmake_minimum_required (VERSION 2.6)
set (CMAKE_Fortran_COMPILER "mpif90")
project (MCCCS Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None DEBUG RELEASE."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the current directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "mpif90*")
  # mpif90 wrapper
  set (CMAKE_Fortran_FLAGS_RELEASE "-free -fpp -O3 -xHost -ipo -no-prec-div")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-free -fpp -O0 -g -traceback -check all -warn unused -warn declarations -warn truncated_source -warn uncalled -warn ignore_loc -D__DEBUG__")
elseif (Fortran_COMPILER_NAME MATCHES "ifort*")
  # Intel Compiler for x86 architecture
  set (CMAKE_Fortran_FLAGS_RELEASE "-free -fpp -O3 -xHost -ipo -no-prec-div")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-free -fpp -O0 -g -traceback -check all -warn unused -warn declarations -warn truncated_source -warn uncalled -warn ignore_loc -D__DEBUG__")
  # SGI Altix: Intel Compiler for ia64 architecture
  #set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -tpp2 -w90 -w95 -cm -f77rtl")
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-check all -traceback -w90 -w95 -cm -CB -CU -p -g")
elseif (Fortran_COMPILER_NAME MATCHES "f95*")
  # NAG Compiler for OSX
  set (CMAKE_Fortran_FLAGS_RELEASE "-O0 -gline -f77 -mismatch_all -C=all -C=undefined -nan")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O2")
elseif (Fortran_COMPILER_NAME MATCHES "g95*")
  # g95
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
elseif (Fortran_COMPILER_NAME MATCHES "gfortran*")
  # GNU Compiler
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
elseif (Fortran_COMPILER_NAME MATCHES "g77*")
  # GNU Compiler
  # set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
elseif (Fortran_COMPILER_NAME MATCHES "pgf95*")
  # Portland Compiler for x86 architecture
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -tp k8-64 -Mvect=sse")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-tp k8-64 -Mbounds -g -Mchkptr -check_bounds")
elseif (Fortran_COMPILER_NAME MATCHES "f77*")
  # SGI Compiler for Origin
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-march=opteron -mips4 -r12000 -check_bounds -g")
elseif (Fortran_COMPILER_NAME MATCHES "xlf*")
  # IBM AIX Compiler for SP
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict -bmaxdata:2000000000")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-check_bounds -g")
else (Fortran_COMPILER_NAME MATCHES "mpif90*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "mpif90*")

#find_package(MPI REQUIRED)

#include_directories(${MPI_Fortran_INCLUDE_PATH})

# build executables
set (NMPROGRAMS "topmon.f" "global_data.f" "var_type.f" "const_phys.f" "const_math.f" "util_math.f" "util_string.f" "util_files.f" "util_timings.f" "sim_cell.F90" "sim_particle.F90" "sim_zeolite.F90" "parser_cssr.F90" "parser_pdb.F90" "parser_cif.F90" "Atom_energy.f" "Atom_traxyz.f" "Intra_energy.f" "alignplanes.f" "anes.f" "boltz.f" "bondlength.f" "calcsolpar.f" "calctor.f" "charge.f" "cleanup.f" "close.f" "cone.f" "coneangle.f" "config.f" "corp.f" "coru.f" "ctrmas.f" "deriv.f" "dipole.f" "dump.f" "ee_energy.f" "ee_index_swap.f" "ee_recip.f" "eemove.f" "eesetup.f" "energy.f" "exgrph.f" "exmuir.f" "expand.f" "explct.f" "exsami.f" "exsix.f" "flucq.f" "garo3body.f" "garofalini.f" "genlj.f" "geometry.f" "inclus.f" "init_vars.f" "initia.f" "ipswsetup.f" "linkcell.f" "ljmuir.f" "ljpsur.f" "ljsami.f" "matops.f" "mimage.f" "mmff.f" "molsetup.f" "monola.f" "monper.f" "newton.f" "ninesix.f" "place.f" "pressure.f" "prvolume.f" "random.f" "readdat.f" "recip.f" "recip_atom.f" "recippress.f" "recipsum.f" "rigfix.f" "rigrot.f" "rosenbluth.f" "rotxyz.f" "safecbmc.f" "safeschedule.f" "schedule.f" "setnn.f" "setpbc.f" "slitpore.f" "sphere.f" "suijtab.f" "sumup.f" "susami.f" "suvibe.f" "swap.f" "swatch.f" "tabulated.f" "torsion.f" "traxyz.f" "update.f" "update1.f" "updnn.f" "v_elect_field.f" "virial.f" "volume.f" "vtorso.f" "writepdb.f" "zeolite.f" "fparser.F90")
set (EXECUTABLES "topmon")
set (SCRIPTS "reconf")

add_executable (topmon ${NMPROGRAMS})
#foreach (p ${NMPROGRAMS})
#  add_executable (${p} "${p}.f")
#endforeach (p)

#target_link_libraries(topmon ${MPI_Fortran_LIBRARIES})

#if(MPI_COMPILE_FLAGS)
#  set_target_properties(topmon PROPERTIES
#    COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
#endif()

#if(MPI_LINK_FLAGS)
#  set_target_properties(topmon PROPERTIES
#    LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#endif()

# install executables and scripts
install (TARGETS ${EXECUTABLES} 
         RUNTIME DESTINATION "bin")
install (PROGRAMS ${SCRIPTS}
         DESTINATION "bin") 
