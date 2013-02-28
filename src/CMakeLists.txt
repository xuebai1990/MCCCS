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
  set (CMAKE_Fortran_FLAGS_RELEASE "-free -fpp -O3 -xHost -ipo -no-prec-div -integer-size 32 -real-size 64 -D__MPI__")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-free -fpp -O0 -integer-size 32 -real-size 64 -g -traceback -check all -warn unused -warn declarations -warn truncated_source -warn uncalled -warn ignore_loc -D__MPI__ -D__DEBUG__")
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
set (NMPROGRAMS "var_type.F90" "const_math.F90" "const_phys.F90" "util_math.F90" "util_search.F90" "util_string.F90" "util_files.F90" "util_timings.F90" "util_memory.F90" "util_runtime.F90" "util_random.F90" "util_prng.F90" "util_mp.F90" "util_mp_aux.F90" "sim_system.F90" "sim_cell.F90" "sim_particle.F90" "sim_zeolite.F90" "fparser.F90" "parser_cif.F90" "parser_pdb.F90" "parser_cssr.F90" "zeolite.F90" "energy_kspace.F90" "energy_intramolecular.F90" "energy_external.F90" "energy_pairwise.F90" "energy_sami.F90" "energy_garofalini.F90" "energy_3body.F90" "energy_4body.F90" "moves_simple.F90" "moves_cbmc.F90" "moves_ee.F90" "transfer_shared.F90" "transfer_swap.F90" "transfer_swatch.F90" "prop_pressure.F90" "topmon.F90" "alignplanes.f" "anes.f" "calcsolpar.f" "charge.f" "ctrmas.f" "deriv.f" "dipole.f" "dump.f" "flucq.f" "inclus.f" "initia.f" "ipswsetup.f" "molsetup.f" "monola.f" "monper.f" "readdat.f" "update.f" "update1.f" "virial.f")
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
