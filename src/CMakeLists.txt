# CMake project file for MCCCS Topmon

cmake_minimum_required (VERSION 2.6)
set(CMAKE_Fortran_COMPILER "mpif90")
project (topmon)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the topmon directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "mpif90")
  # mpif90 wrapper
  #set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # Intel Compiler for x86 architecture
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -w90 -w95 -cm")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-w90 -w95 -cm -CB -CU -p -g -check_bounds -check all -traceback")
  # SGI Altix: Intel Compiler for ia64 architecture
  #set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -tpp2 -w90 -w95 -cm")
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-check all -traceback -w90 -w95 -cm -CB -CU -g")
elseif (Fortran_COMPILER_NAME STREQUAL "f95")
  # NAG Compiler for OSX
  set (CMAKE_Fortran_FLAGS_RELEASE "-O0 -gline -f77 -mismatch_all -C=all -C=undefined -nan")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O2")
elseif (Fortran_COMPILER_NAME STREQUAL "g95")
  # g95
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
elseif (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # GNU Compiler
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # GNU Compiler
  # set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
elseif (Fortran_COMPILER_NAME STREQUAL "pgf95")
  # Portland Compiler for x86 architecture
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -tp k8-64 -Mvect=sse")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-tp k8-64 -Mbounds -g -Mchkptr -check_bounds")
elseif (Fortran_COMPILER_NAME STREQUAL "f77")
  # SGI Compiler for Origin
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-march=opteron -mips4 -r12000 -check_bounds -g")
elseif (Fortran_COMPILER_NAME STREQUAL "xlf")
  # IBM AIX Compiler for SP
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict -bmaxdata:2000000000")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-check_bounds -g")
else (Fortran_COMPILER_NAME STREQUAL "mpif90")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME STREQUAL "mpif90")

# build executables
set (NMPROGRAMS "topmon.f" "Atom_energy.f" "Atom_traxyz.f" "Intra_energy.f" "alignplanes.f" "anes.f" "atomtype.f" "boltz.f" "bondlength.f" "calcsolpar.f" "calctor.f" "charge.f" "chempt.f" "cleanup.f" "close.f" "cone.f" "coneangle.f" "config.f" "corp.f" "coru.f" "coruz.f" "ctrmas.f" "deriv.f" "dipole.f" "dump.f" "ee_energy.f" "ee_index_swap.f" "ee_recip.f" "eemove.f" "eesetup.f" "energy.f" "erfunc.f" "exgrph.f" "exmuir.f" "expand.f" "explct.f" "exsami.f" "exsix.f" "exzeo.f" "exzeof.f" "flucq.f" "garo3body.f" "garofalini.f" "genlj.f" "geometry.f" "grid.f" "inclus.f" "init_vars.f" "initia.f" "ipswsetup.f" "jtrim.f" "linkcell.f" "ljmuir.f" "ljpsur.f" "ljsami.f" "matops.f" "mbessel.f" "mimage.f" "mmff.f" "molsetup.f" "monola.f" "monper.f" "newton.f" "ninesix.f" "nrtab.f" "pgrid.f" "place.f" "pressure.f" "prvolume.f" "random.f" "readdat.f" "recip.f" "recip_atom.f" "recippress.f" "recipsum.f" "rigfix.f" "rigrot.f" "rosenbluth.f" "rotxyz.f" "safecbmc.f" "safeschedule.f" "schedule.f" "setnn.f" "setpbc.f" "slitpore.f" "sphere.f" "suijtab.f" "sumup.f" "susami.f" "suvibe.f" "suzeo.f" "swap.f" "swatch.f" "tabulated.f" "torsion.f" "traxyz.f" "trim_swap.f" "update.f" "update1.f" "updnn.f" "v_elect_field.f" "virial.f" "volume.f" "vtorso.f" "writepdb.f" "zeocoord.f" "ztest.f")
#set (EXECUTABLES "topmon" ${NMPROGRAMS})
set (EXECUTABLES "topmon")
set (SCRIPTS "reconf")

add_executable ("topmon" ${NMPROGRAMS})
#foreach (p ${NMPROGRAMS})
#  add_executable (${p} "${p}.f")
#endforeach (p)

# install executables and scripts
install (TARGETS ${EXECUTABLES} 
         RUNTIME DESTINATION "bin")
install (PROGRAMS ${SCRIPTS}
         DESTINATION "bin") 
