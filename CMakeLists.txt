# CMake project file for MCCCS Topmon

cmake_minimum_required (VERSION 2.6)
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

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  # set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check")
elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # ifort (untested)
  # set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  set (CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -check bounds")
elseif (Fortran_COMPILER_NAME STREQUAL "mpif90")
  # mpif90 wrapper
  set (CMAKE_Fortran_FLAGS_RELEASE "-ipo -O3 -no-prec-div -xHost")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -check bounds")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  # set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  # set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

#########################
### Compiler definitions
#########################

### Intel Compiler for x86 architecture
#COMPFLAGS      =  -w90  -check all -traceback 
#COMPFLAGS     =     -w90 -w95 -cm -check all -traceback#-check_bounds# -g
#COMPFLAGS     =  -O2 -w90 -w95 -cm -p#-check_bounds# -g
## only for debugging
#COMPFLAGS     = -w90 -w95 -cm -CB -CU -g

#COMPFLAGS     =  -O2 -w90 -w95 #-cm -p
#COMPFLAGS     =  -ipo -O3 -no-prec-div -xHost 
#COMPFLAGS     =  -O3 -no-prec-div -check bounds 
#LD            = mpif90 
#LDS           = mpif90 -c $(COMPFLAGS)
#LDFLAGS       =  $(COMPFLAGS)

### NAG Compiler for OSX
#COMPFLAGS     = -O0 -gline -f77 -mismatch_all -C=all -C=undefined -nan
## only for debugging
#COMPFLAGS     = -O2 
#LD            = f95
#LDS           = f95 -c $(COMPFLAGS)
#LDFLAGS       = $(COMPFLAGS)

### SGI Altix: Intel Compiler for ia64 architecture
#COMPFLAGS     = -O2 -tpp2 -w90 -w95 -cm #-check all -traceback
# only for debugging
#COMPFLAGS     =  -check all -traceback -w90 -w95 -cm -CB -CU -g
#LD            = ifort
#LDS           = ifort -c $(COMPFLAGS)
#LDFLAGS       = -static $(COMPFLAGS)

### GNU Compiler
#COMPFLAGS     =  -ftrace=full -pg -fbounds-check
#COMPFLAGS     = -O3 -fomit-frame-pointer -finline-functions -ffast-math 
#LD	      = g95
#LDS 	      = g95 -c $(COMPFLAGS)
#LDFLAGS	      = $(COMPFLAGS)

### Portland Compiler for x86 architecture
#COMPFLAGS     =  -O2 -tp k8-64 -Mvect=sse#-check_bounds# -g
#COMPFLAGS     =  -tp k8-64 -Mbounds -g -Mchkptr#-check_bounds# -g
#LD	      = /opt/pgi/linux86-64/6.1/bin/pgf95
#LDS 	      = /opt/pgi/linux86-64/6.1/bin/pgf95 -c $(COMPFLAGS)
#LDFLAGS	      = $(COMPFLAGS)

### SGI Compiler for Origin
#COMPFLAGS     = -O2 #-march=opteron # -mips4 -r12000 #-check_bounds# -g
#LD            = f77 
#LDS           = f77 -c $(COMPFLAGS)
#LDFLAGS       = $(COMPFLAGS)

### IBM AIX Compiler for SP
#COMPFLAGS     = -O3 -qstrict  -bmaxdata:2000000000 #-check_bounds# -g
#LD            = xlf
#LDS           = xlf -c $(COMPFLAGS)
#LDFLAGS	      = $(COMPFLAGS)

### GNU Compiler
#COMPFLAGS     =  -fbounds-check -g
#COMPFLAGS     = -O3 #-fomit-frame-pointer -finline-functions -ffast-math #-fbounds-check #-mcpu=i686 -march=i686# 
#LD	      = g77
#LDS 	      = g77 -c $(COMPFLAGS)
#LDFLAGS	      = $(COMPFLAGS)

#############################
### End compiler definitions 
#############################

# build executables
set (NMPROGRAMS "topmon.f" "Atom_energy.f" "Atom_traxyz.f" "Intra_energy.f" "alignplanes.f" "analysis.f" "anes.f" "atomtype.f" "boltz.f" "bondlength.f" "calcsolpar.f" "calctor.f" "charge.f" "chempt.f" "cleanup.f" "close.f" "cone.f" "coneangle.f" "config.f" "corp.f" "coru.f" "coruz.f" "ctrmas.f" "deriv.f" "dipole.f" "dump.f" "ee_energy.f" "ee_index_swap.f" "ee_recip.f" "eemove.f" "eesetup.f" "energy.f" "erfunc.f" "exgrph.f" "exmuir.f" "expand.f" "explct.f" "exsami.f" "exsix.f" "exzeo.f" "exzeof.f" "exzeot.f" "flucq.f" "forcefield.f" "garo3body.f" "garofalini.f" "genlj.f" "geometry.f" "grid.f" "inclus.f" "init_vars.f" "initia.f" "ipswsetup.f" "jtrim.f" "linkcell.f" "ljmuir.f" "ljpsur.f" "ljsami.f" "matops.f" "mbessel.f" "mimage.f" "mmff.f" "molsetup.f" "monola.f" "monper.f" "newton.f" "ninesix.f" "nrtab.f" "pgrid.f" "place.f" "pressure.f" "prvolume.f" "qqcheck.f" "random.f" "readRES.f" "readdat.f" "recip.f" "recip_atom.f" "recippress.f" "recipsum.f" "rigfix.f" "rigrot.f" "rosenbluth.f" "rotxyz.f" "safecbmc.f" "safeschedule.f" "schedule.f" "setnn.f" "setpbc.f" "slitpore.f" "sphere.f" "suijtab.f" "sumup.f" "susami.f" "suvibe.f" "suzeo.f" "swap.f" "swatch.f" "tabulated.f" "test1.f" "testupdate.f" "thermopress.f" "torsion.f" "traxyz.f" "trim_swap.f" "update.f" "update1.f" "updnn.f" "v_elect_field.f" "virial.f" "volume.f" "vtorso.f" "writepdb.f" "zeocoord.f" "ztest.f")
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
