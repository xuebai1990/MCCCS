# CMake project file for MCCCS Topmon

INCLUDE(CMakeForceCompiler)

cmake_minimum_required(VERSION 2.6)
project(MCCCS Fortran)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build, options are: None DEBUG RELEASE." FORCE)
endif()

# default installation
get_filename_component(default_prefix "." ABSOLUTE)
set(CMAKE_INSTALL_PREFIX "${default_prefix}" CACHE STRING "Choose the installation directory; by default it installs in the current directory.")

########################################################################
# Check and warn if cache generated on a different host is being reused
########################################################################
if(CMAKE_HOST_UNIX)
    execute_process(COMMAND hostname OUTPUT_VARIABLE TMP_HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(DEFINED BUILD_HOSTNAME AND NOT "${BUILD_HOSTNAME}" STREQUAL "${TMP_HOSTNAME}")
        message(WARNING "
            The CMake cache, probably generated on a different host (${BUILD_HOSTNAME}),
            is being reused! This could lead to inconsitencies; therefore, it is
            recommended to regenerate the cache!")
    endif()
    set(BUILD_HOSTNAME "${TMP_HOSTNAME}" CACHE INTERNAL "Hostname of the machine where the cache was generated." FORCE)
endif()

SET(SHARED_LIBS_DEFAULT ON)
IF(CMAKE_GENERATOR MATCHES "Visual Studio")
    SET(SHARED_LIBS_DEFAULT OFF)
ENDIF()
option(BUILD_SHARED_LIBS "Enable shared libraries (can be problematic with MPI, Windows)" ${SHARED_LIBS_DEFAULT})

if(UNIX AND NOT APPLE)
    option(PREFER_STATIC_LIBS "When finding libraries prefer \".a\" static archives (NOTE: this is enabled only for UNIX (excluding APPLE) platforms but it might not always work!" OFF)
    mark_as_advanced(PREFER_STATIC_LIBS)
    if(PREFER_STATIC_LIBS)
      SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
      if(BUILD_SHARED_LIBS)
        message(WARNING "Static libraries requested, the libraries will also be build static (BUILD_SHARED_LIBS=OFF).")
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "Enable shared libraries (can be problematic with MPI, Windows)." FORCE)
      endif()
    endif()
endif()

# FFLAGS depend on the compiler; or use CMAKE_Fortran_COMPILER_ID
get_filename_component(TMP_Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if(NOT "${Fortran_COMPILER_NAME}" STREQUAL "${TMP_Fortran_COMPILER_NAME}")
  set(Fortran_COMPILER_NAME "${TMP_Fortran_COMPILER_NAME}" CACHE INTERNAL "Compiler name when the cache was generated." FORCE)

  if(Fortran_COMPILER_NAME STREQUAL "gfortran")
    # GNU Compiler
    set(CMAKE_Fortran_FLAGS_RELEASE "-ffree-form -ffree-line-length-none -cpp -O3 -fomit-frame-pointer -finline-functions -ffast-math")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-ffree-form -ffree-line-length-none -cpp -O0 -ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
  elseif(Fortran_COMPILER_NAME STREQUAL "g77")
    # GNU Compiler
    # set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    # set(CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
  elseif(Fortran_COMPILER_NAME STREQUAL "g95")
    # g95
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fomit-frame-pointer -finline-functions -ffast-math")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-ftrace=full -pg -fbounds-check -mcpu=i686 -march=i686")
  elseif(Fortran_COMPILER_NAME STREQUAL "ifort")
    # Intel Compiler for x86 architecture
    set(CMAKE_Fortran_FLAGS_RELEASE "-free -fpp -O3 -xHost -ipo -no-prec-div")# -opt-streaming-stores always -no-unroll-aggressive -scalar-rep -Zp16 -nolib-inline")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-free -fpp -O0 -g -traceback -check all -warn unused -warn declarations -warn truncated_source -warn uncalled -warn ignore_loc -D__DEBUG__")
    # SGI Altix: Intel Compiler for ia64 architecture
    #set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -tpp2 -w90 -w95 -cm -f77rtl")
    #set(CMAKE_Fortran_FLAGS_DEBUG   "-check all -traceback -w90 -w95 -cm -CB -CU -p -g")
  elseif(Fortran_COMPILER_NAME STREQUAL "f95")
    # NAG Compiler for OSX
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2n")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -gline -f77 -mismatch_all -C=all -C=undefined -na")
  elseif(Fortran_COMPILER_NAME STREQUAL "pgf95")
    # Portland Compiler for x86 architecture
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -tp k8-64 -Mvect=sse")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-tp k8-64 -Mbounds -g -Mchkptr -check_bounds")
  elseif(Fortran_COMPILER_NAME STREQUAL "f77")
    # SGI Compiler for Origin
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-march=opteron -mips4 -r12000 -check_bounds -g")
  elseif(Fortran_COMPILER_NAME STREQUAL "xlf")
    # IBM AIX Compiler for SP
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -qstrict -bmaxdata:2000000000")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-check_bounds -g")
  else()
    message("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message("No optimized Fortran compiler flags are known, we just try -O2...")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif()

  set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE} CACHE STRING "Flags used by the compiler during RELEASE builds." FORCE)
  set(CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS_DEBUG} CACHE STRING "Flags used by the compiler during DEBUG builds." FORCE)
endif()

# Set Profiler
set(PROFILER CACHE STRING "Enable a compiler wrapper. Options: mpi, tau, scalasca.")
if(PROFILER STREQUAL "mpi")
  set(COMPILER_WRAPPER "mpif90")
elseif(PROFILER STREQUAL "tau")
  set(COMPILER_WRAPPER "tau_f90.sh")
elseif(PROFILER STREQUAL "scalasca")
  set(COMPILER_WRAPPER "skin_${Fortran_COMPILER_NAME}.sh")
  execute_process(COMMAND echo -e "#!/bin/sh\nskin ${Fortran_COMPILER_NAME} \"$@\"" OUTPUT_FILE ${COMPILER_WRAPPER})
  execute_process(COMMAND chmod +x "${COMPILER_WRAPPER}")
  message("
            NOTE: Please add the absolute path of the current
            directory, where ${COMPILER_WRAPPER} was just generated,
            to the environment variable PATH before running make.
")
endif()

# Set the compiler name to the wrapper, if necessary
if (DEFINED COMPILER_WRAPPER)
  CMAKE_FORCE_Fortran_COMPILER(${COMPILER_WRAPPER} ${CMAKE_Fortran_COMPILER_ID})
endif ()

if(WIN32)
  if (CMAKE_CL_64)
    set(ARCH_IS_64 "TRUE")
  else()
    set(ARCH_IS_64 "FALSE")
  endif()
else()
  get_property(ARCH_IS_64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)
endif()

set(EXTRA_LIBRARIES)

if(APPLE)
   find_library(ACCELERATE_FRAMEWORK Accelerate)
   list(APPEND EXTRA_LIBRARIES ${ACCELERATE_FRAMEWORK})
endif()

set(Fortran_INCLUDE_DIRS)
set(Fortran_LINK_FLAGS)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{FC_FLAGS}" CACHE STRING "Flags for Fortran compiler." FORCE)

option(USE_MPI "Build a parallel version using message passing interface" ON)
option(USE_OPENMP "Build a parallel version using threads (shared memory). Cannot combine with USE_MPI yet" OFF)

if(USE_MPI)
  if(USE_OPENMP)
    message(WARNING "Thread-based parallelization conflicts with MPI.")
    set(USE_OPENMP OFF CACHE BOOL "Thread-based parallelization conflicts with MPI." FORCE)
  endif()
  # fix the problem with ifort not being able to find 64-bit libraries
  if(${ARCH_IS_64} STREQUAL "TRUE" AND Fortran_COMPILER_NAME MATCHES "ifort*")
    list(APPEND CMAKE_LIBRARY_PATH "/usr/lib64")
  endif()

  find_package(MPI) # REQUIRED)
  if(MPI_Fortran_FOUND)
    add_definitions(-D__MPI__)
    list(APPEND EXTRA_LIBRARIES ${MPI_Fortran_LIBRARIES})
    list(APPEND Fortran_INCLUDE_DIRS ${MPI_Fortran_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
    set(Fortran_LINK_FLAGS "${Fortran_LINK_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
  else()
    message(FATAL_ERROR "MPI support requested, but no MPI compiler found.")
  endif()
endif()

if(USE_OPENMP)
  add_definitions(-D__OPENMP__)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp")
endif()

option(DOUBLE_PRECISION "Turn on to use double precision (real*8) floating numbers" ON)
if(DOUBLE_PRECISION)
  add_definitions(-D__DOUBLE_PRECISION__)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
endif()

option(USE_OWN "Use own implementation of certain math functions instead of system-provided ones" ON)
if(USE_OWN)
  add_definitions(-D__USEOWN__)
endif()

file(GLOB_RECURSE Fortran_HEADERS "*.inc")
foreach(_headerFile ${Fortran_HEADERS})
  get_filename_component(_dir ${_headerFile} PATH)
  list(APPEND Fortran_INCLUDE_DIRS ${_dir})
endforeach()

if(Fortran_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES Fortran_INCLUDE_DIRS)
endif()
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${Fortran_INCLUDE_DIRS}")

add_subdirectory(src)
