file(GLOB_RECURSE Fortran_HEADERS "*.inc")
foreach(_headerFile ${Fortran_HEADERS})
  get_filename_component(_dir ${_headerFile} PATH)
  list(APPEND Fortran_INCLUDE_DIRS ${_dir})
endforeach()

if(Fortran_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES Fortran_INCLUDE_DIRS)
endif()
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${Fortran_INCLUDE_DIRS}")

set(ANALYSIS_SOURCES "global.F90" "rdf.F90" "hbond.F90")
add_library(analysis ${ANALYSIS_SOURCES})
target_link_libraries(analysis ${EXTRA_LIBRARIES})

set_target_properties(analysis PROPERTIES LINK_FLAGS "${Fortran_LINK_FLAGS}")
