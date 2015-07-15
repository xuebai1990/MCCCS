set(ANALYSIS_SOURCES "global.F90" "rdf.F90" "hbond.F90")
add_library(analysis ${ANALYSIS_SOURCES})
target_link_libraries(analysis ${EXTRA_LIBRARIES})
set_target_properties(analysis PROPERTIES LINK_FLAGS "${Fortran_LINK_FLAGS}")

set(AVGSYMM_SOURCES "../../../src/var_type.F90" "../../../src/const_math.F90" "../../../src/util_string.F90" "../../../src/util_files.F90" "../../../src/util_runtime.F90" "../../../src/fparser.F90")
add_executable(avgsymm ${AVGSYMM_SOURCES} "avgSymm.F90")
target_link_libraries(avgsymm ${EXTRA_LIBRARIES})
set_target_properties(avgsymm PROPERTIES LINK_FLAGS "${Fortran_LINK_FLAGS}")
