# build executables
set(MCCCS_SOURCES "var_type.F90" "const_math.F90" "const_phys.F90" "util_math.F90" "util_search.F90" "util_string.F90" "util_files.F90" "util_timings.F90" "util_memory.F90" "util_runtime.F90" "util_random.F90" "util_prng.F90" "util_mp.F90" "util_mp_aux.F90" "sim_system.F90" "sim_cell.F90" "sim_particle.F90" "sim_initia.F90" "sim_zeolite.F90" "fparser.F90" "parser_cif.F90" "parser_pdb.F90" "parser_cssr.F90" "zeolite.F90" "energy_kspace.F90" "energy_intramolecular.F90" "energy_external.F90" "energy_pairwise.F90" "energy_sami.F90" "energy_garofalini.F90" "energy_3body.F90" "energy_4body.F90" "moves_simple.F90" "moves_volume.F90" "moves_cbmc.F90" "moves_ee.F90" "transfer_shared.F90" "transfer_swap.F90" "transfer_swatch.F90" "prop_pressure.F90" "topmon_main.F90" "topmon.F90" "anes.f" "calcsolpar.f" "charge.f" "deriv.f" "dipole.f" "dump.F90" "flucq.f" "inclus.F90" "ipswsetup.f" "molsetup.f" "virial.f")
add_executable(topmon ${MCCCS_SOURCES})
target_link_libraries(topmon ${EXTRA_LIBRARIES})
set_target_properties(topmon PROPERTIES LINK_FLAGS "${Fortran_LINK_FLAGS}")

# install executables and scripts
set(EXECUTABLES "topmon")
install(TARGETS ${EXECUTABLES} RUNTIME DESTINATION ".")

set(SCRIPTS "../tests/LONG/sub.script" "../tests/LONG/qbatch.script")
install(PROGRAMS ${SCRIPTS} DESTINATION ".")
