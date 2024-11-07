if (Pism_BUILD_EXTRA_EXECS)

  foreach (TEST IN ITEMS
      ssa_test_i ssa_test_j ssa_test_const ssa_test_linear ssa_test_plug ssa_test_cfbc)
    add_executable (pism_${TEST} ${TEST}.cc)
    target_link_libraries (pism_${TEST} libpism)
    install (TARGETS pism_${TEST} DESTINATION ${CMAKE_INSTALL_BINDIR})
  endforeach()

endif ()
