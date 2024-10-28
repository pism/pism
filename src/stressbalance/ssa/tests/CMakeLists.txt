if (Pism_BUILD_EXTRA_EXECS)

  foreach (TEST IN ITEMS
      ssa_testi ssa_testj ssa_test_const ssa_test_linear ssa_test_plug ssa_test_cfbc)
    add_executable (${TEST} ${TEST}.cc)
    target_link_libraries (${TEST} libpism)
    install (TARGETS ${TEST} DESTINATION ${CMAKE_INSTALL_BINDIR})
  endforeach()

endif ()
