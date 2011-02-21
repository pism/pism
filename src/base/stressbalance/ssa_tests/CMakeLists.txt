include_directories (../../../verif/tests/ ../../util ../../..)

if (Pism_BUILD_EXTRA_EXECS)
    add_executable (ssa_testi ssa_testi.cc ../../../verif/tests/exactTestsIJ.c ../../materials.cc)
    target_link_libraries (ssa_testi pismstressbalance)

    add_executable (ssa_testj ssa_testj.cc ../../../verif/tests/exactTestsIJ.c ../../materials.cc)
    target_link_libraries (ssa_testj pismstressbalance)

    add_executable (ssa_test_const ssa_test_const.cc ../../materials.cc)
    target_link_libraries (ssa_test_const pismstressbalance)

    add_executable (ssa_test_linear ssa_test_linear.cc ../../materials.cc)
    target_link_libraries (ssa_test_linear pismstressbalance)

    add_executable (ssa_test_plug ssa_test_plug.cc ../../materials.cc)
    target_link_libraries (ssa_test_plug pismstressbalance)

    install (TARGETS  ssa_testi ssa_testj ssa_test_const ssa_test_linear ssa_test_plug 
      DESTINATION ${Pism_BIN_DIR})
endif ()