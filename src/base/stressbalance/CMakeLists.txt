include_directories (../../verif/tests/ ../util)

add_library (pismstressbalance
  PISMStressBalance.cc
  PISMStressBalance_diagnostics.cc
  ShallowStressBalance.cc
  SSB_Modifier.cc
  PISMBedSmoother.cc
  SIAFD.cc
  SIAFD_diagnostics.cc
  SSA.cc
  SSA_diagnostics.cc
  SSAFD.cc
  SSAFEM.cc  # uncomment to compile
  SSAFEM_util.cc
  SSATestCase.cc
)
target_link_libraries (pismstressbalance pismutil pismflowlaws)

install (TARGETS pismstressbalance DESTINATION ${Pism_LIB_DIR})

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (ssafd_test ssafd_test.cc ../../verif/tests/exactTestsIJ.c ../materials.cc)
  target_link_libraries (ssafd_test pismstressbalance)

  add_executable (ssa_testi ssa_testi.cc ../../verif/tests/exactTestsIJ.c ../materials.cc)
  target_link_libraries (ssa_testi pismstressbalance)

  add_executable (ssa_testj ssa_testj.cc ../../verif/tests/exactTestsIJ.c ../materials.cc)
  target_link_libraries (ssa_testj pismstressbalance)

  add_executable (siafd_test siafd_test.cc ../../verif/tests/exactTestsFG.c ../materials.cc)
  target_link_libraries (siafd_test pismstressbalance)

  install (TARGETS ssafd_test siafd_test ssa_testi ssa_testj
    DESTINATION ${Pism_BIN_DIR})
endif ()
