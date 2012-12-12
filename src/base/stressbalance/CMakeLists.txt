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
  SSAFD.cc
  SSAFEM.cc
  FETools.cc
  FEEvaluator.cc
  SSATestCase.cc
  SIA_Sliding.cc
#  BlatterStressBalance.cc
#  THI.cc
)
target_link_libraries (pismstressbalance pismutil pismflowlaws pismrevision)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismstressbalance DESTINATION ${Pism_LIB_DIR})
endif()

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (siafd_test siafd_test.cc ../../verif/tests/exactTestsFG.c ../basalstrength/basal_resistance.cc)
  target_link_libraries (siafd_test pismstressbalance)

  # add_executable (blatter_test blatter_test.cc ../basal_strength/basal_resistance.cc)
  # target_link_libraries (blatter_test pismstressbalance)

  install (TARGETS
    siafd_test
#    blatter_test
    DESTINATION ${Pism_BIN_DIR})
endif ()

add_subdirectory (ssa_tests)

