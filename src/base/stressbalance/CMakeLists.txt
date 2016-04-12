add_library (pismstressbalance OBJECT
  PISMStressBalance.cc
  PISMStressBalance_diagnostics.cc
  ShallowStressBalance.cc
  SSB_Modifier.cc
  ssa/FETools.cc
  ssa/node_types.cc
  ssa/SSA.cc
  ssa/SSAFD.cc
  ssa/SSAFEM.cc
  ssa/SSATestCase.cc
  sia/PISMBedSmoother.cc
  sia/SIAFD.cc
  sia/SIAFD_diagnostics.cc
  ../basalstrength/basal_resistance.cc
)

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (siafd_test sia/siafd_test.cc)

  target_link_libraries (siafd_test pism)

  install (TARGETS
    siafd_test
    DESTINATION ${Pism_BIN_DIR})
endif ()

add_subdirectory (ssa/tests)
