add_library (stressbalance OBJECT
  StressBalance.cc
  StressBalance_diagnostics.cc
  ShallowStressBalance.cc
  WeertmanSliding.cc
  SSB_Modifier.cc
  ssa/SSA.cc
  ssa/SSAFD.cc
  ssa/SSAFDBase.cc
  ssa/SSAFD_SNES.cc
  ssa/SSAFEM.cc
  ssa/SSATestCase.cc
  sia/BedSmoother.cc
  sia/SIAFD.cc
  sia/SIAFD_diagnostics.cc
  ../basalstrength/basal_resistance.cc
  factory.cc
)

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (siafd_test sia/siafd_test.cc)

  target_link_libraries (siafd_test pism)

  install (TARGETS
    siafd_test
    DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

add_subdirectory (ssa/tests)

add_subdirectory (blatter)
