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

add_executable (pism_siafd_test sia/siafd_test.cc)
target_link_libraries (pism_siafd_test libpism)

install (TARGETS
  pism_siafd_test
  DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory (ssa/tests)

add_subdirectory (blatter)
