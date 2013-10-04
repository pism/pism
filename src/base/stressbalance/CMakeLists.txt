include_directories (../../verif/tests/ ../util)

add_library (pismstressbalance
  PISMStressBalance.cc
  PISMStressBalance_diagnostics.cc
  ShallowStressBalance.cc
  SSB_Modifier.cc
  ssa/FETools.cc
  ssa/SSA.cc
  ssa/SSAFD.cc
  ssa/SSAFEM.cc
  ssa/SSATestCase.cc
  sia/PISMBedSmoother.cc
  sia/SIAFD.cc
  sia/SIAFD_diagnostics.cc
  sia/SIA_Sliding.cc
)
target_link_libraries (pismstressbalance
  pismutil
  pismflowlaws
  pismrevision)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismstressbalance DESTINATION ${Pism_LIB_DIR})
endif()

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (siafd_test
    sia/siafd_test.cc
    ${PROJECT_SOURCE_DIR}/src/verif/tests/exactTestsFG.c
    ${PROJECT_SOURCE_DIR}/src/base/basalstrength/basal_resistance.cc)

  target_link_libraries (siafd_test pismstressbalance)

  install (TARGETS
    siafd_test
    DESTINATION ${Pism_BIN_DIR})
endif ()

if (Pism_BUILD_EXTRA_EXECS)
  add_executable(blatter_test
    blatter/Blatter_implementation.c
    blatter/BlatterStressBalance.cc
    blatter/blatter_test.cc)
  target_link_libraries (blatter_test pismutil pismbase pismstressbalance)
endif()

add_subdirectory (ssa/tests)
