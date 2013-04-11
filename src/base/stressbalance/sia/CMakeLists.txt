add_library(pism_sia STATIC
  PISMBedSmoother.cc
  PISMBedSmoother.hh
  SIAFD.cc
  SIAFD.hh
  SIAFD_diagnostics.cc
  SIA_Sliding.cc
  SIA_Sliding.hh
)

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (siafd_test
    siafd_test.cc
    ${PROJECT_SOURCE_DIR}/src/verif/tests/exactTestsFG.c
    ${PROJECT_SOURCE_DIR}/src/base/basalstrength/basal_resistance.cc)

  target_link_libraries (siafd_test pismstressbalance)

  install (TARGETS
    siafd_test
    DESTINATION ${Pism_BIN_DIR})
endif ()
