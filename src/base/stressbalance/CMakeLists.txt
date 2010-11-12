include_directories (../../verif/tests/)

add_library (pismstressbalance
  PISMStressBalance.cc
  ShallowStressBalance.cc
  SSB_Modifier.cc
  SIAFD.cc
  SSAFD.cc
)
target_link_libraries (pismstressbalance pismmisc)

install (TARGETS pismstressbalance DESTINATION ${Pism_LIB_DIR})

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (ssafd_test ssafd_test.cc ../../verif/tests/exactTestsIJ.c)
  target_link_libraries (ssafd_test pismstressbalance)

  add_executable (siafd_test siafd_test.cc ../../verif/tests/exactTestsFG.c)
  target_link_libraries (siafd_test pismstressbalance)

  install (TARGETS ssafd_test siafd_test
    DESTINATION ${Pism_BIN_DIR})
endif ()
