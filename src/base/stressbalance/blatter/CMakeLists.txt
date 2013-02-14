if (Pism_BUILD_EXTRA_EXECS)
  add_executable(blatter_test
    Blatter_implementation.c
    BlatterStressBalance.cc
    blatter_test.cc
    IceModelVec3Scaled.cc)
  target_link_libraries (blatter_test pismutil pismbase pismstressbalance)
endif()
