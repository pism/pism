add_library (blatter OBJECT
  Blatter.cc
  residual.cc
  jacobian.cc
  BlatterMod.cc
  util/grid_hierarchy.cc
  verification/BlatterTestXY.cc
  verification/BlatterTestXZ.cc
  verification/BlatterTestCFBC.cc
  verification/BlatterTestHalfar.cc
  verification/BlatterTestvanderVeen.cc
  verification/manufactured_solutions.cc
  ismip-hom/BlatterISMIPHOM.cc
  )

add_executable (pism_blatter_test
  verification/blatter_test.cc
  verification/BlatterTestCase.cc)
target_link_libraries (pism_blatter_test libpism)
install (TARGETS pism_blatter_test DESTINATION ${CMAKE_INSTALL_BINDIR})
