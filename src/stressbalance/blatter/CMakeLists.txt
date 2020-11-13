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
  verification/manufactured_solutions.cc
  ismip-hom/BlatterISMIPHOM.cc
  )
