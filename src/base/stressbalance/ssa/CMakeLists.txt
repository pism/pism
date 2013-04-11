add_library(pism_ssa STATIC
  FEEvaluator.cc
  FETools.cc
  SSA.cc
  SSAFD.cc
  SSAFEM.cc
  SSATestCase.cc
  )

add_subdirectory(ssa_tests)
