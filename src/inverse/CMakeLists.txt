# Classes for inverse problems.

set(PISM_INVERSE_SOURCES
  ./IPDesignVariableParameterization.cc
  ./IP_SSATaucForwardProblem.cc
  ./IP_SSAHardavForwardProblem.cc
  ./IPTwoBlockVec.cc
  ./IP_SSATaucTikhonovGNSolver.cc  
  ./functional/IPFunctional.cc
  ./functional/IP_L2NormFunctional.cc
  ./functional/IP_H1NormFunctional.cc
  ./functional/IPGroundedIceH1NormFunctional.cc
  ./functional/IPLogRatioFunctional.cc
  ./functional/IPLogRelativeFunctional.cc
  ./functional/IPTotalVariationFunctional.cc
  ./functional/IPMeanSquareFunctional.cc
  ./TaoUtil.cc
  ./IP_SSATaucTaoTikhonovProblem.cc
  ./IP_SSATaucTaoTikhonovProblemLCL.cc
  ./IP_SSAHardavTaoTikhonovProblem.cc
  )

add_library (pisminverse OBJECT ${PISM_INVERSE_SOURCES})
