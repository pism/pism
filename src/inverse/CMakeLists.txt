# Classes for inverse problems.

include_directories ( functional )

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
)

if(Pism_USE_TAO)
  list (APPEND PISM_INVERSE_SOURCES
    ./TaoUtil.cc
    ./IP_SSATaucTaoTikhonovProblem.cc
    ./IP_SSATaucTaoTikhonovProblemLCL.cc
    )
endif()

add_library (pisminverse ${PISM_INVERSE_SOURCES})
target_link_libraries (pisminverse pismbase pismutil)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})
endif()
