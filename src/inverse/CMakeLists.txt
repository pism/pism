# Classes for inverse problems.

set(PISM_INVERSE_SOURCES
./IPTaucParameterization.cc
./IP_SSATaucForwardProblem.cc
./IP_SSATaucForwardProblem_dep.cc
./IPTwoBlockVec.cc
./IP_SSATaucTikhonovGNSolver.cc  
./functional/Functional.cc
./functional/L2NormFunctional.cc
./functional/H1NormFunctional.cc
./functional/GroundedIceH1NormFunctional.cc
./functional/LogRatioFunctional.cc
./functional/LogRelativeFunctional.cc
./functional/TotalVariationFunctional.cc
./functional/MeanSquareFunctional.cc
)

if(Pism_USE_TAO)
  list (APPEND PISM_INVERSE_SOURCES
    ./TaoUtil.cc
    ./IP_SSATaucTikhonovProblem.cc
    ./IP_SSATaucTikhonovProblemLCL.cc
    )
endif()

add_library (pisminverse ${PISM_INVERSE_SOURCES})
target_link_libraries (pisminverse pismbase pismutil)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})
endif()
