# Classes for inverse problems.

set(PISM_INVERSE_SOURCES
./InvTaucParameterization.cc
./InvSSAForwardProblem.cc
./InvSSAForwardProblem_dep.cc
./TwoBlockVec.cc
./InvSSATikhonovGN.cc  
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
    ./SSATaucTikhonovProblem.cc
    ./InvSSATikhonovLCL.cc
    )
endif()

add_library (pisminverse ${PISM_INVERSE_SOURCES})
target_link_libraries (pisminverse pismbase pismutil)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})
endif()
