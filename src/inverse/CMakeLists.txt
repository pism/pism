# Classes for inverse problems.

set( PISM_INVERSE_SOURCES 
./InvTaucParameterization.cc
./InvSSAForwardProblem.cc
./InvSSAForwardProblem_dep.cc
./Functional.cc
./L2NormFunctional.cc
./H1NormFunctional.cc
./MeanSquareFunctional.cc
./TwoBlockVec.cc
./InvSSATikhonovGN.cc  
)

if(TAO_FOUND)
  set(PISM_INVERSE_SOURCES ${PISM_INVERSE_SOURCES}
  ./TaoUtil.cc
  ./InvSSATikhonov.cc
  ./InvSSATikhonovLCL.cc  
  )
endif()

add_library (pisminverse
  ${PISM_INVERSE_SOURCES}
)

target_link_libraries (pisminverse pismbase pismutil)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})
endif()
