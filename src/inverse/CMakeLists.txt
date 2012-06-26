# Classes for inverse problems.

set( PISM_INVERSE_SOURCES 
./InvTaucParameterization.cc
./InvSSAForwardProblem.cc
./Functional.cc
./L2NormFunctional.cc
./H1NormFunctional.cc
./MeanSquareFunctional.cc
./TwoBlockVec.cc
)

if(TAO_FOUND)
  set(PISM_INVERSE_SOURCES ${PISM_INVERSE_SOURCES}
  ./TaoUtil.cc
  ./InvSchrodTikhonov.cc
  ./InvSSATikhonov.cc
  ./InvSSABasicTikhonov.cc
  ./InvSSA_LCLTikhonov.cc  
  )
endif()

add_library (pisminverse
  ${PISM_INVERSE_SOURCES}
)

target_link_libraries (pisminverse pismbase pismutil)

if(TAO_FOUND)
  if (Pism_BUILD_EXTRA_EXECS)
    add_executable (taobasicproblem taobasicproblem.cc)
    target_link_libraries (taobasicproblem pisminverse)
  endif()
endif()

install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})
