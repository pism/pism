# Classes for inverse problems.
add_library (pisminverse
  ./TaoUtil.cc
  ./InvTaucParameterization.cc
  ./InvSSAForwardProblem.cc
  ./InvSchrodTikhonov.cc
  ./Functional.cc
  ./H1NormFunctional.cc
  )
target_link_libraries (pisminverse pismbase pismutil)

if (Pism_BUILD_EXTRA_EXECS)
  add_executable (taobasicproblem taobasicproblem.cc)
  target_link_libraries (taobasicproblem pisminverse)
endif()

install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})
