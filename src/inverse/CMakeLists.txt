# Classes for inverse problems.
add_library (pisminverse
  ./InvTaucParameterization.cc
  ./InvSSAForwardProblem.cc
  )
target_link_libraries (pisminverse pismbase pismutil)

install (TARGETS pisminverse DESTINATION ${Pism_LIB_DIR})
