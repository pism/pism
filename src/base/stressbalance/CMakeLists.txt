include_directories (../../verif/tests/ ../util)

add_library (pismstressbalance
  PISMStressBalance.cc
  PISMStressBalance_diagnostics.cc
  ShallowStressBalance.cc
  SSB_Modifier.cc
)
target_link_libraries (pismstressbalance
  pism_sia
  pism_ssa
  pismutil
  pismflowlaws
  pismrevision)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismstressbalance DESTINATION ${Pism_LIB_DIR})
endif()

add_subdirectory (sia)
add_subdirectory (ssa)
add_subdirectory (blatter)
