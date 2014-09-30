add_library(pismfevor
  PISMFEvoR.cc)
target_link_libraries(pismfevor pismutil pismstressbalance)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismfevor DESTINATION ${Pism_LIB_DIR})
endif()

