include_directories(${FEvoR_INCLUDES} ${CGAL_INCLUDE_DIR})

add_library(pismfevor
  PISMFEvoR.cc)

target_link_libraries(pismfevor pismutil pismstressbalance ${FEvoR_LIB} ${CGAL_LIB})

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismfevor DESTINATION ${Pism_LIB_DIR})
endif()
