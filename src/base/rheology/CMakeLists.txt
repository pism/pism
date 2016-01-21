# Flow laws.
add_library (pismflowlaws
  FlowLaw.cc
  FlowLawFactory.cc
  GPBLD.cc
  GoldsbyKohlstedt.cc
  Hooke.cc
  IsothermalGlen.cc
  PatersonBudd.cc
  PatersonBuddCold.cc
  PatersonBuddWarm.cc
  GPBLD3.cc
  approximate/gpbld.cc
  )
target_link_libraries (pismflowlaws pismutil ${Pism_EXTERNAL_LIBS})

# Install libraries only if BUILD_SHARED_LIBS is ON.
if (${BUILD_SHARED_LIBS})
  install (TARGETS
    pismflowlaws
    LIBRARY DESTINATION ${Pism_LIB_DIR}
    ARCHIVE DESTINATION ${Pism_LIB_DIR})
endif ()
