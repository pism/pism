# Boundary models (surface, atmosphere, ocean).
add_library (pismicebin
  IBIceModel.cc
  IBSurfaceModel.cc
  MassEnergyBudget.cc
  NullTransportHydrology.cc
  VecBundleWriter.cc
  )
target_link_libraries (pismicebin pismutil)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismicebin DESTINATION ${Pism_LIB_DIR})
endif()
