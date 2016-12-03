# PISM-side code used by IceBin (https://github.com/citibeth/icebin).

add_library (pismicebin
  IBIceModel.cc
  IBSurfaceModel.cc
  MassEnergyBudget.cc
  NullTransportHydrology.cc
  VecBundleWriter.cc
  )
target_link_libraries (pismicebin pism)

install (TARGETS pismicebin DESTINATION ${Pism_LIB_DIR})
