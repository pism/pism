# PISM-side code used by IceBin (https://github.com/citibeth/icebin).

add_library (pismicebin
  ${EVERYTRACE_cf_mpi_REFADDR}
  IBIceModel.cc
  IBSurfaceModel.cc
  MassEnergyBudget.cc
  VecBundleWriter.cc
  )
target_link_libraries (pismicebin pism)

install (TARGETS pismicebin DESTINATION ${Pism_LIB_DIR})
