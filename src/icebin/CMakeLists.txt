# PISM-side code used by IceBin (https://github.com/citibeth/icebin).

add_library (pismicebin
  ${EVERYTRACE_cf_mpi_REFADDR}
  IBIceModel.cc
  IBSurfaceModel.cc
  MassEnergyBudget.cc
  NullTransportHydrology.cc
  VecBundleWriter.cc
  )
target_link_libraries (pismicebin pism)

install (TARGETS pismicebin DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Generate the pkg-config file
configure_file(pismicebin.pc.in ${CMAKE_CURRENT_BINARY_DIR}/pismicebin.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pismicebin.pc DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
