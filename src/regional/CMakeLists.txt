# Subclasses for regional problems.
add_library (pismregional
  SIAFD_Regional.cc
  SSAFD_Regional.cc
  RegionalDefaultYieldStress.cc
  IceRegionalModel.cc
  )
target_link_libraries (pismregional pismbase)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismregional DESTINATION ${Pism_LIB_DIR})
endif()


# Outlet glacier 
add_executable (pismo pismo.cc)
target_link_libraries (pismo pismbase pismregional)
install (TARGETS pismo DESTINATION ${Pism_BIN_DIR})
