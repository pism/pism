# Subclasses for regional problems.
add_library (pismregional OBJECT
  SIAFD_Regional.cc
  SSAFD_Regional.cc
  RegionalDefaultYieldStress.cc
  IceRegionalModel.cc
  IceGrid_Regional.cc
  )

# Outlet glacier 
add_executable (pismo pismo.cc)
target_link_libraries (pismo pism)
install (TARGETS pismo DESTINATION ${Pism_BIN_DIR})
