# Subclasses for regional problems.
add_library (pismregional
  ./regional.cc
  )
target_link_libraries (pismregional pismbase)
install (TARGETS pismregional DESTINATION ${Pism_LIB_DIR})



# Outlet glacier 
add_executable (pismo pismo.cc)
target_link_libraries (pismo pismbase pismregional )
install (TARGETS pismo DESTINATION ${Pism_BIN_DIR})
