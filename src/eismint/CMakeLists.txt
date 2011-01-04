add_executable (pross pross.cc ../base/materials.cc)
target_link_libraries(pross pismstressbalance)

install (TARGETS pross
  DESTINATION ${Pism_BIN_DIR})
