include_directories (../base/stressbalance/)

add_executable (pross_new pross_new.cc)
target_link_libraries(pross_new pismstressbalance)

install (TARGETS pross_new
  DESTINATION ${Pism_BIN_DIR})
