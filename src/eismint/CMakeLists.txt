include_directories (../base/stressbalance/)

add_executable (pross pross.cc)
target_link_libraries(pross pismstressbalance)

install (TARGETS pross
  DESTINATION ${Pism_BIN_DIR})
