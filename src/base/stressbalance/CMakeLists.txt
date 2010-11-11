add_library (pismstressbalance
  PISMStressBalance.cc
  ShallowStressBalance.cc
  SSB_Modifier.cc
  SIAFD.cc
  SSAFD.cc
)
target_link_libraries (pismstressbalance pismmisc)

add_executable (ssafd_test ssafd_test.cc)
target_link_libraries (ssafd_test pismstressbalance)