# Flow laws.
add_library (flowlaws OBJECT
  FlowLaw.cc
  FlowLawFactory.cc
  GPBLD.cc
  GoldsbyKohlstedt.cc
  Hooke.cc
  IsothermalGlen.cc
  PatersonBudd.cc
  PatersonBuddCold.cc
  PatersonBuddWarm.cc
  GPBLD3.cc
  approximate/gpbld_n.cc
  )

add_executable(gpbld3_test
  gpbld3_test.cc)
target_link_libraries (gpbld3_test pism)
