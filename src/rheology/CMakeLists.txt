# Flow laws.
add_library (pismflowlaws OBJECT
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

add_executable(gpbld_speed_test
  gpbld_speed_test.cc)
target_link_libraries (gpbld_speed_test pism)
