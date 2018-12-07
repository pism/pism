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
  )
