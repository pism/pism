add_library(cubature OBJECT
  cubature.c)
set_target_properties(cubature PROPERTIES COMPILE_FLAGS "-w")
