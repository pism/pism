# Bed deformation models.
add_library(earth OBJECT
  Given.cc
  PointwiseIsostasy.cc
  BedDef.cc
  LingleClark.cc
  Null.cc
  LingleClarkSerial.cc
  greens.cc
  matlablike.cc
  )
