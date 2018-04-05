# Bed deformation models.
add_library(earth OBJECT
  PointwiseIsostasy.cc
  BedDef.cc
  LingleClark.cc
  Null.cc
  deformation.cc
  greens.cc
  matlablike.cc
  )
