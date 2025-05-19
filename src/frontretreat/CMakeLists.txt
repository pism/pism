add_library(frontretreat OBJECT
  FrontRetreat.cc
  PrescribedRetreat.cc
  calving/CalvingAtThickness.cc
  calving/CliffCalvingShear.cc
  calving/CliffCalvingTensile.cc
  calving/EigenCalving.cc
  calving/FloatKill.cc
  calving/HayhurstCalving.cc
  calving/StressCalving.cc
  calving/vonMisesCalving.cc
  util/IcebergRemover.cc
  util/IcebergRemoverFEM.cc
  util/remove_narrow_tongues.cc
  )
