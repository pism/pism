add_library(frontretreat OBJECT
  FrontRetreat.cc
  FrontalMelt.cc
  calving/CalvingAtThickness.cc
  calving/EigenCalving.cc
  calving/FloatKill.cc
  calving/OceanKill.cc
  calving/StressCalving.cc
  calving/vonMisesCalving.cc
  util/IcebergRemover.cc
  util/connected_components.cc
  util/remove_narrow_tongues.cc
  )
