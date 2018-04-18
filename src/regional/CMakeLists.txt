# Subclasses for regional problems.
add_library (regional OBJECT
  SIAFD_Regional.cc
  SSAFD_Regional.cc
  RegionalDefaultYieldStress.cc
  IceRegionalModel.cc
  IceGrid_Regional.cc
  EnthalpyModel_Regional.cc
  )
