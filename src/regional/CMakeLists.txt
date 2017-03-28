# Subclasses for regional problems.
add_library (pismregional OBJECT
  SIAFD_Regional.cc
  SSAFD_Regional.cc
  RegionalDefaultYieldStress.cc
  IceRegionalModel.cc
  IceGrid_Regional.cc
  EnthalpyModel_Regional.cc
  )
