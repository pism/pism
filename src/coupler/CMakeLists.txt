include_directories (./util)

# Boundary models (surface, atmosphere, ocean).
add_library (pismboundary
  ./atmosphere/PAConstantPIK.cc
  ./atmosphere/PAEismintGreenland.cc
  ./atmosphere/PASeariseGreenland.cc
  ./atmosphere/PAYearlyCycle.cc
  ./atmosphere/PACosineYearlyCycle.cc
  ./atmosphere/PAGivenClimate.cc
  ./atmosphere/PALapseRates.cc
  ./atmosphere/PA_delta_T.cc
  ./atmosphere/PA_delta_P.cc
  ./atmosphere/PAAnomaly.cc
  ./ocean/POConstantPIK.cc
  ./ocean/POConstant.cc
  ./ocean/POGivenClimate.cc
  ./ocean/PO_delta_SL.cc
  ./ocean/PO_delta_T.cc
  ./ocean/PO_delta_SMB.cc
  ./surface/localMassBalance.cc
  ./surface/PISMSurface.cc
  ./surface/PSConstantPIK.cc
  ./surface/PSElevation.cc
  ./surface/PSStuffAsAnomaly.cc
  ./surface/PSGivenClimate.cc
  ./surface/PSLapseRates.cc
  ./surface/PS_delta_T.cc
  ./surface/PSSimple.cc
  ./surface/PSTemperatureIndex.cc
  ./surface/PSForceThickness.cc
  ./surface/PSAnomaly.cc
  ./util/PCFactory.cc
  )
target_link_libraries (pismboundary pismutil)

install (TARGETS pismboundary DESTINATION ${Pism_LIB_DIR})
