include_directories (./util)

# Boundary models (surface, atmosphere, ocean).
add_library (pismboundary
  ./atmosphere/PAConstant.cc
  ./atmosphere/PAConstantPIK.cc
  ./atmosphere/PAEismintGreenland.cc
  ./atmosphere/PASeariseGreenland.cc
  ./atmosphere/PAYearlyCycle.cc
  ./atmosphere/PAGivenClimate.cc
  ./atmosphere/PALapseRates.cc
  ./atmosphere/PAdTforcing.cc
  ./atmosphere/PAAnomaly.cc
  ./ocean/POConstant.cc
  ./ocean/POConstantPIK.cc
  ./ocean/POGivenClimate.cc
  ./ocean/POdSLforcing.cc
  ./surface/ebm_driver.cc
  ./surface/localMassBalance.cc
  ./surface/PISMSurface.cc
  ./surface/PSConstantPIK.cc
  ./surface/PSElevation.cc
  ./surface/PSExternal.cc
  ./surface/PSStuffAsAnomaly.cc
  ./surface/PSGivenClimate.cc
  ./surface/PSLapseRates.cc
  ./surface/PSdTforcing.cc
  ./surface/PSSimple.cc
  ./surface/PSTemperatureIndex.cc
  ./surface/PSConstant.cc
  ./surface/PSForceThickness.cc
  ./surface/PSAnomaly.cc
  ./util/PCFactory.cc
  )
target_link_libraries (pismboundary pismutil)

install (TARGETS pismboundary DESTINATION ${Pism_LIB_DIR})
