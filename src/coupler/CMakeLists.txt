# Boundary models (surface, atmosphere, ocean).
add_library (pismboundary
  ./atmosphere/PAConstant.cc
  ./atmosphere/PAConstantPIK.cc
  ./atmosphere/PAEismintGreenland.cc
  ./atmosphere/PAForcing.cc
  ./atmosphere/PASeariseGreenland.cc
  ./atmosphere/PAYearlyCycle.cc
  ./ocean/PISMOcean.cc
  ./ocean/POConstantPIK.cc
  ./ocean/PODirectForcing.cc
  ./surface/ebm_driver.cc
  ./surface/localMassBalance.cc
  ./surface/PISMSurface.cc
  ./surface/PSConstantPIK.cc
  ./surface/PSDirectAnomalies.cc
  ./surface/PSElevation.cc
  ./surface/PSExternal.cc
  ./PASDirectForcing.cc
  ./PASLapseRates.cc
  ./PCFactory.cc
  ./PScalarForcing.cc
  )
target_link_libraries (pismboundary pismutil)

install (TARGETS pismboundary DESTINATION ${Pism_LIB_DIR})
