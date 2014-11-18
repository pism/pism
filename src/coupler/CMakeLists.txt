include_directories (./util)

# Boundary models (surface, atmosphere, ocean).
add_library (pismboundary
  ./atmosphere/PAConstantPIK.cc
  ./atmosphere/PASeariseGreenland.cc
  ./atmosphere/PAYearlyCycle.cc
  ./atmosphere/PACosineYearlyCycle.cc
  ./atmosphere/PAGivenClimate.cc
  ./atmosphere/PALapseRates.cc
  ./atmosphere/PA_delta_T.cc
  ./atmosphere/PA_delta_P.cc
  ./atmosphere/PA_frac_P.cc
  ./atmosphere/PA_paleo_precip.cc
  ./atmosphere/PAAnomaly.cc
  ./atmosphere/PAWeatherStation.cc
  ./ocean/PISMOceanModel.cc
  ./ocean/POCache.cc
  ./ocean/POConstantPIK.cc
  ./ocean/POConstant.cc
  ./ocean/POGivenClimate.cc
  ./ocean/POGivenTH.cc
  ./ocean/PO_delta_SL.cc
  ./ocean/PO_delta_T.cc
  ./ocean/PO_delta_SMB.cc
  ./ocean/PO_delta_MBP.cc
  ./surface/localMassBalance.cc
  ./surface/PISMSurface.cc
  ./surface/PSCache.cc
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
  ./surface/PSTemperatureIndex_Old.cc
  ./surface/localMassBalance_old.cc
  ./surface/PS_FEvoR.cc
  ./surface/PS_EISMINTII.cc
  )
target_link_libraries (pismboundary pismutil)

if (${BUILD_SHARED_LIBS})
  install (TARGETS pismboundary DESTINATION ${Pism_LIB_DIR})
endif()
