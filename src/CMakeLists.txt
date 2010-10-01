add_library (pismmisc
  base/LocalInterpCtx.cc
  base/NCVariable.cc
  base/PISMBedSmoother.cc
  base/PISMComponent.cc
  base/PISMIO.cc
  base/PISMProf.cc
  base/PISMVars.cc
  base/Timeseries.cc
  base/enthalpyConverter.cc
  base/flowlaw_factory.cc
  base/flowlaws.cc
  base/grid.cc
  base/iceModelVec.cc
  base/iceModelVec2.cc
  base/iceModelVec2V.cc
  base/iceModelVec3.cc
  base/materials.cc
  base/nc_util.cc
  base/pism_const.cc
  base/pism_default_config.cc
  base/pism_revision.cc
  base/pism_signal.c
  coupler/PAConstant.cc
  coupler/PAForcing.cc
  coupler/PALapseRates.cc
  coupler/PAYearlyCycle.cc
  coupler/PCFactory.cc
  coupler/PISMOcean.cc
  coupler/PISMSurface.cc
  coupler/iceModelVec2T.cc
  coupler/localMassBalance.cc
)
target_link_libraries (pismmisc pismudunits ${Pism_EXTERNAL_LIBS})
add_dependencies(pismmisc pism_config)

add_library (pismearth
  earth/PBLingleClark.cc
  earth/PBPointwiseIsostasy.cc
  earth/PISMBedDef.cc
  earth/cubature.c
  earth/deformation.cc
  earth/greens.cc
  earth/matlablike.cc
)
target_link_libraries (pismearth pismmisc)

set_source_files_properties (base/pism_revision.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION}\"')
set_source_files_properties (base/pism_default_config.cc
  PROPERTIES COMPILE_FLAGS -DPISM_DEFAULT_CONFIG_FILE='\"${Pism_DEFAULT_CONFIG_FILE}\"')
set_source_files_properties (udunits/utlib.c
  PROPERTIES COMPILE_FLAGS
  "-DUT_INSTALL_PATH='\"${CMAKE_INSTALL_PREFIX}/lib/pism/pismudunits.dat\"' -DUT_SOURCE_PATH='\"${CMAKE_SOURCE_DIR}/src/udunits/pismudunits.dat\"'")

include_directories (base udunits)

add_library (pismbase
  base/materials.cc
  base/bedrockOnlySystem.cc
  base/columnSystem.cc
  base/combinedSystem.cc
  base/iceenthOnlySystem.cc
  base/iMIO.cc
  base/iMadaptive.cc
  base/iMbasal.cc
  base/iMbeddef.cc
  base/iMbootstrap.cc
  base/iMdefaults.cc
  base/iMenthalpy.cc
  base/iMgeometry.cc
  base/iMgrainsize.cc
  base/iMinit.cc
  base/iMmatlab.cc
  base/iMoptions.cc
  base/iMreport.cc
  base/iMsia.cc
  base/iMssa.cc
  base/iMtemp.cc
  base/iMtimeseries.cc
  base/iMutil.cc
  base/iMvelocity.cc
  base/iMviewers.cc
  base/iceModel.cc
)
target_link_libraries (pismbase pismearth pismmisc)


add_library (pismverif
  verif/iCMthermo.cc
  verif/iceCompModel.cc
  verif/iceExactSSAModel.cc
  verif/tests/exactTestH.c
  verif/tests/exactTestK.c
  verif/tests/exactTestL.c
  verif/tests/exactTestM.c
  verif/tests/exactTestsABCDE.c
  verif/tests/exactTestsFG.c
  verif/tests/exactTestsIJ.c
)
target_link_libraries (pismverif pismbase)

add_library (pismudunits
  udunits/udalloc.c
  udunits/utlib.c
  udunits/utparse.c
  udunits/utscan.c
)

# Main executables:
add_executable (pismr pismr.cc)

add_executable (pisms pisms.cc 
  eismint/iceEISModel.cc
  eismint/icePSTexModel.cc
  ismip/iceMISMIPModel.cc)

add_executable (pgrn eismint/pgrn.cc eismint/pgrn_atmosphere.cc)

add_executable (pross eismint/pross.cc eismint/iceROSSModel.cc)

# All of the following are linked against pismbase
foreach (EXEC pismr pisms pgrn pross)
  target_link_libraries (${EXEC} pismbase)
endforeach (EXEC)

add_executable (pismv pismv.cc
  verif/iceCompModel.cc
  verif/iCMthermo.cc
  verif/iceExactSSAModel.cc)
target_link_libraries (pismv pismverif pismbase)

add_executable (pclimate pclimate.cc eismint/pgrn_atmosphere.cc)
target_link_libraries (pclimate pismmisc)

add_custom_command (OUTPUT pism_config.nc
  COMMAND ncgen -o ${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_SOURCE_DIR}/pism_config.cdl
  DEPENDS pism_config.cdl
)
add_custom_target (pism_config DEPENDS pism_config.nc)

install (TARGETS
  pismr pisms pismv pgrn pross pclimate ## executables
  pismmisc pismverif pismbase pismearth pismudunits # libraries
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib/pism
  ARCHIVE DESTINATION lib)

install (FILES
  "udunits/pismudunits.dat"
  "${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pism")

# miscellaneous executables needed by software tests
set (UNIT_EXECS bedrough_test ssa_test flowlaw_test)
foreach (EXEC ${UNIT_EXECS})
  add_executable (${EXEC} software_tests/${EXEC}.cc)
  target_link_libraries (${EXEC} pismmisc)
  install (TARGETS ${EXEC} RUNTIME DESTINATION bin)
endforeach (EXEC)

if (Pism_BUILD_EXTRA_EXECS)
  set (SIMPLE_EXECS simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleL)
  foreach (EXEC ${SIMPLE_EXECS})
    add_executable (${EXEC} verif/tests/${EXEC}.c)
    target_link_libraries (${EXEC} pismverif)
  endforeach (EXEC)

  add_executable (tryLCbd earth/tryLCbd.cc)
  target_link_libraries (tryLCbd pismearth)

  add_executable (pismtests software_tests/pismtests.cc software_tests/iMtests.cc)
  target_link_libraries (pismtests pismbase)

  add_executable (flowTable software_tests/flowTable.cc)
  target_link_libraries (flowTable pismmisc)

  install (TARGETS
    ${SIMPLE_EXECS} tryLCbd pismtests flowTable
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib/pism)
endif (Pism_BUILD_EXTRA_EXECS)
