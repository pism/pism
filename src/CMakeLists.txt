include_directories (
  base
  base/stressbalance
  base/util
  base/energy
  base/rheology
  base/basal_strength
  earth
  coupler
  coupler/atmosphere
  coupler/surface
  coupler/ocean
  udunits
  )

# Computation grid, I/O, profiling, PISM "component" interface, miscellaneous
# utility functions and classes plus the "enthalpy converter" (which is used
# all over the place).
add_library (pismutil
  base/util/LocalInterpCtx.cc
  base/util/NCVariable.cc
  base/util/PISMComponent.cc
  base/util/PISMIO.cc
  base/util/PISMProf.cc
  base/util/PISMVars.cc
  base/util/Timeseries.cc
  base/util/IceGrid.cc
  base/util/iceModelVec.cc
  base/util/iceModelVec2.cc
  base/util/iceModelVec2V.cc
  base/util/iceModelVec3.cc
  base/util/iceModelVec2T.cc
  base/util/NCTool.cc
  base/util/pism_const.cc
  base/util/pism_default_config.cc
  base/util/PISMTime.cc
  base/enthalpyConverter.cc
  base/varcEnthalpyConverter.cc
  )
target_link_libraries (pismutil pismudunits ${Pism_EXTERNAL_LIBS})
add_dependencies (pismutil pism_config)

add_custom_target (touch_pism_revision
  COMMAND ${CMAKE_COMMAND} -E touch pism_revision.cc
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/base
  VERBATIM)

# PISM Revision string
set_source_files_properties (base/pism_revision.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION_TAG}\"')
add_library (pismrevision STATIC base/pism_revision.cc)
add_dependencies (pismrevision touch_pism_revision)

# Add a preprocessor flag if NetCDF-4 output was requested.
if (Pism_NETCDF4_OUTPUT)
  message (STATUS "PISM will write NetCDF-4 files.")
  set_source_files_properties (base/util/NCTool.cc
    PROPERTIES COMPILE_FLAGS -DPISM_WRITE_NETCDF4=1)
else ()
  set_source_files_properties (base/util/NCTool.cc
    PROPERTIES COMPILE_FLAGS -DPISM_WRITE_NETCDF4=0)
endif ()

if (Pism_NETCDF4_COMPRESS)
  message (STATUS "PISM will use NetCDF-4 variable compression.")
  set (NCVariable_flags "-DPISM_WRITE_COMPRESSED_NETCDF4=1")
else ()
  set (NCVariable_flags "-DPISM_WRITE_COMPRESSED_NETCDF4=0")
endif ()

# Set variable storage order in PISM output files.
if     (Pism_VARIABLE_ORDER STREQUAL "xyz")
  message (STATUS "Using the t,x,y,z variable order (fastest).")
  set (fast_write "1")
elseif (Pism_VARIABLE_ORDER STREQUAL "yxz")
  message (STATUS "Using the t,y,x,z variable order.")
  set (fast_write "2")
elseif (Pism_VARIABLE_ORDER STREQUAL "zyx")
  message (STATUS "Using the t,z,y,x variable order (slowest).")
  set (fast_write "3")
else ()
  message (FATAL_ERROR "Pism_VARIABLE_ORDER has to be set to one of 'xyz','yxz','zyx'.")
endif ()

set (NCVariable_flags "${NCVariable_flags} -DPISM_VARIABLE_ORDER=${fast_write}")

set_source_files_properties (base/util/NCVariable.cc
  PROPERTIES COMPILE_FLAGS ${NCVariable_flags})


# Flow laws.
add_library (pismflowlaws
  base/rheology/flowlaw_factory.cc
  base/rheology/flowlaws.cc
)
target_link_libraries (pismflowlaws pismutil pismudunits ${Pism_EXTERNAL_LIBS})

# This library contains PISM code implementing the ice-sheet model itself
# (using other PISM libraries and a good deal of non-trivial code).
add_library (pismbase
  base/pism_signal.c
  base/columnSystem.cc
  base/energy/bedrockThermalUnit.cc
  base/energy/enthSystem.cc
  base/energy/varenthSystem.cc
  base/energy/tempSystem.cc
  base/iMadaptive.cc
  base/iMage.cc
  base/iMbeddef.cc
  base/iMbootstrap.cc
  base/iMcalving.cc
  base/iMenergy.cc
  base/iMenthalpy.cc
  base/iMgeometry.cc
  base/iMhydrology.cc
  base/iMicebergs.cc
  base/iMinit.cc
  base/iMIO.cc
  base/iMoptions.cc
  base/iMpartgrid.cc
  base/iMreport.cc
  base/iMtemp.cc
  base/iMtimeseries.cc
  base/iMutil.cc
  base/iMcell_areas.cc
  base/iMviewers.cc
  base/iceModel.cc
  base/iceModel_diagnostics.cc
  base/basal_strength/PISMYieldStress.cc
  base/basal_strength/PISMConstantYieldStress.cc
  base/basal_strength/basal_resistance.cc
)
target_link_libraries (pismbase pismearth pismboundary pismutil pismstressbalance pismrevision)

# Bed deformation models.
add_library (pismearth
  earth/PBLingleClark.cc
  earth/PBPointwiseIsostasy.cc
  earth/PISMBedDef.cc
  earth/cubature.c
  earth/deformation.cc
  earth/greens.cc
  earth/matlablike.cc
)
target_link_libraries (pismearth pismutil)

# Disable compiler warnings for cubature.c (source file not by PISM authors)
set_source_files_properties (earth/cubature.c
  PROPERTIES COMPILE_FLAGS "-w")

set_source_files_properties (base/util/pism_default_config.cc
  PROPERTIES COMPILE_FLAGS -DPISM_DEFAULT_CONFIG_FILE='\"${Pism_DEFAULT_CONFIG_FILE}\"')
set_source_files_properties (udunits/utlib.c
  PROPERTIES COMPILE_FLAGS
  "-DUT_INSTALL_PATH='\"${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}/pismudunits.dat\"' -DUT_SOURCE_PATH='\"${CMAKE_SOURCE_DIR}/src/udunits/pismudunits.dat\"'")


# Verification tests.
add_library (pismverif
  verif/iCMthermo.cc
  verif/iceCompModel.cc
  verif/tests/exactTestH.c
  verif/tests/exactTestK.c
  verif/tests/exactTestO.c
  verif/tests/exactTestL.c
  verif/tests/exactTestsABCDE.c
  verif/tests/exactTestsFG.c
  verif/tests/exactTestsIJ.c
)
target_link_libraries (pismverif pismbase pismstressbalance)

# A custom version of UDUNITS (version 1).
add_library (pismudunits
  udunits/udalloc.c
  udunits/utlib.c
  udunits/utparse.c
  udunits/utscan.c
)

# Disable compiler warnings for UDUNITS (sources are not by PISM authors)
set_target_properties(pismudunits
  PROPERTIES COMPILE_FLAGS "-w")

# Main executables:
add_executable (pismr pismr.cc)

# Outlet glacier 
add_executable (pismo regional/pismo.cc)

# Simplified geometry
add_executable (pisms pisms.cc 
  eismint/iceEISModel.cc
  eismint/icePSTexModel.cc
  ismip/iceMISMIPModel.cc)

add_executable (pismebm coupler/surface/ebm_driver.cc pismebm.cc)

# All of the following are linked against pismbase
foreach (EXEC pismr pismo pisms pismebm)
  target_link_libraries (${EXEC} pismbase)
endforeach (EXEC)

# Verification mode "driver".
add_executable (pismv pismv.cc)
target_link_libraries (pismv pismverif pismbase)

# Testing boundary models without ice dynamics.
add_executable (pclimate pclimate.cc)
target_link_libraries (pclimate pismboundary pismrevision)

find_program (NCGEN_PROGRAM "ncgen")
mark_as_advanced(NCGEN_PROGRAM)

add_custom_command (OUTPUT pism_config.nc
  COMMAND ${NCGEN_PROGRAM} -o ${PROJECT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_SOURCE_DIR}/pism_config.cdl
  DEPENDS pism_config.cdl
)
add_custom_target (pism_config DEPENDS pism_config.nc)

install (TARGETS
  pismr pismo pisms pismv pclimate pismebm ## executables
  pismutil pismverif pismbase pismflowlaws pismearth pismudunits # libraries
  RUNTIME DESTINATION ${Pism_BIN_DIR}
  LIBRARY DESTINATION ${Pism_LIB_DIR}
  ARCHIVE DESTINATION ${Pism_LIB_DIR})

install (FILES
  "udunits/pismudunits.dat"
  "${PROJECT_BINARY_DIR}/pism_config.nc"
  DESTINATION ${Pism_LIB_DIR})

if (NOT Pism_BUILD_DEBIAN_PACKAGE)
  # miscellaneous executables needed by software tests
  add_executable (flowlaw_test
    software_tests/flowlaw_test.cc)
  target_link_libraries (flowlaw_test pismutil pismflowlaws)
  install (TARGETS flowlaw_test RUNTIME DESTINATION ${Pism_BIN_DIR})
  
  add_executable (bedrough_test
    software_tests/bedrough_test.cc
    base/stressbalance/PISMBedSmoother.cc)
  target_link_libraries (bedrough_test pismutil)
  install (TARGETS bedrough_test RUNTIME DESTINATION ${Pism_BIN_DIR})
endif ()

if (Pism_BUILD_EXTRA_EXECS)
  set (SIMPLE_EXECS simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleL)
  foreach (EXEC ${SIMPLE_EXECS})
    add_executable (${EXEC} verif/tests/${EXEC}.c)
    target_link_libraries (${EXEC} pismverif)
  endforeach (EXEC)

  add_executable (tryLCbd earth/tryLCbd.cc)
  target_link_libraries (tryLCbd pismearth pismrevision)

  add_executable (pismtests software_tests/pismtests.cc software_tests/iMtests.cc)
  target_link_libraries (pismtests pismbase)

  add_executable (flowTable software_tests/flowTable.cc base/rheology/flowlaws.cc)
  target_link_libraries (flowTable pismutil pismrevision)

  add_executable (btutest base/energy/btutest.cc base/energy/bedrockThermalUnit.cc verif/tests/exactTestK.c)
  target_link_libraries (btutest pismutil pismrevision)

  install (TARGETS
    ${SIMPLE_EXECS} tryLCbd pismtests flowTable btutest
    RUNTIME DESTINATION ${Pism_BIN_DIR}
    LIBRARY DESTINATION ${Pism_LIB_DIR}
    ARCHIVE DESTINATION ${Pism_LIB_DIR})
endif (Pism_BUILD_EXTRA_EXECS)

add_subdirectory (coupler)
add_subdirectory (base/stressbalance)
add_subdirectory (eismint)
add_subdirectory (pythonbindings)
