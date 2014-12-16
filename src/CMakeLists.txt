include_directories (
  base
  base/stressbalance
  base/stressbalance/sia
  base/stressbalance/ssa
  base/stressbalance/blatter
  base/calving
  base/util
  base/util/io
  base/energy
  base/rheology
  base/basalstrength
  base/hydrology
  earth
  coupler
  coupler/atmosphere
  coupler/surface
  coupler/ocean
  coupler/util
  calcalcs
  )

# Computation grid, I/O, profiling, PISM "component" interface, miscellaneous
# utility functions and classes plus the "enthalpy converter" (which is used
# all over the place).
set(PISM_PARALLEL_IO_FLAGS "")

# Create a list of files making up libpismutil so that we can add to it later:
set(PISMUTIL_SRC
  base/util/PetscInitializer.cc
  base/util/error_handling.cc
  base/util/Context.cc
  base/enthalpyConverter.cc
  base/util/IceGrid.cc
  base/util/IceGrid.cc
  base/util/ColumnInterpolation.cc
  base/util/Viewer.cc
  base/util/Mask.cc
  base/util/NCVariable.cc
  base/util/PISMConfig.cc
  base/util/PISMComponent.cc
  base/util/PISMTime.cc
  base/util/PISMTime_Calendar.cc
  base/util/PISMVars.cc
  base/util/Timeseries.cc
  base/util/iceModelVec.cc
  base/util/iceModelVec2.cc
  base/util/iceModelVec2T.cc
  base/util/iceModelVec2V.cc
  base/util/iceModelVec3.cc
  base/util/iceModelVec3Custom.cc
  base/util/io/LocalInterpCtx.cc
  base/util/io/PIO.cc
  base/util/io/PISMNC3File.cc
  base/util/io/PISMNC4File.cc
  base/util/io/PISMNC4_Quilt.cc
  base/util/io/PISMNCFile.cc
  base/util/pism_const.cc
  base/util/pism_default_config.cc
  base/util/pism_options.cc
  base/util/TerminationReason.cc
  base/util/PISMUnits.cc
  base/varcEnthalpyConverter.cc
  )

# Check if NetCDF-4 parallel I/O is enabled. If so, set compiler flags and add a source code file.
if (Pism_USE_PARALLEL_NETCDF4)
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PARALLEL_NETCDF4=1")
  list(APPEND PISMUTIL_SRC base/util/io/PISMNC4_Par.cc)
else()
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PARALLEL_NETCDF4=0")
endif()

# Check if PnetCDF parallel I/O is enabled. If so, set compiler flags and add a source code file.
if (Pism_USE_PNETCDF)
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PNETCDF=1")
  list(APPEND PISMUTIL_SRC base/util/io/PISMPNCFile.cc)
else()
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PNETCDF=0")
endif()

# Check if HDF5-based parallel I/O is enabled. If so, set compiler flags and add a source code file.
if (Pism_USE_PARALLEL_HDF5)
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_HDF5=1")
  list(APPEND PISMUTIL_SRC base/util/io/PISMNC4_HDF5.cc)
else()
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_HDF5=0")
endif()


# The list of headers included in PIO.cc depends on parallel I/O settings.
# These flags determine which headers are included.
set_source_files_properties (base/util/io/PIO.cc
  PROPERTIES COMPILE_FLAGS ${PISM_PARALLEL_IO_FLAGS})

add_library (pismutil ${PISMUTIL_SRC})
target_link_libraries (pismutil pismcalcalcs ${Pism_EXTERNAL_LIBS})
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

# Flow laws.
add_library (pismflowlaws
  base/rheology/flowlaw_factory.cc
  base/rheology/flowlaws.cc
)
target_link_libraries (pismflowlaws pismutil ${Pism_EXTERNAL_LIBS})

# This library contains PISM code implementing the ice-sheet model itself
# (using other PISM libraries and a good deal of non-trivial code).
add_library (pismbase
  base/pism_signal.c
  base/columnSystem.cc
  base/energy/bedrockThermalUnit.cc
  base/energy/enthSystem.cc
  base/energy/tempSystem.cc
  base/iMadaptive.cc
  base/iMage.cc
  base/iMbootstrap.cc
  base/iMcalving.cc
  base/iMenergy.cc
  base/iMenthalpy.cc
  base/iMfractures.cc
  base/iMgeometry.cc
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
  base/basalstrength/PISMMohrCoulombYieldStress.cc
  base/basalstrength/PISMConstantYieldStress.cc
  base/hydrology/PISMHydrology.cc
  base/hydrology/PISMNullTransportHydrology.cc
  base/hydrology/PISMRoutingHydrology.cc
  base/hydrology/PISMDistributedHydrology.cc
  base/hydrology/hydrology_diagnostics.cc
  base/calving/connected_components.cc
  base/calving/PISMIcebergRemover.cc
  base/calving/PISMCalvingAtThickness.cc
  base/calving/PISMFloatKill.cc
  base/calving/PISMOceanKill.cc
  base/calving/PISMEigenCalving.cc
)
target_link_libraries (pismbase pismearth pismboundary pismutil pismstressbalance pismrevision)

# Bed deformation models.
set(PISMEARTH_SRC
  earth/PBPointwiseIsostasy.cc
  earth/PISMBedDef.cc
  earth/PBLingleClark.cc
  earth/deformation.cc
  earth/greens.cc
  earth/cubature.c
  earth/matlablike.cc
  )

add_library(pismearth ${PISMEARTH_SRC})
target_link_libraries (pismearth pismutil)

# Disable compiler warnings for cubature.c (source file not by PISM authors)
set_source_files_properties (earth/cubature.c
  PROPERTIES COMPILE_FLAGS "-w")

set_source_files_properties (base/util/pism_default_config.cc
  PROPERTIES COMPILE_FLAGS -DPISM_DEFAULT_CONFIG_FILE='\"${Pism_DEFAULT_CONFIG_FILE}\"')

# Verification tests.
add_library (pismverif
  verif/iCMthermo.cc
  verif/iceCompModel.cc
  verif/PSVerification.cc
  verif/tests/exactTestH.c
  verif/tests/exactTestK.c
  verif/tests/exactTestO.c
  verif/tests/exactTestL.c
  verif/tests/exactTestsABCDE.c
  verif/tests/exactTestsFG.c
  verif/tests/exactTestsIJ.c
)
target_link_libraries (pismverif pismbase pismstressbalance)

# Main executables:
add_executable (pismr pismr.cc)

# Simplified geometry
add_executable (pisms pisms.cc
  eismint/iceEISModel.cc)

# All of the following are linked against pismbase
foreach (EXEC pismr pisms)
  target_link_libraries (${EXEC} pismbase)
endforeach (EXEC)

# Verification mode "driver".
add_executable (pismv pismv.cc)
target_link_libraries (pismv pismverif pismbase)

add_executable (pismmerge
  pismmerge/pismmerge.cc
  pismmerge/metadata.cc
  pismmerge/variables.cc
  pismmerge/util.cc
  base/util/io/PISMNC4_Serial.cc)
target_link_libraries (pismmerge pismutil pismrevision)

find_program (NCGEN_PROGRAM "ncgen" REQUIRED)
mark_as_advanced(NCGEN_PROGRAM)

add_custom_command (OUTPUT pism_config.nc
  COMMAND ${NCGEN_PROGRAM} -o ${PROJECT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_SOURCE_DIR}/pism_config.cdl
  DEPENDS pism_config.cdl
)
add_custom_target (pism_config DEPENDS pism_config.nc)

# Install libraries only if BUILD_SHARED_LIBS is ON.
if (${BUILD_SHARED_LIBS})
  install (TARGETS
    pismutil pismverif pismbase pismflowlaws pismearth # libraries
    LIBRARY DESTINATION ${Pism_LIB_DIR}
    ARCHIVE DESTINATION ${Pism_LIB_DIR})
endif ()

# Always install executables.
install (TARGETS
  pismr pisms pismv pismmerge ## executables
  RUNTIME DESTINATION ${Pism_BIN_DIR})

install (FILES
  "${PROJECT_BINARY_DIR}/pism_config.nc"
  DESTINATION ${Pism_LIB_DIR})

# miscellaneous executables needed by software tests
add_executable (flowlaw_test
  software_tests/flowlaw_test.cc)
target_link_libraries (flowlaw_test pismutil pismflowlaws)
install (TARGETS flowlaw_test RUNTIME DESTINATION ${Pism_BIN_DIR})

add_executable (bedrough_test
  software_tests/bedrough_test.cc
  base/stressbalance/sia/PISMBedSmoother.cc)
target_link_libraries (bedrough_test pismutil)
install (TARGETS bedrough_test RUNTIME DESTINATION ${Pism_BIN_DIR})

if (Pism_BUILD_EXTRA_EXECS)
  set (EXTRA_EXECS simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleL)
  foreach (EXEC ${EXTRA_EXECS})
    add_executable (${EXEC} verif/tests/${EXEC}.c)
    target_link_libraries (${EXEC} pismverif)
  endforeach (EXEC)

  add_executable (tryLCbd earth/tryLCbd.cc)
  target_link_libraries (tryLCbd pismearth pismrevision)
  list (APPEND EXTRA_EXECS tryLCbd)

  add_executable (btutest base/energy/btutest.cc base/energy/bedrockThermalUnit.cc verif/tests/exactTestK.c)
  target_link_libraries (btutest pismutil pismrevision)
  list (APPEND EXTRA_EXECS btutest)

  install (TARGETS
    ${EXTRA_EXECS}
    RUNTIME DESTINATION ${Pism_BIN_DIR}
    LIBRARY DESTINATION ${Pism_LIB_DIR}
    ARCHIVE DESTINATION ${Pism_LIB_DIR})
endif (Pism_BUILD_EXTRA_EXECS)

add_subdirectory (calcalcs)
add_subdirectory (coupler)
add_subdirectory (base/stressbalance)
add_subdirectory (regional)

if (Pism_DEBUG OR Pism_BUILD_PYTHON_BINDINGS)
  add_subdirectory (pythonbindings)
  # Right now the inverse library is not very useful without Python
  # bindings. We can go back to building it unconditionally if someone
  # uses it in C++ code (please let us know).
  add_subdirectory (inverse)
endif()
