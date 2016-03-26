include_directories (.)

# Computation grid, I/O, profiling, PISM "component" interface, miscellaneous
# utility functions and classes plus the "enthalpy converter" (which is used
# all over the place).
set(PISM_PARALLEL_IO_FLAGS "")

# Create a list of files making up libpismutil so that we can add to it later:
set(PISMUTIL_SRC
  base/util/petscwrappers/PetscInitializer.cc
  base/util/petscwrappers/DM.cc
  base/util/petscwrappers/IS.cc
  base/util/petscwrappers/KSP.cc
  base/util/petscwrappers/Mat.cc
  base/util/petscwrappers/SNES.cc
  base/util/petscwrappers/Vec.cc
  base/util/petscwrappers/VecScatter.cc
  base/util/petscwrappers/Viewer.cc
  base/util/Logger.cc
  base/util/error_handling.cc
  base/util/options.cc
  base/util/interpolation.cc
  base/util/Context.cc
  base/util/MaxTimestep.cc
  base/util/PISMDiagnostic.cc
  base/enthalpyConverter.cc
  base/util/IceGrid.cc
  base/util/ColumnInterpolation.cc
  base/util/Mask.cc
  base/util/VariableMetadata.cc
  base/util/PISMConfig.cc
  base/util/PISMConfigInterface.cc
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
  base/util/io/io_helpers.cc
  base/util/io/LocalInterpCtx.cc
  base/util/io/PIO.cc
  base/util/io/PISMNC3File.cc
  base/util/io/PISMNC4File.cc
  base/util/io/PISMNC4_Quilt.cc
  base/util/io/PISMNCFile.cc
  base/util/pism_const.cc
  base/util/pism_utilities.cc
  base/util/Profiling.cc
  base/util/pism_default_config.cc
  base/util/pism_options.cc
  base/util/TerminationReason.cc
  base/util/PISMUnits.cc
  base/util/projection.cc
  base/pism_revision.cc
  )

# PISM Revision string
set_source_files_properties (base/pism_revision.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION_TAG}\"')

if(Pism_USE_TAO)
  list(APPEND PISMUTIL_SRC base/util/petscwrappers/Tao.cc)
endif()

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

add_library (pismutil OBJECT ${PISMUTIL_SRC})

add_custom_target (touch_pism_revision
  COMMAND ${CMAKE_COMMAND} -E touch pism_revision.cc
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/base
  VERBATIM)

# This library contains PISM code implementing the ice-sheet model itself
# (using other PISM libraries and a good deal of non-trivial code).
add_library (pism
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
  base/grounded_cell_fraction.cc
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
  base/basalstrength/PISMYieldStress.cc
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
  $<TARGET_OBJECTS:pismflowlaws>
  $<TARGET_OBJECTS:pismstressbalance>
  $<TARGET_OBJECTS:pismcalcalcs>
  $<TARGET_OBJECTS:pismutil>
  $<TARGET_OBJECTS:pismearth>
  $<TARGET_OBJECTS:pismverif>
  $<TARGET_OBJECTS:pismboundary>
  $<TARGET_OBJECTS:pisminverse>
  $<TARGET_OBJECTS:pismregional>
)
target_link_libraries (pism ${Pism_EXTERNAL_LIBS})
add_dependencies (pism pism_config)
add_dependencies (pism touch_pism_revision)

# Bed deformation models.
add_library(pismearth OBJECT
  earth/PBPointwiseIsostasy.cc
  earth/PISMBedDef.cc
  earth/PBLingleClark.cc
  earth/PBNull.cc
  earth/deformation.cc
  earth/greens.cc
  earth/cubature.c
  earth/matlablike.cc
  )

# Disable compiler warnings for cubature.c (source file not by PISM authors)
set_source_files_properties (earth/cubature.c
  PROPERTIES COMPILE_FLAGS "-w")

set_source_files_properties (base/util/pism_default_config.cc
  PROPERTIES COMPILE_FLAGS -DPISM_DEFAULT_CONFIG_FILE='\"${Pism_DEFAULT_CONFIG_FILE}\"')

# Verification tests.
add_library (pismverif OBJECT
  verif/iCMthermo.cc
  verif/iceCompModel.cc
  verif/PSVerification.cc
  verif/tests/exactTestH.c
  verif/tests/exactTestK.c
  verif/tests/exactTestO.c
  verif/tests/exactTestL.c
  verif/tests/exactTestsABCD.c
  verif/tests/exactTestsFG.c
  verif/tests/exactTestsIJ.c
)

# Main executables:
add_executable (pismr pismr.cc)
target_link_libraries (pismr pism)

# Simplified geometry
add_executable (pisms pisms.cc
  eismint/iceEISModel.cc)
target_link_libraries (pisms pism)

# Verification mode "driver".
add_executable (pismv pismv.cc)
target_link_libraries (pismv pism)

add_executable (pismmerge
  pismmerge/pismmerge.cc
  pismmerge/metadata.cc
  pismmerge/variables.cc
  pismmerge/util.cc
  base/util/io/PISMNC4_Serial.cc)
target_link_libraries (pismmerge pism)

find_program (NCGEN_PROGRAM "ncgen" REQUIRED)
mark_as_advanced(NCGEN_PROGRAM)

add_custom_command (OUTPUT pism_config.nc
  COMMAND ${NCGEN_PROGRAM} -o ${PROJECT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_SOURCE_DIR}/pism_config.cdl
  DEPENDS pism_config.cdl
)
add_custom_target (pism_config DEPENDS pism_config.nc)

# Install the library
install (
  TARGETS pism
  LIBRARY DESTINATION ${Pism_LIB_DIR}
  ARCHIVE DESTINATION ${Pism_LIB_DIR})

# Install executables.
install (TARGETS
  pismr pisms pismv pismmerge ## executables
  RUNTIME DESTINATION ${Pism_BIN_DIR})

install (FILES
  "${PROJECT_BINARY_DIR}/pism_config.nc"
  DESTINATION ${Pism_SHARE_DIR})

# miscellaneous executables needed by software tests
add_executable (flowlaw_test software_tests/flowlaw_test.cc)
target_link_libraries (flowlaw_test pism)
install (TARGETS flowlaw_test RUNTIME DESTINATION ${Pism_BIN_DIR})

add_executable (bedrough_test software_tests/bedrough_test.cc)
target_link_libraries (bedrough_test pism)
install (TARGETS bedrough_test RUNTIME DESTINATION ${Pism_BIN_DIR})

if (Pism_BUILD_EXTRA_EXECS)
  set (EXTRA_EXECS simpleABCD simpleFG simpleH simpleI simpleJ simpleL)
  foreach (EXEC ${EXTRA_EXECS})
    add_executable (${EXEC} verif/tests/${EXEC}.c)
    target_link_libraries (${EXEC} pism)
  endforeach (EXEC)

  add_executable (tryLCbd earth/tryLCbd.cc)
  target_link_libraries (tryLCbd pism)
  list (APPEND EXTRA_EXECS tryLCbd)

  add_executable (btutest base/energy/btutest.cc)
  target_link_libraries (btutest pism)
  list (APPEND EXTRA_EXECS btutest)

  install (TARGETS
    ${EXTRA_EXECS}
    RUNTIME DESTINATION ${Pism_BIN_DIR}
    LIBRARY DESTINATION ${Pism_LIB_DIR}
    ARCHIVE DESTINATION ${Pism_LIB_DIR})
endif (Pism_BUILD_EXTRA_EXECS)


# ------------ Header Files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/ DESTINATION include
          FILES_MATCHING PATTERN "*.hh")

add_subdirectory (calcalcs)
add_subdirectory (base/rheology)
add_subdirectory (coupler)
add_subdirectory (base/stressbalance)
add_subdirectory (regional)

# Right now the inverse library is not very useful without Python
# bindings, but it's easier to build it unconditionally.
add_subdirectory (inverse)

if (Pism_BUILD_ICEBIN)
  add_subdirectory (icebin)
endif()

if (Pism_DEBUG OR Pism_BUILD_PYTHON_BINDINGS)
  add_subdirectory (pythonbindings)
endif()
