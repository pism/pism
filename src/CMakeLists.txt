include_directories (.)

# Computation grid, I/O, profiling, PISM "component" interface, miscellaneous
# utility functions and classes plus the "enthalpy converter" (which is used
# all over the place).
set(PISM_PARALLEL_IO_FLAGS "")

# Create a list of files making up libpismutil so that we can add to it later:
set(PISMUTIL_SRC
  base/geometry/Geometry.cc
  base/geometry/GeometryEvolution.cc
  base/util/FETools.cc
  base/util/node_types.cc
  base/util/petscwrappers/PetscInitializer.cc
  base/util/petscwrappers/DM.cc
  base/util/petscwrappers/IS.cc
  base/util/petscwrappers/KSP.cc
  base/util/petscwrappers/Mat.cc
  base/util/petscwrappers/SNES.cc
  base/util/petscwrappers/Tao.cc
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
  base/util/EnthalpyConverter.cc
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
  base/util/io/PISMNCFile.cc
  base/util/pism_const.cc
  base/util/pism_utilities.cc
  base/util/Profiling.cc
  base/util/pism_default_config.cc
  base/util/pism_options.cc
  base/util/TerminationReason.cc
  base/util/PISMUnits.cc
  base/util/projection.cc
  base/util/pism_revision.cc
  )

# PISM Revision string
set_source_files_properties (base/util/pism_revision.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION_TAG}\"')

if(Pism_USE_JANSSON)
  list(APPEND PISMUTIL_SRC base/util/ConfigJSON.cc)
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

# The list of headers included in PIO.cc depends on parallel I/O settings.
# These flags determine which headers are included.
set_source_files_properties (base/util/io/PIO.cc
  PROPERTIES COMPILE_FLAGS ${PISM_PARALLEL_IO_FLAGS})

add_library (pismutil OBJECT ${PISMUTIL_SRC})

add_custom_target (touch_pism_revision
  COMMAND ${CMAKE_COMMAND} -E touch pism_revision.cc
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/base/util
  VERBATIM)

# This library contains PISM code implementing the ice-sheet model itself
# (using other PISM libraries and a good deal of non-trivial code).
add_library (pism
  ${EVERYTRACE_cf_mpi_REFADDR}
  base/age/AgeColumnSystem.cc
  base/age/AgeModel.cc
  base/basalstrength/PISMConstantYieldStress.cc
  base/basalstrength/PISMMohrCoulombYieldStress.cc
  base/basalstrength/PISMYieldStress.cc
  base/calving/CalvingAtThickness.cc
  base/calving/CalvingFrontRetreat.cc
  base/calving/EigenCalving.cc
  base/calving/FloatKill.cc
  base/calving/FrontalMelt.cc
  base/calving/IcebergRemover.cc
  base/calving/OceanKill.cc
  base/calving/StressCalving.cc
  base/calving/connected_components.cc
  base/calving/remove_narrow_tongues.cc
  base/calving/vonMisesCalving.cc
  base/energy/BTU_Full.cc
  base/energy/BTU_Minimal.cc
  base/energy/BedThermalUnit.cc
  base/energy/DummyEnergyModel.cc
  base/energy/EnergyModel.cc
  base/energy/EnthalpyModel.cc
  base/energy/TemperatureModel.cc
  base/energy/bootstrapping.cc
  base/energy/enthSystem.cc
  base/energy/tempSystem.cc
  base/energy/utilities.cc
  base/geometry/grounded_cell_fraction.cc
  base/geometry/part_grid_threshold_thickness.cc
  base/hydrology/PISMDistributedHydrology.cc
  base/hydrology/PISMHydrology.cc
  base/hydrology/PISMNullTransportHydrology.cc
  base/hydrology/PISMRoutingHydrology.cc
  base/hydrology/hydrology_diagnostics.cc
  base/icemodel/IceModel.cc
  base/icemodel/calving.cc
  base/icemodel/diagnostics.cc
  base/icemodel/diagnostics.cc
  base/icemodel/diagnostics.hh
  base/icemodel/energy.cc
  base/icemodel/flux_balance.hh
  base/icemodel/fracture_density.cc
  base/icemodel/initialization.cc
  base/icemodel/output.cc
  base/icemodel/output_backup.cc
  base/icemodel/output_extra.cc
  base/icemodel/output_save.cc
  base/icemodel/output_ts.cc
  base/icemodel/printout.cc
  base/icemodel/timestepping.cc
  base/icemodel/utilities.cc
  base/icemodel/viewers.cc
  base/stressbalance/timestepping.cc
  base/util/ColumnSystem.cc
  base/util/pism_signal.c
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

if (Pism_USE_JANSSON)
  if (Pism_BUILD_JANSSON)
    add_dependencies (pism pism_jansson)
  endif()
endif()

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
  verif/TemperatureModel_Verification.cc
  verif/iCMthermo.cc
  verif/iceCompModel.cc
  verif/PSVerification.cc
  verif/BTU_Verification.cc
  verif/tests/exactTestH.c
  verif/tests/exactTestK.c
  verif/tests/exactTestO.c
  verif/tests/exactTestL.cc
  verif/tests/exactTestsABCD.c
  verif/tests/exactTestsFG.cc
  verif/tests/exactTestsIJ.c
  verif/tests/exactTestM.c
  verif/tests/exactTestN.c
  verif/tests/exactTestP.cc
)

# Main executables:
add_executable (pismr pismr.cc)
target_link_libraries (pismr pism)

# Simplified geometry
add_executable (pisms pisms.cc
  base/icemodel/IceEISModel.cc)
target_link_libraries (pisms pism)

# Verification mode "driver".
add_executable (pismv pismv.cc)
target_link_libraries (pismv pism)

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
  pismr pisms pismv ## executables
  RUNTIME DESTINATION ${Pism_BIN_DIR})

install (FILES
  "${PROJECT_BINARY_DIR}/pism_config.nc"
  DESTINATION ${Pism_SHARE_DIR})

# miscellaneous executables needed by software tests
if (Pism_BUILD_EXTRA_EXECS)
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
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/
  DESTINATION include/pism
  FILES_MATCHING
  PATTERN "*.hh"
  PATTERN "*.h"
  PATTERN "calcalcs" EXCLUDE
  PATTERN "pythonbindings" EXCLUDE
  PATTERN "vdt" EXCLUDE
  PATTERN "doc" EXCLUDE
  PATTERN "figs" EXCLUDE
  PATTERN "ssa/tests" EXCLUDE
  PATTERN "verif/tests/fortran" EXCLUDE
  )

add_subdirectory (external/calcalcs)
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
