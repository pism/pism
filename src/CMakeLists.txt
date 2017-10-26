# Create a link from build/pism to src so that PISM can use #include "pism/..." in its
# source files. This way other code using PISM needs to add only *one* directory to the
# compiler search path. See #166.

execute_process(
  COMMAND
  ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/pism)

# Alternatively: copy headers to the build directory
# file(COPY .
#   DESTINATION ${CMAKE_BINARY_DIR}/pism/
#   FILES_MATCHING
#   PATTERN "*.hh"
#   PATTERN "*.h")

include_directories (${CMAKE_BINARY_DIR})

# Computation grid, I/O, profiling, PISM "component" interface, miscellaneous
# utility functions and classes plus the "enthalpy converter" (which is used
# all over the place).
set(PISM_PARALLEL_IO_FLAGS "")

# Create a list of files making up libpismutil so that we can add to it later:
set(PISMUTIL_SRC
  util/ColumnInterpolation.cc
  util/Context.cc
  util/EnthalpyConverter.cc
  util/FETools.cc
  util/IceGrid.cc
  util/Logger.cc
  util/Mask.cc
  util/MaxTimestep.cc
  util/Component.cc
  util/Config.cc
  util/ConfigInterface.cc
  util/Diagnostic.cc
  util/Time.cc
  util/Time_Calendar.cc
  util/Units.cc
  util/Vars.cc
  util/Profiling.cc
  util/TerminationReason.cc
  util/Timeseries.cc
  util/VariableMetadata.cc
  util/error_handling.cc
  util/iceModelVec.cc
  util/iceModelVec2.cc
  util/iceModelVec2T.cc
  util/iceModelVec2V.cc
  util/iceModelVec3.cc
  util/iceModelVec3Custom.cc
  util/interpolation.cc
  util/io/LocalInterpCtx.cc
  util/io/PIO.cc
  util/io/NC3File.cc
  util/io/NC4File.cc
  util/io/NCFile.cc
  util/io/io_helpers.cc
  util/node_types.cc
  util/options.cc
  util/petscwrappers/DM.cc
  util/petscwrappers/IS.cc
  util/petscwrappers/KSP.cc
  util/petscwrappers/Mat.cc
  util/petscwrappers/PetscInitializer.cc
  util/petscwrappers/SNES.cc
  util/petscwrappers/Tao.cc
  util/petscwrappers/Vec.cc
  util/petscwrappers/VecScatter.cc
  util/petscwrappers/Viewer.cc
  util/pism_const.cc
  util/pism_default_config.cc
  util/pism_options.cc
  util/pism_revision.cc
  util/pism_utilities.cc
  util/projection.cc
  )

# PISM Revision string
set_source_files_properties (util/pism_revision.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION_TAG}\"')

if(Pism_USE_JANSSON)
  list(APPEND PISMUTIL_SRC util/ConfigJSON.cc)
endif()

# Check if NetCDF-4 parallel I/O is enabled. If so, set compiler flags and add a source code file.
if (Pism_USE_PARALLEL_NETCDF4)
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PARALLEL_NETCDF4=1")
  list(APPEND PISMUTIL_SRC util/io/NC4_Par.cc)
else()
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PARALLEL_NETCDF4=0")
endif()

# Check if PnetCDF parallel I/O is enabled. If so, set compiler flags and add a source code file.
if (Pism_USE_PNETCDF)
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PNETCDF=1")
  list(APPEND PISMUTIL_SRC util/io/PNCFile.cc)
else()
  set(PISM_PARALLEL_IO_FLAGS "${PISM_PARALLEL_IO_FLAGS} -DPISM_USE_PNETCDF=0")
endif()

# The list of headers included in PIO.cc depends on parallel I/O settings.
# These flags determine which headers are included.
set_source_files_properties (util/io/PIO.cc
  PROPERTIES COMPILE_FLAGS ${PISM_PARALLEL_IO_FLAGS})

add_library (pismutil OBJECT ${PISMUTIL_SRC})

add_custom_target (touch_pism_revision
  COMMAND ${CMAKE_COMMAND} -E touch pism_revision.cc
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/util
  VERBATIM)

# This library contains PISM code implementing the ice-sheet model itself
# (using other PISM libraries and a good deal of non-trivial code).
add_library (pism
  ${EVERYTRACE_cf_mpi_REFADDR}
  age/AgeColumnSystem.cc
  age/AgeModel.cc
  basalstrength/ConstantYieldStress.cc
  basalstrength/MohrCoulombYieldStress.cc
  basalstrength/YieldStress.cc
  calving/CalvingAtThickness.cc
  calving/CalvingFrontRetreat.cc
  calving/EigenCalving.cc
  calving/FloatKill.cc
  calving/FrontalMelt.cc
  calving/IcebergRemover.cc
  calving/OceanKill.cc
  calving/StressCalving.cc
  calving/connected_components.cc
  calving/remove_narrow_tongues.cc
  calving/vonMisesCalving.cc
  energy/BTU_Full.cc
  energy/BTU_Minimal.cc
  energy/BedThermalUnit.cc
  energy/DummyEnergyModel.cc
  energy/EnergyModel.cc
  energy/EnthalpyModel.cc
  energy/TemperatureModel.cc
  energy/bootstrapping.cc
  energy/enthSystem.cc
  energy/tempSystem.cc
  energy/utilities.cc
  geometry/Geometry.cc
  geometry/GeometryEvolution.cc
  geometry/grounded_cell_fraction.cc
  geometry/part_grid_threshold_thickness.cc
  hydrology/DistributedHydrology.cc
  hydrology/Hydrology.cc
  hydrology/NullTransportHydrology.cc
  hydrology/RoutingHydrology.cc
  hydrology/hydrology_diagnostics.cc
  icemodel/IceModel.cc
  icemodel/calving.cc
  icemodel/diagnostics.cc
  icemodel/diagnostics.cc
  icemodel/diagnostics.hh
  icemodel/energy.cc
  icemodel/flux_balance.hh
  icemodel/fracture_density.cc
  icemodel/initialization.cc
  icemodel/output.cc
  icemodel/output_backup.cc
  icemodel/output_extra.cc
  icemodel/output_save.cc
  icemodel/output_ts.cc
  icemodel/printout.cc
  icemodel/timestepping.cc
  icemodel/utilities.cc
  icemodel/viewers.cc
  stressbalance/timestepping.cc
  util/ColumnSystem.cc
  util/pism_signal.c
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
  earth/PointwiseIsostasy.cc
  earth/BedDef.cc
  earth/LingleClark.cc
  earth/Null.cc
  earth/deformation.cc
  earth/greens.cc
  earth/matlablike.cc
  external/cubature/cubature.c
  )

# Disable compiler warnings for cubature.c (source file not by PISM authors)
set_source_files_properties (earth/cubature.c
  PROPERTIES COMPILE_FLAGS "-w")

set_source_files_properties (util/pism_default_config.cc
  PROPERTIES COMPILE_FLAGS -DPISM_DEFAULT_CONFIG_FILE='\"${Pism_DEFAULT_CONFIG_FILE}\"')

# Verification tests.
add_library (pismverif OBJECT
  verification/TemperatureModel_Verification.cc
  verification/iCMthermo.cc
  verification/iceCompModel.cc
  verification/PSVerification.cc
  verification/BTU_Verification.cc
  verification/tests/exactTestH.c
  verification/tests/exactTestK.c
  verification/tests/exactTestO.c
  verification/tests/exactTestL.cc
  verification/tests/exactTestsABCD.c
  verification/tests/exactTestsFG.cc
  verification/tests/exactTestsIJ.c
  verification/tests/exactTestM.c
  verification/tests/exactTestN.c
  verification/tests/exactTestP.cc
)

# Main executables:
add_executable (pismr pismr.cc)
target_link_libraries (pismr pism)

# Simplified geometry
add_executable (pisms pisms.cc
  icemodel/IceEISModel.cc)
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
  add_executable (btutest energy/btutest.cc)
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
  PATTERN "external" EXCLUDE
  PATTERN "pythonbindings" EXCLUDE
  PATTERN "doc" EXCLUDE
  PATTERN "figs" EXCLUDE
  PATTERN "ssa/tests" EXCLUDE
  PATTERN "verification/tests/fortran" EXCLUDE
  PATTERN "rheology/approximate" EXCLUDE
  PATTERN "tracer" EXCLUDE
  )

add_subdirectory (external/calcalcs)
add_subdirectory (rheology)
add_subdirectory (coupler)
add_subdirectory (stressbalance)
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
