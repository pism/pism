# Create a list of files making up libpismutil so that we can add to it later:
set(PISMUTIL_SRC
  ColumnInterpolation.cc
  Context.cc
  EnthalpyConverter.cc
  FETools.cc
  IceGrid.cc
  Logger.cc
  Mask.cc
  MaxTimestep.cc
  Component.cc
  Config.cc
  ConfigInterface.cc
  Diagnostic.cc
  Time.cc
  Time_Calendar.cc
  Units.cc
  Vars.cc
  Profiling.cc
  TerminationReason.cc
  Timeseries.cc
  VariableMetadata.cc
  error_handling.cc
  iceModelVec.cc
  iceModelVec2.cc
  iceModelVec2T.cc
  iceModelVec2V.cc
  iceModelVec3.cc
  iceModelVec3Custom.cc
  interpolation.cc
  io/LocalInterpCtx.cc
  io/PIO.cc
  io/NC3File.cc
  io/NC4File.cc
  io/NCFile.cc
  io/io_helpers.cc
  node_types.cc
  options.cc
  petscwrappers/DM.cc
  petscwrappers/IS.cc
  petscwrappers/KSP.cc
  petscwrappers/Mat.cc
  petscwrappers/PetscInitializer.cc
  petscwrappers/SNES.cc
  petscwrappers/Tao.cc
  petscwrappers/Vec.cc
  petscwrappers/VecScatter.cc
  petscwrappers/Viewer.cc
  pism_options.cc
  pism_utilities.cc
  projection.cc
  )

if(Pism_USE_JANSSON)
  list(APPEND PISMUTIL_SRC ConfigJSON.cc)
endif()

# Check if NetCDF-4 parallel I/O is enabled. If so, set compiler flags and add a source code file.
if (Pism_USE_PARALLEL_NETCDF4)
  add_definitions(-DPISM_USE_PARALLEL_NETCDF4=1)
  list(APPEND PISMUTIL_SRC io/NC4_Par.cc)
else()
  add_definitions(-DPISM_USE_PARALLEL_NETCDF4=0)
endif()

# Check if PnetCDF parallel I/O is enabled. If so, set compiler flags and add a source code file.
if (Pism_USE_PNETCDF)
  add_definitions(-DPISM_USE_PNETCDF=1)
  list(APPEND PISMUTIL_SRC io/PNCFile.cc)
else()
  add_definitions(-DPISM_USE_PNETCDF=0)
endif()

add_custom_target (touch_pism_utilities
  COMMAND ${CMAKE_COMMAND} -E touch pism_utilities.cc
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM)

add_library (util OBJECT ${PISMUTIL_SRC})
add_dependencies (util touch_pism_utilities)
