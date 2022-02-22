# Create a list of files making up libpismutil so that we can add to it later:
set(PISMUTIL_SRC
  fem/DirichletData.cc
  fem/Element.cc
  fem/ElementIterator.cc
  fem/FEM.cc
  fem/Quadrature.cc
  ColumnInterpolation.cc
  Context.cc
  EnthalpyConverter.cc
  IceGrid.cc
  Logger.cc
  Mask.cc
  MaxTimestep.cc
  Component.cc
  Config.cc
  ConfigInterface.cc
  Diagnostic.cc
  Time.cc
  Units.cc
  Vars.cc
  Profiling.cc
  TerminationReason.cc
  VariableMetadata.cc
  error_handling.cc
  iceModelVec.cc
  iceModelVec2T.cc
  IceModelVec2CellType.cc
  IceModelVec2Int.cc
  iceModelVec2V.cc
  iceModelVec3.cc
  IceModelVec2S.cc
  IceModelVec2Stag.cc
  interpolation.cc
  io/LocalInterpCtx.cc
  io/File.cc
  io/NC_Serial.cc
  io/NC4_Serial.cc
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
  fftw_utilities.cc
  Poisson.cc
  label_components.cc
  connected_components.cc
  ScalarForcing.cc
  )

if(Pism_USE_JANSSON)
  list(APPEND PISMUTIL_SRC ConfigJSON.cc)
endif()

# Check if NetCDF-4 parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PARALLEL_NETCDF4)
  list(APPEND PISMUTIL_SRC io/NC4_Par.cc)
endif()

# Check if PnetCDF parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PNETCDF)
  list(APPEND PISMUTIL_SRC io/PNCFile.cc)
endif()

# Check if PnetCDF parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PIO)
  list(APPEND PISMUTIL_SRC io/ParallelIO.cc)
endif()

add_library (util OBJECT ${PISMUTIL_SRC})
