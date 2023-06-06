add_library (util OBJECT
  ColumnInterpolation.cc
  Component.cc
  Config.cc
  ConfigInterface.cc
  Context.cc
  Diagnostic.cc
  EnthalpyConverter.cc
  Grid.cc
  Logger.cc
  MaxTimestep.cc
  Profiling.cc
  ScalarForcing.cc
  TerminationReason.cc
  Time.cc
  Units.cc
  VariableMetadata.cc
  Vars.cc
  array/Array.cc
  array/Array3D.cc
  array/CellType.cc
  array/Forcing.cc
  array/Scalar.cc
  array/Staggered.cc
  array/Vector.cc
  GeometryCalculator.cc
  connected_components.cc
  error_handling.cc
  fem/DirichletData.cc
  fem/Element.cc
  fem/ElementIterator.cc
  fem/FEM.cc
  fem/Quadrature.cc
  fftw_utilities.cc
  interpolation.cc
  io/File.cc
  io/LocalInterpCtx.cc
  io/NC4File.cc
  io/NC4_Serial.cc
  io/NCFile.cc
  io/NC_Serial.cc
  io/io_helpers.cc
  label_components.cc
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

if(Pism_DEBUG)
  target_sources(util PRIVATE Poisson.cc)
endif()

if(Pism_USE_JANSSON)
  target_sources(util PRIVATE ConfigJSON.cc)
endif()

# Check if NetCDF-4 parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PARALLEL_NETCDF4)
  target_sources(util PRIVATE io/NC4_Par.cc)
endif()

# Check if PnetCDF parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PNETCDF)
  target_sources(util PRIVATE io/PNCFile.cc)
endif()

# Check if PnetCDF parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PIO)
  target_sources(util PRIVATE io/ParallelIO.cc)
endif()
