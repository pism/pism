add_library (util OBJECT
  fem/DirichletData.cc
  fem/Element.cc
  fem/ElementIterator.cc
  fem/FEM.cc
  fem/Quadrature.cc
  ColumnInterpolation.cc
  Context.cc
  EnthalpyConverter.cc
  Grid.cc
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
  array/CellType.cc
  array/Array.cc
  array/Forcing.cc
  array/Vector.cc
  array/Array3D.cc
  array/Scalar.cc
  array/Staggered.cc
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
  label_components.cc
  connected_components.cc
  ScalarForcing.cc
  Interpolation1D.cc
  InputInterpolation.cc
)

if(Pism_USE_YAC_INTERPOLATION)
  target_sources(util PRIVATE YACInterpolation.cc)
endif()

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
