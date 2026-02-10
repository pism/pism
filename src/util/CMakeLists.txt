add_library (util OBJECT
  ColumnInterpolation.cc
  ColumnSystem.cc
  Component.cc
  Config.cc
  Context.cc
  Diagnostic.cc
  EnthalpyConverter.cc
  Grid.cc
  InputInterpolation.cc
  Interpolation1D.cc
  Logger.cc
  Mask.cc
  MaxTimestep.cc
  NetCDFConfig.cc
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
  connected_components/label_components_parallel.cc
  connected_components/label_components_serial.cc
  error_handling.cc
  fem/DirichletData.cc
  fem/Element.cc
  fem/ElementIterator.cc
  fem/FEM.cc
  fem/Quadrature.cc
  fftw_utilities.cc
  io/File.cc
  io/LocalInterpCtx.cc
  io/NC4File.cc
  io/NC4_Serial.cc
  io/NCFile.cc
  io/NC_Serial.cc
  io/OutputFile.cc
  io/OutputWriter.cc
  io/SynchronousOutputWriter.cc
  io/input_helpers.cc
  io/output_helpers.cc
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
  pism_initialization.cc
  pism_options.cc
  pism_signal.c
  pism_utilities.cc
  projection.cc
)

if(Pism_USE_YAC)
  set_source_files_properties(InputInterpolationYAC.cc
    PROPERTIES
    INCLUDE_DIRECTORIES "${YAC_INCLUDE_DIRS}")
  target_sources(util PRIVATE InputInterpolationYAC.cc)

  set_source_files_properties(pism_initialization.cc
    PROPERTIES
    INCLUDE_DIRECTORIES "${YAC_INCLUDE_DIRS}")

  set_source_files_properties(io/YacOutputWriter.cc
    PROPERTIES
    INCLUDE_DIRECTORIES "${YAC_INCLUDE_DIRS}")
  target_sources(util PRIVATE io/YacOutputWriter.cc)

  set_source_files_properties(yaxt_wrapper.c
    PROPERTIES
    INCLUDE_DIRECTORIES "${YAXT_INCLUDE_DIRS}")
  target_sources(util PRIVATE yaxt_wrapper.c)
endif()

if(Pism_DEBUG)
  target_sources(util PRIVATE Poisson.cc)
endif()

# Check if NetCDF-4 parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PARALLEL_NETCDF4)
  target_sources(util PRIVATE io/NC4_Par.cc)
endif()

# Check if PnetCDF parallel I/O is enabled and add a source code file if necessary.
if (Pism_USE_PNETCDF)
  target_sources(util PRIVATE io/PNCFile.cc)
endif()
