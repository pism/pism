# Create a subdirectory in the build directory containing symlinks to directories in
# PISM's source tree. This way other code using PISM needs to add only *one* directory to
# the compiler search path; see #166. This subdirectory contains a generated header,
# otherwise we could create *one* symlink to ${CMAKE_CURRENT_SOURCE_DIR} instead.

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pism)

subdirlist(directories ${CMAKE_CURRENT_SOURCE_DIR})
foreach(d IN LISTS directories)
  execute_process(
    COMMAND
    ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/${d} ${CMAKE_BINARY_DIR}/pism/${d})
endforeach(d)
unset(directories)

# Generate the header containing information about PISM's build-time configuration.
configure_file(pism_config.hh.in ${CMAKE_BINARY_DIR}/pism/pism_config.hh ESCAPE_QUOTES)
install(FILES ${CMAKE_BINARY_DIR}/pism/pism_config.hh
  DESTINATION include/pism)

# Generate the source file containing constants defined in pism_config.hh Note that we
# don't escape quotes here: PETSc's configure flags should have them escaped already (if
# we're lucky).
configure_file(pism_config.cc.in pism_config.cc)

# PISM's build directory contains the symlink to its source directory (see above). This
# makes it possible to #include PISM's headers using '#include "pism/header.hh"' both in
# PISM's sources and in code linking to PISM as a library.
include_directories (${CMAKE_BINARY_DIR})

# This library contains PISM code implementing the ice-sheet model itself
# (using other PISM libraries and a good deal of non-trivial code).
add_library (pism
  ${EVERYTRACE_cf_mpi_REFADDR}
  ${CMAKE_CURRENT_BINARY_DIR}/pism_config.cc
  age/AgeColumnSystem.cc
  age/AgeModel.cc
  basalstrength/ConstantYieldStress.cc
  basalstrength/MohrCoulombYieldStress.cc
  basalstrength/MohrCoulombPointwise.cc
  basalstrength/OptTillphiYieldStress.cc
  basalstrength/YieldStress.cc
  energy/BedrockColumn.cc
  energy/BTU_Full.cc
  energy/BTU_Minimal.cc
  energy/BedThermalUnit.cc
  energy/DummyEnergyModel.cc
  energy/EnergyModel.cc
  energy/EnthalpyModel.cc
  energy/CHSystem.cc
  energy/TemperatureModel.cc
  energy/bootstrapping.cc
  energy/enthSystem.cc
  energy/tempSystem.cc
  energy/utilities.cc  
  geometry/Geometry.cc
  geometry/GeometryEvolution.cc
  geometry/grounded_cell_fraction.cc
  geometry/flux_limiter.cc
  geometry/part_grid_threshold_thickness.cc
  icemodel/IceModel.cc
  icemodel/IceEISModel.cc
  icemodel/frontretreat.cc
  icemodel/diagnostics.cc
  icemodel/diagnostics.cc
  icemodel/energy.cc
  icemodel/flux_balance.hh
  icemodel/fracture_density.cc
  icemodel/initialization.cc
  icemodel/output.cc
  icemodel/output_checkpoint.cc
  icemodel/output_extra.cc
  icemodel/output_save.cc
  icemodel/output_ts.cc
  icemodel/printout.cc
  icemodel/timestepping.cc
  icemodel/utilities.cc
  icemodel/viewers.cc
  stressbalance/timestepping.cc
  fracturedensity/FractureDensity.cc
  util/ColumnSystem.cc
  util/pism_signal.c
  $<TARGET_OBJECTS:frontretreat>
  $<TARGET_OBJECTS:hydrology>
  $<TARGET_OBJECTS:flowlaws>
  $<TARGET_OBJECTS:stressbalance>
  $<TARGET_OBJECTS:blatter>
  $<TARGET_OBJECTS:util>
  $<TARGET_OBJECTS:earth>
  $<TARGET_OBJECTS:verif>
  $<TARGET_OBJECTS:boundary>
  $<TARGET_OBJECTS:inverse>
  $<TARGET_OBJECTS:regional>
  $<TARGET_OBJECTS:calcalcs>
  $<TARGET_OBJECTS:cubature>
)
target_link_libraries (pism ${Pism_EXTERNAL_LIBS})
add_dependencies (pism pism_config)

if (Pism_USE_JANSSON)
  if (Pism_BUILD_JANSSON)
    add_dependencies (pism pism_jansson)
  endif()
endif()

# Main executable:
add_executable (pismr pismr.cc)
target_link_libraries (pismr pism)

# Verification mode "driver":
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
  pismr pismv # executables
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
  PATTERN "*_impl.hh" EXCLUDE
  PATTERN "external" EXCLUDE
  PATTERN "pythonbindings" EXCLUDE
  PATTERN "doc" EXCLUDE
  PATTERN "figs" EXCLUDE
  PATTERN "ssa/tests" EXCLUDE
  PATTERN "verification/tests/fortran" EXCLUDE
  PATTERN "rheology/approximate" EXCLUDE
  PATTERN "tracer" EXCLUDE
  )

add_subdirectory (coupler)
add_subdirectory (earth)
add_subdirectory (external)
add_subdirectory (frontretreat)
add_subdirectory (hydrology)
add_subdirectory (inverse)
add_subdirectory (regional)
add_subdirectory (rheology)
add_subdirectory (stressbalance)
add_subdirectory (util)
add_subdirectory (verification)

if (Pism_BUILD_ICEBIN)
  add_subdirectory (icebin)
endif()

if (Pism_DEBUG OR Pism_BUILD_PYTHON_BINDINGS)
  add_subdirectory (pythonbindings)
endif()
