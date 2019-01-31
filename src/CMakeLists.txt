add_definitions(
  "-DPISM_REVISION=\"${Pism_REVISION_TAG}\""
  "-DPISM_DEFAULT_CONFIG_FILE=\"${Pism_DEFAULT_CONFIG_FILE}\"")

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
  geometry/part_grid_threshold_thickness.cc
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
  $<TARGET_OBJECTS:hydrology>
  $<TARGET_OBJECTS:flowlaws>
  $<TARGET_OBJECTS:stressbalance>
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

add_subdirectory (coupler)
add_subdirectory (earth)
add_subdirectory (external)
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
