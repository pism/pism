# The terrain-insolation kernels are free of PISM, PETSc, and MPI, so this test compiles
# the kernel directly instead of linking libpism. This keeps it fast and lets it run in
# builds without the Python bindings.
add_executable (pism_test_terrain_insolation_kernel
  terrain_insolation_kernel_test.cc
  ../terrain_insolation_kernel.cc)
