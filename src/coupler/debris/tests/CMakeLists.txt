# The debris column and terminus kernels are free of PISM, PETSc, and MPI, so this test
# compiles them directly instead of linking libpism. This keeps it fast and lets it run in
# builds without the Python bindings.
add_executable (pism_test_debris_kernels
  debris_kernels_test.cc
  ../column_kernels.cc
  ../terminus_kernels.cc)
