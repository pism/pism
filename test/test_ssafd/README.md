This script goes through all the possible geometric configurations of
the 9 cells in a 3x3 neighborhood and runs PISM for 0 years to test
the SSAFD solver.

If the solver fails to converge, PISM will save the model state to a
file.

This experiment confirms that the only 16 out of 256 configurations
that lead to a failure are ones that have a non-zero floating ice cell
(zero basal drag) with all four immediate neighbors ice-free.

See
`src/base/stressbalance/ssa/doc/discretization/ssa_test_configuration.mac`
for more.
