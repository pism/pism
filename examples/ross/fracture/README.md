Ross ice shelf model (prognostic fracture softening mode)
=================

This example demonstrates a model of ice shelf evolution with a calculation of
fracture density and macroscopic softening.  This example uses
`-calving ocean_kill` calving to hold the calving front in a fixed location;
compare the example in `examples/ross/prognostic/`.

The user should probably run the example in `../diagnostic/` before this one,
and read the documentation for diagnostic example in section 12.2 of the PISM
User's Manual.

As in the diagnostic example, start by running `preprocess.py` in the parent
directory.  Then do a default 3000 year run this way:

    $ ./run_frac.sh 4 211 0.6 &> out.Mx211yr3000 &

Note `run_frac.sh` accepts up to four arguments: `run_frac.sh N Mx E Y` does
a run with `N` MPI processes, a `Mx`x`Mx` grid, option `-ssa_e E`, and duration
`-y Y`.

