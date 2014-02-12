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
directory.  Then do

    $ ./run_frac.sh 8 211 0.6 500 7e16 0.5 5.0e-10 0.1 1.0

