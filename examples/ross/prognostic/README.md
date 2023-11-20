Ross ice shelf model (prognostic eigencalving mode)
=================

This example demonstrates prognostic (geometry-evolving) modeling of an
ice shelf using the `-calving eigen_calving,thickness_calving` combination.

The user should probably run the example in `../diagnostic/` before this one,
and read the documentation for diagnostic example in section 12.2 of the PISM
User's Manual.

With the applied calving parameterization we reach a (quasi) steady state after
runs of a few centuries.  Calving parameters, specifically the minimal ice
thickness at the front (`-thickness_calving_threshold`) and the eigencalving
constant (`-eigen_calving_K`) are among the parameters that can be modified
either to get a different steady state or to model changes in the ice shelf.

As in the diagnostic example, start by running `preprocess.py` in the parent
directory.  Then do

    $ ./run_prog.sh 4 211 0.6 100

This 100 model year run on 4 processes and a 5 km grid took about about twenty
minutes.  It starts with a bootstrapping stage which does a `-y 0` run and then
generates `startfile_Mx211.nc`.  It then re-initializes to start the prognostic
run itself.

Note `run_prog.sh` accepts four arguments: `run_prog.sh N Mx E Y` does
a run with `N` MPI processes, a `Mx`x`Mx` grid, option `-ssa_e E`, and duration
`-y Y`.

To view a "movie" of the resulting thicknesses and surface velocities
do

    $ ncview ex-prog_Mx211_yr100.nc

To see a time-series of ice volume or ice area, among other other quantities,
do

    $ ncview ts-prog_Mx211_yr100.nc

To generate figures comparing the final-time modeled velocity to present-day
observed velocity do

    $ ../plot.py prog_Mx211_yr100.nc
