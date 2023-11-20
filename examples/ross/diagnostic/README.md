Ross ice shelf model (diagnostic validation mode)
=================

This example is documented in subsection 12.2 of the PISM User's Manual as an
example of using observed velocities to validate PISM's model for the flow in
ice shelves.

This example is a regional model of an ice shelf using observed, fixed geometry
and observed velocities as boundary conditions.  PISM does a `-y 0` (zero
duration) diagostic run; see subsection 5.5 of the User's Manual.  A python plot
script can compare observed and numerical velocities in the interior of the ice
shelf once the run completes.

To do the example, run `preprocess.py` in the parent directory, and then do the
following in the current directory:

    $ ./run_diag.sh 2 211 0.6

At this resolution the run takes only a few seconds.  To compare to
observations, and generate figures like those in the User's Manual, do

    $ ../plot.py diag_Mx211.nc

Note `run_diag.sh` accepts three arguments: `run_diag.sh N Mx E` does a run with
`N` MPI processes, a `Mx`x`Mx` grid, and option `-ssa_e E`.
