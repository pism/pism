Antarctica example using SeaRISE-Antarctica data
=========

The scripts in this directory apply the PISM model to the present-day Antarctic
ice sheet.  The data are from the SeaRISE collaboration, specifically the file
`Antarctica_5km_dev1.0.nc` (105 MB) described at
<http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica>.

These scripts do not prepare an actual SeaRISE submission.  Some elements of the
modeling strategy in this example is explained in

* M. A. Martin, R. Winkelmann, M. Haseloff, T. Albrecht, E. Bueler, C. Khroulev, A. Levermann (2011).  _The Potsdam Parallel Ice Sheet Model (PISM-PIK) - Part 2: Dynamic equilibrium simulation of the Antarctic ice sheet_, **The Cryosphere** 5, 727--740. <http://www.the-cryosphere.net/5/727/2011/tc-5-727-2011.pdf>

Getting and preprocessing the data
---------

This script downloads `Antarctica_5km_dev1.0.nc` with `wget` and uses NCO to preprocess to PISM-ready condition:

    $ ./preprocess.sh

Running the coarse-grid spinup
---------

Next look at what the run script will do (i.e. a dry-run):

    $ PISM_DO=echo ./antspin-coarse.sh 4

Then actually do the run, saving its `stdout` output in a file; this will take a number of processor-hours:

    $ ./antspin-coarse.sh 4 &> out.ant30km

This is a 30 km grid (the default) run with 4 processes.  The first stage essentially smooths the surface, the second stage improves the enthalpy field (a "no mass continuity" run), and then the third stage uses the "full" physics, which has sliding and PIK calving-front physics.

Adjust the `antspin-coarse.sh` script to use other grids for this coarse grid stage.

Regridding to a fine grid to complete the spinup
---------

This is a 15 km grid run with 4 processes which continues the above "third" stage using the "full" physics, which has sliding and PIK calving-front physics.

    $ PISM_DO=echo ./antspin-regridtofine.sh 4
    $ ./antspin-regridtofine.sh 4 &> out.ant15km

Adjust the `antspin-regridtofine.sh` script to use other grids for the fine grid stage.

SeaRISE experiments
---------

To perform the actual SeaRISE experiments, see scripts in PISM release `stable0.4`.
Look in `examples/searise-antarctica`.  These scripts will require modifications
to run under more recent versions of PISM.
