SeaRISE-Antarctica example
=========

These are minimal scripts which are less than needed for an actual
SeaRISE submission.  The modeling strategy in this example is explained in

* M. A. Martin, R. Winkelmann, M. Haseloff, T. Albrecht, E. Bueler, C. Khroulev, A. Levermann (2011).  _The Potsdam Parallel Ice Sheet Model (PISM-PIK) - Part 2: Dynamic equilibrium simulation of the Antarctic ice sheet_, **The Cryosphere** 5, 727--740. <http://www.the-cryosphere.net/5/727/2011/tc-5-727-2011.pdf>

Getting and preprocessing the data
---------

First download the data and do the preprocessing by calling the `preprocess.sh`
script.  The downloaded files `Antarctica_5km_dev1.0.nc` (105 MB) and
`ANT_climate_forcing_2004_2098_v3.nc` (1.1 GB) are described at
<http://websrv.cs.umt.edu/isis/index.php/Present_Day_Antarctica>.
This script downloads with `wget` and uses NCO to preprocess to PISM-ready
condition:

    $ ./preprocess.sh

Running the example spinup
---------

Next look at the script which would run (i.e. a dry-run):

    $ PISM_DO=echo ./antspinCC.sh N | less 

Then actually do the run in the background, saving its `stdout` output in a
file; this will take a number of processor-hours:

    $ ./antspinCC.sh N &> out.ant30km &

This is a 30 km grid run with N processes.  The first
stage essentially smooths the surface, the second stage improves the enthalpy
field, and then the third stage uses the "full" physics.  (I.e. sliding plus
PIK calving front physics).

Higher-resolution runs can be achieved by modifying the script; reset the `GRID`
and `SKIP` variables.  Different parameter values can be set also.

SeaRISE experiments
---------

To perform the SeaRISE experiments, see scripts in PISM release `stable0.4`.
Look in `examples/searise-antarctica`.  These scripts will require modifications
to run under more recent versions of PISM.

