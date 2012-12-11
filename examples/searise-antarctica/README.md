SeaRISE-Antarctica example
=========

These are minimal scripts which are less than needed for an actual
SeaRISE submission.  Most of the modeling in this example is explained in

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

Finally, do the run.  This is a 15 km grid run with N processes.  The first
stage essentially smooths the surface, the second stage improves the enthalpy
field, and then the third stage uses the "full" physics.  (I.e. sliding plus
PIK calving front physics).  This is a long run which will take many
processor-hours:

    $ ./antspinCC.sh N

Higher-resolution runs can be achieved by modifying the script.

SeaRISE experiments
---------

To perform the SeaRISE experiments, there are scripts that ran under PISM
`stable0.4`; they will require modifications to run under `stable0.5`.

    $ ./experiments.sh N
    $ ./postprocess.sh       # calls postprocess_mask.py

