Antarctica example using ALBMAPv1 data
=========

The scripts in this directory apply PISM to the present-day Antarctic ice sheet. The data
are from the ALBMAP dataset, specifically the file `ALBMAPv1.nc.zip` (16 MB) described in

> A. M. Le Brocq, A. J. Payne, and A. Vieli, “An improved Antarctic dataset for high
> resolution numerical ice sheet models (ALBMAP v1),” Earth System Science Data, vol. 2,
> no. 2, pp. 247–260, Oct. 2010, doi:
> [10.5194/essd-2-247-2010](https://doi.org/10.5194/essd-2-247-2010).

This data set can be downloaded from <https://doi.org/10.1594/PANGAEA.734145>.

Some elements of the modeling strategy in this example are explained in

> M. A. Martin, R. Winkelmann, M. Haseloff, T. Albrecht, E. Bueler, C. Khroulev, and A.
> Levermann. The potsdam parallel ice sheet model (PISM-PIK) – part 2: dynamic equilibrium
> simulation of the antarctic ice sheet. The Cryosphere, 5:727–740, 2011. doi:
> [10.5194/tc-5-727-2011](https://doi.org/10.5194/tc-5-727-2011).

Getting and preprocessing the data
---------

This script downloads `ALBMAPv1.nc` with `wget` and uses NCO to prepare it for use with PISM:

    $ ./preprocess.sh

Running the coarse grid spinup
---------

Next look at what the run script will do (i.e. a dry-run):

    $ PISM_DO=echo ./antspin-coarse.sh 4

Then actually do the run, saving its `stdout` output in a file; this will take a number of processor-hours:

    $ ./antspin-coarse.sh 4 &> out.ant30km

This is a ~30 km grid (the default) run with 4 processes.  The first stage essentially smooths the surface, the second stage improves the enthalpy field (a "no mass continuity" run), and then the third stage uses the "full" physics, which has sliding and PIK calving-front physics.

Adjust the `antspin-coarse.sh` script to use other grids for this coarse grid stage.

Regridding to a fine grid to complete the spinup
---------

This is a ~15 km grid run with 4 processes which continues the above "third" stage using the "full" physics, which has sliding and PIK calving-front physics.

    $ PISM_DO=echo ./antspin-regridtofine.sh 4
    $ ./antspin-regridtofine.sh 4 &> out.ant15km

Adjust the `antspin-regridtofine.sh` script to use other grids for the fine grid stage.

SeaRISE experiments
---------

To perform the actual SeaRISE experiments, see scripts in PISM release `stable0.4`.
Look in `examples/searise-antarctica`.  These scripts will require modifications
to run under more recent versions of PISM.
