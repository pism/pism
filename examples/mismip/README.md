MISMIP in PISM
==============

This directory contains scripts that can be used to run [MISMIP](http://homepages.ulb.ac.be/~fpattyn/mismip/) experiments using PISM.

Older PISM versions included C++ code managing MISMIP experiments. With the addition of more sophisticated reporting code that old code became unnecessary.

Step by step instructions
-------------------------

First of all, you will need to copy or symlink `util/PISMNC.py` into the
current directory to make sure that Python will be able to import it.

The script `run.py` is used to generate a `bash` scrips performing MISMIP
experiments.


Running `run.py --help` produces the following:

    Usage: run.py [options]

    Creates a script running MISMIP experiments.

    Options:
      -h, --help            show this help message and exit
      --initials=INITIALS   Initials (3 letters)
      -e EXPERIMENT, --experiment=EXPERIMENT
                            MISMIP experiments (one of '1a', '1b', '2a', '2b',
                            '3a', '3b')
      -u, --uniform_thickness
                            Use uniform 10 m ice thickness
      -a, --all             Run all experiments
      -m MODE, --mode=MODE  MISMIP grid mode
      --Mx=MX               Custom grid size; use with --mode=3
      --model=MODEL         Models: 1 - SSA only; 2 - SIA+SSA
      --executable=EXECUTABLE
                            Executable to run, e.g. 'mpiexec -n 4 pismr'

For example, to set up MISMIP experiment `1a` using grid mode 1 (150 grid
spaces in the MISMIP modeling domain or 301 grid points in PISM's (doubled)
domain), run

    ./run.py -e 1a --mode=1 > experiment-1a-mode-1.sh

This will create `experiment-1a-mode-1.sh` as well as the bootstrapping file
`MISMIP_boot_1a_M1_A1.nc` and configuration files corresponding to each "step"
in the experiment.

Now you can either run `experiment-1a-mode-1.sh` or copy it (along with
`MISMIP_boot_1a_M1_A1.nc` and `MISMIP_conf_1a_A*.nc`) to a supercomputer to
do the run later.

Note that the `experiment-1a-mode-1.sh` uses environment variables `PISM_DO`,
`PISM_PREFIX` and `PISM_MPIDO`. For example,

    PISM_MPIDO="aprun -n " ./experiment-1a-mode-1.sh 32

will use `aprun` on 32 cores.

Alternatively, you can use

    ./run.py -e 1a --mode=1 --executable="aprun -n 32 pismr"

or similar to skip the "preamble" handling environment variables and get "raw"
commands.

Technical details
-----------------

The script `MISMIP.py` contains MISMIP parameters and the code needed to
compute the semi-analytic grounding line location and the corresponding
thickness profile for each experiment.

The script `prepare.py` contains functions using `MISMIP.py` to generate
PISM-readable NetCDF files with semi-analytic ice thickness profiles, and
the prescribed accumulation map. This script can be imported as a module or run
as a script to generate PISM bootstrapping files.

The script `run.py` generates `bash` scripts performing MISMIP runs using
`MISMIP.py` and `prepare.py`.

Implementation details
----------------------

The only addition to the PISM code necessary to run MISMIP experiment is the
sliding law; see `src/base/basal_strength/MISMIPBasalResistanceLaw.cc`. This
code is turned "on" using the "`-mismip_sliding`" command-line option.

Once selected, `MISMIPBasalResistanceLaw` expects to find configuration
parameters `MISMIP_m` (the the sliding law exponent),
`MISMIP_C` (a multiplicative factor in the sliding law),
and `MISMIP_r` (a regularization parameter) in the configuation database.
We use the `-config_override` option to provide these, along with
MISMIP-specific values of the ice softness, ice density, etc.

Note that PISM does not implement stopping criteria described in the MISMIP specs;
we use maximum run lengths that are provided as an alternative. On the other hand,
PISM's output files contain all the information necessary to compute the rate of change
of the grounding line position and the thickness rate of change during post-processing.

Post-processing
---------------

Converting PISM output files to ASCII files following MISMIP
specifications is left as an exercise.

We do provide the script `showflux.py`, though, which plots ice flux as a
function of the distance from the divide. It shows the discontinuity in
the flux at the drounding line. This is an issue in PISM that needs to be
addressed to improve its handling of the grounding line motion.

Note also that the variable `iareag` in `ts_ABC1_1a_M1_A1.nc` and similar
allows one to see time-dependent changes in the grounding line location
(grounded ice area is proportional to the distance from the divide to the
grounding line).
