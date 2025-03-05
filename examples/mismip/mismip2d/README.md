MISMIP in PISM
==============

This directory contains scripts that can be used to run MISMIP experiments using PISM.  To understand the intent of these experiments, please see the MISMIP website at https://frank.pattyn.web.ulb.be/mismip/welcome.html, and download the intercomparison description PDF from that site.

Older PISM versions included C++ code managing MISMIP experiments. With the addition of more sophisticated reporting code that old code became unnecessary.  Here all MISMIP-specific code is in Python scripts.

Step by step instructions
-------------------------

First of all, you will need to copy or symlink `util/PISMNC.py` into the
current directory to make sure that Python will be able to import it.

The script `run.py` is used to generate a `bash` scripts performing MISMIP
experiments.  Running `run.py --help` produces the following:

    Usage: run.py [options]

    Creates a script running MISMIP experiments.

    Options:
      -h, --help            show this help message and exit
      --initials=INITIALS   Initials (3 letters)
      -e EXPERIMENT, --experiment=EXPERIMENT
                            MISMIP experiments (one of '1a', '1b', '2a', '2b',
                            '3a', '3b')
      -s STEP, --step=STEP  MISMIP step number
      -u, --uniform_thickness
                            Use uniform 10 m ice thickness
      -a, --all             Run all experiments
      -m MODE, --mode=MODE  MISMIP grid mode
      --Mx=MX               Custom grid size; use with --mode=3
      --model=MODEL         Models: 1 - SSA only; 2 - SIA+SSA
      --executable=EXECUTABLE
                            Executable to run, e.g. 'mpiexec -n 4 pism'

For example, to set up MISMIP experiment `1a` using grid mode 1, a 12 km grid, run

    ./run.py -e 1a --mode=1 > experiment-1a-mode-1.sh

This will create `experiment-1a-mode-1.sh` as well as the bootstrapping file
`MISMIP_boot_1a_M1_A1.nc` and configuration files corresponding to each "step"
in the experiment.

Run in the backround with 2 cores and saving output to a text file this way:

    bash experiment-1a-mode-1.sh 2 >& out.1a-mode-1 &

You can also copy the script (along with
`MISMIP_boot_1a_M1_A1.nc` and `MISMIP_conf_1a_A*.nc`) to a supercomputer to
do the run later.  For such application, the script helpfully uses environment variables `PISM_DO`,
`PISM_BIN` and `PISM_MPIDO`. For example, on some Cray machines you might do

    PISM_MPIDO="aprun -n " bash experiment-1a-mode-1.sh 32

will use `aprun` on 32 cores.  Alternatively, you can use

    ./run.py -e 1a --mode=1 --executable="aprun -n 32 pism"

or similar to skip the "preamble" handling environment variables and get "raw"
commands.


Refined grid runs
-----------------

The above "grid mode 1" runs use 150 grid spaces in the MISMIP modeling domain,
which is 301 grid points in PISM's (doubled) domain.  The domain is doubled because
PISM is easiest configure as a whole ice sheet model with ice free ocean at the
edge of the computation domain.  (Compare the example in `examples/jako/`, however.)

To run a higher resolution 3 km grid, with somewhat-improved grounding line 
performance, ask the `run.py` script to put option `-Mx 1201` into the bash
script:

    ./run.py -e 1a --mode=3 --Mx=1201 > experiment-1a-mode-3.sh

Then this is a 4 core run:

    bash experiment-1a-mode-3.sh 4 >& out.1a-mode-3 &


Technical details
-----------------

The script `MISMIP.py` contains MISMIP parameters and the code needed to
compute the semi-analytic grounding line location and the corresponding
thickness profile for each experiment.

The script `prepare.py` contains functions using `MISMIP.py` to generate
PISM-readable NetCDF files with semi-analytic ice thickness profiles, and
the prescribed accumulation map. It can be imported as a module or run
as a script to generate PISM bootstrapping files.

The script `run.py` generates `bash` scripts performing MISMIP runs using
`MISMIP.py` and `prepare.py`.

Implementation details
----------------------

We can turn PISM's default sliding law into MISMIP's power law by setting the
threshold speed to 1 meter per second, which will make it inactive.

The `-pseudo_plastic_uthreshold` command-line option takes an argument in meters per year, so we use `-pseudo_plastic_uthreshold 3.15569259747e7`, where `3.15569259747e7` is the number of seconds in a year.

The MISMIP parameter C corresponds to `tauc` in PISM. It can be set using `-yield_stress constant -tauc C`.

The MISMIP power law exponent `m` corresponds to `-pseudo_plastic_q` in PISM.

We use the `-config_override` option to set other MISMIP-specific parameters, such as ice softness, ice density and others.

Note that PISM does not at this time implement the stopping criteria described in the MISMIP specification.  Instead we use the maximum run lengths that are provided as an alternative. On the other hand, PISM's output files contain all the information necessary to compute the rate of change of the grounding line position and the thickness rate of change during post-processing.

Post-processing
---------------

Converting PISM output files to ASCII files following MISMIP
specifications is left as an exercise.

However, we do provide the `plot.py` script for visualization.

It plots modeled ice flux as a function of the distance from the divide and the geometry profile.

We see a discontinuity in the flux at the grounding
line. This is an issue in PISM that needs to be addressed to improve its
handling of the grounding line motion. For example, try

    ./plot.py ABC1_1a_M1_A1.nc -f -o flux.png

Run `plot.py --help` for a list of command-line options. You can also produce plots for several PISM output files at once. For example,

    ./plot.py ABC*.nc

will create geometry profile *and* ice flux plots for all matched files.

Also, note that the variable `ice_area_glacierized_grounded` in `ts_ABC1_1a_M1_A1.nc` and similar
allows one to see time-dependent changes in the grounding line location
because grounded ice area is proportional to the distance from the divide to the
grounding line.
