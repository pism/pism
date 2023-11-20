MISMIP3D in PISM
==============

This directory contains scripts that can be used to run MISMIP3d experiments using PISM.  See section 10  .3 of the PISM User's Manual.

In outline, running `preprocess.py` creates a config-override file `MISMIP3D_conf.nc`.  Then running `createscript.py` creates a shell script for running the experiments.  This shell script uses  `setup_Stnd.py` or `setup_PXXS.py` for creating the initial setup `.nc`-files which contain the geometry and basal sliding fields, and then it runs PISM for the experiment.

Note that the script `setup_PXXS.py` needs the result of the Standard (Stnd) experiment for creating the setup of the P10S and P75S experiments.  The `PXXR` runs require the output of the corresponding `PXXS` experiment.


Usage
-------

For example, to set up and run the MISMIP3d "Stnd" experiment with SIA+SSA computation (model 2), grounding line interpolation, a resolution of dx = dy = 2 km, and a (shortened) duration of 3000 years, do:

    $ ./preprocess      # generate MISMIP3D_conf.nc, used by all experiments
    $ ./createscript.py -n 2 -r 3 -d 3000 -s -e Stnd > runStnd.sh
    $ bash runStnd.sh > out.Stnd &

This 2 process run (`-n 2`) takes about three minutes on a 2013 laptop.

This `Stnd` run uses `My=3` as would a MISMIP flowline experiment, and it uses `Mx=801` because of the 2 km resolution.  Three files will be output, `Stnd.nc`, `ts_Stnd.nc`, and `ex_Stnd.nc`.

Now do the first experiment, which uses information from `ex_Stnd.nc` and builds a `Mx=801` by `My=51` grid, again with 2 km resolution:

    $ ./createscript.py -n 2 -r 3 -d 100 -s -e P75S > runP75S.sh
    $ bash runP75S.sh > out.P75S &

This 2 process run and the next each take about 25 minutes on a 2013 laptop.  Next we try a reversibility experiment:

    $ ./createscript.py -n 2 -r 3 -d 100 -s -e P75R > runP75R.sh
    $ bash runP75R.sh > out.P75R &

Note that more complete reversibility requires longer runs.  In fact, the PISM submission used `-d 30000` for the `Stnd` experiment and `-d 500` for the other experiments.  Also, the submitted MISMIP3d runs used resolution mode (`-r 2`) with 1 km spacing.

The remaining implemented experiments are `P10S` and `P10R`.  They are run with the obvious modification, namely option `-e P10S` or `-e P10R` to the `createscript.py` script.


Options for `createscript.py`
-------------------------

Run

    $ ./createscript.p -h

to see a usage message.

Regarding the `-s` or `--subgl` option, if it is set then a subgrid grounding line (SGL) interpolation method is used.  This kind of method is called "LI" in Gladstone et al., 2010, "Parameterising the grounding line in flow-line ice sheet models; The Cryosphere 4, p.605--619.  If this option is set then the field `gl_mask` is written into the extra file to determine the SGL position. The SGL interpolation method modifies basal friction in grid cells where the grounding lines position is identified.  By default this option is not set.

The `-r MODE` option sets the resolution mode.  Because the MISMP3d domain has the fixed size of 800 km x 50 km, and for ensuring equally spaced boxes in x- and y-direction, only discrete choices of resolution are offered:

    - MODE = 1: 0.5 km
    - MODE = 2: 1.0 km
    - MODE = 3: 2.0 km
    - MODE = 4: 2.5 km
    - MODE = 5: 5.0 km  [the default]
    - MODE = 6: 10.0 km
    - MODE = 7: 16.6 km

Note that the computational domain is doubled in both directions to 1600 km x 100 km.  This avoids implementation of boundary conditions along symmerty lines, though at additional computational cost.

The initial ice sheet configuration, from which the Stnd experiment is started, has constant thickness of 500 meters.

For all the experiments, the ice shelf is cut off at the distance of x=700 km (option `-calving ocean_kill`) from the center of the computational domain.  At this calving front a stress boundary condition (`-cfbc`) is applied.

For more details read the source code of `createscript.py`.


Implementation details
----------------------

We turn PISM's default sliding law into MISMIP's power law by setting the
threshold speed to 1 meter per second.

The `-pseudo_plastic_uthreshold` command-line option takes an argument in meters per year, so we use `-pseudo_plastic_uthreshold 3.15569259747e7`, where `3.15569259747e7` is the number of seconds in a year.

The MISMIP parameter C corresponds to `tauc` in PISM.  It can be set using `-yield_stress constant -tauc C`.

The MISMIP power law exponent `m` corresponds to `-pseudo_plastic_q` in PISM.

We use the `-config_override` option to set other MISMIP-specific parameters, such as ice softness, ice density and others.

Note that PISM does not at this time implement the stopping criteria described in the MISMIP specification.  Instead we use the maximum run lengths that are provided as an alternative. On the other hand, PISM's output files contain all the information necessary to compute the rate of change of the grounding line position and the thickness rate of change during post-processing.


Post-processing
---------------

Converting PISM output files to ASCII files following MISMIP specifications is left as an exercise.  See the additional variables saved in the extra file for each run.
