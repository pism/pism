MISMIP3D in PISM
==============

This directory contains scripts that can be used to run MISMIP3d experiments using PISM.  See section 10  .3 of the PISM User's Manual.

In outline, running `preprocess.py` creates a config-override file `MISMIP3D_conf.nc`.  Then running `create_runscript.py` creates a shell script for running the experiments.  This shell script uses  `createSetup_Stnd.py` or `createSetup_PXXS.py` for creating the initial setup `.nc`-files which contain the geometry and basal sliding fields, and then it runs PISM for the experiment.

Note that the script `createSetup_PXXS.py` needs the result of the Standard (Stnd) experiment for creating the setup of the P10S and P75S experiments.  The `PXXR` runs require the output of the corresponding `PXXS` experiment.


Usage
-------

For example, to set up and run the MISMIP3d "Stnd" experiment with SIA+SSA computation, grounding line interpolation and accumulation rate of 0.5m/a, and a resolution of dx=dy=1km, do:

    $ ./preprocess      # generate MISMIP3D_conf.nc, used by all experiments
    $ ./create_runscript.py -m 2 -s -a 0.5 -r 2 -e Stnd > runscript.sh
    $ bash runscript.sh > Stnd.out &

FIXME: with adapted paths for the Python and PISM executables
FIXME: add P10S


Options for `create_runscript.py`
-------------------------

The script `create_runscript.py` is used to generate a `bash` script performing the set of MISMIP3d experiments. Following options can be set:

      -m MODEL, --model=MODEL
			    Choose SSA only or SIA+SSA computation. model=1 sets SSA only, model=2 SSA+SIA computation. If no option is set, model 2 is used.

      -s, --subgl
                if this option is set, subgrid groundling line interpolation method is used (as defined as "LI" in Gladstone et al., 2010, "Parameterising the grounding line in flow-line ice sheet models; The Cryosphere 4, p.605--619). Also, the field gl_mask is written into the extra file to determine the subgrid groundling line position. The subgrid groundling line interpolation method modifies basal friction in grid cells where the grounding lines position is identified. If no option is set, subgrid groundling line interpolation is not used.

      -a ACCUMRATE, --accumrate=ACCUMRATE
			    sets the accumulation rate in meters per year. If no option is set, standard accumulation rate of a=0.5 m/a is used.

      -r RESOLUTIONMODE, --resolutionmode=RESOLUTIONMODE
			    sets the grid resolution of the computational domain. Because the MISMP3d domain has the fixed size of 800km x 50km and for ensuring equally spaced boxes in x- and y-direction, only discrete choices of resolution are offered:
			     resolutionmode=1: 0.5km
			     resolutionmode=2: 1.0km
			     resolutionmode=3: 2.5km
			     resolutionmode=4: 5.0km
			     resolutionmode=5: 10.0km
			     resolutionmode=6: 16.6km
			     If no option is set, standard resolution of 16.6 km is used. Note that the computational domain is doubled in both directions to 1600km x 100km to omit the implementation of boundary conditions along symmerty lines, but for additional computational cost.

The initial ice sheet configuration from which the Stnd-experiment is started has the constant thickness of 500 meters. For all experiments, the ice shelf is cut off at the distance of x=700 km (option -ocean_kill) from the center of the computational domain, where the stress boundary condition is applied.


Implementation details
----------------------

We can turn PISM's default sliding law into MISMIP's power law by setting the
threshold speed to 1 meter per second, which will make it inactive.

The `-pseudo_plastic_uthreshold` command-line option takes an argument in meters per year, so we use `-pseudo_plastic_uthreshold 3.15569259747e7`, where `3.15569259747e7` is the number of seconds in a year.

The MISMIP parameter C corresponds to `tauc` in PISM. It can be set using `-hold_tauc -tauc C`.

The MISMIP power law exponent `m` corresponds to `-pseudo_plastic_q` in PISM.

We use the `-config_override` option to set other MISMIP-specific parameters, such as ice softness, ice density and others.

Note that PISM does not at this time implement the stopping criteria described in the MISMIP specification.  Instead we use the maximum run lengths that are provided as an alternative. On the other hand, PISM's output files contain all the information necessary to compute the rate of change of the grounding line position and the thickness rate of change during post-processing.


Post-processing
---------------

Converting PISM output files to ASCII files following MISMIP specifications is left as an exercise.  See the additional variables saved in the extra file for each run.

