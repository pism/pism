MISMIP3D in PISM
==============

This directory contains scripts that can be used to run MISMIP3d experiments using PISM. To understand the intent of these experiments, please see the MISMIP website at http://homepages.ulb.ac.be/~fpattyn/mismip3d/, and download the intercomparison description PDF from that site. Results will be soon published in 
Pattyn et al., 2012; "Grounding-line migration in plan-view marine ice-sheet models: results of the ice2sea MISMIP3d intercomparison" (submitted to JoGl).

To run the complete set of the MISMIP3d experiments, the three files create_runscript.py, createSetup_Stnd.py and createSetup_PXXS.py are needed. By running the file create_runscript.py a very simple structured shell script for running the experiments is created (details and options below), as well as a config-override file MISMIP3D_conf.nc. This shell script then uses the files createSetup_Stnd.py and createSetup_PXXS.py for creating the initial setup nc-files used by PISM and runs the experiments one after another. (Note that the script createSetup_PXXS.py needs the result of the Standard (Stnd) experiment for creating the setup of the P10S and P75S experiments. See the above mentioned description PDF for further details).


Step by step instructions
-------------------------

The script `create_runscript.py` is used to generate a `bash` script performing the set of MISMIP3d experiments. Following options can be set:

      -m MODEL, --model=MODEL
			    Choose SSA only or SIA+SSA computation. model=1 sets SSA only, model=2 SSA+SIA computation. If no option is set, model 2 is used.

      -s, --subgl  	    if this option is set, subgrid groundling line interpolation method is used (as defined as "LI" in Gladstone et al., 2010, "Parameterising the grounding line in flow-line ice sheet models; The Cryosphere 4, p.605--619). Also, the field gl_mask is written into the extra file to determine the subgrid groundling line position. The subgrid groundling line interpolation method is available from PISM in the development version (https://github.com/pism/pism/commit/626b309de4342f379c13ed879da694bbd96bada3). It modifies basal friction in grid cells, where the grounding lines position is identified. If no option is set, subgrid groundling line interpolation is not used.

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
 


Example
-------

For example, to set up a MISMIP3d experiment with SIA+SSA computation, grounding line interpolation and accumulation rate of 0.5m/a for a resolution of dx=dy=1km run

    ./create_runscript.py -m 2 -s -a 0.5 -r 2 > runscript.sh

This will create `runscript.sh` as well as the bootstrapping file
`MISMIP3D_conf.nc`. Running this script with the accordingly adapted path for the Python and PISM executable at the file header will create a full set of MISMIP3d experiments.


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

Converting PISM output files to ASCII files following MISMIP
specifications is left as an exercise. A bunch of additional variables for evaluation is saved to an extra file for each run. For output of the deviatoric stress components, see development version (https://github.com/pism/pism/commit/03f1c4f5676875efb3a4ae855426a7bf79c2fa63).
