Ross prognostic flow model example
=================

This example demonstrates regional evolutionary modeling of ice shelves. According to the diagnostic mode we modify slightly bed topography under the floating ice shelf and the temperatures in the initially open ocean. With applied "eigencalving" parameterization this yields a (quasi) steady state.

It is recommended that you run `preprocess.py` in the parent directory `../`.  Then link the needed data files here to avoid downloading them again:

    $ ln -s ../antarctica_ice_velocity.nc 
    $ ln -s ../ALBMAPv1.nc.zip

Run scripts according to the diagnostic case (preprocess, run, plot).

    $ ./preprocess_prog.py
    $ ./run_prog.sh 2 211 0.6 50 1e17
    $ ./plot_prog.py --pism-output=Mx211_year-00050.nc
