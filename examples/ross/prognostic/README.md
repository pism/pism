Ross prognostic flow model example
=================

This example demonstrates regional evolutionary modeling of ice shelves. According to the diagnostic mode we modify slightly bed topography under the floating ice shelf and the temperatures in the initially open ocean. With applied "eigencalving" parameterization this yields a (quasi) steady state.

Run scripts according to the diagnostic case (preprocess, run, plot).

    ./preprocess_prog.py
    ./run_prog.sh 2 211 0.6 50 1e17
    ./plot_prog.py --pism-output=Mx211_year-00050.nc
