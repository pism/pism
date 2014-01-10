Ross prognostic flow model example
=================

This example demonstrates regional evolutionary modeling of ice shelves.

As in the diagnostic calculation (in the parent directory `examples/ross/`) we modify slightly the bed topography under the floating ice shelf and the temperatures in the initially-open ocean.  With the applied "eigencalving" parameterization this yields a (quasi) steady state.

It is recommended that you run `preprocess.py` in the parent directory `../`.  Then link the needed data files here to avoid downloading them again:

    $ ln -s ../antarctica_ice_velocity.nc 
    $ ln -s ../ALBMAPv1.nc.zip

Run scripts according to the diagnostic case (preprocess, run, plot).

    $ ./preprocess_prog.py
    $ ./run_prog.sh 2 211 0.6 50 1e17
    $ ../plot.py Mx211year50.nc           # generate figures comparing to present-day velocity
    $ ncview ex-Mx211year50.nc            # view "movie" of evolving thk,mask,csurf
