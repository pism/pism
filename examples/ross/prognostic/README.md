Ross ice shelf model (prognostic eigencalving mode)
=================

This example demonstrates regional evolutionary modeling of ice shelves.

As in the diagnostic calculation we modify slightly the bed topography under the floating ice shelf and the temperatures in the initially-open ocean.  With the applied "eigencalving" parameterization this yields a (quasi) steady state. Calving options (minimal ice thickness at the front and eigencalving constant) can be modified as additional options.

It is recommended that you run `./preprocess.py diag` in the parent directory. 

Run scripts according to the diagnostic case (preprocess, run, plot).

    $ ./run_prog.sh 2 211 0.6 50 1e17
    $ ../plot.py Ross_result_prog_Mx211_yr-500.nc          # generate figures comparing to present-day velocity
    $ ncview ex-Ross_result_prog_Mx211_yr-500.nc           # view "movie" of evolving thk,mask,csurf
