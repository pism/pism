Ross diagnostic flow model example
=================

This example demonstrates regional modeling of ice shelves for prescribed geometry (basic diagnostic mode).

You can modify three parameters as options in a row: 

the number of processors used, the resolution and the enhancement factor

    $ ./run_diag.sh 2 211 0.6
    $ ../plot.py Ross_result_diag_Mx211.nc           # generate figures comparing to present-day velocity
