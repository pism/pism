Ross prognostic flow model example with applied fracture softening
=================

This example demonstrates regional modeling of ice shelf evolution with applied calculation of fracture density and respective macroscopic softening.

This example uses `-calving ocean_kill` calving to hold the calving front in a fixed location.  Compare the example in `examples/ross/prognostic/`.

It is recommended that you run `preprocess.py prog` in the parent directory.

Now run the script for simulations. Modify fracture options in the runscript or as additional options:

    $ ./run_frac.sh 8 211 0.6 500 7e16 0.5 5.0e-10 0.1 1.0

    $ ncview ex-Ross_result_frac_Mx211_yr-500.nc        # view "movie" of evolving fracture_density,thk,csurf

