Ross prognostic flow model example with applied fracture softening
=================

This example demonstrates regional modeling of ice shelf evolution with applied calculation of fracture density and respective macroscopic softening.

This example uses `-calving ocean_kill` calving to hold the calving front in a fixed location.  Compare the example in `examples/ross/prognostic/`.

It is recommended that you run `preprocess.py` in the parent directory `../`.  Then link the needed data files here to avoid downloading them again:

    $ ln -s ../antarctica_ice_velocity.nc 
    $ ln -s ../ALBMAPv1.nc.zip

Now run the preprocess script to build computational setup (use nc2cdo.py in util/) and comment and activate the run script for simulations. Modify fracture options in runscript.

    $ ./preprocess_frac.py
    $ ./run_frac.sh 8 211 0.6 500 7e16 0.5 5.0e-10 0.1 1.0

