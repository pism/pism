Ross prognostic flow model example with applied fracture softening
=================

This example demonstrates regional modeling of ice shelf evolution with applied calculation of fracture density and respective macroscopic softening.

It is recommended that you run `preprocess.py` in the parent directory `../`.  Then link the needed data files here to avoid downloading them again:

    $ ln -s ../antarctica_ice_velocity.nc 
    $ ln -s ../ALBMAPv1.nc.zip

Now run the preprocess script to build computational setup (use nc2cdo.py in util/) and comment and activate the run script for simulations. Modify fracture options in runscript.

    $ ./preprocess_prog.py
    $ ./run_frac.sh 8 211 0.6 500 7e16 0.5 5.0e-10 0.1 1.0

