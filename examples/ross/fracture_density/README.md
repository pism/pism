Ross prognostic flow model example with applied fracture softening
=================

This example demonstrates regional modeling of ice shelf evolution with applied calculation of fracture density and respective macroscopic softening.

Run preprocess script to build computational setup (use nc2cdo.py in util/) and comment and activate the run script for simulations. Modify fracture options in runscript.

    ./preprocess_prog.py
    ./run_prog.sh 8 211 0.6 500 7e16 0.5 5.0e-10 0.1 1.0

