## PICO vs PICOP

This directory contains scrips to compare PICO vs PICOP using a
simple Antartica setup.

## Basic usage

First make sure have the necessary input file:

    $ cd ../antarctica
    $ sh preprocess.sh
    $ cd ../picop

Then do

    $ ./pico_vs_picop.sh

Now you can compare `pico_basal_melt_rate` and picop_basal_melt_rate`:

    $ ncview spatial_pico_1s.nc spatial_picop_1s.nc &
