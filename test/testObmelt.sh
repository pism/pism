#!/bin/bash

# test script to evaluate bmelt calculation
 
pism -test O -Mx 4 -My 4 -Mz 41 -Mbz 11 -Lbz 1000 -y 0 -verbose 3 -o foo.nc


# this combination "gets it right", i.e. bmelt = 0.00088055 m/a, but only because it
#   is one timestep; see bmelt var in bar.nc
pism -i foo.nc -no_mass -no_bmr_in_cont -y 1 -verbose 3 -o bar.nc


#gets it wrong; see bmelt var in baz.nc
pism -i foo.nc -no_mass -no_bmr_in_cont -max_dt 1 -y 10 -verbose 3 -o baz.nc
