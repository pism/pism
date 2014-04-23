Re-deriving details of the three-equation model
===============================================

This directory contains Maxima scripts re-deriving all the formulas
necessary to implement the three-equation parameterization of the
thermodynamic ice-ocean interaction at the ice shelf bottom.

There are three cases:
 1. Basal melt.
 2. Basal freeze-on.
 3. Neither melt nor freeze-on at the base.

Run
    make melt
    make freeze-on
    make diffusion

to produce coefficients of the quadratic equation for the salinity at the shelf base.

Once basal salinity is computed, a linearized melting point equation
is used to compute the basal temperature, and the sal flux equation is
used to compute the sub-shelf melt rate.
