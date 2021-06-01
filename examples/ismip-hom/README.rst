.. default-role:: literal

This directory contains scripts that run ISMIP-HOM_ experiments A-D (subdirectory `abcd`)
and experiment E (subdirectory `e-arolla`) and plot resulting ice velocities.

Run

.. code::

   make -C abcd all

to download and unpack inputs from the ISMIP-HOM_ supplement, run PISM, and produce ice
velocity figures for experiments A, B, C, and D.

Run

.. code::

   make -C e-arolla all

to download inputs, run PISM, and produce figures for experiment E.

.. _ISMIP-HOM: https://tc.copernicus.org/articles/2/95/2008/
