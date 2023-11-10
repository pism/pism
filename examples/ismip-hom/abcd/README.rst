To run ISMIP-HOM experiments A, B, C, and D using PISM, install PISM with Python bindings
and run

.. code::

   make

This will download submitted model results from the ISMIP-HOM_ paper supplement, run all
the experiments, and plot results. See the `ISMIP-HOM section of the manual
<ISMIP-HOM-PISM_>`_ for results.

The script `run-ismiphom.py` uses PISM's Python bindings to run the Blatter solver in
PISM. See the top of this script for details.

The script `convert-ismiphom.py` reads submitted model results (see the supplement to
ISMIP-HOM_), samples them along the line (x, 0.25), (in scaled coordinates) and saves to
files. To run this script, place the `ismip_all` directory (or a symlink to it) in this
directory.

The script `plot-ismiphom.py` uses Matplotlib and data processed by `convert-ismiphom.py`.
It produces figures used in the Manual.

.. _ISMIP-HOM: https://tc.copernicus.org/articles/2/95/2008/
.. _ISMIP-HOM-PISM: https://www.pism.io/docs/manual/simplified-geometry/ismip-hom.html
