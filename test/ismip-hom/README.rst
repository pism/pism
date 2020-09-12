.. default-role:: literal

This directory contains scripts that run ISMIP-HOM_ experiments A-D and plot results.

In short, download and unpack the ISMIP-HOM_ supplement, then run:

.. code::

   python3 run-ismiphom.py

   ln -s path/to/ismip_all .

   python3 convert-ismiphom.py
   python3 plot-ismiphom.py

The script `run-ismiphom.py` uses PISM's Python bindings to run the Blatter solver in
PISM. See the top of this script for details.

The script `convert-ismiphom.py` reads submitted model results (see the supplement to
ISMIP-HOM_), samples them along the line (x, 0.25), and saves to files. To run this
script, place the `ismip_all` directory (or a symlink to it) in this directory.

The script `plot-ismiphom.py` uses Bokeh_ and data processed by `convert-ismiphom.py`.

.. _Bokeh: https://bokeh.org/_
.. _ISMIP-HOM: https://tc.copernicus.org/articles/2/95/2008/
