Ross flow model example
=================

This directory contains `preprocess.py`, the script downloading and preprocessing input data, and `run.sh`, using PISM to do a "diagnostic" run.

Input datasets
==========

We use the two datasets:
* *An improved Antarctic dataset for high resolution numerical ice sheet models (ALBMAP v1)*  by A. M. Le Brocq, A. J. Payne, and A. Vieli (http://doi.pangaea.de/10.1594/PANGAEA.734145)
* *MEaSUREs InSAR-Based Antarctica Velocity Map*  by Rignot, E., J. Mouginot, and B. Scheuchl. 2011, (http://nsidc.org/data/nsidc-0484.html)

Fields used
========

* ice thickness (or the upper/lower surface elevation pair)
* ice surface temperature
* ice surface mass balance
* bedrock topography
* u and v components of the surface ice velocity

Notes
====

* periodic boundary conditions at boundaries of the modeling domain
* several fields have jumps at domain boundaries
* need boundary conditions at the ice shelf front (`-cfbc`)







