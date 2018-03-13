PISM, a Parallel Ice Sheet Model
================================

The Parallel Ice Sheet Model is an open source, parallel, high-resolution ice sheet model:

* hierarchy of available stress balances
* marine ice sheet physics, dynamic calving fronts
* polythermal, enthalpy-based conservation of energy scheme
* extensible coupling to atmospheric and ocean models
* verification and validation tools
* complete [documentation](http://www.pism-docs.org/) for users and developers
* uses [MPI](http://www-unix.mcs.anl.gov/mpi/) and [PETSc](http://www-unix.mcs.anl.gov/petsc/petsc-as/) for parallel simulations
* reads and writes [CF-compliant](http://cf-pcmdi.llnl.gov/) [NetCDF](http://www.unidata.ucar.edu/software/netcdf/) files



The **pik/paleo_07dev** branch is based on the development version [e9d2d1f8](https://github.com/pism/pism/commit/e9d2d1f8b5cba9d0fc47d13753d838aa6b49bf01) from March 7th, 2017 
and adds methods that are used in paleo simulations of the Antarctic Ice Sheet

* sub-ice shelf melting calculations with [PICO](https://www.the-cryosphere-discuss.net/tc-2017-70/) 
* iterative optimization of till-friction angle to present-day grounded surface elevation
* air temperature paramterization based on a multiregression fit to [ERA-Interim](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim) reanalysis data
kill ocean
* non-linear lapse rate scaling for precipitation
* reads ocean_kill_mask from file to constrain maximum ice extent
* fix of the calculation of sea-level relevant volume in the output timeseries

You find in the examples/paleo-antarctica folder a working example of a paleo spin-up using all added functionality. 






PISM is jointly developed at the [University of Alaska, Fairbanks (UAF)](http://www.uaf.edu/) and the [Potsdam Institute for Climate Impact Research (PIK)](http://www.pik-potsdam.de/).  UAF developers are based in the [Glaciers Group](http://www.gi.alaska.edu/snowice/glaciers/) at the [Geophysical Institute](http://www.gi.alaska.edu).

PISM development is supported by the [NASA Modeling, Analysis, and Prediction program](http://map.nasa.gov/) (grant #NNX13AM16G) and the [NASA Cryospheric Sciences program](http://ice.nasa.gov/) (grant NNX16AQ40G) and by [NSF Polar Programs](https://nsf.gov/geo/plr/about.jsp) grants PLR-1603799 and PLR-1644277.


Homepage
--------

[www.pism-docs.org](http://www.pism-docs.org/)


Download and Install
--------------------

See [instructions for getting the latest release](http://www.pism-docs.org/wiki/doku.php?id=stable_version).


Generating Documentation
------------------------

See the [INSTALL.md](INSTALL.md) file in this directory.

Contributing
------------

Want to contribute? Great! See [Committing to PISM](http://www.pism-docs.org/wiki/doku.php?id=committing).
