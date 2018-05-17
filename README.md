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

PISM is jointly developed at the [University of Alaska, Fairbanks (UAF)](http://www.uaf.edu/) and the [Potsdam Institute for Climate Impact Research (PIK)](http://www.pik-potsdam.de/).  UAF developers are based in the [Glaciers Group](http://www.gi.alaska.edu/snowice/glaciers/) at the [Geophysical Institute](http://www.gi.alaska.edu).

PISM development is supported by the [NASA Modeling, Analysis, and Prediction program](http://map.nasa.gov/) (grant #NNX13AM16G) and the [NASA Cryospheric Sciences program](http://ice.nasa.gov/) (grant NNX16AQ40G) and by [NSF Polar Programs](https://nsf.gov/geo/plr/about.jsp) grants PLR-1603799 and PLR-1644277.


PICO, a model for sub-shelf melt rates 
================================

The **pik/pico_dev** branch is based on the development version [b659de3](https://github.com/pism/pism/commit/80896b36f7444f78923a12d1c57cea47e30a6b08) from Feb 13th, 2017 (after the v0.7 release but before v1.0) 
and adds an implementation of the Potsdam Ice-Shelf Cavity mOdel 

Code release `pik-pico-07` was used for Reese et al., The Cryosphere (2018), forthcoming.

For further use of PICO, please refer to the latest version on https://github.com/pism/pism.git.
If you make use of PICO, please cite the respective paper and Olbers & Hellmer (2010).


Homepage
--------

[www.pism-docs.org](http://www.pism-docs.org/)


Download and Install
--------------------

See [instructions for getting the latest release](http://www.pism-docs.org/wiki/doku.php?id=stable_version).

License
--------------------

PISM is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License, or (at your option) any later version.

PISM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html) for more details.

You find a copy of the GNU General Public License along with PISM in the file COPYING;



Generating Documentation
------------------------

See the [INSTALL.md](INSTALL.md) file in this directory.

Contributing
------------

Want to contribute? Great! See [Committing to PISM](http://www.pism-docs.org/wiki/doku.php?id=committing).
