PISM, a Parallel Ice Sheet Model
================================

The Parallel Ice Sheet Model is an open source, parallel, high-resolution ice sheet model:

- hierarchy of available stress balances
- marine ice sheet physics, dynamic calving fronts
- polythermal, enthalpy-based conservation of energy scheme
- extensible coupling to atmospheric and ocean models
- verification and validation tools
- `documentation <pism-docs_>`_ for users and developers
- uses MPI_ and PETSc_ for parallel simulations
- reads and writes `CF-compliant <cf_>`_  NetCDF_ files

PISM is jointly developed at the `University of Alaska, Fairbanks (UAF) <uaf_>`_ and the
`Potsdam Institute for Climate Impact Research (PIK) <pik_>`_. UAF developers are based in
the `Glaciers Group <glaciers_>`_ at the `Geophysical Institute <gi_>`_.

PISM development is supported by the `NASA Modeling, Analysis, and Prediction program
<NASA-MAP_>`_ (grant #NNX13AM16G) and the `NASA Cryospheric Sciences program
<NASA-Cryosphere_>`_ (grant NNX16AQ40G) and by `NSF Polar Programs <NSF-Polar_>`_ grants
PLR-1603799 and PLR-1644277.

Homepage
--------

    http://www.pism-docs.org/

Download and Install
--------------------

See `instructions for getting the latest release <pism-stable_>`_.

Generating Documentation
------------------------

See the `INSTALL.rst <INSTALL.rst>`_ file in this directory.

Contributing
------------

Want to contribute? Great! See `Committing to PISM <pism-contribute_>`_.

.. URLs

.. _uaf: http://www.uaf.edu/
.. _pik: http://www.pik-potsdam.de/
.. _pism-docs: http://www.pism-docs.org/
.. _pism-stable: http://www.pism-docs.org/wiki/doku.php?id=stable_version
.. _pism-contribute: http://www.pism-docs.org/wiki/doku.php?id=committing
.. _mpi: http://www.mcs.anl.gov/research/projects/mpi/
.. _petsc: http://www.mcs.anl.gov/petsc/
.. _cf: http://cf-pcmdi.llnl.gov/
.. _netcdf: http://www.unidata.ucar.edu/software/netcdf/
.. _glaciers: http://www.gi.alaska.edu/snowice/glaciers/
.. _gi: http://www.gi.alaska.edu
.. _NASA-MAP: http://map.nasa.gov/
.. _NASA-Cryosphere: http://ice.nasa.gov/
.. _NSF-Polar: https://nsf.gov/geo/plr/about.jsp

..
   Local Variables:
   fill-column: 90
   End:
