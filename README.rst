PISM, a Parallel Ice Sheet Model
================================
|cipism|_

The Parallel Ice Sheet Model is an open source, parallel, high-resolution ice sheet model:

- hierarchy of available stress balances
- marine ice sheet physics, dynamic calving fronts
- polythermal, enthalpy-based conservation of energy scheme
- extensible coupling to atmospheric and ocean models
- verification and validation tools
- `documentation <pism-manual_>`_ for users and developers
- uses MPI_ and PETSc_ for parallel simulations
- reads and writes `CF-compliant <cf_>`_  NetCDF_ files

PISM is jointly developed at the `University of Alaska, Fairbanks (UAF) <uaf_>`_ and the
`Potsdam Institute for Climate Impact Research (PIK) <pik_>`_. UAF developers are based in
the `Glaciers Group <glaciers_>`_ at the `Geophysical Institute <gi_>`_.

Please see ``ACKNOWLEDGE.rst`` and ``doc/funding.csv`` for a list of grants supporting
PISM development.

Homepage
--------

    http://www.pism.io/

Download and Install
--------------------

See the `Installing PISM <pism-installation_>`_ on ``pism.io``.

Support
-------

Please e-mail `uaf-pism@alaska.edu <uaf-pism_>`_ with questions about PISM.

You can also join the PISM workspace on `Slack <Slack-PISM_>`_.

Contributing
------------

Want to contribute? Great! See `Contributing to PISM <pism-contributing_>`_.

.. URLs

.. |cipism| image:: https://circleci.com/gh/pism/pism/tree/master.svg?style=svg
.. _cipism: https://circleci.com/gh/pism/pism/tree/master
.. _uaf: http://www.uaf.edu/
.. _pik: http://www.pik-potsdam.de/
.. _pism-manual: http://www.pism.io/docs
.. _pism-contributing: http://www.pism.io/docs/contributing
.. _pism-installation: http://www.pism.io/docs/installation
.. _mpi: http://www.mcs.anl.gov/research/projects/mpi/
.. _petsc: http://www.mcs.anl.gov/petsc/
.. _cf: http://cf-pcmdi.llnl.gov/
.. _netcdf: http://www.unidata.ucar.edu/software/netcdf/
.. _glaciers: http://www.gi.alaska.edu/snowice/glaciers/
.. _gi: http://www.gi.alaska.edu
.. _NASA-MAP: http://map.nasa.gov/
.. _NASA-Cryosphere: http://ice.nasa.gov/
.. _NSF-Polar: https://nsf.gov/geo/plr/about.jsp
.. _Slack-PISM: https://join.slack.com/t/uaf-pism/shared_invite/enQtODc3Njc1ODg0ODM5LThmOTEyNjEwN2I3ZTU4YTc5OGFhNGMzOWQ1ZmUzMWUwZDAyMzRlMzBhZDg1NDY5MmQ1YWFjNDU4MDZiNTk3YmE
.. _uaf-pism: mailto:uaf-pism@alaska.edu

..
   Local Variables:
   fill-column: 90
   End:
