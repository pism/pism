PISM, a Parallel Ice Sheet Model
================================
|doi| |gpl| |cipism| |fairchecklist|


.. image:: https://github.com/pism/pism/blob/dev/images/Greenland_RCP_85_2008_2300_comp_1080p30.gif


The Parallel Ice Sheet Model is an open source, parallel, high-resolution ice sheet model that includes:

- A hierarchy of available stress balances
- Marine ice sheet physics, dynamic calving fronts
- A polythermal, enthalpy-based conservation of energy scheme
- Extensible coupling to atmospheric and ocean models
- Verification and validation tools
- `Documentation <pism-manual_>`_ for users and developers
- Links to MPI_ and PETSc_ for parallel simulations
- Use of `CF-compliant <cf_>`_  NetCDF_ files for input and output

PISM is jointly developed at the `University of Alaska, Fairbanks (UAF) <uaf_>`_ and the
`Potsdam Institute for Climate Impact Research (PIK) <pik_>`_. UAF developers are based in
the `Snow, Ice, and Permafrost Group <sip_>`_ at the `Geophysical Institute <gi_>`_.

Please see PISM's manual https://www.pism.io/docs/citing for how to acknowledge the use of
PISM and PISM's funding sources. If the manual is not available, please refer to the file
``doc/sphinx/citing/index.rst`` in the PISM repository.

Homepage
--------

    http://www.pism.io/

Download and Install
--------------------

See the section `Installing PISM <pism-installation_>`_ on ``pism.io``.

Support
-------

Please e-mail `uaf-pism@alaska.edu <uaf-pism_>`_ with questions about PISM.

You can also join the PISM workspace on `Slack <Slack-PISM_>`_.

Contributing
------------

Want to contribute? Great! See `Contributing to PISM <pism-contributing_>`_.

.. URLs

.. |fairchecklist| image:: https://fairsoftwarechecklist.net/badge.svg
.. _fairchecklist: https://fairsoftwarechecklist.net/v0.2?f=31&a=32113&i=31311&r=123
.. |doi| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1199019.svg
.. _doi: https://doi.org/10.5281/zenodo.1199019
.. |gpl| image:: https://img.shields.io/badge/License-GPL-green.svg
.. _gpl: https://opensource.org/licenses/GPL-3.0
.. |cipism| image:: https://circleci.com/gh/pism/pism/tree/dev.svg?style=svg
.. _cipism: https://circleci.com/gh/pism/pism/tree/dev
.. _uaf: http://www.uaf.edu/
.. _pik: http://www.pik-potsdam.de/
.. _pism-manual: http://www.pism.io/docs
.. _pism-contributing: http://www.pism.io/docs/contributing
.. _pism-installation: http://www.pism.io/docs/installation
.. _mpi: http://www.mcs.anl.gov/research/projects/mpi/
.. _petsc: https://petsc.org/
.. _cf: http://cf-pcmdi.llnl.gov/
.. _netcdf: http://www.unidata.ucar.edu/software/netcdf/
.. _sip: https://www.gi.alaska.edu/research/snow-ice-and-permafrost
.. _gi: http://www.gi.alaska.edu
.. _Slack-PISM: https://uaf-pism.slack.com/join/shared_invite/enQtODc3Njc1ODg0ODM5LThmOTEyNjEwN2I3ZTU4YTc5OGFhNGMzOWQ1ZmUzMWUwZDAyMzRlMzBhZDg1NDY5MmQ1YWFjNDU4MDZiNTk3YmE
.. _uaf-pism: mailto:uaf-pism@alaska.edu

..
   Local Variables:
   fill-column: 90
   End:
