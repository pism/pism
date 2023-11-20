Changes from 0.5 (May 2012) to 0.6
==================================

Basal strength and basal hydrology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  A significant change to the model physics for till: effective
   pressure in till is now computed by empirically-based (and
   exponential) formula from [``Tulaczyketal2000``\ ] relating amount in
   till to effective pressure.
-  New mass-conserving subglacial hydrology model.

Marine ice sheet modeling
~~~~~~~~~~~~~~~~~~~~~~~~~

-  Implement a parameterization of the melange back pressure as a part
   of the stress boundary condition at the ice shelf front.
-  Implement fracture density model and fracture-induced ice softening
   (see [``AlbrechtLevermann2012``\ ]).
-  Parameterization of the sub-grid position of the grounding line and
   related improvements (see
   [``Feldmannetal2014``,\ ``Gladstoneetal2012``]).
-  MISMIP3D example (see [``MISMIP3d2013``\ ]), grounding line
   reversibility.
-  Fixes related to the handling of the grounding line (i.e. driving
   stress, etc)
-  Use an implementation of a serial two-scan connected-component
   algorithm to identify icebergs.
-  Report cumulative discharge (calving) flux as a 2D field.

Climate inputs and ocean inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``PISMSurfaceModel`` returns climatic mass balance in
   ``[kg m-2 s-1]``
-  Remove ``pclimate``; add the ``-test_climate_models`` option.
-  ``-part_grid`` extended to grounded marine termini.
-  Require time bounds in time-dependent climate forcing.
-  Improve the PDD scheme. (Keeps track of the snow depth during the
   year to choose PDD factors but does not model multi-year snow cover.)
-  Implement the 3-equation sub-shelf melt parameterization (see
   [``Hellmer1998``\ ]).
-  Implement caching surface and ocean models ``-ocean ...,cache`` and
   ``-surface ...,cache``. (Update boundary inputs every X years but
   include interannual variability.)
-  Remove ``EISMINT-Greenland``.

Inverse modeling
~~~~~~~~~~~~~~~~

-  Inverse modeling tools *are* a part of the release now.

Energy
~~~~~~

-  New bootstrapping heuristic filling ice temperature at depth from
   surface mass balance (if available).
-  Allow "regridding" enthalpy from files containing ice temperature and
   liquid water fraction.
-  Allow cold flow laws in the enthalpy mode and GPBLD in the cold mode
   (same as Paterson-Budd).
-  Corrected basal boundary condition in the enthalpy system.

Usability
~~~~~~~~~

-  Implement poor man's parallel I/O, with compression
-  Ensure "continuity" in time of reported cumulative diagnostic fields.
-  Let the user precisely specify the dates corresponding to the run
   (``-time_file``).
-  Many more diagnostic quantities; use ``-list_diagnostics`` to see.
-  PISM keeps track of options and parameters that were set but were not
   used.
-  Use projection info to compute latitude and longitude bounds.
   (Reduces post-processing needed to work with PISM's output.)

Miscellaneous
~~~~~~~~~~~~~

-  Improve mass conservation and mass transport.
-  New validation example based on laboratory experiment with xanthan
   gum. Uses millimeter grid spacing and shows configurability of PISM.
-  Implement a "constant 2D velocity" stress balance object (for testing
   and prescribing sliding)

Under the hood
~~~~~~~~~~~~~~

-  Improve the build system.
-  Switch to PETSC 3.3 or 3.4, stop supporting 3.2.
-  Switch to
   `UDUNITS-2 <http://www.unidata.ucar.edu/software/udunits/>`__.
-  Require `FFTW-3 <http://www.fftw.org>`__ to build PISM.
-  Use `CalCalcs <http://meteora.ucsd.edu/~pierce/calcalcs/>`__ for
   proper calendar support.
-  Remove the Debian [meta-] package.
-  Clean up command-line options selecting sub-models (calving, stress
   balance, energy, basal yield stress).
-  Better ``SSAFD`` recovery logic; save model state on failure.
-  Use non-zero initial guess in the SSAFD KSP solver.
-  Improved basal yield stress code.
