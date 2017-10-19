Using IceModelVec2T to load space- and time-dependent forcing fields
====================================================================

.. contents::

These notes document ``IceModelVec2T`` and implement a minimal ``pism::ocean::OceanModel``
subclass using ``IceModelVec2T`` to load and store space-and-time-dependent sub-shelf melt
rates.

Rationale
---------

The class ``IceModelVec2T`` was created to allow using time-dependent climate forcing
fields. Unlike three-dimensional fields used by PISM (ice enthalpy, bedrock temperature,
velocity components), climate forcing fields

- are not needed "all at once:" only the portion covering the current time step has to be
  available at a given time, and
- may cover very long time periods, making it impossible to hold the whole data set in RAM.
- may be *periodic* with a given period (whole years)

In addition to this, PDD (temperature index) surface process models use *time-series* of
near-surface air temperature and precipitation at each grid point, and ``IceModelVec2T``
supports this use.

``IceModelVec2T`` is derived from ``IceModelVec2S`` and can be used as a "regular" scalar
field in PISM *without* ghosts.

Setup
-----

Allocating storage
^^^^^^^^^^^^^^^^^^

``IceModelVec2T`` contains a "buffer" storing a number of temporal records of a 2D field.
To allocate an instance you need to specify the maximum number of records it can hold.
(Climate forcing code in PISM uses the configuration parameter
``climate_forcing.buffer_size``.)

.. code-block:: c++

   unsigned int N = 100; // hold at most N records
   IceModelVec2T field;
   field.set_n_records(N);
   field.create(grid, "netcdf_variable_name");
   field.set_attrs("intent", "long name", "units", "CF_standard_name");

When using periodic forcing data the buffer has to be big enough to contain *all* the
records covering the period. (Implemented climate forcing code gets the number of records
from an input file.)

The call ``field.create(...)`` allocates storage.

The call ``field.set_attrs()`` sets some metadata. PISM will convert data read from a file
from the units stored in the NetCDF variable attribute "units" into internal units given
in this call. The CF standard name (if present) can be used to find a variable in an input
file.

Initialization
^^^^^^^^^^^^^^

To specify the file to read data from, call ``field.init(filename)``:

.. code-block:: c++

   field.init("input_file.nc");

In the "regular" (non-periodic) case this reads the times stored in the file, preparing to
read array data when ``update()`` is called. In the periodic case, on the other hand, all
records are read in right away.

To set an ``IceModelVec2T`` to a constant, call

.. code-block:: c++

   field.init_constant(value);

This call initializes a "fake" time line, turning subsequent ``update()`` calls into
no-ops.

Actual use
----------

To make sure that records covering a time step from ``t`` to ``t + dt`` are available,
call

.. code-block:: c++

   field.update(t, dt);

This reads some records from a file while avoiding unnecessary I/O.

If ``update(...)`` calls go in sequence every records should be read only once, even if
some time steps overlap previous ones.

Records stored in an ``IceModelVec2T`` are interpreted as *piece-wise constant in time.*
Depending on the application this may require limiting time steps taken by PISM so that no
time step spans more than one of the temporal intervals a field is defined on.

To get the length of such a time step PISM can take at time ``t``, call

.. code-block:: c++

   double max_timestep = field.max_timestep(t);

To get a "snapshot" at a given time, call

.. code-block:: c++

   field.interp(t);

Sometimes (e.g. for precipitation) it makes sense to average over a time step. In this
case, use ``average(t, dt)``.

.. code-block:: c++

   field.average(t, dt);

This call uses the rectangle rule to approximate the average. The interval ``(t, t + dt)``
is split into ``N`` sub-intervals, where ``N`` depends on the length of the time step. The
number of sub-intervals per model year can be set by calling

.. code-block:: c++

   field.set_n_evaluations_per_year(N_per_year);

The configuration parameter ``climate_forcing.evaluations_per_year`` provides the default.

Extracting time-series of a field at a particular grid point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Initializing point-wise access to time-series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given an array ``ts`` of times in the interval ``(t, t + dt)``, the call

.. code-block:: c++

   field.init_interpolation(ts);

will prepare ``field`` for extracting time-series at grid point ``(i, j)``.

``IceModelVec2T`` will use constant extrapolation if some times in the array ``ts`` are
outside the interval given to the last call ``update(t, dt)``.

Accessing point-wise time-series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Just like all other classes derived from ``IceModelVec``, ``IceModelVec2T`` requires that
``begin_access()`` is called before point-wise access and ``end_access()`` after. (Use
``IceModelVec::AccessList`` to avoid doing this "by hand.")

.. code-block:: c++

    IceModelVec::AccessList list(field);

    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      std::vector<double> values;

      field.interp(i, j, values);
      
      // use time-series "values" at times "ts"
    }

Putting it all together to implement an "ocean model"
-----------------------------------------------------

Class declaration
^^^^^^^^^^^^^^^^^

.. only:: html

   Click :download:`here <code/Example.hh>` to download this file.

.. literalinclude:: code/Example.hh
   :language: c++

Class definition
^^^^^^^^^^^^^^^^

.. only:: html

   Click :download:`here <code/Example.cc>` to download this file.

.. literalinclude:: code/Example.cc
   :language: c++

All ocean models need to provide implementatons (``_impl(...)`` methods) corresponding to
the public API of ``pism::ocean::OceanModel``. (See ``src/coupler/PISMOcean.hh`` and note
that some have default implementations.)
