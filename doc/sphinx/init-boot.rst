.. _sec-initboot:

Initialization and bootstrapping
================================

There are three ways to start PISM:

- option :opt:`-i` reads a previously-saved "complete" PISM model state in a NetCDF file, or
- option :opt:`-i` :opt:`-bootstrap` reads an "incomplete" NetCDF file and uses heuristics to fill in needed fields, or
- one of the executables ``pisms`` or ``pismv`` is used to initialize simplified-geometry experiments or verification tests from formulas in the source code, and thus no input file is required.

One of the first two choices is required when using the executable ``pismr``.  Modeling usually starts with the ``-i input.nc -bootstrap`` option because real ice sheet observations are never complete initial conditions.  Runs with multiple stages often use the :opt:`-i` option after the first stage.

Initialization from a saved model state
---------------------------------------

"Initialization" has the specific, simple meaning in PISM that option ":opt:`-i`" was used.  If a previous PISM run has saved a NetCDF file using ":opt:`-o`" then that file will contain complete initial conditions for continuing the run.  The output file from the last run can be loaded with ":opt:`-i`": }

.. code-block:: none

   pisms -eisII A -y 100 -o foo.nc
   pisms -eisII A -i foo.nc -y 100 -o bar.nc

As noted verification tests (section :ref:`sec-verif`) and simplified-geometry experiments (section :ref:`sec-simp`) do not need input files at all because they initialize from formulas in the source code.  They can, however, be continued from saved model states using :opt:`-i`.  Specifying the simplified geometry experiment or verification test *is*, however, necessary if the run is to continue with the climate inputs for that experiment or test.  For example, based on the above ``pisms`` runs, it is valid to do

.. code-block:: none

   pismr -i foo.nc -y 100 -o bar.nc

but the climate and other parameters use PISM default values, and thus are not (necessarily) the values specified in EISMINT II.

As a technical note about saved states, a PISM run with :opt:`-stress_balance ssa` also saves the last SSA velocities to the output file in variables :var:`u_ssa` and :var:`v_ssa`.  The presence of these velocities adds efficiency in restarting because an initial estimate speeds up the solution of the SSA stress balance equations.  If you want to use :opt:`-i` but also ignore these velocities then use option :opt:`-dontreadSSAvels`.

.. _sec-i-format:

``-i`` file format
^^^^^^^^^^^^^^^^^^

PISM produces CF-1.5 compliant NetCDF files.  The easiest way to learn the output format *and* the :opt:`-i` format is to do a simple run and then look at the metadata in the resulting file, like this:

.. code-block:: none

   pisms -eisII A -y 10 -o foo.nc
   ncdump -h foo.nc | less


Note that variables in the output file have a ``pism_intent`` attribute} attribute.  When ``pism_intent`` = ``diagnostic``, the variable can be deleted from the file without affecting whether PISM can use it as a :opt:`-i` input file.  Variables with ``pism_intent`` = ``model_state``, by contrast, must be present when using :opt:`-i`.

The automatically-produced :var:`time` variable has a ``units`` attribute like ``"seconds since 1-1-1"`` because the CF metadata conventions require a reference date.  By default PISM ignores this reference date except when it is used in unit conversions based on a calendar (see below).

.. _sec-bootstrapping:

Bootstrapping
-------------

"Bootstrapping" in PISM means starting a modeling run with less than sufficient data, and letting essentially heuristic models fill in needed fields.  These heuristics are applied before the first time step is taken, so they are part of an initialization process.  Bootstrapping uses the option :opt:`-bootstrap`; see subsection :ref:`sec-runscript` for an example.

The need for an identified stage like "bootstrapping" comes from the fact that initial conditions for the evolution equations describing an ice sheet are not all observable. As a principal example of this problem, these initial conditions include the temperature within the ice. Glaciological observations, specifically remote-sensed observations which cover a large fraction or all of an ice sheet, never include this temperature field in practice. Thus ice sheet modelling often does something like this to get "reasonable" initial fields within the ice:

#. start only with (potentially) observable quantities like surface elevation, ice thickness, ice surface temperature, surface mass balance, and geothermal flux,
#. "bootstrap" as defined here, using heuristics to fill in temperatures at depth and to give a preliminary estimate of the basal sliding condition and the three-dimensional velocity field, and
#. #. *either* do a long run, often holding the current geometry and surface conditions steady, to evolve toward a steady state which has compatible temperature, stress, and velocity fields,
   #. *or* do a long run using an additional (typically spatially-imprecise) historical record from an ice core or a sea bed core (or both), to apply forcing to the surface temperature or sea level (for instance), but with the same functional result of filling in temperature, stress, and velocity fields.
      
When using :opt:`-bootstrap` you will need to specify both grid dimensions (using :opt:`-Mx`, :opt:`-My` and :opt:`-Mz`; see subsection :ref:`sec-grid`) and the height of the computational box for the ice with :opt:`-Lz` (subsection :ref:`sec-coords`). The data read from the file can determine the horizontal extent of the model, if options :opt:`-Lx`, :opt:`-Ly` are not set. The additional required specification of vertical extent by :opt:`-Lz` is reasonably natural because typical data used in "bootstrapping" are two-dimensional. Using :opt:`-bootstrap` without specifying all four options :opt:`-Mx`, :opt:`-My`, :opt:`-Mz`, :opt:`-Lz` is an error.

If :opt:`-Lx` and :opt:`-Ly` specify horizontal grid dimensions smaller than in the bootstrapping file, PISM will cut out the center portion of the domain. Alternatively, options :opt:`-x_range` and :opt:`-y_range` each take a list of two numbers, a list of minimum and maximum :math:`x` and :math:`y` coordinates, respectively (in meters), which makes it possible to select a subset that is not centered in the bootstrapping file's grid.

For the key issue of what heuristic is used to determine the temperatures at depth, there are two methods. The default method uses ice thickness, surface temperature, surface mass balance, and geothermal flux. The temperature is set to the solution of a steady one-dimensional differential equation in which conduction and vertical advection are in balance, and the vertical velocity linearly-interpolates between the surface mass balance rate at the top and zero at the bottom. The non-default method, set with option :opt:`-boot_temperature_heuristic quartic_guess`, was the default in older PISM versions (``stable0.5`` and earlier); it does not use the surface mass balance and instead makes a more-heuristic estimate of the vertical temperature profile based only on the ice thickness, surface temperature, and geothermal flux.

.. _sec-bootstrapping-format:

``-bootstrap`` file format
^^^^^^^^^^^^^^^^^^^^^^^^^^

Allowed formats for a bootstrapping file are relatively simple to describe. 

#. NetCDF variables should have the ``units`` containing a UDUNITS-2-compatible string. If this attribute is missing, PISM will assume that a field uses MKS units. [#]_
#. NetCDF coordinate variables should have ``standard_name`` or ``axis`` attributes. These are used to determine which *spatial* dimension a NetCDF dimension corresponds to; for example see ``ncdump -h`` output from a file produced by PISM. The :var:`x` and :var:`y` dimensions need not be called ":var:`x`" and ":var:`y`".
#. Coordinate variables have to be strictly-increasing.
#. Three-dimensional variables will be ignored in bootstrapping.
#. The ``standard_name`` attribute is used, when available, to identify a variable, so variable names need not match corresponding variables in a PISM output file. See the `PISM Source Code browser <pism-code-browser_>` for a list of CF standard names used in PISM. Specifically, the bed elevation (topography) is read by ``standard_name`` = ``bedrock_altitude`` and the ice thickness by ``standard_name`` = ``land_ice_thickness``.
#. Any two-dimensional variable except bed topography and ice thickness may be missing. For missing variables some heuristic will be applied. See table :numref:`tab-modelhierarchy` for a sketch of the data necessary for bootstrapping; see ``src/base/iMbootstrap.cc`` for all further details.
#. Surface elevation is ignored if present. Users with surface elevation and bed elevation data should compute the ice thickness variable, put it in the bootstrapping file, and set its ``standard_name`` to ``land_ice_thickness``.

.. [#] PISM uses a library called UDUNITS-2 to convert data present in an input file to MKS. This means that having ice thickness in feet or temperature in Fahrenheit *is* allowed.
