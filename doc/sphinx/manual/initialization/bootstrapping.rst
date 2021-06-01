.. include:: ../../global.txt

.. _sec-bootstrapping:

Bootstrapping
-------------

"Bootstrapping" in PISM means starting a modeling run with less than sufficient data, and
then either

- interpolating some of the missing fields from a separate file, and
- letting essentially heuristic models fill in remaining ones.

So, "bootstrapping" is used whenever some fields are missing *or* interpolation is
necessary, for example when going to a finer grid during grid sequencing.

These steps are performed before the first time step is taken, so they are part of an
initialization process. Bootstrapping uses the option :opt:`-bootstrap`; see section
:ref:`sec-runscript` for an example.

The need for an identified stage like "bootstrapping" comes from the fact that initial
conditions for the evolution equations describing an ice sheet are not all observable. As
a principal example of this problem, these initial conditions include the temperature
within the ice. Glaciological observations, specifically remote-sensed observations which
cover a large fraction or all of an ice sheet, never include this temperature field in
practice.

Ice sheet models often need to do something like this to get "reasonable" initial fields
within the ice:

#. start only with (potentially) observable quantities like surface elevation, ice
   thickness, ice surface temperature, surface mass balance, and geothermal flux,
#. "bootstrap" as defined here, using heuristics to fill in temperatures at depth and to
   give a preliminary estimate of the basal sliding condition and the three-dimensional
   velocity field, and
#. #. *either* do a long run, often holding the current geometry and surface conditions
      steady, to evolve toward a steady state which has compatible temperature, stress,
      and velocity fields,
   #. *or* do a long run using an additional (typically spatially-imprecise) historical
      record from an ice core or a sea bed core (or both), to apply forcing to the surface
      temperature or sea level (for instance), but with the same functional result of
      filling in temperature, stress, and velocity fields.
      
When using :opt:`-bootstrap` you will need to specify both grid dimensions (using
:opt:`-Mx`, :opt:`-My` and :opt:`-Mz`; see section :ref:`sec-grid`) and the height of the
computational box for the ice with :opt:`-Lz` (section :ref:`sec-coords`). The data read
from the file can determine the horizontal extent of the model, if options :opt:`-Lx`,
:opt:`-Ly` are not set. The additional required specification of vertical extent by
:opt:`-Lz` is reasonably natural because input data used in "bootstrapping" are
two-dimensional. Using :opt:`-bootstrap` without specifying all four options :opt:`-Mx`,
:opt:`-My`, :opt:`-Mz`, :opt:`-Lz` is an error.

If :opt:`-Lx` and :opt:`-Ly` specify horizontal grid dimensions smaller than in the
bootstrapping file, PISM will cut out the center portion of the domain. In PISM's regional
mode, options :opt:`-x_range` and :opt:`-y_range` each take a list of two numbers, a list
of minimum and maximum `x` and `y` coordinates, respectively (in meters), which makes it
possible to select a subset that is not centered in the bootstrapping file's grid.

For the key issue of what heuristic is used to determine the temperatures at depth, there
are two methods. The default method uses ice thickness, surface temperature, surface mass
balance, and geothermal flux. The temperature is set to the solution of a steady
one-dimensional differential equation in which conduction and vertical advection are in
balance, and the vertical velocity linearly-interpolates between the surface mass balance
rate at the top and zero at the bottom. The non-default method, selected by setting
:config:`bootstrapping.temperature_heuristic` to ``quartic_guess``, was the default in
older PISM versions (``stable0.5`` and earlier); it does not use the surface mass balance
and instead makes a more-heuristic estimate of the vertical temperature profile based only
on the ice thickness, surface temperature, and geothermal flux.

.. _sec-bootstrapping-format:

``-bootstrap`` file format
^^^^^^^^^^^^^^^^^^^^^^^^^^

Allowed formats for a bootstrapping file are relatively simple to describe. 

#. NetCDF variables should have the ``units`` containing a UDUNITS_\-compatible string. If
   this attribute is missing, PISM will assume that a field uses MKS units.\ [#]_
#. NetCDF coordinate variables should have ``standard_name`` or ``axis`` attributes. These
   are used to determine which *spatial* dimension a NetCDF dimension corresponds to; for
   example see ``ncdump -h`` output from a file produced by PISM. The :var:`x` and
   :var:`y` dimensions need not be called ":var:`x`" and ":var:`y`".
#. Coordinate variables have to be strictly-increasing.
#. Three-dimensional variables will be ignored in bootstrapping.
#. The ``standard_name`` attribute is used, when available, to identify a variable, so
   variable names need not match corresponding variables in a PISM output file. See the
   :ref:`sec-cf-standard-names` for a list of CF standard names used in
   PISM.

   For example, the bed elevation (topography) is read by ``standard_name`` =
   ``bedrock_altitude`` and the ice thickness by ``standard_name`` =
   ``land_ice_thickness``.
#. Any two-dimensional variable except bed topography and ice thickness may be missing.
   For missing variables some heuristic will be applied. See :numref:`tab-modelhierarchy`
   for a sketch of the data necessary for bootstrapping.
#. Surface elevation is ignored if present. Users with surface elevation and bed elevation
   data should compute the ice thickness variable, put it in the bootstrapping file, and
   set its ``standard_name`` to ``land_ice_thickness``.

.. [#] PISM automatically converts data present in an input file to MKS. This means that
       having ice thickness in feet or temperature in Fahrenheit *is* allowed.
