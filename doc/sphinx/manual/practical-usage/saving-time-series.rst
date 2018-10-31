.. include:: ../../global.txt

.. _sec-saving-time-series:

Saving time series of scalar diagnostic quantities
--------------------------------------------------

It is also possible to save time-series of certain scalar diagnostic quantities using a
combination of the options ``-ts_file``, ``-ts_times``, and ``-ts_vars``. For example,

.. code-block:: none

   pismr -i foo.nc -y 1e4 -o output.nc -ts_file time-series.nc \
         -ts_times 0:1:1e4 -ts_vars ice_volume_glacierized,ice_area_glacierized_grounded


will run for 10000 years, saving total ice volume and grounded ice area to
``time-series.nc`` *yearly*. See tables :numref:`tab-time-series-opts` for the list of
options and :ref:`sec-ts_vars` for the full list of supported time-series.

Note that, similarly to the snapshot-saving code (section :ref:`sec-snapshots`), this
mechanism does not affect adaptive time-stepping. Here, however, PISM will save exactly
the number of time-series records requested.

Omitting the ``-ts_vars`` option makes PISM save *all* available variables listed in
:ref:`sec-ts_vars`. Because scalar time-series take minimal storage space, compared to
spatially-varying data, this is usually a reasonable choice. Run PISM with the
:opt:`-list_diagnostics` option to see the list of all available time-series.

If the file ``foo.nc``, specified by ``-ts_file foo.nc``, already exists then by default
the existing file will be moved to ``foo.nc~`` and the new time series will go into
``foo.nc``. To append the time series onto the end of the existing file, use option
``-ts_append``.

PISM buffers time-series data and writes it at the end of the run, once 10000 values are
stored, or when an ``-extra_file`` is saved, whichever comes first. Sending an ``USR1``
(or ``USR2``) signal to a PISM process flushes these buffers, making it possible to
monitor the run. (See section :ref:`sec-signal` for more about PISM's signal handling.)

.. list-table:: Command-line options controlling saving scalar time-series
   :name: tab-time-series-opts
   :header-rows: 1
   :widths: 1,3

   * - Option
     - Description

   * - :opt:`-ts_file`
     - Specifies the file to save to.

   * - :opt:`-ts_times`
     - Specifies times to save at as a MATLAB-style range :math:`a:\Delta t:b`, a
       comma-separated list, or a keyword (``hourly``, ``daily``, ``monthly``,
       ``yearly``). See section :ref:`sec-saving-diagnostics`.

   * - :opt:`-ts_vars`
     - Comma-separated list of variables. Omitting this option is equivalent to listing
       the *all* variables.

   * - :opt:`-ts_append`
     - Append time series to file if it already exists. No effect if file does not yet
       exist.

Besides the above information on usage, here are comments on the physical significance of
several scalar diagnostics:

- For each variable named ``..._flux``, positive values mean ice sheet mass gain.

- PISM reports ice volume, ice mass, and several other quantities for "glacierized" areas.
  These quantities do not include contributions from areas where the ice thickness is
  equal to or below the value of the configuration parameter
  ``output.ice_free_thickness_standard`` (in meters). Corresponding quantities without the
  suffix *do* include areas with a thin, "seasonal" ice cover.

- Ice volume and area are computed and then split among floating and grounded portions:
  ``ice_volume_glacierized`` :math:`\mapsto` (``ice_volume_glacierized_shelf``,
  ``ice_volume_glacierized_grounded``) while ``ice_area_glacierized`` :math:`\mapsto`
  (``ice_area_glacierized_shelf``, ``ice_area_glacierized_grounded``). The volumes have units
  :math:`m^3` and the areas have units :math:`m^2`.

- The thermodynamic state of the ice sheet can be assessed, in part, by the amount of cold
  or temperate ("``temp``") ice. Thus there is another splitting: ``ice_volume_glacierized``
  :math:`\mapsto` (``ice_volume_glacierized_cold``, ``ice_volume_glacierized_temperate``) and
  ``ice_area_glacierized`` :math:`\mapsto`
  (``ice_area_glacierized_cold_base``, ``ice_area_glacierized_temperate_base``).

- The sea level rise potential :var:`sea_level_rise_potential` is the increase in sea
  level (in meters) that would result from melting all the grounded ice not displacing sea
  water and distributing the corresponding *fresh water* volume uniformly over the entire
  global ocean (:math:`362.5 \cdot 10^6\, km^2`, see :cite:`Cogley2011`). This follows the
  definition used in the SeaRISE project :cite:`Bindschadler2013SeaRISE`.

- Fields ``max_diffusivity`` and ``max_hor_vel`` relate to PISM time-stepping. These
  quantities appear in per-time-step form in the standard output from PISM (i.e. at
  default verbosity). ``max_diffusivity`` determines the length of the mass continuity
  sub-steps for the SIA stress balance (sub-)model. ``max_hor_vel`` determines the
  CFL-type restriction for mass continuity and conservation of energy contributions of the
  SSA stress balance (i.e. sliding) velocity.
