.. include:: ../../global.txt

.. _sec-ismip6:

ISMIP6 Greenland
----------------

Running ISMIP6-Greenland_ projections required implementing some additional sub-models as
well as several modifications needed to follow ISMIP6 conventions. This section describes
these modifications and explains how to use PISM to run ISMIP6 projections.

Top surface mass balance and temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the ``ismip6`` surface model to implement ISMIP6 surface mass balance forcing.

.. code-block:: bash

   pismr -surface ismip6 \
         -surface_ismip6_file climate_forcing.nc \
         -surface_ismip6_reference_file climate_forcing_reference.nc

Here ``climate_forcing.nc`` should contain time-dependent variables

- :var:`climatic_mass_balance_anomaly` (units: `kg / (m^2 s)`) and
- :var:`ice_surface_temp_anomaly` (units: *Kelvin*).

The file ``climate_forcing_reference.nc`` should contain time-independent (2D) variables

- :var:`climatic_mass_balance_reference` (units: `kg / (m^2 s)`),
- :var:`climatic_mass_balance_gradient` (units: `(kg / (m^2 s)) / m`),
- :var:`ice_surface_temp_reference` (units: *Kelvin*),
- :var:`ice_surface_temp_gradient` (units: *Kelvin / m*),
- :var:`surface_elevation` (units: *m*)

The surface mass balance is computed using the following formula:

.. code-block:: none

   SMB(x,y,t) = SMB_ref(x,y) + aSMB(x,y,t) + dSMBdz(x,y) * [h(x,y,t) - h_ref(x,y)]

.. _sec-ismip6-frontal-melt:

Frontal melt parameterization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the ``discharge_given`` frontal melt model to implement the ISMIP6 frontal melt
parameterization.

.. code-block:: bash

   pismr -frontal_melt discharge_given \
         -frontal_melt_discharge_given_file forcing.nc ...

The file ``forcing.nc`` has to contain variables :var:`theta_ocean` (potential temperature
of adjacent ocean, *degrees Celsius*) and :var:`subglacial_discharge` (water flux per unit
area of submerged ice front, `kg / (m^2 s)`).

These inputs are used in the frontal melt parameterization described in
:cite:`Rignotetal2016`:

.. math::

   q_m = (A\, h\, q_{sg}^{\alpha} + B)\, \theta^{\beta}

Here `q_m` is the frontal melt rate in *m/day*, `h` is the water depth at an ice front,
`\theta` in the *thermal forcing* and `A`, `B`, `\alpha`, `\beta` are model parameters.


Parameterized front retreat
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the `parameterized front retreat mechanism <ismip6-greenland_>`_ use the
:ref:`sec-prescribed-retreat` mechanism.

.. code-block:: bash

   pismr -front_retreat_file retreat_forcing.nc ...

The file ``retreat_forcing.nc`` should contain the variable
``land_ice_area_fraction_retreat`` which defines the maximum ice extent at a given time.

Mass losses resulting from applying this mechanism are reported as a part of
:var:`tendency_of_ice_amount_due_to_discharge` and related diagnostics (i.e. they are
*not* attributed to calving or frontal melt).


Output variables
^^^^^^^^^^^^^^^^

See :numref:`tab-ismip6-variables` for a list of variables requested by ISMIP6. Note that
they have names different from the ones listed in :ref:`sec-extra_vars` and use MKS units.
To reduce the amount of post-processing output files require PISM can follow these
conventions.

Setting :config:`output.ISMIP6` makes PISM save diagnostics using MKS units and recognize
ISMIP6 variable names.

To save *all* the diagnostics requested by ISMIP6 use the short-cut

.. code-block:: bash

   pismr -extra_vars ismip6 ...

The list of variables is stored in the configuration parameter
:config:`output.ISMIP6_extra_variables` and contains variables Greenland projections are
required to provide. (Add ``base,ligroundf`` to this list for Antarctic projections.)

To save all the time series supported by PISM, omit the ``-ts_vars`` option:

.. code-block:: bash

   pismr -ts_times TIMES -ts_file ts.nc

To save all variables requested by ISMIP6, use ``-ts_vars ismip6``:

.. code-block:: bash

   pismr -ts_times TIMES -ts_file ts.nc -ts_vars ismip6

.. list-table:: ISMIP6 variables
   :name: tab-ismip6-variables
   :header-rows: 1
   :widths: 1,1,3

   * - Variable
     - Units
     - Description

   * - ``lithk(x,y,t)``
     - m
     - Ice thickness

   * - ``orog(x,y,t)``
     - m
     - Surface elevation

   * - ``base(x,y,t)``
     - m
     - Base elevation

   * - ``topg(x,y,t)``
     - m
     - Bedrock elevation

   * - ``hfgeoubed(x,y)``
     - W m-2
     - Geothermal heat flux

   * - ``acabf(x,y,t)``
     - kg m-2 s-1
     - Surface mass balance flux

   * - ``libmassbfgr(x,y,t)``
     - kg m-2 s-1
     - Basal mass balance flux beneath grounded ice

   * - ``libmassbffl(x,y,t)``
     - kg m-2 s-1
     - Basal mass balance flux beneath floating ice

   * - ``dlithkdt(x,y,t)``
     - m s-1
     - Ice thickness imbalance

   * - ``xvelsurf(x,y,t)``
     - m s-1
     - Surface velocity in x

   * - ``yvelsurf(x,y,t)``
     - m s-1
     - Surface velocity in y

   * - ``zvelsurf(x,y,t)``
     - m s-1
     - Surface velocity in z

   * - ``xvelbase(x,y,t)``
     - m s-1
     - Basal velocity in x

   * - ``yvelbase(x,y,t)``
     - m s-1
     - Basal velocity in y

   * - ``zvelbase(x,y,t)``
     - m s-1
     - Basal velocity in z

   * - ``xvelmean(x,y,t)``
     - m s-1
     - Mean velocity in x

   * - ``yvelmean(x,y,t)``
     - m s-1
     - Mean velocity in y

   * - ``litemptop(x,y,t)``
     - K
     - Surface temperature

   * - ``litempbotgr(x,y,t)``
     - K
     - Basal temperature beneath grounded ice sheet

   * - ``litempbotfl(x,y,t)``
     - K
     - Basal temperature beneath floating ice shelf

   * - ``strbasemag(x,y,t)``
     - Pa
     - Basal drag

   * - ``licalvf(x,y,t)``
     - kg m-2 s-1
     - Calving flux

   * - ``lifmassbf(x,y,t)``
     - kg m-2 s-1
     - Ice front melt and calving flux

   * - ``ligroundf(x,y,t)``
     - kg m-2 s-1
     - Grounding line flux

   * - ``sftgif(x,y,t)``
     - 1
     - Land ice area fraction

   * - ``sftgrf(x,y,t)``
     - 1
     - Grounded ice sheet area fraction

   * - ``sftflf(x,y,t)``
     - 1
     - Floating ice sheet area fraction

   * - ``lim(t)``
     - kg
     - Total ice mass

   * - ``limnsw(t)``
     - kg
     - Mass above floatation

   * - ``iareagr(t)``
     - m^2
     - Grounded ice area

   * - ``iareafl(t)``
     - m^2
     - Floating ice area

   * - ``tendacabf(t)``
     - kg s-1
     - Total SMB flux

   * - ``tendlibmassbf(t)``
     - kg s-1
     - Total BMB flux

   * - ``tendlibmassbffl(t)``
     - kg s-1
     - Total BMB flux beneath floating ice

   * - ``tendlicalvf(t)``
     - kg s-1
     - Total calving flux

   * - ``tendlifmassbf(t)``
     - kg s-1
     - Total calving and ice front melting flux

   * - ``tendligroundf(t)``
     - kg s-1
     - Total grounding line flux
