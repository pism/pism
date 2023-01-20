.. include:: ../../global.txt

.. default-role:: literal

.. _sec-mass-conservation:

Balancing the books
===================

2D diagnostics
--------------

PISM provides a number of 2D diagnostics to keep track of mass conservation.\ [#f1]_

All of them are computed as time-averaged fluxes over requested reporting intervals.
Positive values correspond to mass gain.

For ice mass, at every grid point we have

.. literalinclude:: conservation/ice_mass_accounting_error.txt

.. only:: html

   Click :download:`here <conservation/ice_mass_accounting_error.txt>` to download this
   `ncap2` script.

All names on the right-hand side correspond to valid PISM diagnostic quantities.

To check that all changes in mass are accounted for, download the script above and run\ [#f2]_

.. code-block:: bash

   ncap2 -v -S ice_mass_accounting_error.txt \
         pism_output.nc mass_accounting_error.nc

The variable `ice_mass_accounting_error` in `mass_accounting_error.nc` will contain ice
mass accounting errors at each point. All values of this variable should be close to or
equal to zero. They are not zero (in general) due to rounding errors, though.

Use a shortcut

.. code-block:: none

   pismr -extra_file ex.nc -extra_times N -extra_vars mass_fluxes,...

to save all fluxes needed to "balance the books" in terms of ice mass.

Alternatively, use fluxes in terms of "ice amount" (mass per unit area):

.. literalinclude:: conservation/ice_amount_accounting_error.txt

.. only:: html

   Click :download:`here <conservation/ice_amount_accounting_error.txt>` to download this
   `ncap2` script.

To save these, use the shortcut

.. code-block:: none

   pismr -extra_file ex.nc -extra_times N -extra_vars amount_fluxes,...

.. rubric:: Comments

- `tendency_of_ice_mass_due_to_flow` is the change in ice mass corresponding to flux
  divergence
- `tendency_of_ice_mass_due_to_conservation_error` is the artificial change in ice mass
  needed to "balance the books". It is uniformly zero in most simulations.
- `tendency_of_ice_mass_due_to_surface_mass_balance` is the change due to the surface
  mass balance; note that this is *not* the same as the provided SMB: in ablation areas
  this is the *effective* mass balance taking into account the amount of ice present.
- `tendency_of_ice_mass_due_to_basal_mass_balance` is the *effective* change due to
  basal (grounded and sub-shelf) melting.
- `tendency_of_ice_mass_due_to_discharge` combines changes due to calving and frontal
  melt.


Scalar diagnostics
------------------

Diagnostics listed above are also available as scalars, integrated over the whole
computational domain. The "integrated" mass accounting error can be computed using the
`ncap2` script below.

.. literalinclude:: conservation/scalar_accounting_error.txt

.. only:: html

   Click :download:`here <conservation/scalar_accounting_error.txt>` to download this
   `ncap2` script.

.. rubric:: Comments

- `tendency_of_ice_mass_due_to_flow` is the integral of :math:`-\nabla \cdot Q` over the
  computational domain. This should be zero (up to the effect of rounding errors) in
  simulations that *do not* use Dirichlet boundary conditions for ice thickness.
  Prescribing ice thickness creates sources and sinks, and this diagnostic describes their
  influence.
- `tendency_of_ice_mass_due_to_conservation_error` should be zero (or close to zero) in
  most simulations

.. _sec-mass-conservation-hydrology:

Mass accounting in subglacial hydrology models
----------------------------------------------

PISM's hydrology models provide all the diagnostic fields needed to keep track of changes
in subglacial water thickness.

.. note::

   We keep track of :math:`W_{\text{till}} + W`, i.e. the sum of the effective thickness
   of subglacial water stored in till *and* the effective thickness of subglacial water in
   the transport layer (if applicable).

At every grid point we have

.. literalinclude:: conservation/water_mass_accounting_error.txt

.. only:: html

   Click :download:`here <conservation/water_mass_accounting_error.txt>` to download this
   `ncap2` script.

All names on the right-hand side correspond to valid PISM diagnostic quantities.

Use a shortcut

.. code-block:: none

   pismr -extra_file ex.nc -extra_times N -extra_vars hydrology_fluxes,...

to save all diagnostics mentioned above.

See :ref:`sec-subhydro` for more information about hydrology models.

Mass accounting in the PDD model
--------------------------------

PISM's PDD model provides diagnostics needed to compare computed accumulation, melt, and
runoff to the effective mass balance. Use diagnostic quantities
`surface_accumulation_flux`, `surface_melt_flux`, and `surface_runoff_flux` (units of mass
per area per time) and `surface_accumulation_rate`, `surface_melt_rate`,
`surface_runoff_rate` (units of mass per time).

To save all these, use `-extra_vars` shortcuts `pdd_fluxes` and `pdd_rates`.

.. _sec-mass-conservation-rough-bed:

Mass conservation and "rough" bed topography
--------------------------------------------

Jarosch and others :cite:`JaroschSchoofAnslow2013` show that Mahaffy's :cite:`Mahaffy` SIA
discretization used by PISM suffers from mass conservation errors near sufficiently abrupt
changes in bed elevation. (It may overestimate ice fluxes through the boundary of a grid
cell and remove more ice than available, producing a negative ice thickness.)

PISM uses a "projection step" to ensure non-negativity of ice thickness :math:`H`:

.. math::
   :label: eq-H-projection-step

   H^{n+1}_{i,j} = \max(\widetilde H^{n+1}_{i,j}, 0),

where :math:`\widetilde H^{n+1}_{i,j}` is the "tentative" ice thickness at a grid point
:math:`(i,j)` and the time step :math:`n+1` computed using an explicit-in-time
finite-volume discretization of the mass continuity equation.

This step is performed *after* computing the change in ice thickness due to flow and
*before* applying top surface and basal mass balance fluxes (i.e. PISM uses operator
splitting in its approximation of the mass continuity equation).

Prior to version `2.0.6` PISM fully relied on :eq:`eq-H-projection-step` to maintain
non-negativity of ice thickness and `tendency_of_ice_mass_due_to_conservation_error`
reported the rate at which mass is created by the projection step.

The current mass transport scheme includes a flux limiter (see section 3 and the appendix
of :cite:`Smolarkiewicz1989`) that ensures non-negativity of :math:`\widetilde
H^{n+1}_{i,j}`, making the projection step :eq:`eq-H-projection-step` unnecessary.

.. note::

   - PISM still performs the projections step to guarantee that :math:`H \ge 0` is true
     even if the flux limiter fails.

   - See `examples/bedrock_step` for PISM's implementation of the "cliff benchmark"
     described in :cite:`JaroschSchoofAnslow2013`.

.. [#f1] See :ref:`sec-diagnostics-list` for the full list of diagnostics.

.. [#f2] `ncap2` is a part of NCO_.
