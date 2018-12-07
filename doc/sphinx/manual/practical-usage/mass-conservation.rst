.. include:: ../../global.txt

.. default-role:: literal

.. _sec-mass-conservation:

Balancing the books
===================

2D diagnostics
--------------

PISM provides a number of 2D diagnostics to keep track of mass conservation.\ [#]_

All of them are computed as time-averaged fluxes over requested reporting intervals.
Positive values correspond to mass gain.

For ice mass, at every grid point we have

.. literalinclude:: conservation/ice_mass_accounting_error.txt

.. only:: html

   Click :download:`here <conservation/ice_mass_accounting_error.txt>` to download this
   `ncap2` script.

All names on the right-hand side correspond to valid PISM diagnostic quantities.

To check that all changes in mass are accounted for, download the script above and run\ [#]_

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

Comments
^^^^^^^^

- `tendency_of_ice_mass_due_to_flow` is the change in ice mass corresponding to flux
  divergence
- `tendency_of_ice_mass_due_to_conservation_error` is the artificial change in ice mass
  needed to preserve non-negativity of ice thickness. While PISM's mass transport scheme
  is not proven to be mass-conserving *in every case* (see
  :cite:`JaroschSchoofAnslow2013`), in most simulations the field is uniformly zero.
- `tendency_of_ice_mass_due_to_surface_mass_balance` is the change due to the surface
  mass balance; note that this is *not* the same as the provided SMB: in ablation areas
  this is the *effective* mass balance taking into account the amount of ice present
- `tendency_of_ice_mass_due_to_basal_mass_balance` is the *effective* change due to
  basal (grounded and sub-shelf) melting
- `tendency_of_ice_mass_due_to_discharge` combines changes due to calving and frontal
  melt

Scalar diagnostics
------------------

Diagnostics listed above are also available as scalars, integrated over the whole
computational domain. The "integrated" mass accounting error can be computed using the
`ncap2` script below.

.. literalinclude:: conservation/scalar_accounting_error.txt

.. only:: html

   Click :download:`here <conservation/scalar_accounting_error.txt>` to download this
   `ncap2` script.

Comments
^^^^^^^^

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
in subglacial water thickness. [#]_

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

.. [#] See :ref:`sec-diagnostics-list` for the full list of diagnostics.

.. [#] `ncap2` is a part of NCO_.

.. [#] We keep track of :math:`W_{\text{till}} + W`, i.e. the sum of the effective
       thickness of subglacial water stored in till *and* the effective thickness of
       subglacial water in the transport layer (if applicable).
