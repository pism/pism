.. include:: ../../global.txt

.. _sec-mass-conservation:

Balancing the books
===================

2D diagnostics
--------------

PISM provides the following 2D diagnostics to keep track of mass conservation.\ [#]_

All these are computed as time-averaged fluxes over reporting intervals. Positive values
correspond to mass gain.

At every grid point (point-wise) we have

.. code-block:: none

   tendency_of_ice_amount = (tendency_of_ice_amount_due_to_flow +
                             tendency_of_ice_amount_due_to_conservation_error +
                             tendency_of_ice_amount_due_to_surface_mass_balance +
                             tendency_of_ice_amount_due_to_basal_mass_balance +
                             tendency_of_ice_amount_due_to_discharge)

Here "ice amount" is "ice mass per unit area", units of `kg / m^2`.

Alternatively one can request the same 2D diagnostics in terms of ice mass:

.. code-block:: none

   tendency_of_ice_mass = (tendency_of_ice_mass_due_to_flow +
                           tendency_of_ice_mass_due_to_conservation_error +
                           tendency_of_ice_mass_due_to_surface_mass_balance +
                           tendency_of_ice_mass_due_to_basal_mass_balance +
                           tendency_of_ice_mass_due_to_discharge)

Use a shortcut

.. code-block:: none

   pismr -extra_file ex.nc -extra_times N -extra_vars amount_fluxes,...

to save fluxes in terms of "ice amount" and

.. code-block:: none

   pismr -extra_file ex.nc -extra_times N -extra_vars mass_fluxes,...

to save fluxes in terms of "ice mass."

Comments
^^^^^^^^

- ``tendency_of_ice_mass_due_to_flow`` is the change in ice mass corresponding to flux
  divergence
- ``tendency_of_ice_mass_due_to_conservation_error`` is the artificial change in ice mass
  needed to preserve non-negativity of ice thickness.
- ``tendency_of_ice_mass_due_to_surface_mass_balance`` is the change due to the surface
  mass balance; note that this is *not* the same as the provided SMB: in ablation areas
  this is the *effective* mass balance taking into account the amount of ice present
- ``tendency_of_ice_mass_due_to_basal_mass_balance`` is the *effective* change due to
  basal (grounded and sub-shelf) melting
- ``tendency_of_ice_mass_due_to_discharge`` combines changes due to calving and frontal
  melt

Scalar diagnostics
------------------

Diagnostics listed above are also available as scalar diagnostics, integrated over the
whole computational domain.

.. code-block:: none

   tendency_of_ice_mass = (tendency_of_ice_mass_due_to_influx +
                           tendency_of_ice_mass_due_to_conservation_error +
                           tendency_of_ice_mass_due_to_basal_mass_balance +
                           tendency_of_ice_mass_due_to_surface_mass_balance +
                           tendency_of_ice_mass_due_to_discharge)

Comments
^^^^^^^^

- ``tendency_of_ice_mass_due_to_influx`` is the integral of `-\nabla \cdot Q` over the
  computational domain. This should be zero (up to the effect of rounding errors) in
  simulations that *do not* use Dirichlet boundary conditions for ice thickness.
  Prescribing ice thickness creates sources and sinks, and this diagnostic describes their
  influence.
- ``tendency_of_ice_mass_due_to_conservation_error`` should be zero (or close to zero) in
  most simulations

.. [#] See :ref:`sec-diagnostics-list` for the full list of diagnostics.

.. [#] While PISM's mass transport scheme is not proven to be mass-conserving *in every
  case* (see :cite:`JaroschSchoofAnslow2013`), in most simulations this field is uniformly
  zero.
