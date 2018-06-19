.. include:: ../../../global.txt

.. contents::

.. _sec-regional:

Modeling individual outlet glaciers
===================================

PISM was created to model ice sheets in entirety; in this context we can assume that ice
does not extend to the edge of the computational domain and there is no need to provide
lateral boundary conditions.

However, in some interesting cases the ice *does* extend to the edge of the domain, but we
can assume that changes near the boundary do not affect the behavior in the region we want
to model. Two examples come to mind.

1. Modeling an individual alpine glacier: there may be ice near a domain boundary, but we
   can select the domain so that corresponding ice masses are not connected to the glacier
   we're modeling. Note that it is not always possible to "remove" glaciers we don't care
   about because they may re-appear due to the mass balance forcing.
2. Modeling an outlet glacier in an ice sheet: we cut out a region from an ice sheet, so
   ice *will* extend to the edge of the domain, but we can select the domain so that it
   contains the whole *drainage basin* of the outlet glacier of interest.

   In some ways this case is harder to model: the shape and size of drainage basins of an
   ice sheet changes with its geometry, so we have to assume that our simulation is short
   enough to ensure that the drainage basin we're modeling remains roughly the same in
   shape and extent.

   See :ref:`sec-jako` for an example.

Use

.. code:: bash

   pismr -regional ...

to enable PISM's "regional mode."

Ideally, modeling a region containing an ice mass extending to the edge of the domain
would use the following time-dependent lateral boundary conditions

- energy balance: enthalpy flux across the boundary,
- mass continuity: mass flux across the boundary,
- stress balance: sliding speed at the boundary.

PISM's regional mode uses a special mask :var:`no_model_mask` (zeros in the interior of
the modeling domain, ones at the edge of the domain or in other areas that are *not
modeled*) to implement modifications at domain boundaries. This mask is saved to output
files and read back in when the model is re-started. Set :config:`regional.no_model_strip`
during bootstrapping to create a "non-modeled" strip of a given width along the domain
boundary.

Energy
------

PISM assumes that ice enthalpy and the basal melt rate (i.e. parts of the model state that
capture the energy state) near the boundary of the domain *remain constant*: at the end of
each time step updated enthalpy and basal melt rate are re-set to values read from an
input file or computed during bootstrapping at all grid points where :var:`no_model_mask`
is `1`.

Stress balance
--------------

When prescribing the sliding velocity, the :var:`no_model_mask` overrides the basal
sliding B.C. mask: all :var:`no_model_mask` locations are *also* the Dirichlet B.C.
locations for the sliding velocity. This makes it possible to prescribe the sliding
velocity of the ice across the domain boundary. Set
:config:`stress_balance.ssa.dirichlet_bc` to ``true`` to enable this feature.

In many cases it makes sense to *disable* sliding at the boundary. When the sliding
velocity near the boundary is not prescribed, PISM sets the basal yield stress to a high
value (see :config:`regional.no_model_yield_stress`).

The domain in PISM "wraps around", which means that we can not accurately compute
gradients near the boundary in the non-periodic case.

Note, though, that updating the velocity field requires computing the gravitational
driving stress, which depends on gradients of the ice thickness and surface elevation.

To avoid using finite differences across the domain boundary when computing these
gradients, PISM stores ice thickness and surface elevation near the edge of the domain and
uses them to modify surface elevation and thickness gradients.

.. note::

   In the SIA stress balance model, prescribing ice thickness and surface elevation near
   the edge of the domain is equivalent to prescribing the *flux* across the domain
   boundary.

To use *zero* surface elevation and thickness gradients, set
:config:`regional.zero_gradient`. (This disables SIA flow across the boundary.)

.. warning::

   High surface elevation and ice thickness gradients near the domain boundary *will*
   affect time-stepping even if they do now affect model evolution.

   The resulting high SIA diffusivity will force PISM to take unreasonably short time
   steps, wasting computational time.

   Consider setting :config:`regional.zero_gradient` if you see high SIA diffusivities
   near domain boundaries (save :var:`diffusivity_staggered` to check).

Mass continuity
---------------

PISM uses the SSA Dirichlet B.C. mask as the ice thickness Dirichlet B.C. mask, i.e. ice
thickness is fixed wherever the sliding velocity is fixed. (In other words, PISM allows
prescribing the *ice flux* at a given location.)

This means that the *ice thickness does not evolve* in the :var:`no_model_mask` area.

Mass balance adjustment
^^^^^^^^^^^^^^^^^^^^^^^

Prescribing the ice thickness near the boundary when the ice in the interior of the domain
thins would lead to high thickness and surface elevation gradients at the inner boundary
of the "non-modeled" strip. Use :ref:`sec-surface-forcing` to keep the ice geometry from
deviating from the target *without* sharp transitions at the boundary from fixed to
evolving ice thickness.

Calving
^^^^^^^

Set :config:`calving.front_retreat.wrap_around` to ``true`` to allow calving front retreat
due to calving to "wrap around" the computational domain. This may be necessary in some
regional synthetic-geometry setups.
