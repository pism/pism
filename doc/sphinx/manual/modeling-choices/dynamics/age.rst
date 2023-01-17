.. include:: ../../../global.txt

.. _sec-age:

Computing ice age
-----------------

By default, PISM does not compute the age of the ice because it does not directly impact
ice flow when using the default flow laws. (But see :ref:`sec-sia-age-coupling` for SIA
parameterizations that use it.)

Set :config:`age.enabled` to turn it on. A 3D variable :var:`age` will appear in output
files. It is read during model initialization if :config:`age.enabled` is set and ignored
otherwise. If :config:`age.enabled` is set and the variable :var:`age` is absent in the
input file then the initial age is set to :config:`age.initial_value`.

The first order upwinding method used to approximate the evolution of :var:`age` is
diffusive, which leads to the loss of detail, especially closer to the base of the ice
where age varies more quickly with depth. Increasing vertical grid resolution and using
quadratic grid spacing (finer near the bed, coarser by a factor of :config:`grid.lambda`
near the surface; see :ref:`sec-grid`) would reduce the effect of numerical diffusion but
cannot eliminate it.

.. _sec-isochronal-layers:

Isochronal layer tracing
========================

To model closely-spaced isochrones PISM implements the *isochronal layer tracing scheme*
described in :cite:`Born2016` and :cite:`Born2021`; see section 2.4 of the latter for
details. This method uses a "Lagrangian" approximation in the vertical direction and a
first-order upwinding method in horizontal directions. This eliminates numerical diffusion
in the critical, *vertical* direction.

Set :config:`isochrones.deposition_times` to enable this mechanism and see the
:var:`isochrone_depth` diagnostic it provides.

.. note::

   - The parameter :config:`isochrones.deposition_times` takes an argument in the format
     similar to :config:`output.extra.times`; see :ref:`sec-saving-diagnostics`.

   - This approximation assumes that the age of the ice increases monotonically with
     depth, so it cannot be used to model overturning folds :cite:`Bons2016` and plume
     formation :cite:`Vieli2018`.
