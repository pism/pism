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
diffusive, which leads to the loss of detail, especially near the base of the ice where
the variation of age with depth is more pronounced. Increasing vertical grid resolution
and using quadratic grid spacing (finer near the bed, coarser by a factor of
:config:`grid.lambda` near the surface; see :ref:`sec-grid`) would reduce the effect of
numerical diffusion but cannot eliminate it.

.. _sec-isochronal-layers:

Isochronal layer tracing
========================

To model closely-spaced isochrones PISM implements the *isochronal layer tracing scheme*
described in :cite:`Born2016` and :cite:`Born2021`; see section 2.4 of the latter for
details. This method uses a "Lagrangian" approximation in the vertical direction and a
first-order upwinding method in horizontal directions. This eliminates numerical diffusion
in the critical, *vertical* direction.

.. note::

   This approximation assumes that the age of the ice increases monotonically with depth,
   so it cannot be used to model overturning folds :cite:`Bons2016` and plume formation
   :cite:`Vieli2018`.

.. rubric:: Summary

Ice masses are interpreted as "stacks" of layers of varying thickness. Isochrones are
represented by boundaries between these layers; the depth of an isochrone is the sum of
thicknesses of all the layers above it (some may have zero thickness, e.g. in ablation
areas). The surface mass balance is applied to the topmost layer; if a layer is depleted
by a negative mass balance, the remaining mass loss is used to reduce the thickness of the
next layer below it. Similarly, the basal melt rate is applied to the bottom layer, then
the one above it, etc. Within each layer the mass is transported according to the
horizontal components of the 3D modeled ice velocity by sampling it at the depth
corresponding the middle of a layer. There is no mass transport between layers. A new
layer is added to the top of the stack each time PISM reaches a time in
:config:`isochrones.deposition_times`.

.. figure:: figures/isochrones.png
   :name: fig-isochrones

   Modeled isochrones in the 20000 year simulation using the isothermal SIA stress balance
   and the ``examples/isochrones``

Set :config:`isochrones.deposition_times` to enable this mechanism and see the
:var:`isochrone_depth` diagnostic it provides. The parameter
:config:`isochrones.deposition_times` takes an argument in the format similar to
:config:`output.extra.times`; see :ref:`sec-saving-diagnostics`.

.. rubric:: Bootstrapping

When starting a simulation that does not have the :var:`isochronal_layer_thickness`
information available PISM has to interpret existing ice thickness.

.. rubric:: Diagnostics

The isochronal tracing scheme provides the :var:`isochrone_depth` diagnostic (depth below
the surface for isochrones corresponding to all the requested
:config:`isochrones.deposition_times`) and the :var:`isochronal_layer_thickness` variable
(part of the state of this model).

.. note::

   A restarted simulation may not be able to write diagnostics to the same file because of
   a difference in the number of deposition times handled by the model during the first
   and the second run. In other words: using :config:`output.extra.append` is not
   recommended.

.. rubric:: Parameters

.. pism-parameters::
   :prefix: isochrones.
