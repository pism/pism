.. include:: ../../global.txt

.. _sec-adapt:

Understanding adaptive time-stepping
------------------------------------

It is helpful to keep in mind this fundamental fact:

    *length of time steps taken by a model affects results of a simulation.*

This applies to all evolutionary models and PISM is no different.

We expect model results to converge to the solution of the continuum problem corresponding
to a model as `\dt` goes to zero. Also, results using different `\dt` should be "close" to
this solution and to each other, *but they need not be the same*.

One important consequence is that changes in PISM settings that may not seem to be related
to physical choices may affect results *if* they affect time stepping.

Here is a possibly-incomplete list of PISM components and settings that may affect time
step length.

#. Numerical stability criteria.

   #. Diffusivity-based time step restriction for the mass continuity (mass transport)
      step using SIA diffusivity (or its estimate when the Blatter solver is used).
   #. The value of :config:`time_stepping.adaptive_ratio` adjusting the diffusivity-based
      time step restriction (see :eq:`eq-sia-max-dt`).
   #. CFL time step restriction for the mass continuity step using sliding velocity, (or
      vertically-averaged horizontal velocity with the Blatter solver).
   #. CFL time step restriction for horizontal advection within the ice volume within
      energy balance and age models. Uses horizontal (`u,v`) components of the ice
      velocity within the 3D volume of the ice.

#. Reporting.

   #. If  :config:`time_stepping.hit_ts_times` is set, PISM will adjust time step lengths
      to "hit" times requested with :config:`output.timeseries.times`.
   #. If :config:`time_stepping.hit_extra_times` is set (the default), PISM will adjust
      time step lengths to "hit" times requested with :config:`output.extra.times`.
   #. If :config:`time_stepping.hit_save_times` is set, PISM will adjust time step lengths
      to "hit" times requested with :config:`output.snapshot.times`.

#. Time-step "skipping" to reduce computational costs:
   :config:`time_stepping.skip.enabled` and :config:`time_stepping.skip.max`.

   This mechanism enables PISM to take a number of "cheap" mass-balance steps (including
   SIA diffusivity updates) before more expensive temperature, age, and SSA stress balance
   computations are done.

   Time step "skipping" roughly corresponds to asynchronous coupling between

   - ice flow by shear in planes parallel to the geoid and
   - membrane-type ice flow and advection of energy and tracers (such as age).

   This is only effective if the time step is being limited by the diffusivity time step
   restriction associated to mass continuity using the SIA.

   PISM computes time step restrictions `\{\dt_0, \dt_1, \dots, \dt_n \}` from *all* of
   PISM's sub-modules and sorts them from from smallest to largest. Then the maximum
   allowed time step is

   .. math::

      \dt_{\text{max}} = \dt_0.

   If :config:`time_stepping.skip.enabled` is set *and* the most severe restriction comes
   from the SIA-diffusivity-based stability criterion for mass continuity, it skips

   .. math::
      :label: eq-dt-skip

      N = \min\left(\left\lfloor 0.95\, \frac{\dt_1}{\dt_0} \right\rfloor,
                    N_{\text{max}}  \right)

   energy, age, and 3D velocity updates, where `N_{\text{max}}` is set using
   :config:`time_stepping.skip.max`.

   .. warning::

      The effects of this mechanism are not well understood. Please use with caution.

      The maximum recommended value for :config:`time_stepping.skip.max` depends on the
      context. The temperature field should be updated when the surface changes
      significantly, and likewise the basal sliding velocity if it comes from the SSA
      calculation.

#. Atmosphere, surface process, ocean, and sea level forcing components.

#. The Lingle-Clark bed deformation model (see :ref:`sec-bed-def-lc` and
   :config:`bed_deformation.lc.update_interval`).

#. If :config:`geometry.front_retreat.use_cfl` is set, PISM adjusts time step lengths to
   satisfy the CFL condition that uses the total front retreat rate coming from calving
   and frontal melt models.

#. The time step length never exceeds :config:`time_stepping.maximum_time_step`.

#. If :config:`time_stepping.hit_multiples` is set to a positive number, PISM will "hit"
   multiples of the number of model years specified. For example, if stability criteria
   require a time-step of 11 years and the ``-timestep_hit_multiples 3`` option is set,
   PISM will take a 9 model year long time step. This can be useful to enforce consistent
   sampling of periodic climate data.

#. If the value `R` set using :config:`time_stepping.resolution` is positive PISM
   reduces the time step length so that

   .. math::
      :label: eq-dt-rounding-down

      \dt = N\cdot R

   for some integer `N`.

   The default `R` (`1` second) allows PISM to represent model time more accurately,
   reducing the influence of rounding errors.

   .. note::

      This is an intermediate-term solution for an issue reported by Thomas Kleiner:
      some simulations produced different results with identical inputs but *different*
      start years.

      We tracked it down to the fact that these simulations ended up using slightly
      different time step lengths. This, in turn, was caused by differences in the
      *absolute* precision of the C++ type ``double`` for numbers of different magnitudes.

#. The time step length never exceeds the total length of the run.

At each time step the PISM standard output includes "flags" and then a summary of the
model state using a few numbers. A typical example is

.. code-block:: none

   v$Eh  diffusivity (dt=0.83945 in 2 substeps; av dt_sub_mass_cont=0.41972)
   S -124791.571:  3.11640   2.25720      3.62041    18099.93737
   y  SSA:     3 outer iterations, ~17.0 KSP iterations each

The characters "``v$Eh``" at the beginning of the flags line, the first line in the above
example, give a very terse description of which physical processes were modeled in that
time step. Here "``v``" means that a stress balance was solved to compute the velocity.
Then the enthalpy was updated ("``E``") and the ice thickness and surface elevation were
updated ("``h``"). The rest of the line looks like

.. code-block:: none

   diffusivity (dt=0.83945 in 2 substeps; av dt_sub_mass_cont=0.41972)

Recall that the PISM time step is determined by an adaptive mechanism. Stable mass
conservation and conservation of energy solutions require such an adaptive time-stepping
scheme :cite:`BBL`. The first word we see here, namely "``diffusivity``", is the
adaptive-timestepping "reason". See :numref:`tab-dt-reason`. We also see that
there was a major time step of `0.83945` model years divided into `2` substeps of about
`0.42` years. The parameter :config:`time_stepping.skip.enabled` enables this mechanism,
while :config:`time_stepping.skip.max` sets the maximum number of such substeps. The
adaptive mechanism may choose to take fewer substeps than :config:`time_stepping.skip.max`
so as to satisfy certain numerical stability criteria, however.

The second line in the above, the line which starts with "``S``", is the summary. Its
format, and the units for these numbers, is simple and is given by a couple of lines
printed near the beginning of the standard output for the run:

.. code-block:: none

   P       YEAR:       ivol      iarea  max_diffusivity  max_sliding_vel
   U      years   10^6_km^3  10^6_km^2         m^2 s^-1           m/year

That is, in each summary we have the total ice volume, total ice area, maximum diffusivity
(of the SIA mass conservation equation), and "maximum sliding velocity". Specifically, the
last number is `\max(\max(|u|), \max(|v|))`, i.e. the number that appears in the CFL time
step restriction for the "advective" part of the flow in the mass continuity equation.

.. note::

   ``max_sliding_vel`` reported here is not the same as the maximum sliding speed,
   `\max(\sqrt{u^2 + v^2})`.

The third line of the above example shows that the SSA stress balance was solved.
Information on the number of nonlinear (outer) and linear (inner) iterations is provided
:cite:`BBssasliding`.

.. list-table:: Meaning of some adaptive time-stepping "reasons" in the standard
                output line
   :header-rows: 1
   :name: tab-dt-reason
   :widths: 1,2

   * - PISM output
     - Active adaptive constraint or PISM sub-system that limited time-step size

   * - ``2D CFL``
     - CFL condition for the "advective" part of the mass continuity equation. Uses (`u`,
       `v`) components of the vertically-averaged horizontal ice velocity with the Blatter
       stress balance model. Uses sliding velocity with all other stress balance choices
       :cite:`BBssasliding`.

   * - ``diffusivity``
     - SIA-diffusivity-based time step restriction for the mass continuity equation
       :cite:`BBL`, :cite:`HindmarshPayne`

   * - ``energy``, ``age model``
     - CFL condition using horizontal (`u`, `v`) components of the ice velocity within the
       ice volume. Restricts the time step of enthalpy, temperature, or age advection
       :cite:`BBL`. (This CFL condition does not use the vertical (`w`) component of ice
       velocity: PISM's 3D advection scheme is explicit in `x` and `y` and implicit in
       `z`.)

   * - ``end of the run``
     - end of prescribed run time

   * - ``max``
     - maximum allowed `\dt` set with :config:`time_stepping.maximum_time_step`

   * - ``-ts_... reporting``
     - the ``-ts_times`` option and the configuration parameter
       :config:`time_stepping.hit_ts_times`; see section :ref:`sec-saving-time-series`

   * - ``-extra_... reporting``
     - the ``-extra_times`` option and the configuration parameter
       :config:`time_stepping.hit_extra_times`; see section :ref:`sec-saving-diagnostics`

   * - ``surface``
     - a surface or an atmosphere model

   * - ``ocean``
     - an ocean model

   * - ``hydrology``
     - a hydrology model stability criterion, see section :ref:`sec-subhydro`

   * - ``front_retreat``
     - CFL condition using the 2D horizontal retreat rate such as the eigen-calving model,
       see section :ref:`sec-calving`

.. _sec-time-stepping-parameters:

Parameters
==========

Prefix: ``time_stepping.``

.. pism-parameters::
   :prefix: time_stepping.
