.. include:: ../../global.txt

.. _sec-adapt:

Understanding adaptive time-stepping
------------------------------------

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
updated ("``h``"). The rest of the flags line looks like

.. code-block:: none

   diffusivity (dt=0.83945 in 2 substeps; av dt_sub_mass_cont=0.41972)

Recall that the PISM time step is determined by an adaptive mechanism. Stable mass
conservation and conservation of energy solutions require such an adaptive time-stepping
scheme :cite:`BBL`. The first character we see here, namely "``diffusivity``", is the
adaptive-timestepping "reason" flag. See :numref:`tab-adaptiveflag`. We also see
that there was a major time step of :math:`0.83945` model years divided into :math:`2`
substeps of about :math:`0.42` years. The :opt:`-skip` option enables this mechanism,
while :opt:`-skip_max` sets the maximum number of such substeps. The adaptive mechanism
may choose to take fewer substeps than ``-skip_max`` so as to satisfy certain numerical
stability criteria, however.

The second line in the above, the line which starts with "``S``", is the summary. Its
format, and the units for these numbers, is simple and is given by a couple of lines
printed near the beginning of the standard output for the run:

.. code-block:: none

   P       YEAR:       ivol      iarea  max_diffusivity  max_hor_vel
   U      years   10^6_km^3  10^6_km^2         m^2 s^-1       m/year

That is, in each summary we have the total ice volume, total ice area, maximum diffusivity
(of the SIA mass conservation equation), and maximum horizontal velocity (i.e.
:math:`\max(\max(|u|), \max(|v|))`).

The third line of the above example shows that the SSA stress balance was solved.
Information on the number of nonlinear (outer) and linear (inner) iterations is provided
:cite:`BBssasliding`.

.. list-table:: Meaning of the adaptive time-stepping "reason" flag in the standard output
                flag line.
   :header-rows: 1
   :name: tab-adaptiveflag
   :widths: 1,2

   * - PISM output
     - Active adaptive constraint or PISM sub-system that limited time-step size

   * - ``3D CFL``
     - three-dimensional CFL for temperature/age advection :cite:`BBL`

   * - ``diffusivity``
     - diffusivity for SIA mass conservation :cite:`BBL`, :cite:`HindmarshPayne`

   * - ``end of the run``
     - end of prescribed run time

   * - ``max``
     - maximum allowed :math:`\Delta t` applies; set with ``-max_dt``

   * - ``internal (derived class)``
     - maximum :math:`\Delta t` was temporarily set by a derived class

   * - ``2D CFL``
     - 2D CFL for mass conservation in SSA regions (upwinded; :cite:`BBssasliding`)

   * - ``-ts_... reporting``
     - the ``-ts_times`` option and the configuration flag
       :config:`time_stepping.hit_ts_times`; see section :ref:`sec-saving-time-series`

   * - ``-extra_... reporting``
     - the ``-extra_times`` option; see section :ref:`sec-saving-diagnostics`

   * - ``surface``
     - a surface or an atmosphere model

   * - ``ocean``
     - an ocean model

   * - ``hydrology``
     - a hydrology model stability criterion, see section :ref:`sec-subhydro`

   * - ``BTU``
     - time-the bedrock thermal layer model, see section :ref:`sec-energy`

   * - ``eigencalving``
     - the eigen-calving model, see section :ref:`sec-calving`

.. list-table:: Options controlling time-stepping
   :header-rows: 1
   :name: tab-time-stepping
   :widths: 3,5

   * - Option
     - Description

   * - :opt:`-adapt_ratio`
     - Adaptive time stepping ratio for the explicit scheme for the mass balance equation.

   * - :opt:`-max_dt` (years)
     - The maximum time-step in years. The adaptive time-stepping scheme will make the
       time-step shorter than this as needed for stability, but not longer.

   * - :opt:`-skip`
     - Enables time-step skipping, see below.

   * - :opt:`-skip_max`
     - Number of mass-balance steps, including SIA diffusivity updates, to perform before
       temperature, age, and SSA stress balance computations are done. This is only
       effective if the time step is being limited by the diffusivity time step
       restriction associated to mass continuity using the SIA. The maximum recommended
       value for ``-skip_max`` is, unfortunately, dependent on the context. The
       temperature field should be updated when the surface changes significantly, and
       likewise the basal sliding velocity if it comes (as it should) from the SSA
       calculation.

   * - :opt:`-timestep_hit_multiples` (years)
     - Hit multiples of the number of model years specified. For example, if stability
       criteria require a time-step of 11 years and the ``-timestep_hit_multiples 3``
       option is set, PISM will take a 9 model year long time step. This can be useful to
       enforce consistent sampling of periodic climate data.
