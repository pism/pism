.. include:: ../../../global.txt

.. _sec-stressbalance:

Choosing the stress balance
---------------------------

The basic stress balance used for all grounded ice in PISM is the non-sliding,
thermomechanically-coupled SIA :cite:`BBL`. For the vast majority of most ice sheets, as
measured by area or volume, this is an appropriate model, which is an `O(\epsilon^2)`
approximation to the Stokes model if `\epsilon` is the depth-to-length ratio of the ice
sheet :cite:`Fowler`.

The shallow shelf approximation (SSA) stress balance applies to floating ice. See the Ross
ice shelf example in section :ref:`sec-ross` for an example in which the SSA is only
applied to floating ice.

In PISM the SSA is also used to describe the sliding of grounded ice and the formation of
ice streams :cite:`BBssasliding`. Specifically for the SSA with "plastic" (Coulomb
friction) basal resistance, the locations of ice streams are determined as part of a free
boundary problem of Schoof :cite:`SchoofStream`, a model for emergent ice streams within a
ice sheet and ice shelf system. This model explains ice streams through a combination of
plastic till failure and SSA stress balance.

This SSA description of ice streams is the preferred "sliding law" for the SIA
:cite:`BBssasliding`, :cite:`Winkelmannetal2011`. The SSA should be combined with the SIA,
in this way, in preference to classical SIA sliding laws which make the sliding velocity
of ice a local function of the basal value of the driving stress. The resulting
combination of SIA and SSA is a "hybrid" approximation of the Stokes model
:cite:`Winkelmannetal2011`. Option ``-stress_balance ssa+sia`` turns on this "hybrid"
model. In this use of the SSA as a sliding law, floating ice is also subject to the SSA.

In addition to this, PISM includes an implementation of the first order approximation of
Stokes equations due to Blatter (:cite:`Blatter`, :cite:`Pattyn03`),

:numref:`tab-stress-balance-choice` describes the basic choice of stress balance.

.. list-table:: The basic choice of stress balance
   :name: tab-stress-balance-choice
   :header-rows: 1
   :widths: 1,2

   * - Option
     - Description

   * - :opt:`-stress_balance none`
     - Turn off ice flow completely.

   * - :opt:`-stress_balance sia` (default)
     - Grounded ice flows by the non-sliding SIA. Floating ice does not flow, so this model
       is not recommended for marine ice sheets.

   * - :opt:`-stress_balance ssa`
     - Use the SSA model exclusively. Horizontal ice velocity is constant throughout ice
       columns.

   * - :opt:`-stress_balance prescribed_sliding`
     - Use the constant-in-time prescribed sliding velocity field read from a file set
       using :opt:`-prescribed_sliding_file`, variables ``ubar`` and ``vbar``.
       Horizontal ice velocity is constant throughout ice columns.

   * - :opt:`-stress_balance ssa+sia`
     - The recommended sliding law, which gives the SIA+SSA hybrid stress balance.
       Combines SSA-computed velocity, using pseudo-plastic till, with SIA-computed
       velocity according to the combination in :cite:`Winkelmannetal2011`; similar to
       :cite:`BBssasliding`. Floating ice uses SSA only.

   * - :opt:`-stress_balance prescribed_sliding+sia`
     - Use the constant-in-time prescribed sliding velocity in combination with the
       non-sliding SIA.

   * - :opt:`-stress_balance blatter`
     - Use the first-order approximation of Stokes equations due to Blatter
       (:cite:`Blatter` :cite:`Pattyn03`).

.. _sec-blatter:

Controlling Blatter's stress balance model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Unlike the rest of PISM, the Blatter solver uses a geometry-following vertical grid (see
:numref:`fig-grid-vertical-sigma`) to approximate horizontal components of ice velocity.
The number of vertical "levels" in this grid is controlled by
:config:`stress_balance.blatter.Mz`.

The non-linear system resulting from the discretization of PDEs corresponding to Blatter's
stress balance model\ [#bp-fem]_ is much harder to solve than the one corresponding to the
SSA system (:cite:`BrownSmithAhmadia2013`, :cite:`Tuminaro2016`) and (at this point)
experimentation with preconditioner choices seems inevitable. We use PETSc's command-line
options to control these choices.

.. note::

   The Blatter solver uses the ``-bp_`` command-line option prefix.

   Run PISM like this

   .. code-block:: bash

      pismr -stress_balance blatter [other options] -help | grep "-bp_"

   to see the complete list of PETSc option controlling this solver.

The multigrid (MG) preconditioner using semi-coarsening in the vertical direction followed
by further (horizontal) coarsening using algebraic multigrid methods appears to be
effective :cite:`Tuminaro2016`. The option combination

.. code-block:: bash

   -bp_pc_type mg \
   -bp_pc_mg_levels N \
   -bp_mg_levels_ksp_type gmres \
   -bp_mg_coarse_pc_type gamg

roughly corresponds to this approach (see :ref:`sec-blatter-preconditioners` for more).

Unlike :cite:`Tuminaro2016`, who used a purely
algebraic approach, these options select a combination of geometric and algebraic
multigrid preconditioners.

To use a multigrid preconditioner the user has to specify

- the number of MG levels `N` using ``-bp_pc_mg_levels N``,
- the coarsening factor `C` by setting :config:`stress_balance.blatter.coarsening_factor`, and
- the vertical grid size `M_z` (:config:`stress_balance.blatter.Mz`).

The values of these parameters have to be compatible. Specifically, `M_z` has to have the
form

.. math::
   :label: eq-bp-vertical-grid-size

   M_z = A\cdot C^{N - 1} + 1

for some positive integer `A`.

.. note::

   PISM stops with an error message if :eq:`eq-bp-vertical-grid-size` is not satisfied.

To set up a multigrid preconditioner PISM needs to build a hierarchy of vertical grids\
[#semi-coarsening]_ with `M_z` points on the finest grid.. Starting with this grid, PISM
creates the next one by dividing the number of vertical *spaces* by the coarsening factor
`C`. Then the newly-created grid is coarsened and this process is continued, stopping when
the desired number `N` of grids (MG levels) reached.

Overall, the number of points `M_{z}^k` in the vertical grid number `k` in the hierarchy
is

.. math::

   M_{z}^1 &= M_z,

   M_{z}^k &= (M_{k-1} - 1)\, /\, C + 1.

This process explains the compatibility condition :eq:`eq-bp-vertical-grid-size`: the
number of **spaces** in all vertical grids in the hierarchy *except for the coarsest one*
has to be divisible by `C`.

.. list-table:: Some vertical grid hierarchies
   :name: tab-blatter-mg-levels
   :header-rows: 1
   :widths: 1,3

   * - Coarsening factor `C`
     - Possible sizes of vertical grids in a hierarchy

   * - `2`
     - 2, 3, 5, 9, 17, 33, **65**, 129, 257, 513, 1025, `\dots`

   * - `3`
     - 2, 4, 10, 28, 82, 244, 730, `\dots`

   * - `4`
     - 2, 5, 17, **65**, 257, 1025, `\dots`

   * - `5`
     - 2, 6, 26, 126, 626, 3126, `\dots`

   * - `6`
     - 2, 7, 37, 217, 1297, `\dots`

   * - `7`
     - 2, 8, 50, 344, 2402, `\dots`

   * - `8`
     - 2, 9, **65**, 513, 4097, `\dots`

By default `C = 2`, but *aggressive coarsening* (i.e. larger values of `C`, up to around
`8`) has been observed to work. As highlighted in :numref:`tab-blatter-mg-levels`,
sometimes the same number of vertical grid levels can be achieved using more than one
combination of the coarsening factor and the number of MG levels.

For example, we can set up a solver using `65` vertical levels and `3` MG levels with the
coarsening factor of `8`, or `4` MG levels and the factor of `4`, or `7` MG levels and the
coarsening factor of `2`. In general, the computational cost of an MG preconditioner
application increases with the number of MG levels, so the first hierarchy (`2, 9, 65`,
`C=8`) *may* be the best choice. *However,* coarsening that is too aggressive may make a
less effective preconditioner, requiring more Krylov iterations and increasing the
computational cost. Again, one may have to experiment to find settings that work best in a
particular setup.

The coarsest grid in a hierarchy should be as small as possible (corresponding to `A = 1`
in :eq:`eq-bp-vertical-grid-size`); two levels is the minimum achievable in the context of
the finite element method used to discretize the system (this corresponds to a mesh that
is just one element thick).

.. FIXME: I should document the way PISM computes the maximum allowed time step when the
   Blatter solver is "on" and compare to SIA's diffusivity-driven stability condition
   below in :ref:`sec-sia`. Possibly mention that :config:`time_stepping.adaptive_ratio`
   affects runs with this solver...

.. _sec-blatter-gradient:

Surface gradient computation
############################

Some synthetic geometry experiments with grounded margins show "checkerboard" artifacts in
computed ice velocity near steep margins. A similar issue and an attempt to address it are
described in :cite:`Lipscomb2019`.

This implementation takes a different approach: instead of using an "upwinded" finite
difference approximation of the surface gradient we allow using the `\eta` transformation
described in :ref:`sec-gradient`. Set :config:`stress_balance.blatter.use_eta_transform`
to enable it.

.. _sec-blatter-preconditioners:

Practical preconditioners choices
=================================

The option combination

.. code-block:: bash

   -bp_pc_type mg \
   -bp_pc_mg_levels N \
   -bp_mg_levels_ksp_type gmres \
   -bp_mg_coarse_pc_type gamg

sets up the *kind* is a multigrid preconditioner known to be effective, but it is not the
only one, and most likely not the best one.

Our experiments suggest that

.. code-block:: bash

   -bp_pc_type mg \
   -bp_pc_mg_levels N \
   -bp_snes_ksp_ew  \
   -bp_snes_ksp_ew_version 3 \
   -bp_mg_levels_ksp_type gmres \
   -bp_mg_levels_pc_type bjacobi \
   -bp_mg_levels_sub_pc_type ilu \
   -bp_mg_coarse_ksp_type gmres \
   -bp_mg_coarse_pc_type hypre \
   -bp_mg_coarse_pc_hypre_type boomeramg

may to work better\ [#bp-pc-settings]_, but requires PETSc built with hypre_.

Here ``-bp_snes_ksp_ew -bp_snes_ksp_ew_version 3`` enables Luis Chacón’s variant of the
Eisenstat-Walker :cite:`Eisenstat1996` method of adjusting linear solver tolerances to
avoid oversolving and ``-bp_mg_coarse_pc_type hypre -bp_mg_coarse_pc_hypre_type
boomeramg`` selects the BoomerAMG algebraic MG preconditioner from hypre_ for the coarse
MG level.

Some simulations may benefit from using a direct solver on the coarse MG level. For
example, the following would use MUMPS_ on the coarse grid:

.. code-block:: bash

   -bp_pc_type mg \
   -bp_pc_mg_levels N \
   -bp_snes_ksp_ew  \
   -bp_snes_ksp_ew_version 3 \
   -bp_mg_levels_ksp_type gmres \
   -bp_mg_levels_pc_type bjacobi \
   -bp_mg_levels_sub_pc_type ilu \
   -bp_mg_coarse_ksp_type preonly \
   -bp_mg_coarse_pc_type lu

*if* PETSc is built with MUMPS_.

Note, though, that the multigrid preconditioner, even if it is effective in terms of
reducing the number of Krylov iterations, may not be the cheapest one :cite:`Tezaur2015b`:
there is a trade off between the number of iterations and the cost of a single iteration.
Other preconditioner options may be worth considering as well.

List of parameters
##################

Below is the complete list of configuration parameters controlling this solver (prefix:
``stress_balance.blatter.``):

.. pism-parameters::
   :prefix: stress_balance.blatter.

.. _sec-ssa:

Controlling the SSA stress balance model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the SSA stress balance is used, a choice of two solvers is available, namely
``-ssa_method fd`` (default) or ``-ssa_method fem``. See :numref:`tab-ssa-usage`, which
describes additional controls on the numerical solution of the stress balance equations.
If option ``-ssa_method fd`` is chosen then several more controls on numerics are
available; see :numref:`tab-ssafd-controls`. If the ice sheet being modeled has any
floating ice then the user is advised to read section :ref:`sec-pism-pik` on modeling
marine ice sheets.

When using SSA as a "sliding law" one also needs to model the yield stress, or a
pseudo-yield-stress in the case of power law sliding (section :ref:`sec-basestrength`).

The basal yield stress is normally a function of the amount of water stored in the till
and a (generally) spatially-varying till strength. The amount of stored basal water is
modeled by the subglacial hydrology model (section :ref:`sec-subhydro`) based on the basal
melt rate which is, primarily, thermodynamically-determined (see :ref:`sec-energy`).


.. list-table:: Choice of, and controls on, the numerical SSA stress balance.
   :name: tab-ssa-usage
   :header-rows: 1
   :widths: 1,2

   * - Option
     - Description

   * - :opt:`-ssa_method` [ ``fd | fem`` ]
     - Both finite difference (``fd``; the default) and finite element (``fem``) versions
       of the SSA numerical solver are implemented in PISM. The ``fd`` solver is the only
       one which allows PIK options (section :ref:`sec-pism-pik`). ``fd`` uses Picard
       iteration :cite:`BBssasliding`, while ``fem`` uses a Newton method. The ``fem`` solver
       has surface velocity inversion capability :cite:`Habermannetal2013`.

   * - :opt:`-ssa_eps` (`10^{13}`)
     - The numerical schemes for the SSA compute an effective viscosity `\nu` which
       depends on strain rates and ice hardness (thus temperature). The minimum value of
       the effective viscosity times the thickness (i.e. `\nu H`) largely determines the
       difficulty of solving the numerical SSA. This constant is added to keep `\nu H`
       bounded away from zero: `\nu H \to \nu H + \epsilon_{\text{SSA}}`, where
       `\epsilon_{\text{SSA}}` is set using this option. Units of :opt:`ssa_eps` are
       `\text{Pa}\,\text{m}\,\text{s}`. Set to zero to turn off this lower bound.

   * - :opt:`-ssa_view_nuh`
     - View the product `\nu H` for your simulation as a runtime viewer (section
       :ref:`sec-diagnostic-viewers`). In a typical Greenland run we see a wide range of
       values for `\nu H` from `\sim 10^{14}` to `\sim 10^{20}`
       `\text{Pa}\,\text{m}\,\text{s}`.

.. list-table:: Controls on the numerical iteration of the ``-ssa_method fd`` solver
   :name: tab-ssafd-controls
   :header-rows: 1
   :widths: 1,2

   * - Option
     - Description

   * - :opt:`-ssafd_picard_maxi` (300)
     - Set the maximum allowed number of Picard (nonlinear) iterations in solving the
       shallow shelf approximation.

   * - :opt:`-ssafd_picard_rtol` (`10^{-4}`)
     - The Picard iteration computes a vertically-averaged effective viscosity which is
       used to solve the equations for horizontal velocity. Then the new velocities are
       used to recompute an effective viscosity, and so on. This option sets the relative
       change tolerance for the effective viscosity. The Picard iteration stops when
       successive values `\nu^{(k)}` of the vertically-averaged effective viscosity
       satisfy

       .. math::

          \|(\nu^{(k)} - \nu^{(k-1)}) H\|_1 \le Z \|\nu^{(k)} H\|_1

       where `Z=` ``ssafd_picard_rtol``.

   * - :opt:`-ssafd_ksp_rtol` (`10^{-5}`)
     - Set the relative change tolerance for the iteration inside the Krylov linear solver
       used at each Picard iteration.

   * - :opt:`-ssafd_max_speed` (`50 km/yr`)
     - Limits computed SSA velocities: ice speed is capped at this limit after each Picard
       iteration of the SSAFD solver. This may allow PISM to take longer time steps by
       ignoring high velocities at a few troublesome locations.

.. _sec-sia:

Controlling the SIA stress balance model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The section :ref:`sec-age` describes coupling SIA stress balance to the age of the ice.

The section :ref:`sec-gradient` covers available methods of computing the surface gradient.

.. note::

   The explicit time stepping of the mass continuity equation in the case of the SIA flow
   comes with a severe restriction on time step length:

   .. math::

      \Delta t \le \frac{2 R}{D\left( 1/\Delta x^2 + 1/\Delta y^2 \right)}

   Here `D` is the maximum diffusivity of the SIA flow and `R` is
   :config:`time_stepping.adaptive_ratio`, a tuning parameter that further reduces the
   maximum allowed time step length.

   The maximum diffusivity `D` may be achieved at an isolated grid point near the ice
   margin. In this case it might make sense to limit the diffusivity of the SIA flow,
   sacrificing accuracy at a few grid points to increase time step length and reduce the
   computational cost. Set :config:`stress_balance.sia.limit_diffusivity` to enable this
   mechanism.

   When :config:`stress_balance.sia.limit_diffusivity` is ``false`` PISM stops as soon as
   the SIA diffusivity at any grid point exceeds
   :config:`stress_balance.sia.max_diffusivity`. We do this to make it easier to detect
   problematic model configurations: in many cases it does not make sense to continue a
   simulation if `D` is very large.

.. _sec-weertman:

Weertman-style sliding
^^^^^^^^^^^^^^^^^^^^^^

.. warning::

   This kind of sliding is, in general, a bad idea. We implement it to simplify
   comparisons of the "hybrid" model mentioned above to older studies using this
   parameterization.

PISM implements equation 5 from :cite:`Tomkin2007`:

.. math::
   :label: eq-weertman-sliding

   \mathbf{u}_s = \frac{2 A_s \beta_c (\rho g H)^{n}}{N - P} |\nabla h|^{n-1} \nabla h.

.. list-table:: Notation used in :eq:`eq-weertman-sliding`
   :name: tab-weertman-notation
   :header-rows: 1
   :widths: 1,9

   * - Variable
     - Meaning

   * - `H`
     - ice thickness

   * - `h`
     - ice surface elevation

   * - `n`
     - flow law exponent

   * - `g`
     - acceleration due to gravity

   * - `\rho`
     - ice density

   * - `N`
     - ice overburden pressure, `N = \rho g H`

   * - `P`
     - basal water pressure

   * - `A_s`
     - sliding parameter

   * - `\beta_c`
     - "constriction parameter" capturing the effect of valley walls on the flow;
       set to `1` in this implementation

We assume that the basal water pressure is a given constant fraction of the overburden
pressure: `P = k N`. This simplifies :eq:`eq-weertman-sliding` to

.. math::

   \mathbf{u}_s = \frac{2 A_s}{1 - k} ( \rho g H\, |\nabla h| )^{n-1} \nabla h.

This parameterization is used for grounded ice *where the base of the ice is temperate*.

To enable, use :opt:`-stress_balance weertman_sliding` (this results in constant-in-depth
ice velocity) or :opt:`-stress_balance weertman_sliding+sia` to use this parameterization
as a sliding law with the deformational flow modeled using the SIA model.

Use configuration parameters :config:`stress_balance.weertman_sliding.k` and
:config:`stress_balance.weertman_sliding.A` tot set `k` and `A_s`, respectively. Default
values come from :cite:`Tomkin2007`.

.. rubric:: Footnotes

.. [#semi-coarsening] Horizontal coordinates of grid points are the same in all grids in a
                      hierarchy, i.e. each grid is "extruded" from PISM's 2D grid with
                      uniform spacing in `x` and `y` directions.
.. [#bp-fem] See :ref:`sec-blatter-details`.

.. [#bp-pc-settings] These settings are inspired by :cite:`BrownSmithAhmadia2013`.
