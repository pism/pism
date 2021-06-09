.. include:: ../../../global.txt

.. _sec-blatter:

Blatter's model
^^^^^^^^^^^^^^^

Unlike the rest of PISM, the Blatter solver uses a geometry-following vertical grid (see
:numref:`fig-grid-vertical-sigma`) to approximate horizontal components of ice velocity.
The number of vertical "levels" in this grid is controlled by
:config:`stress_balance.blatter.Mz`.

The non-linear system resulting from the discretization of PDEs corresponding to Blatter's
stress balance model is much harder to solve than the one corresponding to the
SSA system (:cite:`BrownSmithAhmadia2013`, :cite:`Tuminaro2016`) and (at this point)
experimentation with preconditioner choices seems inevitable. We use PETSc's command-line
options to control these choices.

.. note::

   The Blatter solver uses the ``-bp_`` command-line option prefix.

   Run PISM like this

   .. code-block:: bash

      pismr -stress_balance blatter \
            [other options] -help | grep "-bp_"

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

   M_{z}^0 &= M_z,

   M_{z}^k &= (M_{z}^{k-1} - 1)\, /\, C + 1.

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
described in :ref:`sec-sia-gradient`. Set
:config:`stress_balance.blatter.use_eta_transform` to enable it.

.. _sec-blatter-time-stepping:

Adaptive time stepping
######################

PISM's explicit in time mass continuity code is *conditionally stable*. When used with the
SSA + SIA hybrid, the maximum allowed time step is computed using a combination of the CFL
criterion :cite:`MortonMayers` and the maximum diffusivity of the SIA flow
:cite:`BBssasliding`. This time step restriction does not disappear when the same mass
continuity code is used with a stress balance model that does not explicitly compute
"advective" and "diffusive" parts of the flow. We need a work-around.

.. note::

   Very little is known about stability of explicit time stepping methods of the mass
   continuity equation coupled to a "generic" stress balance model.

   We don't have a rigorous justification for the approach described below.

When this BP solver is coupled to PISM, the vertically-averaged ice velocity is used in
place of the "advective" ("sliding") velocity from the SSA. As a result, the CFL-based
time step restriction is applied by existing PISM code.

However, it is almost always the case that the diffusivity-driven time step restriction is
more severe and so we need a replacement: CFL alone does not appear to be sufficient for
stability.

We compute an estimate of the "SIA-like" maximum diffusivity by observing that for the SIA
the vertically-averaged ice flux `Q` satisfies

.. math::

   Q = -D \nabla s.

We solve this for the diffusivity `D`:

.. math::
   :label: eq-bp-max-diffusivity

   D = \frac{H\, |\bar{\uu}|}{|\nabla s| + \epsilon}

.. FIXME: talk about the choice of \epsilon.

and use the maximum of this quantity to determine the maximum allowed time step using
:eq:`eq-sia-max-dt`.

.. note::

   Other models supporting this stress balance model and using an explicit in time
   geometry evolution method (:cite:`Lipscomb2019`, :cite:`Hoffman2018`) report that the
   CFL condition appears to be sufficient in practice.

   Given the lack of a theory describing the maximum time step necessary for stability it
   may make sense to experiment with *increasing* :config:`time_stepping.adaptive_ratio`.

   Setting it to a very large value would *completely disable* the diffusivity-based time
   step restriction.

.. note::

   The "time step skipping mechanism" enabled using :config:`time_stepping.skip.enabled`
   (see :ref:`sec-adapt`) has a different effect when the Blatter stress balance model is
   used: the full 3D ice velocity is updated during every sub-step and only the energy
   balance and age models takes the "long" time step.

   Since the Blatter solver is likely to dominate the computational cost, setting
   :config:`time_stepping.skip.enabled` to "true" is not likely to be beneficial.

.. _sec-blatter-preconditioners:

Practical preconditioners choices
#################################

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
   -bp_mg_levels_ksp_type richardson \
   -bp_mg_levels_pc_type sor \
   -bp_mg_coarse_ksp_type gmres \
   -bp_mg_coarse_pc_type hypre \
   -bp_mg_coarse_pc_hypre_type boomeramg

may work better\ [#bp-pc-settings]_, but requires PETSc built with hypre_.

Here ``-bp_snes_ksp_ew -bp_snes_ksp_ew_version 3`` enables Luis Chacón’s variant of the
Eisenstat-Walker :cite:`Eisenstat1996` method of adjusting linear solver tolerances to
avoid oversolving and ``-bp_mg_coarse_pc_type hypre -bp_mg_coarse_pc_hypre_type
boomeramg`` selects the BoomerAMG algebraic MG preconditioner from hypre_ for the coarse
MG level.

.. note::

   The Eisenstat-Walker adjustment of linear solver tolerances saves time when a
   low-accuracy estimate of the Newton step is sufficient but may lead to solver failures,
   especially when the initial guess is of poor quality. In an attempt to reduce
   computational costs while maintaining robustness PISM disables ``-bp_snes_ksp_ew`` if
   the initial guess is zero (beginning of a simulation) or if the solver fails with
   ``-bp_snes_ksp_ew``.

.. Maybe mention that setting -bp_snes_ksw_ew_rtol0 to a smaller value may make the solver
   more robust.

Some simulations may benefit from using a direct solver on the coarse MG level. For
example, the following would use MUMPS_ on the coarse grid:

.. code-block:: bash

   -bp_pc_type mg \
   -bp_pc_mg_levels N \
   -bp_snes_ksp_ew  \
   -bp_snes_ksp_ew_version 3 \
   -bp_mg_levels_ksp_type richardson \
   -bp_mg_levels_pc_type sor \
   -bp_mg_coarse_ksp_type preonly \
   -bp_mg_coarse_pc_type lu

*if* PETSc is built with MUMPS_.

.. note::

   Parallel direct solvers such as MUMPS really benefit from using optimized BLAS and
   LAPACK libraries.

   Please see section 3.5.3 of :cite:`petsc-user-ref` for instructions. At the time of
   writing

   .. code-block:: bash

      --download-f2cblaslapack --download-blis

   is recommended as a portable high-performance option. However, it makes sense to try
   other freely-available libraries (Intel MKL, OpenBLAS) as well.

Note, though, that the multigrid preconditioner, even if it is effective in terms of
reducing the number of Krylov iterations, may not be the cheapest one :cite:`Tezaur2015b`:
there is a trade off between the number of iterations and the cost of a single iteration.
Other preconditioner options may be worth considering as well.

In some cases node ordering and the way the domain is split among processes in a parallel
run may affect solver performance (see :cite:`BrownSmithAhmadia2013`, :cite:`Tezaur2015b`,
:cite:`Tuminaro2016`). These references mention staggering the unknowns so that `u` and
`v` components at the same node correspond to adjacent equations in the system and using
contiguous ordering of unknowns in the same ice column. This allows the solver to capture
vertical coupling *locally* using incomplete factorization.

In addition to this, :cite:`Tezaur2015b` mention that parallel domain distribution
partitioning ice columns among multiple processes *sometimes* leads to convergence issues.
Following this advice, PISM does not partition the domain in the `z` direction, but some
of our experiments show that if the solver struggles, switching to a *one-dimensional*
domain decomposition along the `y` direction may help (see
:ref:`sec-domain-distribution`).

Run PISM as follows to give this a try:

.. code-block:: bash

   mpiexec -n M pismr -Nx 1 -Ny M ...

This forces PISM to split the domain into `M` parts in the `y` direction instead of the
default (approximately `\sqrt{M}` in both `x` and `y`).

Please see :ref:`sec-blatter-details` for more.

Parameters
##########

Below is the complete list of configuration parameters controlling this solver (prefix:
``stress_balance.blatter.``):

.. pism-parameters::
   :prefix: stress_balance.blatter.

.. rubric:: Footnotes

.. [#semi-coarsening] Horizontal coordinates of grid points are the same in all grids in a
                      hierarchy, i.e. each grid is "extruded" from PISM's 2D grid with
                      uniform spacing in `x` and `y` directions.

.. [#bp-pc-settings] These settings are inspired by :cite:`BrownSmithAhmadia2013`.
