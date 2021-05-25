.. include:: ../../global.txt

.. _sec-verif:

Verification
============

  Two types of errors may be distinguished: modeling errors and numerical errors. Modeling
  errors arise from not solving the right equations. Numerical errors result from not
  solving the equations right. The assessment of modeling errors is *validation*, whereas
  the assessment of numerical errors is called *verification*... Validation makes sense
  only after verification, otherwise agreement between measured and computed results may
  well be fortuitous.

  P. Wesseling, (2001) *Principles of Computational Fluid Dynamics*, pp. 560--561
  :cite:`Wesseling`

Verification is the essentially mathematical task of checking that the predictions of the
numerical code are close to the predictions of a continuum model, the one which the
numerical code claims to approximate. It is a crucial task for a code as complicated as
PISM. In verification there is no comparison between model output and observations of
nature. Instead, one compares exact solutions of the continuum model, in circumstances in
which they are available, to their numerical approximations.

Reference :cite:`Roache` gives a broad discussion of verification and validation in
computational fluid dynamics. See :cite:`BLKCB` and :cite:`BBL` for discussion of verification
issues for the isothermal and thermomechanically coupled shallow ice approximation (SIA),
respectively, and for exact solutions to these models, and :cite:`BBssasliding`,
:cite:`SchoofStream` for verification using an exact solution to the SSA equations for ice
streams.

In PISM there is a separate executable ``pismv`` which is used for SIA-related
verification, and there are additional scripts for SSA-related verification. The source
codes which are verified by ``pismv`` are, however, exactly the same source files as are
run by the normal PISM executable ``pismr``. In technical terms, ``pismv`` runs a derived
class of the PISM base class.

.. list-table:: Exact solutions for verification
   :header-rows: 1
   :name: tab-tests
   :widths: 1,4,2,4

   * - Test
     - Continuum model tested
     - Reference
     - Comments
   * - A
     - isothermal SIA, steady,  flat bed, constant accumulation
     - :cite:`BLKCB`
     -
   * - B
     - isothermal SIA, flat bed, zero accumulation
     - :cite:`BLKCB`, :cite:`Halfar83`
     - similarity solution
   * - C
     - isothermal SIA, flat bed, growing accumulation
     - :cite:`BLKCB`
     - similarity solution
   * - D
     - isothermal SIA, flat bed, oscillating accumulation
     - :cite:`BLKCB`
     - uses compensatory accumulation
   * - E
     - isothermal SIA; as A, but with sliding in a sector
     - :cite:`BLKCB`
     - uses compensatory accumulation
   * - F
     - thermomechanically coupled SIA (mass and energy conservation), steady, flat bed
     - :cite:`BB`, :cite:`BBL`
     - uses compensatory accumulation and heating
   * - G
     - thermomechanically coupled SIA; as F  but with oscillating accumulation
     - :cite:`BB`, :cite:`BBL`
     - ditto
   * - H
     - bed deformation coupled with isothermal SIA
     - :cite:`BLKfastearth`
     - joined similarity solution
   * - I
     - stream velocity computation using SSA (plastic till)
     - :cite:`SchoofStream`, :cite:`BBssasliding`
     -
   * - J
     - shelf velocity computation using SSA
     - (source code)
     -
   * - K
     - pure conduction in ice and bedrock
     - :cite:`BuelerTestK`
     -
   * - L
     - isothermal SIA, steady, non-flat bed
     - (source code)
     - numerical ODE solution

.. csv-table:: Canonical PISM verification runs using the exact solutions listed in
               :numref:`tab-tests`.
   :header: Test, Example invocation
   :name: tab-tests-exec
   :widths: auto

   A, ``pismv -test A -Mx 61 -My 61 -Mz 11 -y 25000``
   B, ``pismv -test B -Mx 61 -My 61 -Mz 11 -ys 422.45 -y 25000``
   C, ``pismv -test C -Mx 61 -My 61 -Mz 11 -y 15208.0``
   D, ``pismv -test D -Mx 61 -My 61 -Mz 11 -y 25000``
   E, ``pismv -test E -Mx 61 -My 61 -Mz 11 -y 25000``
   F, ``pismv -test F -Mx 61 -My 61 -Mz 61 -y 25000``
   G, ``pismv -test G -Mx 61 -My 61 -Mz 61 -y 25000``
   H, ``pismv -test H -Mx 61 -My 61 -Mz 11 -y 40034 -bed_def iso``
   I, ``ssa_testi -ssa_method fd -Mx 5 -My 500 -ssafd_picard_rtol 1e-6 -ssafd_ksp_rtol 1e-11``
   J, ``ssa_testj -ssa_method fd -Mx 60 -My 60 -ssafd_ksp_rtol 1e-12``
   K, ``pismv -test K -Mx 6 -My 6 -Mz 401 -Mbz 101 -y 130000``
   L, ``pismv -test L -Mx 61 -My 61 -Mz 31 -y 25000``

.. csv-table:: ``pismv`` command-line options
   :header: Option, Description
   :name: tab-pismv-options

   ``-test``,        Choose verification test by single character name; see :numref:`tab-tests`.
   ``-no_report``,   Do not report errors at the end of a verification run.
   ``-eo``,          Only evaluate the exact solution; no numerical approximation at all.
   ``-report_file``, Save error report to a netCDF file.
   ``-append``,      Append to a report file.

:numref:`tab-tests` summarizes the many exact solutions currently available in PISM. Most
of these exact solutions are solutions of *free boundary problems* for partial
differential equations; only Tests A, E, J, K are fixed boundary value problems.

:numref:`tab-tests-exec` shows how to run each of them on a coarse grids. Note that tests
I and J require special executables ``ssa_testi,ssa_testj`` which are built with
configuration flag ``Pism_BUILD_EXTRA_EXECS`` equal to ``ON``. :numref:`tab-pismv-options`
gives the special verification-related options of the ``pismv`` executable.

Numerical errors are not, however, the dominant reasons why ice sheet models give
imperfect results. The largest sources of errors include those from using the wrong (e.g.
over-simplified or incorrectly-parameterized) continuum model, and from observational or
pre-processing errors present in input data. Our focus here on numerical errors has a
model-maintenance goal. It is *easier* to maintain code by quantitatively confirming that
it produces small errors in cases where those can be measured, rather than "eyeballing"
results to see that they are "right" according to human judgment.

The goal of verification is not generally to see that the error is zero at any particular
resolution, or even to show that the error is small in a predetermined absolute sense.
Rather the goals are

- to see that the error *is* decreasing,
- to measure the rate at which it decreases, and
- to develop a sense of the magnitude of numerical error before doing realistic ice sheet
  model runs.

Knowing the error decay rate may give a prediction of how fine a grid is necessary to
achieve a desired smallness for the numerical error.

Therefore one must "go down" a grid refinement "path" and measure numerical error for each
grid :cite:`Roache`. The refinement path is defined by a sequence of spatial grid cell sizes
which decrease toward the refinement limit of zero size :cite:`MortonMayers`. In PISM the
timestep `\dt` is determined adaptively by a stability criterion (see
section :ref:`sec-adapt`). In PISM one specifies the number of grid points, thus the
grid cell sizes because the overall dimensions of the computational box are normally
fixed; see section :ref:`sec-coords`. By "measuring the error for each grid" we mean
computing a norm (or norms) of the difference between the numerical solution and the exact
solution.

For a grid refinement path example, in tests of the thermomechanically-coupled SIA model
one refines in three dimensions, and these runs produced Figures 13, 14, and 15 of :cite:`BBL`:

.. code-block:: none

   pismv -test G -max_dt 10.0 -y 25000 -Mx 61 -My 61 -Mz 61 -z_spacing equal
   pismv -test G -max_dt 10.0 -y 25000 -Mx 91 -My 91 -Mz 91 -z_spacing equal
   pismv -test G -max_dt 10.0 -y 25000 -Mx 121 -My 121 -Mz 121 -z_spacing equal
   pismv -test G -max_dt 10.0 -y 25000 -Mx 181 -My 181 -Mz 181 -z_spacing equal
   pismv -test G -max_dt 10.0 -y 25000 -Mx 241 -My 241 -Mz 241 -z_spacing equal
   pismv -test G -max_dt 10.0 -y 25000 -Mx 361 -My 361 -Mz 361 -z_spacing equal

The last two runs require a supercomputer! In fact the :math:`361\times 361\times 361` run
involves more than :math:`100` million unknowns, updated at each of millions of time
steps. Appropriate use of parallelism (``mpiexec -n NN pismv``) and of the ``-skip``
modification to adaptive timestepping accelerates such fine-grid runs; see section
:ref:`sec-adapt`.

Figures :numref:`fig-thickerrsB` through :numref:`fig-velerrsI` in
:ref:`sec-convergence-plots` show a sampling of the results of verifying PISM using the
tests described above. These figures were produced automatically using Python scripts
``test/vfnow.py`` and ``test/vnreport.py``. See section :ref:`sec-scripts`.

These figures *do not* show outstanding rates of convergence, relative to textbook partial
differential equation examples. For the errors in tests B and G, see the discussion of
free margin shape in :cite:`BLKCB`. For the errors in test I, the exact continuum solution is
not very smooth at the free boundary :cite:`SchoofStream`.

.. toctree::

   convergence-figures.rst
