.. include:: ../../../global.txt

.. _sec-ssa:

Shallow shelf approximation (SSA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Parameters
##########

.. pism-parameters::
   :prefix: stress_balance.ssa.

.. _sec-ssa-remarks:

Technical remarks
#################

Ice stress balance models ignoring inertia are "diagnostic" models that do not have a
"state" evolving through time: the ice velocity is fully determined by ice geometry, basal
boundary conditions, and the ice viscosity.

In addition to this, shallow stress balance models (other than the SIA) correspond to
*nonlinear* systems of equations, which means that computing an estimate of the ice
velocity is an iterative process starting with a particular *initial guess*.

The quality of this guess may affect the number of iterations needed and even whether the
solver succeeds at all. It also has an influence on the computed solution: if an equation
has a unique solution, estimates produced using different initial guesses should be
*close*, but they need not be *identical*. If an equation has multiple solutions, even a
"small" change in the initial guess may give a completely different result.

In the context of a evolutionary model we can usually assume that the change in the state
of the model from one step to the next is "small" and we use (`u, v`) estimates from one
time step as the initial guess for the next. To ensure that stopping and re-starting a
simulation does not affect results we save these to the output file as variables
:var:`u_ssa` and :var:`v_ssa` and read them in when re-starting a stopped simulation.

One could say that the continuum SSA model does not have a state, but its implementation
does. Set :config:`stress_balance.ssa.read_initial_guess` to "false" to ignore it during
initialization and use the zero initial guess instead.
