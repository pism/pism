.. include:: ../../global.rst

.. _sec-modeling-dynamics:

Modeling choices: Ice dynamics and thermodynamics
=================================================

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

The SSA is also used in PISM to describe the sliding of grounded ice and the formation of
ice streams :cite:`BBssasliding`. Specifically for the SSA with "plastic" (Coulomb friction)
basal resistance, the locations of ice streams are determined as part of a free boundary
problem of Schoof :cite:`SchoofStream`, a model for emergent ice streams within a ice sheet and
ice shelf system. This model explains ice streams through a combination of plastic till
failure and SSA stress balance.

This SSA description of ice streams is, however, also the preferred "sliding law" for the
SIA :cite:`BBssasliding`, :cite:`Winkelmannetal2011`. The SSA should be combined with the SIA, in
this way, in preference to classical SIA sliding laws which make ice basal velocity a
local function of the basal value of the driving stress. The resulting combination of SIA
and SSA is a "hybrid" approximation of the Stokes model :cite:`Winkelmannetal2011`. Option
``-stress_balance ssa+sia`` turns on this "hybrid" model. In this use of the SSA as a
sliding law, floating ice is also subject to the SSA.

Of course there is more to the use of a stress balance than just turning it on! At all
grounded points a yield stress, or a pseudo-yield-stress in the case of power law sliding
(section :ref:`sec-basestrength`), is computed from the amount of stored basal water
and from a (generally) spatially-varying till strength. The amount of stored basal water
is modeled by the subglacial hydrology mode choice (section :ref:`sec-subhydro`) based
on the basal melt rate which is, primarily, thermodynamically-determined (subsection
:ref:`sec-basestrength`).

:numref:`tab-stress-balance-choice` describes the basic choice of stress balance. If the
SSA stress balance is used, a choice of two solvers is available, namely ``-ssa_method
fd`` (default) or ``-ssa_method fem``. See :numref:`tab-ssa-usage`, which describes
additional controls on the numerical solution of the stress balance equations. If option
``-ssa_method fd`` is chosen then several more controls on numerics are available; see
:numref:`tab-ssafd-controls`. If the ice sheet being modeled has any floating ice then
the user is advised to read section :ref:`sec-pism-pik` on modeling marine ice sheets.

.. list-table:: The basic choice of stress balance
   :name: tab-stress-balance-choice
   :header-rows: 1

   * - Option
     - Description

   * - :opt:`-stress_balance none`
     - Turn off ice flow completely.

   * - :opt:`-stress_balance sia` (default)
     - Grounded ice flows by the non-sliding SIA. Floating ice essentially doesn't flow,
       so this model is not recommended for marine ice sheets.

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

.. list-table:: Choice of, and controls on, the numerical SSA stress balance.
   :name: tab-ssa-usage
   :header-rows: 1

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
   :widths: 20, 80

   * - Option
     - Description
     
   * - :opt:`-ssa_maxi` (300)
     - Set the maximum allowed number of Picard (nonlinear) iterations in solving the
       shallow shelf approximation.

       .. FIXME: this should be "ssafd_picard_maxi"?

   * - :opt:`-ssa_rtol` (`10^{-4}`)
     - The Picard iteration computes a vertically-averaged effective viscosity which is
       used to solve the equations for horizontal velocity. Then the new velocities are
       used to recompute an effective viscosity, and so on. This option sets the relative
       change tolerance for the effective viscosity. The Picard iteration stops when
       successive values `\nu^{(k)}` of the vertically-averaged effective viscosity
       satisfy

       .. FIXME: this should be "ssafd_picard_rtol"?

       .. math::

          \|(\nu^{(k)} - \nu^{(k-1)}) H\|_1 \le Z \|\nu^{(k)} H\|_1

       where `Z=` ``ssa_rtol``. 

   * - :opt:`-ssafd_ksp_rtol` (`10^{-5}`)
     - Set the relative change tolerance for the iteration inside the Krylov linear solver
       used at each Picard iteration.



.. _sec-rheology:

Ice rheology
------------


The "rheology" of a viscous fluid refers to the relation between the applied stress and the resulting deformation, the strain rate.  The models of ice rheology available in PISM are all isotropic :cite:`Paterson`.   A rheology in this class is described by a "flow law", which is, in the most general case in PISM, a function `F(\sigma,T,\omega,P,d)` in the "constitutive relation" form

.. math::
   :name: eq-constitutive

   D_{ij} = F(\sigma,T,\omega,P,d)\, \sigma_{ij}'.

Here `D_{ij}` is the strain rate tensor, `\sigma_{ij}'` is the stress deviator tensor, `T` is the ice temperature, `\omega` is the liquid water fraction, `P` is the pressure, `d` is the grain size, and `\sigma^2 = \frac{1}{2} \|\sigma_{ij}'\|_F = \frac{1}{2} \sigma_{ij}' \sigma_{ij}'` defines the second invariant `\sigma` of the stress deviator tensor.

Form :eq:`eq-constitutive` of the flow law is used in the SIA, but the "viscosity" form of a flow law, found by inverting the constitutive relation :eq:`eq-constitutive`, is needed for ice shelf and ice stream (SSA) flow :cite:`BBssasliding`:

.. math::
   :name: eq-viscosityform

   \sigma_{ij}' = 2 \nu(D,T,\omega,P,d)\,D_{ij}

Here `\nu(D,T,\omega,P,d)` is the "effective viscosity" and `D^2 = \frac{1}{2}
D_{ij} D_{ij}`.

Most of the flow laws in PISM are of Glen-Nye single-power type.  For example,

.. math::
   :name: eq-glen

   F(\sigma,T) = A(T) \sigma^{n-1}

is the common temperature-dependent Glen law :cite:`PatersonBudd`, :cite:`BBL` (which has no
dependence on liquid water fraction, pressure, or grain size). If the ice softness
`A(T)=A_0` is constant then the law is isothermal, whereas if there is dependence on
temperature then `A(T)` is usually a generalization of "Arrhenius" form

.. math::

   A(T) = A \exp(-Q/(R T)).

The more elaborate Goldsby-Kohlstedt law :cite:`GoldsbyKohlstedt` is a function
`F(\sigma,T,P,d)`, but in this case the function `F` cannot be factored into a
product of a function of `T,P,d` and a single power of `\sigma`, as in form
:eq:`eq-glen`.

There is only one choice for the flow law which takes full advantage of the enthalpy mode
of PISM, which is the thermodynamical modeling (i.e. conservation of energy) default.
Namely the Glen-Paterson-Budd-Lliboutry-Duval flow law
:cite:`AschwandenBuelerKhroulevBlatter`, :cite:`LliboutryDuval1985`, :cite:`PatersonBudd`,
which is a function `F(\sigma,T,\omega,P)`. This law is the only one in the literature
where the ice softness depends on both the temperature and the liquid water fraction, so
it parameterizes the (observed) softening of pressure-melting-temperature ice as its
liquid fraction increases. One can use this default polythermal law or one may choose
among a number of "cold ice" laws listed in :numref:`tab-flowlaw` which do not use the
liquid water fraction.

All flow law parameters can be changed using configuration parameters; see section
:ref:`sec-pism-defaults` and the implementation of flow laws in the \emph{Source Code
Browser}. Note that different flow laws have different numbers of parameters, but all have
at least two parameters (e.g. `A_0` and `n` in ``isothermal_glen``). One can
create a new, and reasonably arbitrarily, scalar function `F` by modifying source
code; see source files ``flowlaws.hh``, ``flowlaws.cc`` in ``src/base/rheology/``. To
assist such modifications, note that :numref:`tab-flowlaw` below also lists the C++
classes declared in ``flowlaw.hh``.

Choosing the flow laws for SIA and SSA stress balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Command-line options :opt:`-sia_flow_law` and :opt:`-ssa_flow_law` choose which flow law
is used by the SIA and SSA stress balances, respectively. Allowed arguments are listed in
Tables :numref:`tab-flowlaw` and :numref:`tab-flowlawgk` below. Viscosity form
:eq:`eq-viscosityform` is not known for the Goldsby-Kohlstedt law :cite:`GoldsbyKohlstedt`,
so option "``-ssa_flow_law gk``" is an error.

.. list-table:: Single-power flow laws. Choose the ice rheology using ``-sia_flow_law``
                and ``-ssa_flow_law`` and one of the names in this table. Flow law choices
                other than ``gpbld`` do not use the liquid water fraction `\omega`
                but only the temperature `T`.
   :name: tab-flowlaw
   :header-rows: 1

   * - Name
     - C++ class
     - Comments and References

   * - ``gpbld``
     - :class:`rheology::GPBLD`
     - Glen-Paterson-Budd-Lliboutry-Duval law :cite:`LliboutryDuval1985`, the enthalpy-based
       default in PISM :cite:`AschwandenBuelerKhroulevBlatter`. Extends the Paterson-Budd law
       (below) to positive liquid water fraction. If `A_{c}(T)` is from Paterson-Budd then
       this law returns

       .. math::
       
          A(T,\omega) = A_{c}(T) (1 + C \omega),

       where `\omega` is the liquid water fraction, `C` is a configuration parameter
       :config:`flow_law.gpbld.water_frac_coeff` [default `C=181.25`\], and `\omega` is
       capped at level :config:`flow_law.gpbld.water_frac_observed_limit`.
       
   * - ``pb``
     - :class:`rheology::PatersonBudd`
     - Paterson-Budd law, the cold-mode default. Fixed Glen exponent `n=3`. Has a split
       "Arrhenius" term `A(T) = A \exp(-Q/RT^*)` where

       .. math::

          A &= 3.615 \times 10^{-13}\, \text{s}^{-1}\, \text{Pa}^{-3},

          Q &= 6.0 \times 10^4\, \text{J}\, \text{mol}^{-1}

       if `T^* < 263` K and

       .. math::

          A &= 1.733 \times 10^{3}\, \text{s}^{-1}\, \text{Pa}^{-3},

          Q &= 13.9 \times 10^4\, \text{J}\, \text{mol}^{-1}

       if `T^* > 263` K;

       here `T^*` is pressure-adjusted temperature :cite:`PatersonBudd`.
 
   * - ``arr``
     - :class:`rheology::PatersonBuddCold`
     - *Cold* part of Paterson-Budd. Regardless of temperature, the `A` and `Q` values for
       `T^*<263` K in the Paterson-Budd law apply. This is the flow law used in the
       thermomechanically-coupled exact solutions run by ``pismv -test F`` and
       ``pismv -test G`` :cite:`BBL`, :cite:`BB`.
       
   * - ``arrwarm``
     - :class:`rheology::PatersonBuddWarm`
     - *Warm* part of Paterson-Budd. Regardless of temperature, the `A` and `Q` values for
       `T^*>263` K in Paterson-Budd apply.
  
   * - ``hooke``
     - :class:`rheology::Hooke`
     - Hooke law with

       .. math::

          A(T) = A \exp(-Q/(RT^*) + 3C (T_r - T^*)^\kappa).

       Fixed Glen exponent `n=3` and constants as in :cite:`Hooke`, :cite:`PayneBaldwin`.
       
   * - ``isothermal_glen``
     - :class:`rheology::IsothermalGlen`
     - The isothermal Glen flow law. Here `F(\sigma) = A_0 \sigma^{n-1}` with inverse
       `\nu(D) = \frac{1}{2} B_0 D^{(1-n)/(2n)}` where `A_0` is the ice softness and
       `B_0=A_0^{-1/n}` is the ice hardness.


.. list-table:: The Goldsby-Kohlstedt flow law. Use option ``-sia_flow_law gk``
   :name: tab-flowlawgk
   :header-rows: 1

   * - Name
     - C++ class
     - Comments and References
   * - ``gk``
     - :class:`rheology::GoldsbyKohlstedt`
     - This law has a combination of exponents from `n=1.8` to `n=4`
       :cite:`GoldsbyKohlstedt`. It can only be used by the SIA stress balance. Because it has
       more than one power, option ``-sia_n`` has no effect, though ``-sia_e`` works as
       expected. This law does not use the liquid water fraction, but only the
       temperature.

Choose enhancement factor and exponent
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An enhancement factor can be added to any flow law through a runtime option. Single-power
laws also permit control of the flow law exponent through a runtime option.

Options :opt:`-sia_e` and :opt:`-ssa_e` set flow enhancement factors for the SIA and SSA
respectively. Option ``-sia_e`` sets "`e`" in `D_{ij} = e\, F(\sigma,T,\omega,P,d)\,
\sigma_{ij}',` in equation :eq:`eq-constitutive`. Option ``-ssa_e`` sets "`e`" in the
viscosity form so that `\sigma_{ij}' = e^{-1/n}\, 2\, \nu(D,T,\omega,P,d)\, D_{ij}.`

Options :opt:`-sia_n` and :opt:`-ssa_n` set the exponent when a single-power flow law is
used (see :numref:`tab-flowlaw`). Simply changing to a different value from the default
`n=3` is not recommended without a corresponding change to the enhancement factor,
however. This is because the coefficient and the power are non-trivially linked when a
power law is fit to experimental data :cite:`CuffeyPaterson`, :cite:`PatersonBudd`.

Here is a possible approach to adjusting both the enhancement factor and the exponent.
Suppose `\sigma_0` is preferred as a scale (reference) for the driving stress that
appears in both SIA and SSA models. Typically this is on the order of one bar or
`10^5` Pa. Suppose one wants the same amount of deformation `D_0` at this
reference driving stress as one changes from the old exponent `n_{old}` to the new
exponent `n_{new}`. That is, suppose one wants both

.. math::

   D_0 = E_{old}\, A\, \sigma_0^{n_{old}} \qquad \text{and} \qquad D_0
   = E_{new}\, A\, \sigma_0^{n_{new}}

to be true with a new enhancement factor `E_{new}`. Eliminating `D_0` and
solving for the new enhancement factor gives

.. math::
   :name: eq-renewexponent

   E_{new} = E_{old}\, \sigma_0^{n_{old} - n_{new}}.

It follows, for example, that if one has a run with values

.. code-block:: none

   -sia_e 3.0 -sia_n 3.0

then a new run with exponent `n=6.0` and the same deformation at the reference
driving stress of `10^5` Pa will use

.. code-block:: none

   -sia_e 3.0e-15 -sia_n 6.0

because `E_{new} = 3.0 \sigma_0^{3-6} = 3.0 \times (10^5)^{-3}` from equation
:eq:`eq-renewexponent`.

A corresponding formula applies to ``-ssa_e`` if the ``-ssa_n`` value changes.

.. list-table:: For all flow laws, an enhancement factor can be added by a runtime option.
                For the single-power flow laws in :numref:`tab-flowlaw`, the (Glen)
                exponent can be controlled by a runtime option.
   :name: tab-enhancementandexponent
   :header-rows: 1

   * - Option
     - Configuration parameter
     - Comments

   * - :opt:`-sia_e` (1.0)
     - ``stress_balance.sia.enhancement_factor``
     - Note (see the supplement of :cite:`AschwandenAdalgeirsdottirKhroulev`) used `3.0`
       for Greenland ice sheet simulations while :cite:`Martinetal2011` used `4.5` for
       simulations of the Antarctic ice sheet with PISM-PIK.

   * - :opt:`-sia_n` (3.0)
     - ``stress_balance.sia.Glen_exponent``
     - See text and eqn :eq:`eq-renewexponent` to also set ``-sia_e`` if ``-sia_n`` changes.

   * - :opt:`-ssa_e` (1.0)
     - ``stress_balance.ssa.enhancement_factor``
     - Note :cite:`Martinetal2011` used `0.512` for simulations of the Antarctic ice sheet with
       PISM-PIK.

   * - :opt:`-ssa_n` (3.0)
     - ``stress_balance.ssa.Glen_exponent``
     - See text and eqn :eq:`eq-renewexponent` to also set ``-ssa_e`` if ``-ssa_n``
       changes.

.. _sec-gradient:

Surface gradient method
-----------------------


PISM computes surface gradients to determine the "driving stress"

.. math::

   (\tau_{d,x},\tau_{d,y}) = - \rho g H \nabla h,

where `H` is the ice thickness, and `h = H+b` is the ice surface elevation.
The driving stress enters into both the SIA and SSA stress balances, but in the former the
driving stress is needed on a staggered grid, while in the latter the driving stress is
needed on the regular grid.

Surface gradients are computed by finite differences in several slightly-different ways.
There are options for choosing which method to use, but to the best of our knowledge there
is no theoretical advice on the best, most robust mechanism. There are three
:opt:`-gradient` methods in PISM:

.. list-table:: Options controlling the surface gradient computation in the SIA code
   :name: tab-sia-gradient
   :header-rows: 1

   * - Option
     - Description

   * - :opt:`-gradient mahaffy`
     - This most "standard" way computes the surface slope onto the staggered grid for the
       SIA :cite:`Mahaffy`. It makes `O(\Delta x^2,\Delta y^2)` errors. For computations
       of driving stress on the regular grid, centered differencing is used instead.

   * - :opt:`-gradient haseloff`
     - This is the default method. It only differs from ``mahaffy`` at ice-margin
       locations, where it alters the formula for the slope in cases where an adjacent
       ice-free bedrock surface elevation is above the ice elevation.

   * - :opt:`-gradient eta`
     - In this method we first transform the thickness `H` by `\eta =
       H^{(2n+2)/n}` and then differentiate the sum of the thickness and the bed using
       centered differences:

       .. math::

          \nabla h = \nabla H + \nabla b = \frac{n}{(2n+2)}
          \eta^{(-n-2)/(2n+2)} \nabla \eta + \nabla b.

       Here `b` is the bed elevation and `h` is the surface elevation. This
       transformation sometimes has the benefits that the surface values of the horizontal
       velocity and vertical velocity, and the driving stress, are better behaved near the
       margin. See :cite:`BLKCB` for technical explanation of this transformation and compare
       :cite:`SaitoMargin`. The actual finite difference schemes applied to compute the surface
       slope are similar to option ``mahaffy``.

.. _sec-energy:

Modeling conservation of energy
-------------------------------

In normal use PISM solves the conservation of energy problem within the ice, the thin
subglacial layer, and a layer of thermal bedrock. For the ice and the subglacial layer it
uses an enthalpy-based scheme :cite:`AschwandenBuelerKhroulevBlatter` which allows the energy
to be conserved even when the temperature is at the pressure-melting point.

Ice at the melting point is called "temperate" ice. Part of the thermal energy of
temperate ice is in the latent heat of the liquid water stored between the crystals of the
temperate ice. Part of the thermal energy of the whole glacier is in the latent heat of
the liquid water under the glacier. The enthalpy scheme correctly models these storehouses
of thermal energy, and thus it allows polythermal and fully-temperate glaciers to be
modeled :cite:`AschwandenBlatter`.

The state of the full conservation of energy model includes the 3D ``enthalpy`` variable
plus the 2D ``bwat`` and ``tillwat`` subglacial hydrology state variables (subsection
:ref:`sec-subhydro`), all of which are seen in output files. The important basal melt rate
computation involves all of these energy state variables, because the basal melt rate
(``bmelt`` in output files) comes from conserving energy across the ice-bedrock layer
:cite:`AschwandenBuelerKhroulevBlatter`. Fields ``temp``, ``liqfrac``, and ``temp_pa`` seen in
output files are all actually diagnostic outputs because all of these can be recovered
from the enthalpy and the ice geometry.

Because this part of PISM is just a conservation law, there is little need for the user to
worry about controlling it. If desired, however, conservation of energy can be turned off
entirely with :opt:`-energy none`. The default enthalpy-based conservation of energy model
(i.e. ``-energy enthalpy``) can be replaced by the temperature-based (i.e. "cold ice")
method used in :cite:`BBssasliding` and verified in :cite:`BBL` by setting option :opt:`-energy
cold`.

The thermal bedrock layer model is turned off by setting ``-Mbz 1`` (i.e. zero spaces)
while it is turned on by choosing a depth and number of points, as in ``-Lbz 1000 -Mbz
21``, for example, which gives a layer depth of 1000 m and grid spaces of 50 m (=
1000/20). The input geothermal flux (``bheatflx`` in output files) is applied at the
bottom of the bedrock thermal layer if such a layer is present and otherwise it is applied
at the base of the ice.

.. _sec-age:

Computing ice age
-----------------

By default, PISM does not compute the age of the ice because it does not directly impact
ice flow when using the default flow laws. It is very easy to turn on. Just set
:opt:`-age`. A 3D variable ``age`` will appear in output files. It is read at input if
``-age`` is set and otherwise it is ignored even if present in the input file. If ``-age``
is set and the variable ``age`` is absent in the input file then the initial age is set to
zero.

