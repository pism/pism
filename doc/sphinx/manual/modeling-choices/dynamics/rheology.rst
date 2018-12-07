.. include:: ../../../global.txt

.. _sec-rheology:

Ice rheology
------------

.. contents::

The "rheology" of a viscous fluid refers to the relation between the applied stress and
the resulting deformation, the strain rate. The models of ice rheology available in PISM
are all isotropic :cite:`Paterson`. A rheology in this class is described by a "flow law",
which is, in the most general case in PISM, a function `F(\sigma,T,\omega,P,d)` in the
"constitutive relation" form

.. math::
   :label: eq-constitutive

   D_{ij} = F(\sigma,T,\omega,P,d)\, \sigma_{ij}'.

Here `D_{ij}` is the strain rate tensor, `\sigma_{ij}'` is the stress deviator tensor, `T`
is the ice temperature, `\omega` is the liquid water fraction, `P` is the pressure, `d` is
the grain size, and `\sigma^2 = \frac{1}{2} \|\sigma_{ij}'\|_F = \frac{1}{2} \sigma_{ij}'
\sigma_{ij}'` defines the second invariant `\sigma` of the stress deviator tensor.

Form :eq:`eq-constitutive` of the flow law is used in the SIA, but the "viscosity" form of
a flow law, found by inverting the constitutive relation :eq:`eq-constitutive`, is needed
for ice shelf and ice stream (SSA) flow :cite:`BBssasliding`:

.. math::
   :label: eq-viscosityform

   \sigma_{ij}' = 2 \nu(D,T,\omega,P,d)\,D_{ij}

Here `\nu(D,T,\omega,P,d)` is the "effective viscosity" and `D^2 = \frac{1}{2}
D_{ij} D_{ij}`.

Most of the flow laws in PISM are of Glen-Nye single-power type.  For example,

.. math::
   :label: eq-glen

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
:ref:`sec-pism-defaults` and the implementation of flow laws in the `Source Code Browser
<pism-browser_>`_. Note that different flow laws have different numbers of parameters, but
all have at least two parameters (e.g. `A_0` and `n` in ``isothermal_glen``). One can
create a new, and reasonably arbitrarily, scalar function `F` by modifying source code;
see source files in ``src/base/rheology/``.

Choosing the flow laws for SIA and SSA stress balances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Command-line options :opt:`-sia_flow_law` and :opt:`-ssa_flow_law` choose which flow law
is used by the SIA and SSA stress balances, respectively. Allowed arguments are listed in
:numref:`tab-flowlaw` below. Viscosity form :eq:`eq-viscosityform` is not known for the
Goldsby-Kohlstedt law :cite:`GoldsbyKohlstedt`, so option "``-ssa_flow_law gk``" is an
error.

.. list-table:: Single-power flow laws. Choose the ice rheology using ``-sia_flow_law``
                and ``-ssa_flow_law`` and one of the names in this table. Flow law choices
                other than ``gpbld`` do not use the liquid water fraction `\omega`
                but only the temperature `T`.
   :name: tab-flowlaw
   :header-rows: 1
   :widths: 1,3

   * - Name
     - Comments and References

   * - ``gpbld``
     - Glen-Paterson-Budd-Lliboutry-Duval law :cite:`LliboutryDuval1985`, the enthalpy-based
       default in PISM :cite:`AschwandenBuelerKhroulevBlatter`. Extends the Paterson-Budd law
       (below) to positive liquid water fraction. If `A_{c}(T)` is from Paterson-Budd then
       this law returns

          `A(T,\omega) = A_{c}(T) (1 + C \omega),`

       where `\omega` is the liquid water fraction, `C` is a configuration parameter
       :config:`flow_law.gpbld.water_frac_coeff` [default `C=181.25`\], and `\omega` is
       capped at level :config:`flow_law.gpbld.water_frac_observed_limit`.

   * - ``pb``
     - Paterson-Budd law, the cold-mode default. Fixed Glen exponent `n=3`. Has a split
       "Arrhenius" term `A(T) = A \exp(-Q/RT^*)` where

          `A = 3.615 \times 10^{-13}\, \text{s}^{-1}\, \text{Pa}^{-3},`
          `Q = 6.0 \times 10^4\, \text{J}\, \text{mol}^{-1}`

       if `T^* < 263` K and
       
          `A = 1.733 \times 10^{3}\, \text{s}^{-1}\, \text{Pa}^{-3},`
          `Q = 13.9 \times 10^4\, \text{J}\, \text{mol}^{-1}`

       if `T^* > 263` K.

       Here `T^*` is pressure-adjusted temperature :cite:`PatersonBudd`.

   * - ``arr``
     - *Cold* part of Paterson-Budd. Regardless of temperature, the `A` and `Q` values for
       `T^*<263` K in the Paterson-Budd law apply. This is the flow law used in the
       thermomechanically-coupled exact solutions run by ``pismv -test F`` and
       ``pismv -test G`` :cite:`BBL`, :cite:`BB`.

   * - ``arrwarm``
     - *Warm* part of Paterson-Budd. Regardless of temperature, the `A` and `Q` values for
       `T^*>263` K in Paterson-Budd apply.

   * - ``hooke``
     - Hooke law with

          `A(T) = A \exp(-Q/(RT^*) + 3C (T_r - T^*)^\kappa).`

       Fixed Glen exponent `n=3` and constants as in :cite:`Hooke`, :cite:`PayneBaldwin`.

   * - ``isothermal_glen``
     - The isothermal Glen flow law. Here `F(\sigma) = A_0 \sigma^{n-1}` with inverse
       `\nu(D) = \frac{1}{2} B_0 D^{(1-n)/(2n)}` where `A_0` is the ice softness and
       `B_0=A_0^{-1/n}` is the ice hardness.

   * - ``gk``
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
   :label: eq-renewexponent

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
   :widths: 1,2,2

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
