.. include:: ../global.txt

.. _sec-ocean-th-details:

Three-equation ocean model (implementation details)
===================================================
This model is based on :cite:`HellmerOlbers1989` and :cite:`HollandJenkins1999`.

We use equations for the heat and salt flux balance at the base of the shelf to compute
the temperature at the base of the shelf and the sub-shelf melt rate.

Following :cite:`HellmerOlbers1989`, let `Q_T` be the total heat flux crossing the
interface between the shelf base and the ocean, `Q_T^B` be the amount of heat lost by the
ocean due to melting of glacial ice, and `Q_T^I` be the conductive flux into the ice
column.

`Q_{T}` is parameterized by (see :cite:`HellmerOlbers1989`, equation 10):

.. math::

   Q_{T} = \rho_W\, c_{pW}\, \gamma_{T}\, (T^B - T^W),

where `\rho_{W}` is the sea water density, `c_{pW}` is the heat capacity of sea water, and
`\gamma_{T}` is a turbulent heat exchange coefficient.

We assume that the difference between the basal temperature and adjacent ocean temperature
`T^B - T^W` is well approximated by `\Theta_B - \Theta_W,` where `\Theta_{\cdot}` is the
corresponding potential temperature.

`Q_T^B` is (see :cite:`HellmerOlbers1989`, equation 11):

.. math::

   Q_T^B = \rho_I\, L\, \frac{\partial h}{\partial t},

where `\rho_I` is the ice density, `L` is the latent heat of fusion, and `\partial
h\, /\, \partial t` is the ice thickening rate (equal to minus the melt rate).

The conductive flux into the ice column is (:cite:`Hellmeretal1998`, equation 7):

.. math::

   Q_T^I = \rho_I\, c_{pI}\, \kappa\, T_{\text{grad}},

where `\rho_I` is the ice density, `c_{pI}` is the heat capacity of ice, `\kappa` is the
ice thermal diffusivity, and `T_{\text{grad}}` is the vertical temperature gradient at the
base of a column of ice.

Now, the heat flux balance implies

.. math::

   Q_T = Q_T^B + Q_T^I.

For the salt flux balance, we have

.. math::

   Q_S = Q_S^B + Q_S^I,

where `Q_S` is the total salt flux across the interface, `Q_S^B` is the basal salt flux
(negative for melting), `Q_S^I = 0` is the salt flux due to molecular diffusion of salt
through ice.

`Q_S` is parameterized by (:cite:`Hellmeretal1998`, equation 13)

.. math::

   Q_S = \rho_W\, \gamma_S\, (S^B - S^W),

where `\gamma_S` is a turbulent salt exchange coefficient, `S^B` is salinity at the shelf
base, and `S^W` is the salinity of adjacent ocean.

The basal salt flux `Q_S^B` is (:cite:`Hellmeretal1998`, equation 10)

.. math::

   Q_S^B = \rho_I\, S^B\, {\frac{\partial h}{\partial t}}.

To avoid converting shelf base temperature to shelf base potential temperature and back,
we use two linearizations of the freezing point equation for sea water for in-situ and for
potential temperature, respectively:

.. math::

   T^{B}(S,h) &= a_0\cdot S + a_1 + a_2\cdot h,

   \Theta^{B}(S,h) &= b_0\cdot S + b_1 + b_2\cdot h,

where `S` is salinity and `h` is ice shelf thickness.

.. note::

   The linearized equation for the freezing point of seawater as a function of salinity
   and pressure (ice thickness) is only valid for salinity ranges from 4 to 40 psu (see
   :cite:`HollandJenkins1999`).

The linearization coefficients for the basal temperature `T^B(S,h)` are taken from
:cite:`Hellmeretal1998`, going back to :cite:`FoldvikKvinge1974`.

Given `T^B(S,h)` and a function `\Theta_T^B(T)` one can define
`\Theta^B_{*}(S,h) = \Theta_T^B\left(T^B(S,h)\right)`.

The parameterization `\Theta^B(S,h)` used here was produced by linearizing
`\Theta^B_{*}(S,h)` near the melting point. (The definition of `\Theta_T^B(T)`, converting
in situ temperature into potential temperature, was adopted from FESOM
:cite:`Wangetal2013`).

Treating ice thickness, sea water salinity, and sea water potential temperature as "known"
and choosing an approximation of the temperature gradient at the base `T_{\text{grad}}`
(see below) we can write down a system of equations

.. math::

   Q_T &= Q_T^B + Q_T^I,

   Q_S &= Q_S^B + Q_S^I,

   T^{B}(S,h) &= a_0\cdot S + a_1 + a_2\cdot h,

   \Theta^{B}(S,h) &= b_0\cdot S + b_1 + b_2\cdot h

and simplify it to produce a quadratic equation for the salinity at the shelf base, `S^B`:

.. math::

   A\cdot (S^B)^2 + B\cdot S^B + C = 0.

The coefficients `A,` `B`, and `C` depend on the basal temperature gradient approximation
for the sub-shelf melt, sub-shelf freeze-on, and diffusion-only cases.

- Melt at the base:

  .. math::

     T_{\text{grad}} = -\Delta T\, \frac{\partial h\, /\, \partial t}{\kappa}.

  See equation 13 in :cite:`Hellmeretal1998`.
- Freeze on at the base: we assume that

  .. math::

     T_{\text{grad}} = 0.
- No melt and no freeze on:

  .. math::

     T_{\text{grad}} = \frac{\Delta T}{h}.

  See :cite:`HollandJenkins1999`, equation 21.

One remaining problem is that we cannot compute the basal melt rate without making an
assumption about whether there is basal melt or not, and cannot pick one of the three
cases without computing the basal melt rate first.

Our implementation tries to compute basal salinity that is consistent with the
corresponding basal melt rate. See the code for details.

Once `S_B` is found by solving this quadratic equation, we can compute the basal
temperature using the parameterization for `T^{B}(S,h)`.

To find the basal melt rate, we solve the salt flux balance equation for `{\frac{\partial
h}{\partial t}},` obtaining

.. math::

   w_b = -\frac{\partial h}{\partial t} = \frac{\gamma_S\, \rho_W\, (S^W - S^B)}{\rho_I\, S^B}.


See :ref:`sec-ocean-th` for the user's documentation of this model.
