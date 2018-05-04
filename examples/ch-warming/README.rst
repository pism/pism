.. default-role:: math

One-column cryo-hydrologic warming setup
========================================

This directory contains a set up of the one-column cryo-hydrologic warming experiment
inspired by **Phillips, T. and Rajaram, H. and Steffen, K.**, *Cryo-hydrologic warming: A
potential mechanism for rapid thermal response of ice sheets*, Geophy. Res. Lett., 2010.

In summary: in ablation areas the water generated during the melt season drains to the
base through the channels in the cryo-hydrologic system (CHS). The presence of water
in the CHS leads to warming of the ice adjacent to the channels.

The authors assume that the additional heat flux into the ice due to this mechanism is
proportional to the difference in temperature between the CHS and the ice next to it:

.. math::

   Q = \frac{k_{i}}{R^2} (T_{CH} - T_i)

Here `k_i` is the thermal conductivity of ice and `R` is the spacing between the
channels in the CHS (a poorly-constrained model parameter).

*Phillips et al* model `T_i` using a familiar temperature-based energy balance equation
with a parameterization of cooling due to horizontal advection. The flux `Q` defined above
appears as an additional source term:

.. math::

   \frac{\partial \rho_i C_i T_i}{\partial t} + \frac{\partial (\rho_i C_i w T_i)}{\partial z}
   - \frac{\partial}{\partial z}\left(k_i \frac{\partial T_i}{\partial z}\right)
   = - \frac{\partial \rho_i C_i u T_i}{\partial x} + Q_{s} + Q

Here `\rho_i` is the density of ice, `u` and `w` are horizontal and vertical velocity
components, `C_i` is the specific heat capacity of ice, and `Q_s` is the strain heating.

The temperature in the CHS is modeled using an enthalpy-based energy equation omitting
advection terms:

.. math::

   \frac{\partial \bar{\rho H}}{\partial t}
   - \frac{\partial}{\partial z}\left(\bar{k} \frac{\partial T_{CH}}{\partial z}\right)
   = - Q

Here the enthalpy is defined by

.. math::

   \rho H = (1 - \phi_w) \rho_i C_i T_{CH} + \phi_w \rho_w C_w T_{CH} + L,

where `\phi_w` is the water fraction, `L` is the latent heat of fusion, `\rho_w`, is the
water density, and `C_w` is the specific heat capacity of fresh water.

The thermal conductivity of the mixture, `\bar k`, is defined by

.. math::

   \bar k = (1 - \phi_w) k_i + \phi_w k_w.

The authors use Dirichlet (temperature) boundary conditions at the top boundary (assuming
that the temperature at the top of the snow or ice is equal to the near-surface air
temperature) and Dirichlet or Neumann (mentioned but not not specified in the paper) at the
bottom boundary.

It is also mentioned that they account *for snow cover, which insulates the ice surface
from cold winter temperatures* and *temperature dependence of the thermal properties was
incorporated in the model*, but no details are provided.

During the melt season (when the water is present) the temperature in the CHS is set to
the pressure-melting and the water fraction is set to a pre-determined constant (they use
`0.005`, i.e. one half of a percent) while during the winter the CHS is allowed to cool.

Comparing to PISM's energy conservation model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PISM's energy conservation model includes all the mechanisms modeled by equations above
*except* for the parameterization of advection-driven cooling not necessary in a 3D model.

PISM's enthalpy-based energy conservation model can be used to model *both* the ice and
the CHS columns.

PISM assumes that the specific heat capacity of ice is constant. The thermal conductivity
can be a function of temperature, but does not depend on the water fraction (unlike `\bar
k` above). Note, though, that for water fractions at or below `0.005` the value of `\bar
k` is very close to `k_i`.

By default PISM "ignores" thermal conductivity of temperate ice *without* completely
removing it: this can be thought of as a *regularization*. (The thermal conductivity of
temperate ice is set to a given fraction of `k_i` (see
``energy.temperate_ice_enthalpy_conductivity_ratio``). In this setup the ice column
remains below the pressure-melting point, so setting this parameter to `1` gets us close
to the equation modeling CHS enthalpy evolution in *Phillips et al*.

Technical details
^^^^^^^^^^^^^^^^^

FIXME

Results
^^^^^^^

FIXME
