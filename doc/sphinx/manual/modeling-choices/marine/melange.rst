.. include:: ../../../global.txt

.. _sec-model-melange-pressure:

Modeling melange back-pressure
------------------------------

Equation :eq:`eq-cfbc` above, describing the stress boundary condition for ice shelves,
can be written in terms of velocity components:

.. math::
   :label: eq-cfbc-uv

   2 \nu H (2u_x + u_y) \nx + 2 \nu H (u_y + v_x)  \ny &= \displaystyle \int_{b}^{h}(\pice(z) - \psw(z)) dz\, \nx,

   2 \nu H (u_y + v_x)  \nx + 2 \nu H (2v_y + u_x) \ny &= \displaystyle \int_{b}^{h}(\pice(z) - \psw(z)) dz\, \ny.

Here `\nu` is the vertically-averaged ice viscosity, `H` is the ice thickness, `b` is the
elevation of the bottom and `h` of the top ice surface, `\psw` and `\pice` are pressures
of the column of ice and water, respectively:

.. math::
   :label: eq-cfbc-pressures

   \pice &= \rhoi\, g (h - z),

   \psw &= \rhow\, g\, \max(\zs - z,\, 0).

We call the integral on the right hand side of :eq:`eq-cfbc-uv` the *pressure difference
term*.

It can be re-written as

.. math::
   :label: eq-cfbc-pressure-difference

   \int_b^h \pice(z) - \psw(z) dz &= H (\bar p_{\text{ice}} - \bar p_{\text{water}}),\, \text{where}

   \bar p_{\text{ice}} &= \frac12\, \rho_{\text{ice}}\, g\, H,

   \bar p_{\text{water}} &= \frac12\, \rho_{\text{water}}\, g\, \frac{\max(z_s - b, 0)^2}{H}.

PISM's ocean model components provide `\bar p_{\text{water}}`, the vertically-averaged
pressure of the water column adjacent to an ice margin.

To model the effect of melange :cite:`Amundsonetal2010` on the stress boundary condition
we modify the pressure difference term in :eq:`eq-cfbc-uv`, adding `\pmelange`, the
vertically-averaged melange back pressure:

.. math::
   :label: eq-cfbc-3

   \int_{b}^{h}(\pice - (\psw + \pmelange))\, dz.

By default, `\pmelange` is zero, but PISM implements two ocean model components to support
scalar time-dependent melange pressure forcing. Please see the :ref:`Climate Forcing
Manual <sec-ocean-frac-mbp>` for details.
