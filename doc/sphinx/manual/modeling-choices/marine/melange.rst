.. include:: ../../../global.rst

.. _sec-model-melange-pressure:

Modeling melange back-pressure
------------------------------


Equation :eq:`eq-cfbc` above, describing the stress boundary condition for ice shelves,
can be written in terms of velocity components:

.. math::
   :name: eq-cfbc-uv

   \newcommand{\psw}{p_{\text{ocean}}}
   \newcommand{\pice}{p_{\text{ice}}}
   \newcommand{\pmelange}{p_{\text{melange}}}
   \newcommand{\n}{\mathbf{n}}
   \newcommand{\nx}{\n_{x}}
   \newcommand{\ny}{\n_{y}}
   
   \begin{array}{lclcl}
     2 \nu H (2u_x + u_y) \nx &+& 2 \nu H (u_y + v_x)  \ny &=& \displaystyle \int_{b}^{h}(\pice - \psw) dz\, \nx,\\
     2 \nu H (u_y + v_x)  \nx &+& 2 \nu H (2v_y + u_x) \ny &=& \displaystyle \int_{b}^{h}(\pice - \psw) dz\, \ny.
   \end{array}

Here `\nu` is the vertically-averaged ice viscosity, `b` is the ice base elevation, `h` is
the ice top surface elevation, and `\psw` and `\pice` are pressures of the column of sea
water and ice, respectively.

We call the integral on the right hand side of :eq:`eq-cfbc-uv` the "pressure imbalance
term". To model the effect of melange :cite:`Amundsonetal2010` on the stress boundary
condition, we assume that the melange back-pressure `\pmelange` does not exceed `\pice -
\psw`. Therefore we introduce `\lambda \in [0,1]` (the melange back pressure fraction)
such that

.. math::

   \pmelange = \lambda (\pice - \psw).

Then melange pressure is added to the ordinary ocean pressure so that the pressure
imbalance term scales with `\lambda`:

.. math::
   :name: eq-cfbc-3

   \int_{b}^{h}(\pice - (\psw + \pmelange))\, dz &= \int_{b}^{h}(\pice - (\psw + \lambda(\pice - \psw)))\, dz

   &= (1 - \lambda) \int_{b}^{h} (\pice - \psw)\, dz.

This formula replaces the right hand side of :eq:`eq-cfbc-uv`.

By default, `\lambda` is set to zero, but PISM implements a scalar time-dependent "melange
back pressure fraction offset" forcing in which `\lambda` can be read from a file. Please
see the :ref:`Climate Forcing Manual <sec-ocean-frac-mbp>` for details.