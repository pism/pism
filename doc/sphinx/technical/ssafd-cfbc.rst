.. include:: ../global.txt

Calving front stress boundary condition
=======================================

.. contents::

.. only:: html

   .. math::

      \newcommand{\diff}[2]{\frac{\partial #1}{\partial #2}}
      \newcommand{\n}{\mathbf{n}}
      \newcommand{\nx}{\n_{x}}
      \newcommand{\ny}{\n_{y}}
      \newcommand{\nz}{\n_{z}}
      \newcommand{\psw}{p_{\text{ocean}}}
      \newcommand{\pmelange}{p_{\text{melange}}}
      \newcommand{\pice}{p_{\text{ice}}}
      \newcommand{\rhoi}{\rho_{\text{ice}}}
      \newcommand{\rhosw}{\rho_{\text{ocean}}}
      \newcommand{\zs}{z_{\text{s}}}
      \newcommand{\td}[1]{t^{D}_{#1}}
      \newcommand{\D}{\displaystyle}
      \newcommand{\dx}{\Delta x}
      \newcommand{\dy}{\Delta y}
      \newcommand{\viscosity}{\nu}
      \newcommand{\partI}{(2\tilde N_{xx} + \tilde N_{yy})}
      \newcommand{\partII}{(\tilde N_{xy})}

.. raw:: latex

   \providecommand{\diff}[2]{\frac{\partial #1}{\partial #2}}
   \providecommand{\n}{\mathbf{n}}
   \providecommand{\nx}{\n_{x}}
   \providecommand{\ny}{\n_{y}}
   \providecommand{\nz}{\n_{z}}
   \providecommand{\psw}{p_{\text{ocean}}}
   \providecommand{\pmelange}{p_{\text{melange}}}
   \providecommand{\pice}{p_{\text{ice}}}
   \providecommand{\rhoi}{\rho_{\text{ice}}}
   \providecommand{\rhosw}{\rho_{\text{ocean}}}
   \providecommand{\zs}{z_{\text{s}}}
   \providecommand{\td}[1]{t^{D}_{#1}}
   \providecommand{\D}{\displaystyle}
   \providecommand{\dx}{\Delta x}
   \providecommand{\dy}{\Delta y}
   \providecommand{\viscosity}{\nu}
   \providecommand{\partI}{(2\tilde N_{xx} + \tilde N_{yy})}
   \providecommand{\partII}{(\tilde N_{xy})}

Notation
--------

.. csv-table::
   :name: tab-cfbc-notation
   :header: Variable, Meaning

   `h`,          ice top surface elevation
   `b`,          ice bottom surface elevation
   `H = h - b`,  ice thickness
   `\zs`,        sea level elevation
   `g`,          acceleration due to gravity
   `\rhoi`,      ice density
   `\rhosw`,     sea water density
   `\viscosity`, vertically-averaged viscosity of ice
   `\n`,         normal vector
   `B(T)`,       ice hardness
   `\pice`,      pressure due to the weight of a column of ice
   `\psw`,       pressure due to the weight of a column of seawater
   `D`,          strain rate tensor
   `d_{e}`,      effective strain rate
   `t`,          Cauchy stress tensor
   `t^{D}`,      deviatoric stress tensor; note `\td{ij} = t_{ij} + p \delta_{ij}`

Calving front stress boundary condition
---------------------------------------

In the 3D case the calving front stress boundary condition (:cite:`GreveBlatter2009`, equation
(6.19)) reads

.. math::

   \left.t\right|_{\text{cf}} \cdot \n = -\psw \n.

Expanded in component form, and evaluating the left-hand side at the calving front and
assuming that the calving front face is vertical (`\nz = 0`), this gives

.. math::

   \begin{array}{rcrcl}
     (\td{xx} - p) \nx &+& \td{xy} \ny &=& -\psw \nx,\\
     \td{xy} \nx &+& (\td{yy} - p) \ny &=& -\psw  \ny,\\
     \td{xz} \nx &+& \td{yz} \ny &=& 0.
   \end{array}


Because we see boundary conditions for the SSA stress balance, in which the
vertically-integrated forms of the stresses `\td{xx},\td{xy},\td{yy}` are balanced
separately from the shear stresses `\td{xz},\td{yz}`, the third of the above equations can
be omitted from the remaining analysis.

Let `\pice=\rhoi g (h-z)`. In the hydrostatic approximation, `t_{zz}=-\pice`
(:cite:`GreveBlatter2009`, equation (5.59)). Next, since `\td{}` has trace zero,

.. math::

   p &= p - \td{xx} - \td{yy} - \td{zz}

   &= - t_{zz} - \td{xx} - \td{yy}

   &= \pice - \td{xx} - \td{yy}.

Thus

.. math::
   :label: ssafd-cfbc-1

   \begin{array}{rcrcl}
     (2\td{xx} + \td{yy}) \nx &+& \td{xy} \ny &=& (\pice - \psw) \nx,\\
     \td{xy} \nx &+& (2\td{yy} + \td{xx}) \ny &=& (\pice - \psw) \ny.\\
   \end{array}

Now, using the "viscosity form" of the flow law

.. math::

   \td{} &= 2\viscosity\, D,

   \viscosity &= \frac 12 B(T) d_{e}^{1/n-1}

and the fact that in the shallow shelf approximation `u` and `v` are depth-independent, we
can define the membrane stress `N` as the vertically-integrated deviatoric stress

.. math::

   N_{ij} = \int_{b}^{h} t^{D}_{ij} dz = 2\, \bar \viscosity\, H\, D_{ij}.

Here `\bar \viscosity` is the vertically-averaged effective viscosity.

Integrating :eq:`ssafd-cfbc-1` with respect to `z`, we get

.. math::
   :label: ssafd-cfbc-3

   \begin{array}{rcrcl}
     (2N_{xx} + N_{yy}) \nx &+& N_{xy} \ny &=& \displaystyle \int_{b}^{h}(\pice - \psw) dz\, \nx,\\
     N_{xy} \nx &+& (2N_{yy} + N_{xx}) \ny &=& \displaystyle \int_{b}^{h}(\pice - \psw) dz\, \ny.
   \end{array}


Shallow shelf approximation
---------------------------

The shallow shelf approximation written in terms of `N_{ij}` is

.. math::
   :label: ssafd-cfbc-2

   \begin{array}{rcrcl}
     \D \diff{}{x}(2N_{xx} + N_{yy}) &+& \D \diff{N_{xy}}{y} &=& \D \rho g H \diff{h}{x},\\
     \D \diff{N_{xy}}{x} &+& \D \diff{}{y}(2N_{yy} + N_{xx}) &=& \D \rho g H \diff{h}{y}.
   \end{array}

Implementing the calving front boundary condition
-------------------------------------------------

We use centered finite difference approximations in the discretization of the SSA
:eq:`ssafd-cfbc-2`. Consider the first equation:

.. math::
   :label: ssafd-cfbc-4

   \diff{}{x}(2N_{xx} + N_{yy}) + \D \diff{N_{xy}}{y} = \D \rho g H \diff{h}{x}.

Let `\tilde F` be an approximation of `F` using a finite-difference scheme. Then the first
SSA equation is approximated by

.. math::

   \frac1{\dx}\left(\partI_{i+\frac12,j} - \partI_{i-\frac12,j}\right) +
   \frac1{\dy}\left(\partII_{i,j+\frac12} - \partII_{i,j-\frac12}\right) =
   \rho g H \frac{h_{i+\frac12,j} - h_{i-\frac12,j}}{\dx}.


Now, assume that the cell boundary (face) at `i+\frac12,j` is at the
calving front. Then `\n = (1,0)` and from :eq:`ssafd-cfbc-3` we have

.. math::
   :label: ssafd-cfbc-vertintbdry

   2N_{xx} + N_{yy} = \int_{b}^{h}(\pice - \psw) dz.

We call the right-hand side of :eq:`ssafd-cfbc-vertintbdry` the "pressure difference term."

In forming the matrix approximation of the SSA :cite:`BBssasliding`,
:cite:`Winkelmannetal2011`, instead of assembling a matrix row corresponding to the
generic equation :eq:`ssafd-cfbc-4` we use

.. math::

   \frac1{\dx}\left(\left[\int_{b}^{h}(\pice - \psw) dz\right]_{i+\frac12,j} - \partI_{i-\frac12,j}\right) +
   \frac1{\dy}\left(\partII_{i,j+\frac12} - \partII_{i,j-\frac12}\right) =
   \rho g H \frac{h_{i+\frac12,j} - h_{i-\frac12,j}}{\dx}.

Moving terms that do not depend on the velocity field to the
right-hand side, we get

.. math::

   \frac1{\dx}\left(- \partI_{i-\frac12,j}\right) +
   \frac1{\dy}\left(\partII_{i,j+\frac12} - \partII_{i,j-\frac12}\right) =
   \rho g H \frac{h_{i+\frac12,j} - h_{i-\frac12,j}}{\dx} + \left[\frac{\int_{b}^{h}(\psw - \pice) dz}{\dx}\right]_{i+\frac12,j}.


The second equation and other cases (`\n = (-1,0)`, etc) are treated similarly.

Evaluating the "pressure difference term"
-----------------------------------------

For `z \in [b, h]` the modeled pressures are

.. math::

   \begin{array}{ll}
   \psw &=
     \begin{cases}
       0, & z > \zs,\\
       \rhosw\, g (\zs - z), & z \le \zs,
     \end{cases}\\
   \pice &= \rhoi\, g (h - z).
   \end{array}

Depending on the local geometry `b` is either prescribed (grounded case) or is a function
of the ice thickness and sea level elevation.

Floating case
^^^^^^^^^^^^^

Using the flotation thickness relation `\rhoi H = \rhosw (\zs - b)`, which applies to
floating ice given that `b` here denotes the base elevation of the floating ice, we have

.. math::

   \int_{b}^{h}(\pice - \psw) dz =
   \frac{1}{2}\, g\, \rhoi \left(H^2 - \frac{\rhoi}{\rhosw} H^2\right)

Grounded below sea level
^^^^^^^^^^^^^^^^^^^^^^^^

Because `b` denotes both the base elevation of the grounded ice, and the bedrock
elevation, here `\rhoi H \ge \rhosw (\zs - b)`. The integral simplifies to

.. math::

   \int_{b}^{h}(\pice - \psw) dz =
   \frac{1}{2}\, g\, \rhoi \left(H^2-\frac{\rhosw}{\rhoi}\, \left(\zs-b\right)^2\right)


This, in fact, applies in both floating and grounded-below-sea-level cases.

Grounded above sea level
^^^^^^^^^^^^^^^^^^^^^^^^

In this case `\psw = 0`, so

.. math::

   \int_{b}^{h}(\pice - \psw) dz = \frac{1}{2}\, g\, \rhoi\, H^2


Modeling melange back-pressure
------------------------------

Let `\pmelange` be the additional melange back-pressure. Then `\psw \le \psw + \pmelange
\le \pice`. Put another way,

.. math::
   :label: ssafd-cfbc-13

   0 \le \pmelange \le \pice - \psw.

Let `\lambda` be the "melange back-pressure fraction" (or "relative melange pressure")
ranging from `0` to `1`, so that

.. math::
   :label: ssafd-cfbc-14

   \pmelange = \lambda \cdot (\pice - \psw).

Then the modified pressure difference term is

.. math::
   :label: ssafd-cfbc-15

   \int_{b}^{h}(\pice - (\psw + \pmelange)) dz &= \int_{b}^{h}(\pice - (\psw + \lambda(\pice - \psw)))\, dz

   &= (1 - \lambda) \int_{b}^{h} (\pice - \psw)\, dz.
