.. include:: ../global.txt

Calving front stress boundary condition
=======================================

.. contents::

Notation
--------

.. csv-table::
   :name: tab-cfbc-notation
   :header: Variable, Meaning

   `h`,          ice top surface elevation
   `b`,          ice bottom surface elevation
   `H = h - b`,  ice thickness
   `g`,          acceleration due to gravity
   `\viscosity`, vertically-averaged viscosity of ice
   `\n`,         normal vector
   `B(T)`,       ice hardness
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


Because we seek boundary conditions for the SSA stress balance, in which the
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
