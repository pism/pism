.. include:: ../../global.txt

.. _sec-basicmodes:

Evolutionary versus diagnostic modeling
---------------------------------------

The main goal of a numerical ice sheet model like PISM is to be a dynamical system which
evolves as similarly as possible to the modeled ice sheet. Such a goal assumes one has the
"right" climate inputs and parameter choices at each time step. It also assumes one has
the "right" initial conditions, such as an adequate description of the present state of
the ice sheet, but this assumption is rarely satisfied. Instead a variety of heuristics
must be used to minimally-initialize an ice sheet model. For options associated to
establishing mathematical initial conditions when first starting PISM, see section
:ref:`sec-initboot`.

Inside PISM are evolution-in-time partial differential equations which are solved by
taking small time steps. "Small" may vary from thousandths to tens of model years, in
practice, depending primarily on grid resolution, but also on modeled ice geometry and
flow speed. Time steps are chosen adaptively in PISM, according to the stability criteria
of the combined numerical methods :cite:`BBssasliding`, :cite:`BBL`.

However, especially for ice streams and shelves, non-time-stepping "diagnostic" solution
of the stress balance partial differential equations might be the desired computation, and
PISM can also produce such "diagnostic" velocity fields. Such computations necessarily
assume that the ice geometry, viscosity, and boundary stresses are known. Because of the
slowness of the ice, in the sense that inertia can be neglected in the stress balance
:cite:`Fowler`, such computations can determine the ice velocity.

Sections :ref:`sec-start` and :ref:`sec-ross` give examples illustrating evolutionary and
diagnostic modes of PISM, respectively. The first describes time-stepping evolution models
for the Greenland ice sheet, while the second describes a diagnostic SSA model for the
Ross ice shelf.

