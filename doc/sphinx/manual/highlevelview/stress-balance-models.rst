.. _sec-stress-balance-models:

Stress balance models: SIA, SSA, and the First Order Approximation
------------------------------------------------------------------

At each time-step of a typical PISM run, the geometry, temperature, and basal strength of
the ice sheet are included into stress (momentum) balance equations to determine the
velocity of the flowing ice. The "full" stress balance equations for flowing ice form a
non-Newtonian Stokes model :cite:`Fowler`. PISM does not attempt to solve the Stokes equations
themselves, however. Instead it can numerically solve, in parallel, three different shallow
approximations which are well-suited to ice sheet and ice shelf systems:

- the non-sliding shallow ice approximation (SIA) :cite:`Hutter`, also called the "lubrication
  approximation" :cite:`Fowler`, which describes ice as flowing by shear in planes parallel to
  the geoid, with a strong connection of the ice base to the bedrock, and
- the shallow shelf approximation (SSA) :cite:`WeisGreveHutter`, which describes a
  membrane-type flow of floating ice :cite:`Morland`, or of grounded ice which is sliding over
  a weak base :cite:`MacAyeal`, :cite:`SchoofStream`.
- a first order approximation to the Stokes equations due to Blatter (:cite:`Blatter`,
  :cite:`Pattyn03`). In the remainder, we refer to is as the "Blatter's model."

The SIA equations are easier to solve numerically than the SSA and Blatter's model, and
easier to parallelize, because they are local in each column of ice. Specifically, they
describe the vertical shear stress as a local function of the driving stress
:cite:`Paterson`. They can confidently be applied to those grounded parts of ice sheets
for which the basal ice is frozen to the bedrock, or which is minimally sliding, and where
the bed topography is relatively slowly-varying in the map-plane :cite:`Fowler`. These
characteristics apply to the majority (by area) of the Greenland and Antarctic ice sheets.

We solve the SIA with a non-sliding base because the traditional :cite:`Greve`,
:cite:`HuybrechtsdeWolde`, :cite:`PayneBaldwin` additions of ad hoc "sliding laws" into
the SIA stress balance, and especially schemes which "switch on" at the pressure-melting
temperature :cite:`EISMINT00`, have bad continuum :cite:`Fowler01` and numerical (see
:cite:`BBssasliding`, appendix B) modeling consequences.

The SSA equations can confidently be applied to large floating ice shelves, which have
small depth-to-width ratio and negligible basal resistance :cite:`Morland`,
:cite:`MorlandZainuddin`. The flow speeds in ice shelves are frequently an order-of-magnitude
higher than in the non-sliding, grounded parts of ice sheets.

Terrestrial ice sheets also have fast-flowing grounded parts, however, called "ice
streams" or "outlet glaciers" :cite:`TrufferEchelmeyer`. Such features appear at the
margin of, and sometimes well into the interior of, the Greenland :cite:`Joughinetal2001`
and Antarctic :cite:`BamberVaughanJoughin` ice sheets. Describing these faster-flowing
grounded parts of ice sheets requires something more than the non-sliding SIA. This is
because adjacent columns of ice which have different amounts of basal resistance exert
strong "longitudinal" or "membrane" stresses :cite:`SchoofStream` on each other.

In PISM the SSA may be used as a "sliding law" for grounded ice which is already modeled
everywhere by the non-sliding SIA :cite:`BBssasliding`, :cite:`Winkelmannetal2011`. For
grounded ice, in addition to including shear in planes parallel to the geoid, we must
balance the membrane stresses where there is sliding. This inclusion of a membrane stress
balance is especially important when there are spatial and/or temporal changes in basal
strength. This "sliding law" role for the SSA is in addition to its more obvious role in
ice shelf modeling. The SSA plays both roles in a PISM whole ice sheet model in which
there are large floating ice shelves (e.g. as in Antarctica :cite:`Golledgeetal2012ant`,
:cite:`Martinetal2011`, :cite:`Winkelmannetal2011`; see also :ref:`sec-ross`).

The "SIA+SSA hybrid" model is recommended for most whole ice sheet modeling purposes
because it seems to be a good compromise given currently-available data and computational
power. A related hybrid model described by Pollard and deConto :cite:`PollardDeConto` adds
the shear to the SSA solution in a slightly-different manner, but it confirms the success
of the hybrid concept.

By default, however, PISM does not turn on (activate) the SSA solver. This is because a
decision to solve the SSA must go with a conscious user choice about basal strength. The
user must both use a command-line option to turn on the SSA (e.g. option ``-stress_balance
ssa``; see section :ref:`sec-stressbalance`) and also make choices in input files and
runtime options about basal strength (see section :ref:`sec-basestrength`). Indeed,
uncertainties in basal strength boundary conditions usually dominate the modeling error
made by not including higher-order stresses in the balance.

When the SSA model is applied a parameterized sliding relation must be chosen. A
well-known SSA model with a linear basal resistance relation is the Siple Coast
(Antarctica) ice stream model by MacAyeal :cite:`MacAyeal`. The linear sliding law choice is
explained by supposing the saturated till is a linearly-viscous fluid. A free boundary
problem with the same SSA balance equations but a different sliding law is the Schoof
:cite:`SchoofStream` model of ice streams, using a plastic (Coulomb) sliding relation. In this
model ice streams appear where there is "till failure" :cite:`Paterson`, i.e. where the basal
shear stress exceeds the yield stress. In this model the location of ice streams is not
imposed in advance.

As noted, both the SIA and SSA models are *shallow* approximations. These equations are
derived from the Stokes equations by distinct small-parameter arguments, both based on a
small depth-to-width ratio for the ice sheet. For the small-parameter argument in the SIA
case see :cite:`Fowler`. For the corresponding SSA argument, see :cite:`WeisGreveHutter`
or the appendices of :cite:`SchoofStream`. Schoof and Hindmarsh :cite:`SchoofHindmarsh`
have analyzed the connections between these shallowest models and higher-order models,
while :cite:`GreveBlatter2009` discusses ice dynamics and stress balances comprehensively.
Note that SIA, SSA, and higher-order models all approximate the pressure as hydrostatic.

Instead of a SIA+SSA hybrid model implemented in PISM one might use the Stokes equations,
or a "higher-order" model (e.g. Blatter's model :cite:`Blatter`, :cite:`Pattyn03`), but
this immediately leads to a resolution-versus-stress-inclusion tradeoff. The amount of
computation per map-plane grid location is much higher in higher-order models, although
careful numerical analysis can generate large performance improvements for such equations
:cite:`BrownSmithAhmadia2013`.

Time-stepping solutions of the mass conservation and energy conservation equations, which
use the ice velocity for advection, can use any of the SIA or SSA or SIA+SSA hybrid stress
balances. No user action is required to turn on these conservation models. They can be
turned off by user options ``-no_mass`` (ice geometry does not evolve) or ``-energy none``
(ice enthalpy and temperature does not evolve), respectively.
